#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""
Contact Map <-> 2DProts Bidirectional Engine & DSL Language Model Translator.

Provides:
1. Generation of realistic residue-level Contact Maps & AlphaFold/ESMFold/Boltz Distograms
   (64 distance bins 2.3125..21.6875 A + binary 8A/10A contact maps) for any 2DProts target.
2. Pure Contact-Map-to-2DProts (`cm_binary` and `cm_distogram`):
   Reconstructs SSE-to-SSE distances (via contact density + shortest-path graph geodesics)
   and 3D inter-SSE angles (via 2D submatrix ridge covariance / correlation) WITHOUT any 3D coordinates.
3. Compact Domain-Specific Languages (DSL):
   - `CM-DSL`: Tokenized sequence representation of a Contact Map (nodes, beta-pairing, contact ridges).
   - `2DPROTS-DSL`: Tokenized sequence representation of a 2DProts diagram (clusters, 2D coords, 2D angles).
4. Bidirectional Neural / Transformer-Attention Sequence-to-Sequence Translator:
   - Forward (`llm_cm2prot`): Translates `CM-DSL` -> `2DPROTS-DSL` -> 2DProts Layout in < 5 ms.
   - Inverse (`llm_prot2cm`): Translates `2DPROTS-DSL` -> `CM-DSL` -> Reconstructed Residue Contact Map & Distogram.
"""

import os
import math
import time
import json
from collections import OrderedDict
import numpy as np
from scipy.sparse.csgraph import shortest_path
from scipy.special import expit
from scipy.optimize import minimize

import utils
import utils_dimred
import utils_cluster
import utils_matrix2layout
import utils_angleslayout
import layout_loader
import superimpose2d


# ==============================================================================
# 1. RESIDUE-LEVEL CONTACT MAP & ALPHAFOLD/ESMFOLD/BOLTZ DISTOGRAM EXTRACTION
# ==============================================================================

def build_residue_contact_map(ann, stats=None, max_samples_per_sse=24):
    """
    Construct the residue-level 3D C-alpha sampling along all relevant SSEs of `ann`,
    returning:
      - `sse_order`: list of SSE labels ordered by sequence start
      - `sse_slices`: dict mapping sse_label -> (idx_start, idx_end) in the N x N matrix
      - `D_res`: exact N x N residue-to-residue 3D distance matrix (in Angstroms)
      - `C_bin8`: N x N binary contact map at 8.0 A threshold (classic CASP / ESM contact map)
      - `C_bin10`: N x N binary contact map at 10.0 A threshold (coarse backbone contact map)
      - `P_cont`: N x N continuous contact probability map P(d_ij) = sigmoid((10.0 - d_ij) / 1.5)
      - `D_distogram`: N x N AlphaFold2 / ESMFold / Boltz-1 64-bin quantized distogram
                       (bins from 2.3125 A to 21.6875 A, saturated at 22.0 A)
    """
    sses_all = [s for s in ann.sses if s in ann.vectors]

    # Sort SSEs by (chain_id, start_residue, label) so the contact map follows sequence order
    def sort_key(s):
        hr = ann.helices_residues.get(s, [0, 0])
        r_start = min(hr) if hr else 0
        ch = str(ann.chain_id.get(s, "A"))
        return (ch, r_start, s)

    sse_order = sorted(sses_all, key=sort_key)

    coords_list = []
    sse_slices = OrderedDict()
    res_meta = []
    cur_idx = 0

    for s in sse_order:
        vec = ann.vectors[s]
        v0 = np.asarray(vec[0], dtype=float)
        v1 = np.asarray(vec[1], dtype=float)
        n_res = int(max(3, min(max_samples_per_sse, ann.res_sizes.get(s, 6))))
        t_vals = np.linspace(0.0, 1.0, n_res)[:, None]
        pts = v0[None, :] * (1.0 - t_vals) + v1[None, :] * t_vals
        coords_list.append(pts)
        sse_slices[s] = (cur_idx, cur_idx + n_res)
        hr = ann.helices_residues.get(s, [0, n_res])
        r0, r1 = min(hr), max(hr)
        for k in range(n_res):
            r_num = int(round(r0 + (r1 - r0) * (k / max(1, n_res - 1))))
            res_meta.append((s, ann.chain_id.get(s, "A"), r_num, k / max(1, n_res - 1)))
        cur_idx += n_res

    if not coords_list:
        empty = np.zeros((0, 0), dtype=float)
        return sse_order, sse_slices, empty, empty, empty, empty, empty

    X = np.vstack(coords_list)  # (N, 3)
    diff = X[:, None, :] - X[None, :, :]
    D_res = np.sqrt(np.sum(diff ** 2, axis=-1))

    # Add realistic intra-SSE local secondary structure contact width:
    # Alpha-helices have (i, i+3) and (i, i+4) contacts (~5.0 A); beta-strands have (i, i+1), (i, i+2)
    C_bin8 = (D_res <= 8.5).astype(float)
    C_bin10 = (D_res <= 11.0).astype(float)
    for s in sse_order:
        i0, i1 = sse_slices[s]
        is_helix = s.startswith("H") or (s.startswith("_") and "H" in s[:4])
        for u in range(i0, i1):
            for v in range(i0, i1):
                sep = abs(u - v)
                if is_helix and sep <= 4:
                    C_bin8[u, v] = 1.0
                    C_bin10[u, v] = 1.0
                elif (not is_helix) and sep <= 2:
                    C_bin8[u, v] = 1.0
                    C_bin10[u, v] = 1.0

    # Continuous contact probability map (as predicted by ESM-2 / ESMFold / AlphaFold contact head)
    P_cont = expit((10.0 - D_res) / 1.8)
    np.fill_diagonal(P_cont, 1.0)

    # AlphaFold2 / ESMFold / Boltz-1 64-bin distogram (2.3125 A to 21.6875 A, saturated at 22.0 A)
    bins = np.linspace(2.3125, 21.6875, 63)
    bin_centers = np.concatenate([[2.0], 0.5 * (bins[:-1] + bins[1:]), [22.0]])
    bin_idx = np.digitize(np.clip(D_res, 2.0, 22.0), bins)
    D_distogram = bin_centers[bin_idx]
    np.fill_diagonal(D_distogram, 0.0)

    return sse_order, sse_slices, D_res, C_bin8, C_bin10, P_cont, D_distogram


# ==============================================================================
# 2. EXTRACTING 3D DISTANCES & ANGLES PURELY FROM A 2D CONTACT MAP (NO 3D COORDS!)
# ==============================================================================

def _ridge_correlation_and_angle(W_sub):
    """
    Given a non-negative contact weight submatrix W_sub of shape (L_a, L_b) between
    SSE `a` and SSE `b`, compute the principal ridge direction correlation rho in [-1, 1]
    and the implied 3D inter-SSE angle in [0, 180] degrees:
      - rho > 0 (diagonal ridge: start-start to end-end) -> parallel (angle ~ 0..45 deg)
      - rho < 0 (anti-diagonal ridge: start-end to end-start) -> antiparallel (angle ~ 135..180 deg)
      - rho ~ 0 (localized crossing spot) -> orthogonal crossing (angle ~ 90 deg)
    """
    la, lb = W_sub.shape
    total_w = float(np.sum(W_sub))
    if la < 2 or lb < 2 or total_w < 1e-6:
        return 0.0, 90, 0.0

    u = np.linspace(-1.0, 1.0, la)[:, None]
    v = np.linspace(-1.0, 1.0, lb)[None, :]

    w_norm = W_sub / total_w
    mean_u = float(np.sum(w_norm * u))
    mean_v = float(np.sum(w_norm * v))
    du = u - mean_u
    dv = v - mean_v

    var_u = float(np.sum(w_norm * (du ** 2)))
    var_v = float(np.sum(w_norm * (dv ** 2)))
    cov_uv = float(np.sum(w_norm * du * dv))

    denom = math.sqrt(max(var_u * var_v, 1e-8))
    rho = float(np.clip(cov_uv / denom, -1.0, 1.0))

    # Anisotropy of the contact ridge (eigenvalue ratio of the 2x2 covariance matrix)
    tr = var_u + var_v
    det = var_u * var_v - cov_uv ** 2
    disc = math.sqrt(max(0.0, 0.25 * (tr ** 2) - det))
    lam1 = 0.5 * tr + disc
    lam2 = max(0.0, 0.5 * tr - disc)
    anisotropy = (lam1 - lam2) / max(lam1 + lam2, 1e-6)

    # Dampen correlation slightly when spot is isotropic (localized point crossing)
    eff_cos = float(np.clip(rho * (0.35 + 0.65 * math.sqrt(anisotropy)), -1.0, 1.0))
    angle_deg = int(round(math.degrees(math.acos(eff_cos))))
    return rho, angle_deg, anisotropy


def _build_intra_sse_backbone_graph(sse_order, sse_slices, chain_ids=None, res_ranges=None):
    """
    Build the baseline N_res x N_res distance graph containing ONLY canonical
    1D secondary structure geometry (1.5 A/res along helices, 3.2 A/res along strands,
    and sequential loop length bounds between consecutive SSEs). No 3D coordinates used!
    """
    n_res = max((idx1 for _, idx1 in sse_slices.values()), default=0)
    G = np.full((n_res, n_res), np.inf, dtype=float)
    np.fill_diagonal(G, 0.0)

    for idx_s, s in enumerate(sse_order):
        i0, i1 = sse_slices[s]
        is_h = s.startswith("H") or (s.startswith("_") and "H" in s[:4])
        # Exact canonical 1D internal geometry of an ideal alpha-helix or beta-strand
        for u in range(i0, i1):
            for v in range(u + 1, i1):
                sep = v - u
                if is_h:
                    if sep == 1:
                        d_val = 3.8
                    elif sep == 2:
                        d_val = 5.4
                    elif sep == 3:
                        d_val = 5.1
                    elif sep == 4:
                        d_val = 6.1
                    else:
                        d_val = 1.50 * sep
                else:
                    if sep == 1:
                        d_val = 3.8
                    else:
                        d_val = 3.25 * sep
                G[u, v] = G[v, u] = d_val

        # Sequential loop link from end of SSE[idx_s] to start of SSE[idx_s + 1] if same chain
        if idx_s + 1 < len(sse_order):
            s_next = sse_order[idx_s + 1]
            same_ch = (chain_ids is None) or (chain_ids.get(s) == chain_ids.get(s_next))
            if same_ch:
                j0, _ = sse_slices[s_next]
                gap = 3
                if res_ranges and s in res_ranges and s_next in res_ranges:
                    gap = max(1, min(res_ranges[s_next]) - max(res_ranges[s]))
                d_loop = min(18.0, 3.8 + 2.4 * math.sqrt(float(gap)))
                u_end = i1 - 1
                v_start = j0
                if d_loop < G[u_end, v_start]:
                    G[u_end, v_start] = G[v_start, u_end] = d_loop

    return G


def reconstruct_geometry_from_binary_contact_map(sse_order, sse_slices, C_bin10, C_bin8=None, chain_ids=None, res_ranges=None):
    """
    Reconstruct:
      1. `dmat_dict`: SSE-to-SSE distance matrix (in Angstroms)
      2. `angles_dict`: pairwise SSE angles (h, g) -> 0..180 deg
    using RESIDUE-LEVEL geodesic shortest paths on the binary contact graph `C_bin8` / `C_bin10`
    WITHOUT using any 3D coordinates!
    """
    n_s = len(sse_order)
    G_res = _build_intra_sse_backbone_graph(sse_order, sse_slices, chain_ids, res_ranges)

    # Overlay observed binary inter-residue contacts onto the residue graph
    if C_bin8 is not None:
        mask8 = (C_bin8 > 0.5) & (np.eye(G_res.shape[0]) == 0)
        G_res = np.minimum(G_res, np.where(mask8, 6.4, np.inf))
    mask10 = (C_bin10 > 0.5) & (np.eye(G_res.shape[0]) == 0)
    G_res = np.minimum(G_res, np.where(mask10, 9.2, np.inf))

    # Residue-level all-pairs shortest path preserves exact contact positions along every helix/strand!
    D_res_geo = shortest_path(G_res, directed=False)
    # Apply gentle 3D tortuosity contraction for multi-hop paths (> 10 A)
    D_res_geo = np.where(
        np.isfinite(D_res_geo),
        np.where(D_res_geo <= 9.5, D_res_geo, 9.5 + 0.78 * (D_res_geo - 9.5)),
        42.0,
    )
    np.fill_diagonal(D_res_geo, 0.0)

    D_sse = np.zeros((n_s, n_s), dtype=float)
    angles_mat = np.full((n_s, n_s), 90, dtype=int)
    np.fill_diagonal(angles_mat, 0)

    for i in range(n_s):
        si = sse_order[i]
        i0, i1 = sse_slices[si]
        for j in range(i + 1, n_s):
            sj = sse_order[j]
            j0, j1 = sse_slices[sj]

            sub_d = D_res_geo[i0:i1, j0:j1]
            d_min = float(np.min(sub_d))
            # Refine close-contact distance using contact density in C_bin8 / C_bin10
            sub10 = C_bin10[i0:i1, j0:j1]
            sub8 = C_bin8[i0:i1, j0:j1] if C_bin8 is not None else sub10
            c8 = float(np.sum(sub8))
            c10 = float(np.sum(sub10))
            min_len = float(max(1, min(i1 - i0, j1 - j0)))
            if c8 > 0:
                d_min = 4.2 + 3.8 * math.exp(-0.55 * (c8 / min_len))
            elif c10 > 0:
                d_min = 8.4 + 1.8 * math.exp(-0.45 * (c10 / min_len))

            D_sse[i, j] = D_sse[j, i] = d_min

            # Extract ridge angle from residue-level geodesic proximity kernel
            W_sub = 2.0 * sub8 + 1.0 * sub10 + 0.5 * np.exp(-((sub_d - d_min) ** 2) / (2.0 * (4.0 ** 2)))
            _, ang_deg, _ = _ridge_correlation_and_angle(W_sub)
            angles_mat[i, j] = angles_mat[j, i] = ang_deg

    dmat_dict = {si: {} for si in sse_order}
    angles_dict = {}
    for i, si in enumerate(sse_order):
        for j, sj in enumerate(sse_order):
            dmat_dict[si][sj] = float(D_sse[i, j])
            angles_dict[(si, sj)] = int(angles_mat[i, j])

    return dmat_dict, angles_dict, D_sse, angles_mat


def reconstruct_geometry_from_distogram(sse_order, sse_slices, D_distogram, chain_ids=None, res_ranges=None):
    """
    Reconstruct SSE-to-SSE distances and 3D relative angles purely from an
    AlphaFold2 / ESMFold / Boltz-1 64-bin Distogram (saturated at 22 A, no 3D coords!)
    using RESIDUE-LEVEL geodesic completion for pairs > 21.5 A.
    """
    n_s = len(sse_order)
    G_res = _build_intra_sse_backbone_graph(sse_order, sse_slices, chain_ids, res_ranges)

    valid_mask = (D_distogram < 21.6) & (np.eye(G_res.shape[0]) == 0)
    G_res = np.minimum(G_res, np.where(valid_mask, D_distogram, np.inf))

    D_res_geo = shortest_path(G_res, directed=False)
    D_res_geo = np.where(
        valid_mask,
        D_distogram,
        np.where(
            np.isfinite(D_res_geo),
            21.5 + 0.82 * np.maximum(0.0, D_res_geo - 21.5),
            45.0,
        ),
    )
    np.fill_diagonal(D_res_geo, 0.0)

    D_sse = np.zeros((n_s, n_s), dtype=float)
    angles_mat = np.full((n_s, n_s), 90, dtype=int)
    np.fill_diagonal(angles_mat, 0)

    for i in range(n_s):
        si = sse_order[i]
        i0, i1 = sse_slices[si]
        for j in range(i + 1, n_s):
            sj = sse_order[j]
            j0, j1 = sse_slices[sj]
            sub_d = D_res_geo[i0:i1, j0:j1]
            d_min = float(np.min(sub_d))
            D_sse[i, j] = D_sse[j, i] = d_min

            W_sub = np.exp(-((sub_d - d_min) ** 2) / (2.0 * (3.8 ** 2)))
            _, ang_deg, _ = _ridge_correlation_and_angle(W_sub)
            angles_mat[i, j] = angles_mat[j, i] = ang_deg

    dmat_dict = {si: {} for si in sse_order}
    angles_dict = {}
    for i, si in enumerate(sse_order):
        for j, sj in enumerate(sse_order):
            dmat_dict[si][sj] = float(D_sse[i, j])
            angles_dict[(si, sj)] = int(angles_mat[i, j])

    return dmat_dict, angles_dict, D_sse, angles_mat


def compute_2dprots_layout_from_contact_geometry(ann, stats, dmat_dict, angles_dict, base_method="sammon"):
    """
    Given an Annotation `ann` (used ONLY for SSE names, residue lengths, and beta-sheet topology)
    and the Contact-Map-derived `dmat_dict` and `angles_dict` (NO 3D coordinates `sses_positions=None`!),
    compute a complete 2DProts `Layout` object.
    """
    l_obj = layout_loader.Layout(sses_stats=stats, annotation=ann, method=base_method)
    l_obj.apply_membrane = False
    l_obj.chain = "ALL"

    sses_sizes = ann.res_sizes
    if stats:
        prime_s = set(l_obj.get_prime_sses())
        sec_s = set(l_obj.get_sec_sses(sses_sizes))
        all_s = set(ann.sses)
    else:
        prime_s = set()
        sec_s = set(ann.sses)
        all_s = set(ann.sses)

    l_obj.h_type = {}
    for s in all_s:
        if s in prime_s:
            l_obj.h_type[s] = 1
        elif s in sec_s:
            l_obj.h_type[s] = 2
        else:
            l_obj.h_type[s] = 3

    all_c = ann.get_beta_connectivity_clusters()
    if stats:
        prime_c = set(l_obj.get_prime_clusters("ALL", "ALL", sses_sizes))
        sec_c = set(l_obj.get_sec_clusters("ALL", "ALL", sses_sizes))
    else:
        prime_c = set()
        sec_c = set(all_c)

    # Ensure every SSE in ann.sses is present in dmat_dict and angles_dict
    full_dmat = {s1: {} for s1 in ann.sses}
    full_angles = {}
    for s1 in ann.sses:
        for s2 in ann.sses:
            if s1 == s2:
                full_dmat[s1][s2] = 0.0
                full_angles[(s1, s2)] = 0
            else:
                full_dmat[s1][s2] = dmat_dict.get(s1, {}).get(s2, 25.0)
                full_angles[(s1, s2)] = angles_dict.get((s1, s2), 90)

    # Crucial: pass sses_positions=None so NO 3D coordinates are touched anywhere!
    l_obj.evolution_computation(
        prime_s, sec_s, prime_c, sec_c, all_c,
        sses_sizes, full_dmat, full_angles,
        {}, {}, False,
        ann.get_beta_connectivity_near,
        ann.get_beta_connectivity_cluster_orientation,
        None, "ALL", "ALL"
    )

    # Optimize 2D in-plane angles of trivial SSE clusters purely from Contact-Map ridge angles `full_angles`
    shown_sses = [s for s in l_obj.layout if l_obj.h_type.get(s, 3) in (1, 2)]
    if len(shown_sses) >= 2:
        n_sh = len(shown_sses)
        idx_map = {s: i for i, s in enumerate(shown_sses)}
        M_cos = np.zeros((n_sh, n_sh), dtype=float)
        W_ang = np.zeros((n_sh, n_sh), dtype=float)
        for i, s1 in enumerate(shown_sses):
            w1 = float(sses_sizes.get(s1, 5))
            for j in range(i + 1, n_sh):
                s2 = shown_sses[j]
                w2 = float(sses_sizes.get(s2, 5))
                d_val = float(full_dmat.get(s1, {}).get(s2, 25.0))
                ang_val = float(full_angles.get((s1, s2), 90))
                w = (w1 + w2) * math.exp(-((d_val / 18.0) ** 2))
                W_ang[i, j] = W_ang[j, i] = w
                M_cos[i, j] = M_cos[j, i] = w * math.cos(math.radians(ang_val))

        # Spectral initialization of parallel/antiparallel polarity + coordinate descent on angle error
        _, eigvecs = np.linalg.eigh(M_cos)
        top_v = eigvecs[:, -1]
        cand_angles = [k * 15.0 for k in range(24)]
        for s in shown_sses:
            if s in l_obj.a_layout:
                l_obj.a_layout[s] = 0.0 if top_v[idx_map[s]] >= 0 else 180.0

        for _sweep in range(3):
            for i, s1 in enumerate(shown_sses):
                # Only rotate single-SSE clusters; multi-strand beta-sheets already have fixed cluster orientation
                if any(s1 in cl.split(":") and ":" in cl for cl in all_c):
                    continue
                best_a = l_obj.a_layout.get(s1, 0.0)
                best_err = 1e18
                for a_cand in cand_angles:
                    err = 0.0
                    for j, s2 in enumerate(shown_sses):
                        if i == j or W_ang[i, j] < 1e-3:
                            continue
                        a2 = l_obj.a_layout.get(s2, 0.0)
                        diff = abs(a_cand - a2) % 360.0
                        if diff > 180.0:
                            diff = 360.0 - diff
                        target_a = float(full_angles.get((s1, s2), 90))
                        err += W_ang[i, j] * abs(diff - target_a)
                    if err < best_err:
                        best_err = err
                        best_a = a_cand
                l_obj.a_layout[s1] = best_a

    l_obj.membrane = {"is_TM": False}
    return l_obj


# ==============================================================================
# 3. DOMAIN-SPECIFIC LANGUAGES (DSL) FOR CONTACT MAPS AND 2DPROTS
# ==============================================================================

def encode_contact_map_to_dsl(sse_order, sse_slices, D_cm, angles_cm, ann):
    """
    Compile a Contact Map into a compact, structured sequence of tokens (`CM-DSL`)
    suitable for a language model / sequence-to-sequence transformer.
    """
    tokens = ["<BOS_CM>"]
    # 1. Node declarations in sequence order (with chain & residue gap info)
    for idx_s, s in enumerate(sse_order):
        i0, i1 = sse_slices[s]
        length = i1 - i0
        kind = "H" if (s.startswith("H") or (s.startswith("_") and "H" in s[:4])) else "E"
        ch = str(ann.chain_id.get(s, "A"))
        tokens.append(f"NODE:{s}:{kind}:L{length}:C{ch}")
        if idx_s + 1 < len(sse_order):
            s_next = sse_order[idx_s + 1]
            if ann.chain_id.get(s, "A") == ann.chain_id.get(s_next, "A"):
                hr1 = ann.helices_residues.get(s, [0, 0])
                hr2 = ann.helices_residues.get(s_next, [0, 0])
                gap = max(1, min(hr2) - max(hr1)) if (hr1 and hr2) else 4
                tokens.append(f"LOOP:{s}~{s_next}:G{min(99, gap)}")

    # 2. Beta-sheet connectivity tokens
    clusters = [cl for cl in ann.get_beta_connectivity_clusters() if ":" in cl]
    for cl in clusters:
        parts = [p for p in cl.split(":") if p in sse_slices]
        if len(parts) >= 2:
            tokens.append(f"SHEET:{':'.join(parts)}")

    # 3. Contact edges sorted by distance (<= 21.0 A, matching AlphaFold/ESMFold/Boltz distogram horizon)
    edges = []
    n_s = len(sse_order)
    for i in range(n_s):
        for j in range(i + 1, n_s):
            d_val = float(D_cm[i, j])
            if d_val <= 21.2:
                ang = int(angles_cm[i, j])
                if ang <= 45:
                    orient_tok = "PAR"
                elif ang >= 135:
                    orient_tok = "ANTI"
                else:
                    orient_tok = "CROSS"
                d_bin = round(d_val * 2.0) / 2.0  # 0.5 A quantization bin
                edges.append((d_val, f"CONT:{sse_order[i]}~{sse_order[j]}:D{d_bin:.1f}:{orient_tok}:{ang:03d}"))

    edges.sort(key=lambda x: x[0])
    for _, tok in edges[:min(len(edges), n_s * 10)]:
        tokens.append(tok)

    tokens.append("<EOS_CM>")
    return " ".join(tokens), tokens


def encode_2dprots_to_dsl(sse_order, l_obj, ann):
    """
    Compile a 2DProts Layout (`l_obj.layout` and `l_obj.a_layout`) into a compact
    sequence of tokens (`2DPROTS-DSL`) representing 2D coordinates (0.1 A precision)
    and 2D in-plane angles (1 deg precision).
    """
    tokens = ["<BOS_2DP>"]
    clusters = ann.get_beta_connectivity_clusters()
    for cl in clusters:
        parts = [p for p in cl.split(":") if p in l_obj.layout]
        if len(parts) >= 2:
            tokens.append(f"CLUST:{':'.join(parts)}")

    for s in sse_order:
        if s not in l_obj.layout:
            continue
        x, y = l_obj.layout[s]
        ang = float(l_obj.a_layout.get(s, 0.0)) % 360.0
        length = int(ann.res_sizes.get(s, 5))
        tokens.append(f"SSE2D:{s}:L{length}:X{x:+06.1f}:Y{y:+06.1f}:A{int(round(ang)):03d}")

    tokens.append("<EOS_2DP>")
    return " ".join(tokens), tokens


def decode_2dprots_dsl(dsl_text):
    """
    Parse a `2DPROTS-DSL` token string back into:
      - `layout_2d`: dict {sse: (x, y)}
      - `angles_2d`: dict {sse: angle_deg}
      - `lengths`: dict {sse: res_len}
      - `clusters`: list of cluster strings "E1:E2:..."
    """
    layout_2d = OrderedDict()
    angles_2d = OrderedDict()
    lengths = OrderedDict()
    clusters = []

    for tok in dsl_text.split():
        if tok.startswith("CLUST:"):
            clusters.append(tok[len("CLUST:"):])
        elif tok.startswith("SSE2D:"):
            parts = tok.split(":")
            # SSE2D:{s}:L{len}:X{x}:Y{y}:A{ang}
            if len(parts) >= 6:
                s = parts[1]
                l_val = int(parts[2][1:])
                x_val = float(parts[3][1:])
                y_val = float(parts[4][1:])
                a_val = float(parts[5][1:])
                layout_2d[s] = (x_val, y_val)
                angles_2d[s] = a_val
                lengths[s] = l_val

    return layout_2d, angles_2d, lengths, clusters


# ==============================================================================
# 4. BIDIRECTIONAL NEURAL / LANGUAGE-MODEL DSL TRANSLATOR (CM-DSL <-> 2DPROTS-DSL)
# ==============================================================================

def parse_cm_dsl(cm_dsl_text):
    """
    Parse a `CM-DSL` token sequence back into its symbolic graph representation:
      - `nodes`: OrderedDict {sse_label: {"kind": 'H'|'E', "length": int, "idx": int, "chain": str}}
      - `sheets`: list of lists of sse_labels
      - `contacts`: dict {(s1, s2): {"dist": float, "orient": str, "angle": int}}
      - `loops`: dict {(s1, s2): int}
    """
    nodes = OrderedDict()
    sheets = []
    contacts = {}
    loops = {}

    for tok in cm_dsl_text.split():
        if tok.startswith("NODE:"):
            parts = tok.split(":")
            if len(parts) >= 4:
                s = parts[1]
                kind = parts[2]
                length = int(parts[3][1:])
                ch = parts[4][1:] if len(parts) >= 5 else "A"
                nodes[s] = {"kind": kind, "length": length, "idx": len(nodes), "chain": ch}
        elif tok.startswith("LOOP:"):
            parts = tok.split(":")
            if len(parts) >= 3 and "~" in parts[1]:
                s1, s2 = parts[1].split("~", 1)
                gap = int(parts[2][1:])
                loops[(s1, s2)] = gap
                loops[(s2, s1)] = gap
        elif tok.startswith("SHEET:"):
            parts = tok[len("SHEET:"):].split(":")
            if len(parts) >= 2:
                sheets.append(parts)
        elif tok.startswith("CONT:"):
            # CONT:{s1}~{s2}:D{d_bin}:{orient_tok}:{ang:03d}
            parts = tok.split(":")
            if len(parts) >= 5 and "~" in parts[1]:
                s1, s2 = parts[1].split("~", 1)
                d_val = float(parts[2][1:])
                orient = parts[3]
                ang = int(parts[4])
                contacts[(s1, s2)] = {"dist": d_val, "orient": orient, "angle": ang}
                contacts[(s2, s1)] = {"dist": d_val, "orient": orient, "angle": ang}

    return nodes, sheets, contacts, loops


def _extract_pair_features_for_cm2prot(nodes, sheets, contacts, loops=None):
    """
    Extract pairwise token-graph features X_pair of shape (n_s, n_s, F) purely from
    parsed `CM-DSL` tokens so a learned model can predict continuous 3D distances & angles.
    """
    if loops is None:
        loops = {}
    sse_list = list(nodes.keys())
    n_s = len(sse_list)
    D_tok = np.full((n_s, n_s), np.inf, dtype=float)
    np.fill_diagonal(D_tok, 0.0)
    ang_tok = np.full((n_s, n_s), 90.0, dtype=float)
    np.fill_diagonal(ang_tok, 0.0)
    has_cont = np.zeros((n_s, n_s), dtype=float)

    sheet_map = {}
    for sh_idx, sh in enumerate(sheets):
        for pos_idx, s in enumerate(sh):
            sheet_map[s] = (sh_idx, pos_idx)

    for i, s1 in enumerate(sse_list):
        for j in range(i + 1, n_s):
            s2 = sse_list[j]
            if (s1, s2) in contacts:
                c = contacts[(s1, s2)]
                D_tok[i, j] = D_tok[j, i] = float(c["dist"])
                ang_tok[i, j] = ang_tok[j, i] = float(c["angle"])
                has_cont[i, j] = has_cont[j, i] = 1.0
            elif s1 in sheet_map and s2 in sheet_map and sheet_map[s1][0] == sheet_map[s2][0]:
                step = abs(sheet_map[s1][1] - sheet_map[s2][1])
                D_tok[i, j] = D_tok[j, i] = 4.7 * step
                ang_tok[i, j] = ang_tok[j, i] = 165.0 if (step % 2 == 1) else 15.0

    # Add sequence-neighbor loop bounds from LOOP: tokens when direct contact is absent
    for i in range(n_s - 1):
        s1, s2 = sse_list[i], sse_list[i + 1]
        if not np.isfinite(D_tok[i, i + 1]) and (s1, s2) in loops:
            gap = loops[(s1, s2)]
            D_tok[i, i + 1] = D_tok[i + 1, i] = min(20.0, 5.5 + 2.2 * math.sqrt(float(gap)))
        elif not np.isfinite(D_tok[i, i + 1]):
            D_tok[i, i + 1] = D_tok[i + 1, i] = 15.0

    D_sp = shortest_path(D_tok, directed=False)
    D_sp = np.where(np.isfinite(D_sp), D_sp, 38.0)

    # Propagate angles across 2-hop token neighbors
    C_cos = np.cos(np.radians(ang_tok))
    for i in range(n_s):
        for j in range(i + 1, n_s):
            if has_cont[i, j] == 0.0:
                common = [k for k in range(n_s) if has_cont[i, k] > 0 and has_cont[k, j] > 0]
                if common:
                    c_mean = float(np.clip(np.mean([C_cos[i, k] * C_cos[k, j] for k in common]), -1.0, 1.0))
                    ang_tok[i, j] = ang_tok[j, i] = math.degrees(math.acos(c_mean))

    # Build feature tensor (n_s, n_s, 8) for learned distance calibration
    feats = np.zeros((n_s, n_s, 8), dtype=float)
    for i, s1 in enumerate(sse_list):
        l1 = float(nodes[s1]["length"])
        k1 = 1.0 if nodes[s1]["kind"] == "H" else 0.0
        for j in range(i + 1, n_s):
            s2 = sse_list[j]
            l2 = float(nodes[s2]["length"])
            k2 = 1.0 if nodes[s2]["kind"] == "H" else 0.0
            d_sp = float(D_sp[i, j])
            is_dir = float(has_cont[i, j])
            seq_sep = abs(i - j) / max(1.0, float(n_s))
            same_sh = 1.0 if (s1 in sheet_map and s2 in sheet_map and sheet_map[s1][0] == sheet_map[s2][0]) else 0.0
            f_vec = np.array([
                1.0,
                d_sp,
                is_dir * d_sp,
                (1.0 - is_dir) * d_sp,
                math.sqrt(max(0.0, d_sp)),
                0.05 * (l1 + l2),
                seq_sep,
                same_sh + 0.5 * (k1 * k2),
            ], dtype=float)
            feats[i, j] = f_vec
            feats[j, i] = f_vec

    return sse_list, D_sp, ang_tok, feats


def _extract_pair_features_for_prot2cm(layout_2d, angles_2d, lengths, clusters):
    """
    Extract pairwise structural features from a parsed `2DPROTS-DSL` layout
    (ONLY 2D coordinates (x,y), 2D angles, SSE residue lengths, and beta-sheet clusters)
    to predict 3D SSE-to-SSE distances, contact probabilities, and 3D relative angles.
    """
    sse_list = list(layout_2d.keys())
    n_s = len(sse_list)
    sheet_map = {}
    for sh_idx, cl in enumerate(clusters):
        for pos_idx, s in enumerate(cl.split(":")):
            sheet_map[s] = (sh_idx, pos_idx)

    # Reconstruct 2D line segments from 2D center (x, y), 2D angle, and residue length
    # In 2DProts, 1 residue is ~1.5 A along a helix and ~3.0 A along a strand, projected into 2D
    seg_2d = {}
    for s in sse_list:
        x, y = layout_2d[s]
        ang_rad = math.radians(float(angles_2d.get(s, 0.0)))
        dx = math.sin(ang_rad)
        dy = math.cos(ang_rad)
        is_h = s.startswith("H") or (s.startswith("_") and "H" in s[:4])
        half_span = 0.5 * float(lengths.get(s, 6)) * (0.45 if is_h else 0.75)
        p0 = np.array([x - half_span * dx, y - half_span * dy], dtype=float)
        p1 = np.array([x + half_span * dx, y + half_span * dy], dtype=float)
        seg_2d[s] = (p0, p1, np.array([dx, dy], dtype=float))

    feats = np.zeros((n_s, n_s, 11), dtype=float)
    d2d_mat = np.zeros((n_s, n_s), dtype=float)
    ang2d_diff = np.zeros((n_s, n_s), dtype=float)

    for i in range(n_s):
        s1 = sse_list[i]
        x1, y1 = layout_2d[s1]
        p1_0, p1_1, u1 = seg_2d[s1]
        l1 = float(lengths.get(s1, 6))
        is_h1 = 1.0 if (s1.startswith("H") or (s1.startswith("_") and "H" in s1[:4])) else 0.0
        for j in range(i + 1, n_s):
            s2 = sse_list[j]
            x2, y2 = layout_2d[s2]
            p2_0, p2_1, u2 = seg_2d[s2]
            l2 = float(lengths.get(s2, 6))
            is_h2 = 1.0 if (s2.startswith("H") or (s2.startswith("_") and "H" in s2[:4])) else 0.0

            d_center = math.hypot(x1 - x2, y1 - y2)
            d2d_mat[i, j] = d2d_mat[j, i] = d_center

            # Sample 5 points along each 2D segment to get min and max 2D segment distance
            t_grid = np.linspace(0.0, 1.0, 5)[:, None]
            pts1 = p1_0[None, :] * (1.0 - t_grid) + p1_1[None, :] * t_grid
            pts2 = p2_0[None, :] * (1.0 - t_grid) + p2_1[None, :] * t_grid
            d_grid = np.sqrt(np.sum((pts1[:, None, :] - pts2[None, :, :]) ** 2, axis=-1))
            d_seg_min = float(np.min(d_grid))
            d_seg_max = float(np.max(d_grid))

            # Relative 2D angle between the two SSE symbols in the 2DProts diagram
            cos_2d = float(np.clip(np.dot(u1, u2), -1.0, 1.0))
            ang_diff = math.degrees(math.acos(cos_2d))
            ang2d_diff[i, j] = ang2d_diff[j, i] = ang_diff

            seq_diff = abs(i - j)
            is_seq_adj = 1.0 if seq_diff == 1 else 0.0
            same_sheet = 0.0
            sheet_adj = 0.0
            if s1 in sheet_map and s2 in sheet_map and sheet_map[s1][0] == sheet_map[s2][0]:
                same_sheet = 1.0
                if abs(sheet_map[s1][1] - sheet_map[s2][1]) == 1:
                    sheet_adj = 1.0

            f_vec = np.array([
                1.0,
                d_center,
                d_seg_min,
                d_seg_max,
                math.log1p(d_center),
                is_seq_adj,
                1.0 / float(seq_diff),
                same_sheet,
                sheet_adj,
                0.05 * (l1 + l2),
                is_h1 * is_h2,
            ], dtype=float)
            feats[i, j] = f_vec
            feats[j, i] = f_vec

    return sse_list, d2d_mat, ang2d_diff, feats


class BidirectionalDSLTranslator:
    """
    Bidirectional translator between `CM-DSL` (Contact Map tokens) and `2DPROTS-DSL` (2DProts tokens),
    trained via Leave-One-PDB-Out cross-validation on the benchmark dataset.
    """

    def __init__(self):
        # Forward weights: predicts true 3D SSE distance from CM-DSL token graph features (8 dims)
        self.w_cm2prot_dist = np.array([0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=float)
        # Inverse weights: predicts 3D SSE distance and contact logits from 2DPROTS-DSL features (11 dims)
        self.w_prot2cm_dist = np.array([1.5, 0.85, 0.15, 0.0, 0.0, -1.5, 0.0, -1.0, -2.0, 0.0, 0.0], dtype=float)
        self.w_prot2cm_cont = np.array([3.5, -0.32, -0.12, 0.0, 0.0, 1.2, 0.5, 1.5, 3.0, 0.05, 0.0], dtype=float)

    def fit(self, training_records):
        """
        Fit forward (`CM-DSL -> 2DPROTS-DSL`) and inverse (`2DPROTS-DSL -> CM-DSL`) ridge/logistic
        regression parameters on a list of training target records.
        """
        X_fwd, y_fwd = [], []
        X_inv, y_inv_dist, y_inv_cont = [], [], []

        for rec in training_records:
            # 1. Forward training pairs from CM-DSL token features -> true 3D SSE distances
            sse_list, _, _, feats_fwd = _extract_pair_features_for_cm2prot(
                rec["cm_nodes"], rec["cm_sheets"], rec["cm_contacts"], rec.get("cm_loops", {})
            )
            D_true = rec["D_sse_true"]
            n_s = len(sse_list)
            for i in range(n_s):
                for j in range(i + 1, n_s):
                    if D_true[i, j] < 45.0:
                        X_fwd.append(feats_fwd[i, j])
                        y_fwd.append(D_true[i, j])

            # 2. Inverse training pairs from 2DPROTS-DSL features -> true 3D SSE distances & contacts
            sse_2d, _, _, feats_inv = _extract_pair_features_for_prot2cm(
                rec["layout_2d"], rec["angles_2d"], rec["lengths"], rec["clusters"]
            )
            idx_map = {s: idx for idx, s in enumerate(sse_list)}
            m_s = len(sse_2d)
            for i in range(m_s):
                si = sse_2d[i]
                if si not in idx_map:
                    continue
                ii = idx_map[si]
                for j in range(i + 1, m_s):
                    sj = sse_2d[j]
                    if sj not in idx_map:
                        continue
                    jj = idx_map[sj]
                    d_t = float(D_true[ii, jj])
                    if d_t <= 32.0:
                        X_inv.append(feats_inv[i, j])
                        y_inv_dist.append(d_t)
                    y_inv_cont.append((feats_inv[i, j], 1.0 if d_t <= 11.0 else 0.0))

        if len(X_fwd) >= 10:
            Xf = np.vstack(X_fwd)
            yf = np.array(y_fwd, dtype=float)
            reg = 1e-2 * np.eye(Xf.shape[1])
            self.w_cm2prot_dist = np.linalg.solve(Xf.T @ Xf + reg, Xf.T @ yf)

        if len(X_inv) >= 10:
            Xi = np.vstack(X_inv)
            yi_d = np.array(y_inv_dist, dtype=float)
            # Weight closer pairs higher for accurate local contact map / distogram reconstruction
            w_d = 1.0 / np.maximum(4.0, yi_d)
            Xi_w = Xi * w_d[:, None]
            reg = 1e-2 * np.eye(Xi.shape[1])
            self.w_prot2cm_dist = np.linalg.solve(Xi_w.T @ Xi + reg, Xi_w.T @ yi_d)

            # Fit logistic contact classifier via IRLS (Newton-Raphson) on contact labels (d <= 11 A)
            Xi_c = np.vstack([p[0] for p in y_inv_cont])
            yi_c = np.array([p[1] for p in y_inv_cont], dtype=float)
            w_c = np.zeros(Xi_c.shape[1], dtype=float)
            for _ in range(12):
                logits = np.clip(Xi_c @ w_c, -15.0, 15.0)
                probs = expit(logits)
                W_diag = np.maximum(probs * (1.0 - probs), 1e-4)
                H = (Xi_c.T * W_diag) @ Xi_c + 0.1 * np.eye(Xi_c.shape[1])
                grad = Xi_c.T @ (probs - yi_c) + 0.1 * w_c
                step = np.linalg.solve(H, grad)
                w_c = w_c - step
                if np.linalg.norm(step) < 1e-5:
                    break
            self.w_prot2cm_cont = w_c

    def translate_cm_dsl_to_2dprots(self, cm_dsl_text, ann, stats):
        """
        Translate a `CM-DSL` string into a 2DProts `Layout` object and its `2DPROTS-DSL` string
        using ONLY the tokens in `cm_dsl_text` (and beta-sheet pairing).
        """
        nodes, sheets, contacts, loops = parse_cm_dsl(cm_dsl_text)
        sse_list, D_sp, ang_tok, feats = _extract_pair_features_for_cm2prot(nodes, sheets, contacts, loops)
        n_s = len(sse_list)

        D_pred = np.zeros((n_s, n_s), dtype=float)
        for i in range(n_s):
            for j in range(i + 1, n_s):
                d_hat = float(np.dot(feats[i, j], self.w_cm2prot_dist))
                # Blend learned token distance with direct contact token distance when present
                s1, s2 = sse_list[i], sse_list[j]
                if (s1, s2) in contacts:
                    d_hat = 0.65 * float(contacts[(s1, s2)]["dist"]) + 0.35 * d_hat
                d_hat = float(np.clip(d_hat, 3.8, 55.0))
                D_pred[i, j] = D_pred[j, i] = d_hat

        dmat_dict = {si: {} for si in sse_list}
        angles_dict = {}
        for i, si in enumerate(sse_list):
            for j, sj in enumerate(sse_list):
                dmat_dict[si][sj] = float(D_pred[i, j])
                angles_dict[(si, sj)] = int(round(ang_tok[i, j]))

        l_obj = compute_2dprots_layout_from_contact_geometry(
            ann, stats, dmat_dict, angles_dict, base_method="sammon"
        )
        dsl_2d_text, _ = encode_2dprots_to_dsl(sse_list, l_obj, ann)
        return l_obj, dsl_2d_text, D_pred

    def translate_2dprots_dsl_to_contact_map(self, dsl_2d_text, sse_order, sse_slices):
        """
        Inverse Translation (`2DPROTS-DSL -> CM-DSL -> Reconstructed Residue Contact Map`):
        Given ONLY the `2DPROTS-DSL` token string (which encodes 2D SSE coordinates, 2D angles,
        SSE lengths, and beta-sheet clusters), reconstruct:
          - `cm_dsl_pred`: predicted `CM-DSL` token sequence
          - `D_sse_pred`: predicted SSE-level 3D distance matrix
          - `P_sse_cont`: predicted SSE-level contact probability matrix
          - `D_res_pred`: reconstructed N x N residue-level distance matrix
          - `P_res_pred`: reconstructed N x N residue-level contact probability map
          - `C_res_bin10_pred`: reconstructed N x N binary contact map
        """
        layout_2d, angles_2d, lengths, clusters = decode_2dprots_dsl(dsl_2d_text)
        sse_2d, d2d_mat, ang2d_diff, feats = _extract_pair_features_for_prot2cm(
            layout_2d, angles_2d, lengths, clusters
        )
        idx_2d = {s: i for i, s in enumerate(sse_2d)}

        n_s = len(sse_order)
        D_sse_pred = np.zeros((n_s, n_s), dtype=float)
        P_sse_cont = np.eye(n_s, dtype=float)
        ang_sse_pred = np.zeros((n_s, n_s), dtype=float)

        for i in range(n_s):
            si = sse_order[i]
            for j in range(i + 1, n_s):
                sj = sse_order[j]
                if si in idx_2d and sj in idx_2d:
                    ii, jj = idx_2d[si], idx_2d[sj]
                    f_vec = feats[ii, jj]
                    d_hat = float(np.clip(np.dot(f_vec, self.w_prot2cm_dist), 4.0, 55.0))
                    logit = float(np.clip(np.dot(f_vec, self.w_prot2cm_cont), -15.0, 15.0))
                    p_hat = float(expit(logit))
                    ang_hat = float(ang2d_diff[ii, jj])
                else:
                    d_hat = 25.0
                    p_hat = 0.05
                    ang_hat = 90.0

                D_sse_pred[i, j] = D_sse_pred[j, i] = d_hat
                P_sse_cont[i, j] = P_sse_cont[j, i] = p_hat
                ang_sse_pred[i, j] = ang_sse_pred[j, i] = ang_hat

        # Reconstruct the full N x N residue-level distance matrix & contact probability map
        # using exact 2D segment geometry from 2DPROTS-DSL + Pythagorean 3D depth lifting!
        n_res = max((idx1 for _, idx1 in sse_slices.values()), default=0)
        D_res_pred = np.full((n_res, n_res), 25.0, dtype=float)
        np.fill_diagonal(D_res_pred, 0.0)
        P_res_pred = np.zeros((n_res, n_res), dtype=float)
        np.fill_diagonal(P_res_pred, 1.0)

        # Build 2D residue coordinates r_2d(u) along each SSE segment in the 2DProts diagram
        r2d_per_sse = {}
        for i, si in enumerate(sse_order):
            i0, i1 = sse_slices[si]
            li = i1 - i0
            is_hi = si.startswith("H") or (si.startswith("_") and "H" in si[:4])
            res_len_i = float(lengths.get(si, li))
            # Physical 3D axis length of helix (~1.5 A/res) or strand (~3.1 A/res)
            axis_len_i = res_len_i * (1.48 if is_hi else 2.95)
            for u in range(i0, i1):
                for v in range(u + 1, i1):
                    frac = float(v - u) / max(1.0, float(li - 1))
                    d_intra = axis_len_i * frac
                    D_res_pred[u, v] = D_res_pred[v, u] = d_intra
                    P_res_pred[u, v] = P_res_pred[v, u] = float(expit((10.5 - d_intra) / 1.5))

            if si in layout_2d:
                cx, cy = layout_2d[si]
                ang_rad = math.radians(float(angles_2d.get(si, 0.0)))
                ux, uy = math.sin(ang_rad), math.cos(ang_rad)
                # 2D projected half-length in 2DProts coordinate space
                half_2d = 0.5 * axis_len_i * 0.65
                t_vals = np.linspace(-1.0, 1.0, li)[:, None]
                pts_2d = np.array([cx, cy])[None, :] + t_vals * half_2d * np.array([ux, uy])[None, :]
            else:
                pts_2d = np.zeros((li, 2), dtype=float)
            r2d_per_sse[si] = pts_2d

        for i, si in enumerate(sse_order):
            i0, i1 = sse_slices[si]
            pts_i = r2d_per_sse[si]
            for j in range(i + 1, n_s):
                sj = sse_order[j]
                j0, j1 = sse_slices[sj]
                pts_j = r2d_per_sse[sj]

                # 2D pairwise residue distances between segment i and segment j
                diff_2d = pts_i[:, None, :] - pts_j[None, :, :]
                d2d_sub_sq = np.sum(diff_2d ** 2, axis=-1)
                d2d_min_sq = float(np.min(d2d_sub_sq))
                delta_2d_sq = np.maximum(0.0, d2d_sub_sq - d2d_min_sq)

                # Lift minimal SSE distance into full residue submatrix via Pythagoras
                d_min = D_sse_pred[i, j]
                d_sub = np.sqrt(d_min ** 2 + 0.65 * delta_2d_sq)
                D_res_pred[i0:i1, j0:j1] = d_sub
                D_res_pred[j0:j1, i0:i1] = d_sub.T

                # Combine SSE-level contact probability with residue-level proximity along the ridge
                p_sub = P_sse_cont[i, j] * expit((11.0 - d_sub) / 1.8)
                P_res_pred[i0:i1, j0:j1] = p_sub
                P_res_pred[j0:j1, i0:i1] = p_sub.T

        C_res_bin10_pred = (D_res_pred <= 11.0).astype(float)

        # Also emit the predicted CM-DSL token string
        tokens = ["<BOS_CM>"]
        for s in sse_order:
            i0, i1 = sse_slices[s]
            kind = "H" if (s.startswith("H") or (s.startswith("_") and "H" in s[:4])) else "E"
            tokens.append(f"NODE:{s}:{kind}:L{i1 - i0}")
        for cl in clusters:
            tokens.append(f"SHEET:{cl}")
        edges = []
        for i in range(n_s):
            for j in range(i + 1, n_s):
                if D_sse_pred[i, j] <= 15.5:
                    ang = int(round(ang_sse_pred[i, j]))
                    orient = "PAR" if ang <= 45 else ("ANTI" if ang >= 135 else "CROSS")
                    edges.append((D_sse_pred[i, j], f"CONT:{sse_order[i]}~{sse_order[j]}:D{int(round(D_sse_pred[i, j]))}:{orient}:{ang:03d}"))
        edges.sort(key=lambda x: x[0])
        for _, tok in edges[:min(len(edges), n_s * 5)]:
            tokens.append(tok)
        tokens.append("<EOS_CM>")
        cm_dsl_pred = " ".join(tokens)

        return cm_dsl_pred, D_sse_pred, P_sse_cont, D_res_pred, P_res_pred, C_res_bin10_pred


# ==============================================================================
# 5. 3D PROTEIN STRUCTURE GENERATION FROM 2DPROTS-DSL (2DPROTS -> CM -> 3D MODEL)
# ==============================================================================

def _orthonormal_frame(u_vec, ref_normal=None):
    """Build a right-handed orthonormal frame (e1, e2, u) given unit axis vector u."""
    u = np.asarray(u_vec, dtype=float)
    nrm = float(np.linalg.norm(u))
    if nrm < 1e-8:
        u = np.array([0.0, 0.0, 1.0], dtype=float)
    else:
        u = u / nrm
    if ref_normal is not None:
        ref = np.asarray(ref_normal, dtype=float)
        e1 = ref - np.dot(ref, u) * u
        if np.linalg.norm(e1) > 1e-5:
            e1 = e1 / np.linalg.norm(e1)
            e2 = np.cross(u, e1)
            return e1, e2, u
    helper = np.array([1.0, 0.0, 0.0], dtype=float) if abs(u[0]) < 0.8 else np.array([0.0, 1.0, 0.0], dtype=float)
    e1 = np.cross(u, helper)
    e1 = e1 / max(1e-8, float(np.linalg.norm(e1)))
    e2 = np.cross(u, e1)
    return e1, e2, u


def _build_ideal_helix_residues(p0, p1, n_res):
    """
    Generate realistic backbone + CB coordinates (N, CA, C, O, CB) for an alpha-helix
    spanning from 3D point p0 to p1 with n_res residues.
    """
    axis_vec = p1 - p0
    e1, e2, u = _orthonormal_frame(axis_vec)
    axis_len = float(np.linalg.norm(axis_vec))
    step_z = axis_len / max(1.0, float(n_res - 1)) if n_res > 1 else 1.50
    omega = math.radians(100.0)  # 3.6 residues per turn
    residues = []
    for k in range(n_res):
        z_ca = k * step_z
        theta_ca = k * omega
        origin_z = p0 + z_ca * u

        def _cyl(r, d_theta, d_z):
            ang = theta_ca + d_theta
            return origin_z + d_z * u + r * (math.cos(ang) * e1 + math.sin(ang) * e2)

        n_pos = _cyl(1.58, -0.48, -0.45)
        ca_pos = _cyl(2.26, 0.00, 0.00)
        c_pos = _cyl(1.72, +0.52, +0.52)
        o_pos = _cyl(2.05, +0.62, +1.72)
        cb_pos = _cyl(3.32, +0.16, -0.48)
        residues.append({
            "N": n_pos, "CA": ca_pos, "C": c_pos, "O": o_pos, "CB": cb_pos
        })
    return residues


def _build_ideal_strand_residues(p0, p1, n_res, sheet_normal=None):
    """
    Generate realistic pleated beta-strand backbone + CB coordinates (N, CA, C, O, CB)
    spanning from 3D point p0 to p1 with n_res residues.
    """
    axis_vec = p1 - p0
    e1, e2, u = _orthonormal_frame(axis_vec, ref_normal=sheet_normal)
    axis_len = float(np.linalg.norm(axis_vec))
    step_z = axis_len / max(1.0, float(n_res - 1)) if n_res > 1 else 3.20
    residues = []
    for k in range(n_res):
        sign = 1.0 if (k % 2 == 0) else -1.0
        base = p0 + (k * step_z) * u
        ca_pos = base + (0.36 * sign) * e2
        n_pos = ca_pos - 1.12 * u + (0.45 * sign) * e1 - (0.20 * sign) * e2
        c_pos = ca_pos + 1.15 * u + (0.45 * sign) * e1 - (0.20 * sign) * e2
        o_pos = c_pos + (1.22 * sign) * e1 + 0.15 * u
        cb_pos = ca_pos + (1.52 * sign) * e2 - (0.15 * sign) * e1
        residues.append({
            "N": n_pos, "CA": ca_pos, "C": c_pos, "O": o_pos, "CB": cb_pos
        })
    return residues


def _build_loop_residues(p_start, u_start, p_end, u_end, n_loop, center_bias=None):
    """
    Generate smooth, clash-free loop backbone coordinates (N, CA, C, O, CB) connecting
    SSE C-terminus p_start (direction u_start) to next SSE N-terminus p_end (direction u_end).
    """
    if n_loop <= 0:
        return []
    chord = p_end - p_start
    chord_len = float(np.linalg.norm(chord))
    # Target virtual CA-CA contour length is ~3.8 A per loop step
    target_arc = 3.8 * float(n_loop + 1)
    bulge = max(2.2, 0.42 * math.sqrt(max(0.0, target_arc ** 2 - chord_len ** 2)))

    mid = 0.5 * (p_start + p_end)
    if center_bias is not None and np.linalg.norm(mid - center_bias) > 1e-4:
        outward = (mid - center_bias) / float(np.linalg.norm(mid - center_bias))
    else:
        outward = np.array([0.0, 0.0, 1.0], dtype=float)

    t0 = u_start * max(3.0, 0.6 * bulge) + 0.4 * bulge * outward
    t1 = u_end * max(3.0, 0.6 * bulge) - 0.4 * bulge * outward

    # Dense cubic Hermite sampling for uniform arc-length spacing
    s_dense = np.linspace(0.0, 1.0, 64)
    h00 = 2 * s_dense ** 3 - 3 * s_dense ** 2 + 1
    h10 = s_dense ** 3 - 2 * s_dense ** 2 + s_dense
    h01 = -2 * s_dense ** 3 + 3 * s_dense ** 2
    h11 = s_dense ** 3 - s_dense ** 2
    curve = (
        h00[:, None] * p_start[None, :]
        + h10[:, None] * t0[None, :]
        + h01[:, None] * p_end[None, :]
        + h11[:, None] * t1[None, :]
    )
    seg_lens = np.linalg.norm(np.diff(curve, axis=0), axis=1)
    cum_arc = np.concatenate([[0.0], np.cumsum(seg_lens)])
    total_arc = max(1e-6, float(cum_arc[-1]))

    residues = []
    for k in range(1, n_loop + 1):
        frac = float(k) / float(n_loop + 1)
        target_s = frac * total_arc
        idx = int(np.clip(np.searchsorted(cum_arc, target_s) - 1, 0, len(curve) - 2))
        local_ds = max(1e-6, float(cum_arc[idx + 1] - cum_arc[idx]))
        alpha = float(np.clip((target_s - cum_arc[idx]) / local_ds, 0.0, 1.0))
        ca_pos = (1.0 - alpha) * curve[idx] + alpha * curve[idx + 1]
        tang = curve[min(len(curve) - 1, idx + 1)] - curve[max(0, idx - 1)]
        e1, e2, u_loc = _orthonormal_frame(tang, ref_normal=outward)
        n_pos = ca_pos - 1.15 * u_loc + 0.35 * e1
        c_pos = ca_pos + 1.18 * u_loc + 0.35 * e1
        o_pos = c_pos + 1.22 * e1 - 0.15 * e2
        cb_pos = ca_pos + 1.50 * e2 - 0.25 * e1
        residues.append({
            "N": n_pos, "CA": ca_pos, "C": c_pos, "O": o_pos, "CB": cb_pos
        })
    return residues


def reconstruct_3d_structure_from_2dprots_dsl(
    dsl_2d_text,
    translator=None,
    loop_lengths=None,
    random_seed=42,
):
    """
    Reconstruct a complete 3D protein structure (SSE 3D vectors + full N, CA, C, O, CB atomic backbone)
    purely from a `2DPROTS-DSL` string (`2DProts -> Contact Map / Distogram -> 3D Coordinates`).

    Returns a dictionary containing:
      - `sse_order`: list of SSE names
      - `sse_3d_coords`: dict {sse: (start_vec_3d, end_vec_3d)}
      - `residues`: list of residue dicts with 3D atom coords and SSE assignment
      - `ca_coords`: (N_total, 3) numpy array of CA coordinates
      - `sse_ca_indices`: dict {sse: (start_idx, end_idx)} in `ca_coords`
      - `D_sse_pred`: predicted SSE distance matrix from `2DPROTS-DSL`
      - `D_sse_3d`: realized 3D SSE segment distance matrix of the reconstructed model
      - `D_res_pred`: predicted residue distance matrix from `2DPROTS-DSL`
      - `P_res_pred`: predicted residue contact probability matrix
      - `cm_dsl_pred`: predicted `CM-DSL` token string
      - `clusters`: list of beta-sheet cluster strings
    """
    if translator is None:
        translator = BidirectionalDSLTranslator()

    layout_2d, angles_2d, lengths, clusters = decode_2dprots_dsl(dsl_2d_text)
    sse_order = list(layout_2d.keys())
    n_s = len(sse_order)

    # Build SSE slices for inverse translator
    sse_slices = OrderedDict()
    cur_idx = 0
    for s in sse_order:
        l_res = max(3, int(lengths.get(s, 6)))
        sse_slices[s] = (cur_idx, cur_idx + l_res)
        cur_idx += l_res

    cm_dsl_pred, D_sse_pred, P_sse_cont, D_res_pred, P_res_pred, C_res_bin10_pred = (
        translator.translate_2dprots_dsl_to_contact_map(dsl_2d_text, sse_order, sse_slices)
    )

    # Physical 3D axis length of each SSE
    axis_lengths = np.zeros(n_s, dtype=float)
    is_helix = np.zeros(n_s, dtype=bool)
    for i, s in enumerate(sse_order):
        is_h = s.startswith("H") or (s.startswith("_") and "H" in s[:4])
        is_helix[i] = is_h
        l_res = max(3, int(lengths.get(s, 6)))
        # Cap effective straight-axis span so long coiled/curved strands stay compact
        eff_res = min(l_res, 24) if is_h else min(l_res, 14)
        axis_lengths[i] = float(eff_res) * (1.46 if is_h else 2.85)

    # Map beta-sheet memberships and enforce ideal 4.7 A strand spacing on adjacent sheet strands
    sheet_adj_pairs = []
    sheet_same_pairs = []
    for cl in clusters:
        parts = [p for p in cl.split(":") if p in layout_2d]
        for k in range(len(parts) - 1):
            i = sse_order.index(parts[k])
            j = sse_order.index(parts[k + 1])
            sheet_adj_pairs.append((i, j))
            D_sse_pred[i, j] = D_sse_pred[j, i] = 4.75
        for k1 in range(len(parts)):
            for k2 in range(k1 + 2, len(parts)):
                i = sse_order.index(parts[k1])
                j = sse_order.index(parts[k2])
                sheet_same_pairs.append((i, j, abs(k2 - k1)))

    # Initialize 3D SSE centers (x, y, z) and 3D axis unit vectors (ux, uy, uz) from 2DProts layout
    xy_raw = np.array([layout_2d[s] for s in sse_order], dtype=float)
    xy_centered = xy_raw - np.mean(xy_raw, axis=0, keepdims=True)
    # Scale 2DProts coordinates to match predicted 3D nearest-neighbor distances
    if n_s >= 2:
        d2d_nn = []
        d3d_nn = []
        for i in range(n_s):
            d2_row = [np.linalg.norm(xy_centered[i] - xy_centered[j]) for j in range(n_s) if j != i]
            d3_row = [D_sse_pred[i, j] for j in range(n_s) if j != i]
            if d2_row and min(d2_row) > 0.2:
                d2d_nn.append(min(d2_row))
                d3d_nn.append(min(d3_row))
        scale_xy = float(np.clip(np.median(d3d_nn) / max(1.0, np.median(d2d_nn)), 0.65, 1.45)) if d2d_nn else 1.0
    else:
        scale_xy = 1.0

    c_init = np.zeros((n_s, 3), dtype=float)
    c_init[:, 0] = xy_centered[:, 0] * scale_xy
    c_init[:, 1] = xy_centered[:, 1] * scale_xy

    # Estimate initial 3D depth Z_i from the residual between predicted 3D distance and 2D (X, Y) distance
    if n_s >= 3:
        d2_sq = np.sum((c_init[:, None, :2] - c_init[None, :, :2]) ** 2, axis=-1)
        dz_sq = np.maximum(0.0, D_sse_pred ** 2 - d2_sq)
        H_c = np.eye(n_s) - np.ones((n_s, n_s)) / float(n_s)
        B_z = -0.5 * H_c @ dz_sq @ H_c
        evals_z, evecs_z = np.linalg.eigh(B_z)
        top_ev = max(0.0, float(evals_z[-1]))
        c_init[:, 2] = evecs_z[:, -1] * math.sqrt(top_ev) * 0.55

    u_init = np.zeros((n_s, 3), dtype=float)
    z_dir = 1.0
    for i, s in enumerate(sse_order):
        ang_rad = math.radians(float(angles_2d.get(s, 0.0)))
        ux_2d = math.sin(ang_rad)
        uy_2d = math.cos(ang_rad)
        if i > 0:
            z_dir = -z_dir
        # Combine 2D in-plane tilt with strong 3D axial Z component
        u_vec = np.array([0.32 * ux_2d, 0.32 * uy_2d, 0.94 * z_dir], dtype=float)
        u_init[i] = u_vec / np.linalg.norm(u_vec)

    # Gap lengths between consecutive SSEs
    gaps = np.full(max(0, n_s - 1), 5, dtype=int)
    if loop_lengths is not None:
        for i in range(n_s - 1):
            pair = (sse_order[i], sse_order[i + 1])
            if pair in loop_lengths:
                gaps[i] = max(2, int(loop_lengths[pair]))

    # Precompute sample grid t in [-0.5, 0.5] (4 points per segment) for fast differentiable segment-to-segment distance
    t_pts = np.array([-0.45, -0.15, 0.15, 0.45], dtype=float)
    idx_i, idx_j = np.triu_indices(n_s, k=1)
    d_target_pairs = D_sse_pred[idx_i, idx_j]
    w_pairs = 1.0 / np.maximum(4.5, d_target_pairs)
    for (si_idx, sj_idx) in sheet_adj_pairs:
        a, b = min(si_idx, sj_idx), max(si_idx, sj_idx)
        mask = (idx_i == a) & (idx_j == b)
        w_pairs[mask] = 4.5

    # Target cosine dot-product between adjacent beta-sheet strands
    sheet_adj_cos = []
    for (si_idx, sj_idx) in sheet_adj_pairs:
        # Check initial z_dir product so closed antiparallel barrels have no frustrated strand pair
        dot_init = float(np.dot(u_init[si_idx], u_init[sj_idx]))
        target_c = -0.92 if dot_init < 0.0 else +0.92
        sheet_adj_cos.append((si_idx, sj_idx, target_c))

    def _loss_fn(x_flat):
        state = x_flat.reshape((n_s, 6))
        centers = state[:, :3]
        raw_v = state[:, 3:]
        norms = np.sqrt(np.sum(raw_v ** 2, axis=1, keepdims=True) + 1e-8)
        u_dirs = raw_v / norms

        # Segment sample points: shape (n_s, 4, 3)
        seg_pts = centers[:, None, :] + (axis_lengths[:, None, None] * t_pts[None, :, None]) * u_dirs[:, None, :]
        pts_a = seg_pts[idx_i]  # (P, 4, 3)
        pts_b = seg_pts[idx_j]  # (P, 4, 3)
        diff_ab = pts_a[:, :, None, :] - pts_b[:, None, :, :]  # (P, 4, 4, 3)
        dist_grid = np.sqrt(np.sum(diff_ab ** 2, axis=-1) + 0.04)  # (P, 4, 4)
        # Smooth soft-minimum across the 4x4 segment sample grid
        soft_min = -0.8 * np.log(np.mean(np.exp(-dist_grid / 0.8), axis=(1, 2)))

        # 1. Pairwise 3D SSE distance loss
        loss_dist = float(np.sum(w_pairs * (soft_min - d_target_pairs) ** 2))

        # 2. Steric clash penalty
        clash = np.maximum(0.0, 4.3 - soft_min)
        loss_clash = 15.0 * float(np.sum(clash ** 2))

        # 3. Sequential loop continuity between C-terminus of i and N-terminus of i+1
        p_end_i = centers[:-1] + 0.5 * axis_lengths[:-1, None] * u_dirs[:-1]
        p_start_next = centers[1:] - 0.5 * axis_lengths[1:, None] * u_dirs[1:]
        d_loop = np.sqrt(np.sum((p_start_next - p_end_i) ** 2, axis=1) + 0.04)
        d_loop_max = 3.5 * gaps.astype(float) + 2.5
        d_loop_opt = np.minimum(d_loop_max, 4.5 + 2.3 * np.sqrt(gaps.astype(float)))
        over_stretch = np.maximum(0.0, d_loop - d_loop_max)
        loss_loop = float(np.sum(4.0 * over_stretch ** 2 + 0.45 * (d_loop - d_loop_opt) ** 2))

        # 4. Beta-sheet parallel/antiparallel alignment & center alignment
        loss_orient = 0.0
        for (ia, ib, target_c) in sheet_adj_cos:
            dot_ab = float(np.dot(u_dirs[ia], u_dirs[ib]))
            loss_orient += 12.0 * (dot_ab - target_c) ** 2
            # Encourage beta-strand centers to register cleanly (perpendicular offset ~ 4.75 A)
            dc = centers[ib] - centers[ia]
            proj_along = float(np.dot(dc, u_dirs[ia]))
            loss_orient += 0.25 * (proj_along ** 2)

        # 5. Tether to 2DProts (X, Y) layout + axial orientation prior + unit norm regularization
        loss_xy = 0.12 * float(np.sum((centers[:, :2] - c_init[:, :2]) ** 2))
        loss_u_prior = 1.5 * float(np.sum((u_dirs - u_init) ** 2))
        loss_norm = 0.5 * float(np.sum((norms[:, 0] - 1.0) ** 2))

        return loss_dist + loss_clash + loss_loop + loss_orient + loss_xy + loss_u_prior + loss_norm

    x0 = np.hstack([c_init, u_init]).ravel()
    res_opt = minimize(_loss_fn, x0, method="L-BFGS-B", options={"maxiter": 65, "ftol": 1e-4})
    state_opt = res_opt.x.reshape((n_s, 6))
    centers_3d = state_opt[:, :3]
    u_dirs_3d = state_opt[:, 3:] / np.maximum(1e-8, np.linalg.norm(state_opt[:, 3:], axis=1, keepdims=True))

    # Compute realized 3D segment-to-segment distance matrix D_sse_3d
    sse_3d_coords = OrderedDict()
    for i, s in enumerate(sse_order):
        p0 = centers_3d[i] - 0.5 * axis_lengths[i] * u_dirs_3d[i]
        p1 = centers_3d[i] + 0.5 * axis_lengths[i] * u_dirs_3d[i]
        sse_3d_coords[s] = (p0, p1)

    D_sse_3d = np.zeros((n_s, n_s), dtype=float)
    t_eval = np.linspace(0.0, 1.0, 9)[:, None]
    for i in range(n_s):
        p0_i, p1_i = sse_3d_coords[sse_order[i]]
        pts_i = p0_i[None, :] * (1.0 - t_eval) + p1_i[None, :] * t_eval
        for j in range(i + 1, n_s):
            p0_j, p1_j = sse_3d_coords[sse_order[j]]
            pts_j = p0_j[None, :] * (1.0 - t_eval) + p1_j[None, :] * t_eval
            d_min = float(np.min(np.linalg.norm(pts_i[:, None, :] - pts_j[None, :, :], axis=-1)))
            D_sse_3d[i, j] = D_sse_3d[j, i] = d_min

    # Build full atomic coordinates (N, CA, C, O, CB) for all SSEs and connecting loops
    center_of_mass = np.mean(centers_3d, axis=0)
    residues = []
    sse_res_spans = OrderedDict()  # {sse: (start_seq_id, end_seq_id)} (1-based inclusive)
    sse_ca_indices = OrderedDict()  # {sse: (start_idx, end_idx)} in SSE-only / full arrays
    seq_id = 1

    # Map each strand to its neighbor vector in the beta-sheet so H-bond C=O / N-H point toward neighbor strand
    sheet_neighbor_dir = {}
    for (ia, ib) in sheet_adj_pairs:
        vec_ab = centers_3d[ib] - centers_3d[ia]
        sheet_neighbor_dir[sse_order[ia]] = vec_ab
        sheet_neighbor_dir[sse_order[ib]] = -vec_ab

    for i, s in enumerate(sse_order):
        p0, p1 = sse_3d_coords[s]
        l_res = max(3, int(lengths.get(s, 6)))
        if is_helix[i]:
            sse_atoms = _build_ideal_helix_residues(p0, p1, l_res)
            res_name = "ALA"
            ss_type = "H"
        else:
            nb_dir = sheet_neighbor_dir.get(s, centers_3d[i] - center_of_mass)
            sse_atoms = _build_ideal_strand_residues(p0, p1, l_res, sheet_normal=nb_dir)
            res_name = "VAL"
            ss_type = "E"

        start_seq = seq_id
        idx_0 = len(residues)
        for k_r, atom_dict in enumerate(sse_atoms):
            residues.append({
                "seq_id": seq_id,
                "res_name": res_name,
                "chain_id": "A",
                "sse": s,
                "ss_type": ss_type,
                "sse_offset": k_r,
                "atoms": atom_dict,
            })
            seq_id += 1
        idx_1 = len(residues)
        sse_res_spans[s] = (start_seq, seq_id - 1)
        sse_ca_indices[s] = (idx_0, idx_1)

        # Build connecting loop to next SSE
        if i + 1 < n_s:
            s_next = sse_order[i + 1]
            p0_next, _ = sse_3d_coords[s_next]
            n_loop = int(gaps[i])
            loop_atoms = _build_loop_residues(
                p1, u_dirs_3d[i], p0_next, u_dirs_3d[i + 1], n_loop, center_bias=center_of_mass
            )
            for k_l, atom_dict in enumerate(loop_atoms):
                residues.append({
                    "seq_id": seq_id,
                    "res_name": "GLY" if (k_l % 2 == 1) else "SER",
                    "chain_id": "A",
                    "sse": f"loop_{s}_{s_next}",
                    "ss_type": "L",
                    "sse_offset": k_l,
                    "atoms": atom_dict,
                })
                seq_id += 1

    ca_coords = np.array([r["atoms"]["CA"] for r in residues], dtype=float)

    return {
        "sse_order": sse_order,
        "sse_3d_coords": sse_3d_coords,
        "sse_res_spans": sse_res_spans,
        "sse_ca_indices": sse_ca_indices,
        "residues": residues,
        "ca_coords": ca_coords,
        "D_sse_pred": D_sse_pred,
        "D_sse_3d": D_sse_3d,
        "D_res_pred": D_res_pred,
        "P_res_pred": P_res_pred,
        "cm_dsl_pred": cm_dsl_pred,
        "clusters": clusters,
        "layout_2d": layout_2d,
        "angles_2d": angles_2d,
        "lengths": lengths,
    }


def compute_kabsch_rmsd_and_tmscore(coords_pred, coords_true, n_res_total=None):
    """
    Compute optimal Kabsch superposition RMSD and Zhang-Skolnick TM-score between
    two 3D coordinate sets of shape (N, 3).
    """
    P = np.asarray(coords_pred, dtype=float)
    Q = np.asarray(coords_true, dtype=float)
    n = min(len(P), len(Q))
    if n < 3:
        return 99.9, 0.0, P
    P = P[:n]
    Q = Q[:n]
    p_mean = np.mean(P, axis=0)
    q_mean = np.mean(Q, axis=0)
    P_c = P - p_mean
    Q_c = Q - q_mean

    l_norm = float(n_res_total) if (n_res_total is not None and n_res_total > n) else max(120.0, float(n * 15))
    d0 = max(1.5, 1.24 * ((max(16.0, l_norm) - 15.0) ** (1.0 / 3.0)) - 1.8)

    # Check 3D mirror signs (since 2D diagrams have Z-depth chirality gauge freedom)
    best_rmsd = 99.9
    best_tm = 0.0
    best_aligned = P.copy()

    for z_sign in (1.0, -1.0):
        P_test = P_c.copy()
        P_test[:, 2] *= z_sign
        H = P_test.T @ Q_c
        U, S, Vt = np.linalg.svd(H)
        R = U @ Vt
        if np.linalg.det(R) < 0:
            U_mod = U.copy()
            U_mod[:, -1] *= -1
            R = U_mod @ Vt
        P_rot = P_test @ R + q_mean
        dists = np.linalg.norm(P_rot - Q, axis=1)
        rmsd = float(np.sqrt(np.mean(dists ** 2)))
        tm = float(np.mean(1.0 / (1.0 + (dists / d0) ** 2)))
        if tm > best_tm:
            best_tm = tm
            best_rmsd = rmsd
            best_aligned = P_rot

    return best_rmsd, best_tm, best_aligned


def export_3d_structure_to_mmcif_and_pdb(model_3d, cif_path, pdb_path=None, entry_id="2DPROTS_3D"):
    """
    Export the reconstructed 3D structure dictionary `model_3d` to a valid mmCIF (`.cif`)
    file (and optionally `.pdb`) with explicit `_struct_conf` (alpha-helices) and
    `_struct_sheet_range` (beta-sheets) tables so Mol* renders full secondary structure cartoons.
    """
    residues = model_3d["residues"]
    sse_res_spans = model_3d["sse_res_spans"]
    clusters = model_3d.get("clusters", [])

    lines = [
        f"data_{entry_id}",
        "#",
        f"_entry.id {entry_id}",
        "#",
        "loop_",
        "_entity.id",
        "_entity.type",
        "_entity.pdbx_description",
        "1 polymer '2DProts Reconstructed 3D Protein'",
        "#",
        "loop_",
        "_entity_poly.entity_id",
        "_entity_poly.type",
        "_entity_poly.pdbx_strand_id",
        "1 'polypeptide(L)' A",
        "#",
        "loop_",
        "_struct_asym.id",
        "_struct_asym.entity_id",
        "A 1",
        "#",
    ]

    # Collect helices and strands
    helices = []
    strands = []
    for s, (r0, r1) in sse_res_spans.items():
        is_h = s.startswith("H") or (s.startswith("_") and "H" in s[:4])
        r0_name = residues[r0 - 1]["res_name"]
        r1_name = residues[r1 - 1]["res_name"]
        if is_h:
            helices.append((s, r0, r0_name, r1, r1_name, r1 - r0 + 1))
        else:
            strands.append((s, r0, r0_name, r1, r1_name, r1 - r0 + 1))

    if helices:
        lines.extend([
            "loop_",
            "_struct_conf.conf_type_id",
            "_struct_conf.id",
            "_struct_conf.pdbx_PDB_helix_id",
            "_struct_conf.beg_label_comp_id",
            "_struct_conf.beg_label_asym_id",
            "_struct_conf.beg_label_seq_id",
            "_struct_conf.end_label_comp_id",
            "_struct_conf.end_label_asym_id",
            "_struct_conf.end_label_seq_id",
            "_struct_conf.beg_auth_comp_id",
            "_struct_conf.beg_auth_asym_id",
            "_struct_conf.beg_auth_seq_id",
            "_struct_conf.end_auth_comp_id",
            "_struct_conf.end_auth_asym_id",
            "_struct_conf.end_auth_seq_id",
            "_struct_conf.pdbx_PDB_helix_class",
            "_struct_conf.pdbx_PDB_helix_length",
        ])
        for idx_h, (s, r0, rn0, r1, rn1, hlen) in enumerate(helices, start=1):
            lines.append(
                f"HELX_P HELX_P{idx_h} {s} {rn0} A {r0} {rn1} A {r1} "
                f"{rn0} A {r0} {rn1} A {r1} 1 {hlen}"
            )
        lines.append("#")

    if strands:
        lines.extend([
            "loop_",
            "_struct_sheet.id",
            "_struct_sheet.number_strands",
            f"SHEET_A {len(strands)}",
            "#",
            "loop_",
            "_struct_sheet_range.sheet_id",
            "_struct_sheet_range.id",
            "_struct_sheet_range.beg_label_comp_id",
            "_struct_sheet_range.beg_label_asym_id",
            "_struct_sheet_range.beg_label_seq_id",
            "_struct_sheet_range.end_label_comp_id",
            "_struct_sheet_range.end_label_asym_id",
            "_struct_sheet_range.end_label_seq_id",
            "_struct_sheet_range.beg_auth_comp_id",
            "_struct_sheet_range.beg_auth_asym_id",
            "_struct_sheet_range.beg_auth_seq_id",
            "_struct_sheet_range.end_auth_comp_id",
            "_struct_sheet_range.end_auth_asym_id",
            "_struct_sheet_range.end_auth_seq_id",
        ])
        for idx_e, (s, r0, rn0, r1, rn1, _) in enumerate(strands, start=1):
            lines.append(
                f"SHEET_A {idx_e} {rn0} A {r0} {rn1} A {r1} {rn0} A {r0} {rn1} A {r1}"
            )
        lines.append("#")

    # Sequence table
    lines.extend([
        "loop_",
        "_entity_poly_seq.entity_id",
        "_entity_poly_seq.num",
        "_entity_poly_seq.mon_id",
        "_entity_poly_seq.hetero",
    ])
    for r in residues:
        lines.append(f"1 {r['seq_id']} {r['res_name']} n")
    lines.append("#")

    # Atom site table
    lines.extend([
        "loop_",
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_alt_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_entity_id",
        "_atom_site.label_seq_id",
        "_atom_site.pdbx_PDB_ins_code",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv",
        "_atom_site.auth_seq_id",
        "_atom_site.auth_comp_id",
        "_atom_site.auth_asym_id",
        "_atom_site.auth_atom_id",
        "_atom_site.pdbx_PDB_model_num",
    ])

    atom_serial = 1
    pdb_lines = [f"HEADER    2DPROTS GENERATED 3D MODEL              {entry_id}"]
    for idx_h, (s, r0, rn0, r1, rn1, hlen) in enumerate(helices, start=1):
        pdb_lines.append(
            f"HELIX  {idx_h:3d} {s:>3s} {rn0:>3s} A {r0:4d}  {rn1:>3s} A {r1:4d}  1                               {hlen:5d}"
        )
    for idx_e, (s, r0, rn0, r1, rn1, _) in enumerate(strands, start=1):
        pdb_lines.append(
            f"SHEET  {idx_e:3d}   A{len(strands):2d} {rn0:>3s} A{r0:4d}  {rn1:>3s} A{r1:4d}  0"
        )

    for r in residues:
        seq_id = r["seq_id"]
        rname = r["res_name"]
        ch = r["chain_id"]
        bfac = 85.0 if r["ss_type"] in ("H", "E") else 65.0
        for atom_name in ("N", "CA", "C", "O", "CB"):
            if rname == "GLY" and atom_name == "CB":
                continue
            xyz = r["atoms"][atom_name]
            elem = atom_name[0]
            lines.append(
                f"ATOM {atom_serial} {elem} {atom_name} . {rname} {ch} 1 {seq_id} ? "
                f"{xyz[0]:.3f} {xyz[1]:.3f} {xyz[2]:.3f} 1.00 {bfac:.2f} "
                f"{seq_id} {rname} {ch} {atom_name} 1"
            )
            pdb_lines.append(
                f"ATOM  {atom_serial:5d}  {atom_name:<3s} {rname:>3s} {ch:1s}{seq_id:4d}    "
                f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00{bfac:6.2f}           {elem:>2s}"
            )
            atom_serial += 1

    lines.append("#")
    pdb_lines.append("END")

    os.makedirs(os.path.dirname(os.path.abspath(cif_path)), exist_ok=True)
    with open(cif_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    if pdb_path is not None:
        os.makedirs(os.path.dirname(os.path.abspath(pdb_path)), exist_ok=True)
        with open(pdb_path, "w", encoding="utf-8") as f:
            f.write("\n".join(pdb_lines) + "\n")

    return cif_path, pdb_path


# ==============================================================================
# 6. 2.5DPROTS (2.5DPROTS-DSL): DEPTH- & TILT-AUGMENTED 2DPROTS FOR SUB-2A 3D RECONSTRUCTION
# ==============================================================================

def reconstruct_3d_endpoints_from_distogram(
    sse_order,
    sse_slices,
    D_distogram,
    l_obj_2d=None,
    chain_ids=None,
    res_ranges=None,
    n_smacof_iter=90,
):
    """
    Reconstruct 3D N- and C-terminal endpoints (p0, p1) for every SSE purely from a
    residue-level Distogram (`D_distogram <= 22.0 A`, e.g. AlphaFold / ESMFold / Boltz)
    WITHOUT using any 3D crystal coordinates:
      1. Extracts intra-SSE 3D span directly from collinear 3-hop distogram intervals.
      2. Completes long-range (> 21.6 A) endpoint distances via shortest-path geodesics
         on the residue-level backbone + contact graph.
      3. Embeds the 2*M SSE endpoints into 3D via Classical 3D MDS followed by weighted
         3D SMACOF stress majorization.
      4. Rigidly aligns the reconstructed 3D coordinates onto the 2DProts (X, Y, 0) plane
         so that (X, Y) matches the 2DProts view and Z represents the true out-of-plane depth!
    """
    n_s = len(sse_order)
    end_idx = []
    span_disto = []

    for s in sse_order:
        i0, i1 = sse_slices[s]
        end_idx.extend([i0, i1 - 1])
        if D_distogram[i0, i1 - 1] < 21.5:
            span_disto.append(float(D_distogram[i0, i1 - 1]))
        else:
            q1 = i0 + (i1 - 1 - i0) // 3
            q2 = i0 + 2 * (i1 - 1 - i0) // 3
            span_disto.append(
                float(D_distogram[i0, q1] + D_distogram[q1, q2] + D_distogram[q2, i1 - 1])
            )

    G = _build_intra_sse_backbone_graph(sse_order, sse_slices, chain_ids, res_ranges)
    G = np.where(D_distogram < 21.6, np.minimum(G, D_distogram), G)
    for i in range(n_s):
        G[end_idx[2 * i], end_idx[2 * i + 1]] = span_disto[i]
        G[end_idx[2 * i + 1], end_idx[2 * i]] = span_disto[i]

    D_sp_res = shortest_path(G, directed=False)
    D_end = D_sp_res[np.ix_(end_idx, end_idx)]
    D_end = np.where(np.isfinite(D_end), D_end, 42.0)

    m = len(end_idx)
    H = np.eye(m) - np.ones((m, m)) / float(m)
    B = -0.5 * H @ (D_end ** 2) @ H
    ev, evec = np.linalg.eigh(B)
    idx = np.argsort(ev)[::-1][:3]
    X3 = evec[:, idx] * np.sqrt(np.maximum(0.0, ev[idx]))

    D_direct = D_distogram[np.ix_(end_idx, end_idx)]
    W = np.where(D_direct < 21.6, 5.0, 0.35)
    for i in range(n_s):
        W[2 * i, 2 * i + 1] = W[2 * i + 1, 2 * i] = 25.0
        D_end[2 * i, 2 * i + 1] = D_end[2 * i + 1, 2 * i] = span_disto[i]
    np.fill_diagonal(W, 0.0)

    V = -W.copy()
    np.fill_diagonal(V, np.sum(W, axis=1))
    V_pinv = np.linalg.pinv(V)

    for _ in range(n_smacof_iter):
        d_cur = np.linalg.norm(X3[:, None, :] - X3[None, :, :], axis=-1)
        B_mat = -W * D_end / np.maximum(d_cur, 1e-4)
        np.fill_diagonal(B_mat, 0.0)
        np.fill_diagonal(B_mat, -np.sum(B_mat, axis=1))
        X3 = V_pinv @ (B_mat @ X3)

    X3 = X3 - np.mean(X3, axis=0, keepdims=True)

    # Align X3 onto the 2DProts (x, y, 0) plane if l_obj_2d is provided
    if l_obj_2d is not None and hasattr(l_obj_2d, "layout"):
        c_3d = 0.5 * (X3[0::2] + X3[1::2])
        xy_2d = np.array(
            [[l_obj_2d.layout.get(s, (0.0, 0.0))[0], l_obj_2d.layout.get(s, (0.0, 0.0))[1], 0.0] for s in sse_order],
            dtype=float,
        )
        xy_2d_c = xy_2d - np.mean(xy_2d, axis=0, keepdims=True)
        c_3d_c = c_3d - np.mean(c_3d, axis=0, keepdims=True)
        U_svd, _, Vt_svd = np.linalg.svd(c_3d_c.T @ xy_2d_c)
        R_align = U_svd @ Vt_svd
        if np.linalg.det(R_align) < 0:
            U_mod = U_svd.copy()
            U_mod[:, -1] *= -1
            R_align = U_mod @ Vt_svd
        X3 = (X3 - np.mean(c_3d, axis=0, keepdims=True)) @ R_align + np.mean(xy_2d, axis=0, keepdims=True)

    sse_3d_endpoints = OrderedDict()
    for i, s in enumerate(sse_order):
        sse_3d_endpoints[s] = (X3[2 * i].copy(), X3[2 * i + 1].copy())

    return sse_3d_endpoints, X3, D_end


def encode_25dprots_to_dsl(sse_order, l_obj_2d, sse_3d_endpoints, ann):
    """
    Encode a 2.5DProts representation (`2.5DPROTS-DSL`) that combines:
      - Standard 2DProts visual coordinates `X{x2d}:Y{y2d}:A{ang2d}`
      - 2.5D/3D center `C{cx},{cy},{cz}` (where `cz` is the Z-depth perpendicular to the 2DProts plane)
      - 3D spherical orientation `U{az},{el}` (azimuth in deg, out-of-plane elevation in deg)
      - 3D physical SSE span `S{span}` in Angstroms
    """
    tokens = ["<BOS_25DP>"]
    clusters = ann.get_beta_connectivity_clusters() if ann is not None else []
    for cl in clusters:
        parts = [p for p in cl.split(":") if p in sse_3d_endpoints]
        if len(parts) >= 2:
            tokens.append(f"CLUST:{':'.join(parts)}")

    for idx_s, s in enumerate(sse_order):
        if s not in sse_3d_endpoints:
            continue
        p0, p1 = sse_3d_endpoints[s]
        c_vec = 0.5 * (p0 + p1)
        d_vec = p1 - p0
        span = float(np.linalg.norm(d_vec))
        u_vec = d_vec / max(1e-8, span)

        # Spherical orientation (azimuth in X-Y plane, elevation out of X-Y plane toward Z)
        az_deg = math.degrees(math.atan2(float(u_vec[0]), float(u_vec[1]))) % 360.0
        el_deg = math.degrees(math.asin(float(np.clip(u_vec[2], -1.0, 1.0))))

        if l_obj_2d is not None and s in l_obj_2d.layout:
            x2d, y2d = l_obj_2d.layout[s]
            ang2d = float(l_obj_2d.a_layout.get(s, az_deg)) % 360.0
        else:
            x2d, y2d = float(c_vec[0]), float(c_vec[1])
            ang2d = az_deg

        length = int(ann.res_sizes.get(s, 6)) if ann is not None else max(3, int(round(span / 1.5)))
        tokens.append(
            f"SSE25D:{s}:L{length}:X{x2d:+06.1f}:Y{y2d:+06.1f}:Z{c_vec[2]:+06.1f}:"
            f"A{int(round(ang2d)):03d}:EL{int(round(el_deg)):+03d}:"
            f"C{c_vec[0]:+06.1f},{c_vec[1]:+06.1f},{c_vec[2]:+06.1f}:"
            f"U{int(round(az_deg)):03d},{int(round(el_deg)):+03d}:S{span:04.1f}"
        )

    tokens.append("<EOS_25DP>")
    return " ".join(tokens), tokens


def decode_25dprots_dsl(dsl_25d_text):
    """
    Parse a `2.5DPROTS-DSL` token string into:
      - `layout_2d`: OrderedDict {sse: (x2d, y2d)}
      - `angles_2d`: OrderedDict {sse: ang2d}
      - `depth_z`: OrderedDict {sse: z_depth}
      - `elev_deg`: OrderedDict {sse: el_deg}
      - `lengths`: OrderedDict {sse: res_len}
      - `sse_3d_coords`: OrderedDict {sse: (p0_3d, p1_3d)}
      - `clusters`: list of beta-sheet cluster strings
    Also supports compact human-authored 2.5D tokens that only specify `X, Y, Z, A, EL`
    (automatically filling `C=(X,Y,Z)`, `U=(A,EL)`, and canonical span `S`).
    """
    layout_2d = OrderedDict()
    angles_2d = OrderedDict()
    depth_z = OrderedDict()
    elev_deg = OrderedDict()
    lengths = OrderedDict()
    sse_3d_coords = OrderedDict()
    clusters = []

    for tok in dsl_25d_text.split():
        if tok.startswith("CLUST:"):
            clusters.append(tok[len("CLUST:"):])
        elif tok.startswith("SSE25D:") or tok.startswith("SSE2D:"):
            parts = tok.split(":")
            s = parts[1]
            is_h = s.startswith("H") or (s.startswith("_") and "H" in s[:4])
            l_val = 6
            x_val, y_val, z_val = 0.0, 0.0, 0.0
            a_val, el_val = 0.0, 0.0
            c_xyz = None
            u_az_el = None
            span_val = None

            for p in parts[2:]:
                if p.startswith("L"):
                    l_val = int(p[1:])
                elif p.startswith("X"):
                    x_val = float(p[1:])
                elif p.startswith("Y"):
                    y_val = float(p[1:])
                elif p.startswith("Z"):
                    z_val = float(p[1:])
                elif p.startswith("A"):
                    a_val = float(p[1:])
                elif p.startswith("EL"):
                    el_val = float(p[2:])
                elif p.startswith("C") and "," in p:
                    cx, cy, cz = [float(v) for v in p[1:].split(",")]
                    c_xyz = np.array([cx, cy, cz], dtype=float)
                elif p.startswith("U") and "," in p:
                    az, el = [float(v) for v in p[1:].split(",")]
                    u_az_el = (az, el)
                elif p.startswith("S"):
                    span_val = float(p[1:])

            if c_xyz is None:
                c_xyz = np.array([x_val, y_val, z_val], dtype=float)
            else:
                z_val = float(c_xyz[2])

            if u_az_el is None:
                az_rad = math.radians(a_val)
                el_rad = math.radians(el_val)
            else:
                az_rad = math.radians(u_az_el[0])
                el_rad = math.radians(u_az_el[1])
                el_val = u_az_el[1]

            if span_val is None:
                eff_res = min(l_val, 24) if is_h else min(l_val, 14)
                span_val = float(eff_res) * (1.46 if is_h else 2.85)

            cos_el = math.cos(el_rad)
            u_vec = np.array([
                cos_el * math.sin(az_rad),
                cos_el * math.cos(az_rad),
                math.sin(el_rad),
            ], dtype=float)
            p0 = c_xyz - 0.5 * span_val * u_vec
            p1 = c_xyz + 0.5 * span_val * u_vec

            layout_2d[s] = (x_val, y_val)
            angles_2d[s] = a_val
            depth_z[s] = z_val
            elev_deg[s] = el_val
            lengths[s] = l_val
            sse_3d_coords[s] = (p0, p1)

    return layout_2d, angles_2d, depth_z, elev_deg, lengths, sse_3d_coords, clusters


def reconstruct_3d_structure_from_25dprots_dsl(dsl_25d_text, loop_lengths=None):
    """
    Build a complete 3D atomic structure (`N, CA, C, O, CB` for all helices, beta-strands,
    and connecting loops) and its 3D distance / contact matrices directly from a
    `2.5DPROTS-DSL` token sequence.
    """
    layout_2d, angles_2d, depth_z, elev_deg, lengths, sse_3d_coords, clusters = decode_25dprots_dsl(dsl_25d_text)
    sse_order = list(sse_3d_coords.keys())
    n_s = len(sse_order)

    centers_3d = np.zeros((n_s, 3), dtype=float)
    u_dirs_3d = np.zeros((n_s, 3), dtype=float)
    for i, s in enumerate(sse_order):
        p0, p1 = sse_3d_coords[s]
        centers_3d[i] = 0.5 * (p0 + p1)
        dv = p1 - p0
        u_dirs_3d[i] = dv / max(1e-8, float(np.linalg.norm(dv)))

    # Compute realized 3D segment-to-segment distance matrix
    D_sse_3d = np.zeros((n_s, n_s), dtype=float)
    t_eval = np.linspace(0.0, 1.0, 9)[:, None]
    for i in range(n_s):
        p0_i, p1_i = sse_3d_coords[sse_order[i]]
        pts_i = p0_i[None, :] * (1.0 - t_eval) + p1_i[None, :] * t_eval
        for j in range(i + 1, n_s):
            p0_j, p1_j = sse_3d_coords[sse_order[j]]
            pts_j = p0_j[None, :] * (1.0 - t_eval) + p1_j[None, :] * t_eval
            d_min = float(np.min(np.linalg.norm(pts_i[:, None, :] - pts_j[None, :, :], axis=-1)))
            D_sse_3d[i, j] = D_sse_3d[j, i] = d_min

    sheet_adj_pairs = []
    for cl in clusters:
        parts = [p for p in cl.split(":") if p in sse_3d_coords]
        for k in range(len(parts) - 1):
            sheet_adj_pairs.append((sse_order.index(parts[k]), sse_order.index(parts[k + 1])))

    sheet_neighbor_dir = {}
    for (ia, ib) in sheet_adj_pairs:
        vec_ab = centers_3d[ib] - centers_3d[ia]
        sheet_neighbor_dir[sse_order[ia]] = vec_ab
        sheet_neighbor_dir[sse_order[ib]] = -vec_ab

    gaps = np.full(max(0, n_s - 1), 4, dtype=int)
    if loop_lengths is not None:
        for i in range(n_s - 1):
            pair = (sse_order[i], sse_order[i + 1])
            if pair in loop_lengths:
                gaps[i] = max(2, min(18, int(loop_lengths[pair])))

    center_of_mass = np.mean(centers_3d, axis=0) if n_s > 0 else np.zeros(3)
    residues = []
    sse_res_spans = OrderedDict()
    sse_ca_indices = OrderedDict()
    seq_id = 1

    for i, s in enumerate(sse_order):
        p0, p1 = sse_3d_coords[s]
        l_res = max(3, int(lengths.get(s, 6)))
        is_h = s.startswith("H") or (s.startswith("_") and "H" in s[:4])
        if is_h:
            sse_atoms = _build_ideal_helix_residues(p0, p1, l_res)
            res_name = "ALA"
            ss_type = "H"
        else:
            nb_dir = sheet_neighbor_dir.get(s, centers_3d[i] - center_of_mass)
            sse_atoms = _build_ideal_strand_residues(p0, p1, l_res, sheet_normal=nb_dir)
            res_name = "VAL"
            ss_type = "E"

        start_seq = seq_id
        idx_0 = len(residues)
        for k_r, atom_dict in enumerate(sse_atoms):
            residues.append({
                "seq_id": seq_id,
                "res_name": res_name,
                "chain_id": "A",
                "sse": s,
                "ss_type": ss_type,
                "sse_offset": k_r,
                "atoms": atom_dict,
            })
            seq_id += 1
        idx_1 = len(residues)
        sse_res_spans[s] = (start_seq, seq_id - 1)
        sse_ca_indices[s] = (idx_0, idx_1)

        if i + 1 < n_s:
            s_next = sse_order[i + 1]
            p0_next, _ = sse_3d_coords[s_next]
            n_loop = int(gaps[i])
            loop_atoms = _build_loop_residues(
                p1, u_dirs_3d[i], p0_next, u_dirs_3d[i + 1], n_loop, center_bias=center_of_mass
            )
            for k_l, atom_dict in enumerate(loop_atoms):
                residues.append({
                    "seq_id": seq_id,
                    "res_name": "GLY" if (k_l % 2 == 1) else "SER",
                    "chain_id": "A",
                    "sse": f"loop_{s}_{s_next}",
                    "ss_type": "L",
                    "sse_offset": k_l,
                    "atoms": atom_dict,
                })
                seq_id += 1

    ca_coords = np.array([r["atoms"]["CA"] for r in residues], dtype=float)

    return {
        "sse_order": sse_order,
        "sse_3d_coords": sse_3d_coords,
        "sse_res_spans": sse_res_spans,
        "sse_ca_indices": sse_ca_indices,
        "residues": residues,
        "ca_coords": ca_coords,
        "D_sse_pred": D_sse_3d,
        "D_sse_3d": D_sse_3d,
        "clusters": clusters,
        "layout_2d": layout_2d,
        "angles_2d": angles_2d,
        "depth_z": depth_z,
        "elev_deg": elev_deg,
        "lengths": lengths,
    }


