#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""
Comprehensive Benchmark of Amino-Acid (C-alpha Residue) Level Dimensionality Reduction
in 2DProts (`benchmark_dimred_aa.py`).

Evaluates the exact same 12 dimensionality reduction algorithms and 38 targets
(single domains & whole multi-chain/multi-domain complexes) as `benchmark_dimred.py`,
but computes the 3D -> 2D projection directly on the positions of individual amino acids
(C-alpha atoms) instead of secondary structure elements (SSEs).
"""

import os
import sys
import json
import time
import math
import logging
import warnings
from collections import OrderedDict

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

warnings.filterwarnings("ignore")

import numpy as np
import matplotlib
matplotlib.use("Agg")

import benchmark_dimred
import utils_dimred_aa
import superimpose2d

logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE_DIR, "benchmark_dimred_aa_results")
CACHE_DIR = os.path.join(BASE_DIR, "benchmark_dimred_results", "cache")
SVG_DIR = os.path.join(OUT_DIR, "svgs")
JSON_DIR = os.path.join(OUT_DIR, "layouts")

os.makedirs(CACHE_DIR, exist_ok=True)
os.makedirs(SVG_DIR, exist_ok=True)
os.makedirs(JSON_DIR, exist_ok=True)

METHODS = [
    ("evolution", "Evolution-AA (Stochastická)", "Víceúrovňová stochastická mutace pozic Cα + relaxace páteře"),
    ("pca", "PCA (Ortogonální 3D→2D)", "Lineární projekce 3D souřadnic Cα aminokyselin na 2 hlavní osy"),
    ("pcoa", "PCoA / Classical MDS", "Spektrální rozklad centrované Gramovy matice vzdáleností Cα"),
    ("mds", "Metric MDS (SMACOF)", "Metrické škálování matice 3D vzdáleností Cα (s landmark akcelerací)"),
    ("nmmds", "Non-Metric MDS", "Monotónní regrese pořadí vzdáleností aminokyselin (Kruskal)"),
    ("sammon", "Sammon / Local wMDS", "Kvazi-Newtonova (L-BFGS-B) optimalizace lokální chyby Cα <40 Å + vazeb páteře"),
    ("isomap", "Isomap (Geodetický)", "Nejkratší cesty po 3D k-NN grafu kontaktů a polypeptidové páteři + PCoA"),
    ("spectral", "Spectral (Laplacian)", "Vlastní vektory normalizovaného Laplaciánu kontaktního grafu Cα"),
    ("lle", "LLE (Locally Linear)", "Lokální lineární barycentrické rekonstrukční váhy sousedství Cα"),
    ("kamada_kawai", "Kamada-Kawai (Pružinový)", "Fyzikální model soustavy pružin na grafu vzdáleností aminokyselin"),
    ("tsne", "t-SNE", "Student-t stochastické vnoření sousedních aminokyselin (perplexity)"),
    ("umap", "UMAP", "Uniform Manifold Approximation nad vzdálenostmi jednotlivých aminokyselin"),
]

BENCHMARK_TARGETS = benchmark_dimred.BENCHMARK_TARGETS


def extract_ca_map_best_numbering(chain_obj, dom_data):
    """
    Extract C-alpha coordinates and 3-letter names from `chain_obj`, automatically
    choosing between `res.label_seq` (_atom_site.label_seq_id) and `res.seqid.num`
    (_atom_site.auth_seq_id) based on which numbering matches `dom_data['helices_residues']`
    and `dom_data['3dlayout']` with lower Kabsch RMSD.
    """
    import gemmi

    hr = dom_data.get("helices_residues", {})
    l3d = dom_data.get("3dlayout", {})

    candidates = []
    for mode in ("label_seq", "seqid"):
        ca_map = {}
        for res in chain_obj:
            if "CA" in res:
                key = res.label_seq if mode == "label_seq" else res.seqid.num
                if key is not None:
                    pos = res["CA"][0].pos
                    ca_map[int(key)] = (np.array([pos.x, pos.y, pos.z], dtype=float), res.name)
        if not ca_map:
            continue

        src, dst = [], []
        for s, pair in hr.items():
            if s not in l3d:
                continue
            r1, r2 = min(pair), max(pair)
            if r1 in ca_map and r2 in ca_map:
                src.extend([l3d[s][0], l3d[s][1]])
                dst.extend([ca_map[r1][0].tolist(), ca_map[r2][0].tolist()])
        if len(src) >= 6:
            sup = gemmi.superpose_positions(
                [gemmi.Position(*p) for p in src],
                [gemmi.Position(*p) for p in dst],
            )
            score = (len(src), -float(sup.rmsd))
        else:
            score = (len(src), -999.0)
        candidates.append((score, ca_map))

    if not candidates:
        return {}
    candidates.sort(key=lambda x: x[0], reverse=True)
    return candidates[0][1]


def extract_target_amino_acids(target):
    """
    Extract ordered amino-acid (C-alpha) records, 3D coordinates, pairwise distance matrix D,
    and covalent backbone bonds for a benchmark target (single domain or multi-domain complex).
    """
    pdb_id = target["pdb"]
    use_au = target.get("use_au", True)
    model = benchmark_dimred.fetch_pdb_structure(pdb_id, use_au=use_au)
    poly_chains = benchmark_dimred.get_polymer_chains(model)

    ndom = len(target["domains"])
    aa_records = []

    for d_idx, (ch_id, dom_id, dom_code) in enumerate(target["domains"]):
        dom_data = benchmark_dimred.fetch_domain_json(dom_code)
        hr = dom_data["helices_residues"]
        sizes = dom_data["size"]
        h_type = dom_data["h_type"]
        colors = dom_data.get("color", {})

        if ch_id in poly_chains:
            ca_full = extract_ca_map_best_numbering(poly_chains[ch_id], dom_data)
            ca_xyz = {r: v[0] for r, v in ca_full.items()}
            l3d = benchmark_dimred.transform_domain_3d_to_chain(dom_data, ca_xyz)
        else:
            ca_full = {}
            l3d = dom_data["3dlayout"]

        prefix = f"{d_idx + 1:02d}" if ndom >= 10 else f"{d_idx + 1}"
        dom_color = f"s{int(820.0 * d_idx / max(1, ndom))}"

        def rename_sse(s):
            if target["kind"] == "domain":
                return s
            if s.startswith("_"):
                return f"{s[:2]}{prefix}_{s[2:]}"
            return f"{s[0]}{prefix}_{s[1:]}"

        # Collect valid SSEs sorted by sequence start residue
        sse_spans = []
        for sse in l3d:
            if sse not in hr or sse not in sizes:
                continue
            ht_val = h_type.get(sse, 3)
            if ndom >= 4 and ht_val == 3:
                continue
            r1, r2 = min(hr[sse]), max(hr[sse])
            sse_spans.append((r1, r2, sse, ht_val))
        sse_spans.sort(key=lambda x: (x[0], x[1]))

        # Map each residue inside SSEs and intra-domain loops
        dom_res_map = OrderedDict()

        # 1. Include up to 3 N-terminal tail residues before the first SSE
        if sse_spans:
            first_r = sse_spans[0][0]
            first_col = dom_color if target["kind"] == "complex" else colors.get(sse_spans[0][2], "s100")
            for r in range(first_r - 3, first_r):
                if r in ca_full:
                    xyz, rname = ca_full[r]
                    dom_res_map[r] = {
                        "res_num": r,
                        "res_name": rname,
                        "chain_id": ch_id,
                        "dom_idx": d_idx,
                        "sse_label": None,
                        "ss_type": "L",
                        "h_type": 3,
                        "color": first_col,
                        "coord_3d_raw": xyz,
                        "coord_3d": xyz,
                    }

        # 2. Include SSE residues and intra-domain loops between consecutive SSEs
        for s_i, (r1, r2, sse_orig, ht_val) in enumerate(sse_spans):
            s_lbl = rename_sse(sse_orig)
            is_helix = sse_orig.startswith("H") or sse_orig.startswith("_H")
            ss_type = "H" if is_helix else "E"
            s_col = dom_color if target["kind"] == "complex" else colors.get(sse_orig, f"s{int(100 + 780 * s_i / max(1, len(sse_spans)))}")
            v_start = np.asarray(l3d[sse_orig][0], dtype=float)
            v_end = np.asarray(l3d[sse_orig][1], dtype=float)
            span_len = max(1, r2 - r1)

            for r in range(r1, r2 + 1):
                frac = float(r - r1) / float(span_len)
                axis_pt = (1.0 - frac) * v_start + frac * v_end
                if r in ca_full:
                    xyz_raw, rname = ca_full[r]
                    # For alpha-helices, blend 35% raw C-alpha with 65% helix axis so 3.6-res turns
                    # are smooth for manifold learners while preserving helical wave in 2D
                    if is_helix:
                        xyz_smooth = 0.35 * xyz_raw + 0.65 * axis_pt
                    else:
                        xyz_smooth = 0.80 * xyz_raw + 0.20 * axis_pt
                else:
                    xyz_raw = axis_pt
                    xyz_smooth = axis_pt
                    rname = "UNK"

                dom_res_map[r] = {
                    "res_num": r,
                    "res_name": rname,
                    "chain_id": ch_id,
                    "dom_idx": d_idx,
                    "sse_label": s_lbl,
                    "ss_type": ss_type,
                    "h_type": 1 if ht_val == 1 else (2 if ht_val == 2 else 3),
                    "color": s_col,
                    "coord_3d_raw": xyz_raw,
                    "coord_3d": xyz_smooth,
                }

            # Intra-domain loop to next SSE (if gap <= 30 residues and not large complex)
            max_loop_gap = 12 if ndom >= 6 else 30
            if s_i + 1 < len(sse_spans):
                next_r1 = sse_spans[s_i + 1][0]
                gap = next_r1 - r2 - 1
                if 0 < gap <= max_loop_gap:
                    for r in range(r2 + 1, next_r1):
                        if r in ca_full and r not in dom_res_map:
                            xyz_raw, rname = ca_full[r]
                            dom_res_map[r] = {
                                "res_num": r,
                                "res_name": rname,
                                "chain_id": ch_id,
                                "dom_idx": d_idx,
                                "sse_label": None,
                                "ss_type": "L",
                                "h_type": 3,
                                "color": s_col,
                                "coord_3d_raw": xyz_raw,
                                "coord_3d": xyz_raw,
                            }

        # 3. Include up to 3 C-terminal tail residues after the last SSE
        if sse_spans:
            last_r = sse_spans[-1][1]
            last_col = dom_color if target["kind"] == "complex" else colors.get(sse_spans[-1][2], "s880")
            for r in range(last_r + 1, last_r + 4):
                if r in ca_full and r not in dom_res_map:
                    xyz, rname = ca_full[r]
                    dom_res_map[r] = {
                        "res_num": r,
                        "res_name": rname,
                        "chain_id": ch_id,
                        "dom_idx": d_idx,
                        "sse_label": None,
                        "ss_type": "L",
                        "h_type": 3,
                        "color": last_col,
                        "coord_3d_raw": xyz,
                        "coord_3d": xyz,
                    }

        for r in sorted(dom_res_map.keys()):
            aa_records.append(dom_res_map[r])

    n_aa = len(aa_records)
    X3d = np.array([rec["coord_3d"] for rec in aa_records], dtype=float)
    X3d_raw = np.array([rec["coord_3d_raw"] for rec in aa_records], dtype=float)

    diff = X3d[:, None, :] - X3d[None, :, :]
    D = np.sqrt(np.sum(diff ** 2, axis=-1))

    # Build covalent backbone bonds (consecutive residues in same chain/domain with 3D distance < 10 Å)
    backbone_bonds = []
    for i in range(n_aa - 1):
        r_i = aa_records[i]
        r_j = aa_records[i + 1]
        if r_i["chain_id"] == r_j["chain_id"] and r_i["dom_idx"] == r_j["dom_idx"]:
            seq_gap = abs(r_j["res_num"] - r_i["res_num"])
            d3 = float(D[i, i + 1])
            if seq_gap == 1 and d3 < 5.5:
                backbone_bonds.append((i, i + 1, d3))
            elif seq_gap <= 6 and d3 < 16.0:
                backbone_bonds.append((i, i + 1, d3))

    return aa_records, X3d, X3d_raw, D, backbone_bonds


def evaluate_aa_metrics(D3, X2d, backbone_bonds):
    """
    Compute quantitative quality metrics directly at the Amino-Acid (C-alpha) level:
      1. aa_local_err : Mean absolute distance error [Å] for local amino-acid pairs (<40 Å)
      2. aa_stress_1  : Kruskal normalized global metric Stress-1 across all AA pairs
      3. aa_knn_pres  : Non-sequential 3D spatial k-NN contact preservation (k=8, |i-j|>2) [%]
      4. backbone_err : Mean absolute deviation of sequential peptide bond lengths in 2D vs 3D [Å]
    """
    n = D3.shape[0]
    if n < 4:
        return {"aa_local_err": 0.0, "aa_stress_1": 0.0, "aa_knn_pres": 100.0, "backbone_err": 0.0}

    # Subsample up to 500 residues for O(N^2) pairwise evaluation if N is huge
    if n > 500:
        idx = np.unique(np.round(np.linspace(0, n - 1, 500)).astype(int))
        D3_eval = D3[np.ix_(idx, idx)]
        X2_eval = X2d[idx]
    else:
        D3_eval = D3
        X2_eval = X2d

    m = D3_eval.shape[0]
    diff2 = X2_eval[:, None, :] - X2_eval[None, :, :]
    D2_eval = np.sqrt(np.sum(diff2 ** 2, axis=-1))

    triu = np.triu(np.ones((m, m), dtype=bool), k=1)
    local_mask = triu & (D3_eval < 40.0)
    if np.any(local_mask):
        aa_local_err = float(np.mean(np.abs(D2_eval[local_mask] - D3_eval[local_mask])))
    else:
        aa_local_err = float(np.mean(np.abs(D2_eval[triu] - D3_eval[triu])))

    denom = float(np.sum(D3_eval[triu] ** 2))
    aa_stress_1 = float(np.sqrt(np.sum((D2_eval[triu] - D3_eval[triu]) ** 2) / max(denom, 1e-9)))

    # Non-trivial 3D spatial k-NN preservation (k=6, excluding immediate sequence neighbors |i-j|<=2)
    k = min(6, max(1, m - 6))
    overlap_sum = 0.0
    valid_cnt = 0
    for i in range(m):
        r3 = D3_eval[i].copy()
        r2 = D2_eval[i].copy()
        for off in range(-2, 3):
            if 0 <= i + off < m:
                r3[i + off] = np.inf
                r2[i + off] = np.inf
        nn3 = set(np.argsort(r3)[:k])
        nn2 = set(np.argsort(r2)[:k])
        overlap_sum += len(nn3 & nn2) / float(k)
        valid_cnt += 1
    aa_knn_pres = float(overlap_sum / max(1, valid_cnt)) * 100.0

    # Covalent peptide bond preservation
    covalent = [(i, j, d0) for (i, j, d0) in (backbone_bonds or []) if d0 < 5.0]
    if covalent:
        errs = [
            abs(float(np.linalg.norm(X2d[i] - X2d[j])) - d0)
            for (i, j, d0) in covalent
        ]
        backbone_err = float(np.mean(errs))
    else:
        backbone_err = 0.0

    return {
        "aa_local_err": round(aa_local_err, 2),
        "aa_stress_1": round(aa_stress_1, 3),
        "aa_knn_pres": round(aa_knn_pres, 1),
        "backbone_err": round(backbone_err, 2),
    }


def _rgb_to_hex(rgb):
    r, g, b = [int(np.clip(round(c * 255), 0, 255)) for c in rgb[:3]]
    return f"#{r:02x}{g:02x}{b:02x}"


def _sse_hex(color_code):
    if not color_code or not isinstance(color_code, str):
        return "#64748b"
    c = color_code.strip()
    if c.startswith("s") and c[1:].isdigit():
        val = float(c[1:]) / 900.0
        cmap = matplotlib.colormaps["turbo"]
        return _rgb_to_hex(cmap(np.clip(val, 0.05, 0.95)))
    try:
        return _rgb_to_hex(matplotlib.colors.to_rgb(c))
    except Exception:
        return "#3b82f6"


def render_aa_svg(aa_records, X2d_aa, D3_aa, sse_layout, sse_angles, svg_path):
    """
    Render a rich, high-precision standalone SVG of the 2D Amino-Acid protein representation:
      - Polypeptide backbone trace (solid for continuous peptide bonds, dashed for loop gaps)
      - Secondary structure ribbons/arrows underneath the constituent amino acids
      - Inter-strand beta-sheet hydrogen-bond contacts (dashed lines for E-E pairs < 5.6 Å in 3D)
      - Individual C-alpha amino acid circles colored by SSE / rainbow / subunit, with tooltips
    """
    n = len(aa_records)
    if n == 0:
        return

    xs = X2d_aa[:, 0]
    ys = -X2d_aa[:, 1]  # Invert Y for SVG screen coordinates so orientation matches matplotlib
    pad = 8.0
    min_x, max_x = float(np.min(xs)) - pad, float(np.max(xs)) + pad
    min_y, max_y = float(np.min(ys)) - pad, float(np.max(ys)) + pad
    width = max(20.0, max_x - min_x)
    height = max(20.0, max_y - min_y)

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{min_x:.2f} {min_y:.2f} {width:.2f} {height:.2f}" width="100%" height="100%">',
        f'<rect x="{min_x:.2f}" y="{min_y:.2f}" width="{width:.2f}" height="{height:.2f}" fill="#ffffff"/>',
    ]

    # 1. Beta-sheet inter-strand H-bond contacts (E-E pairs from different SSEs with 3D distance < 5.6 Å)
    if n <= 1200:
        e_indices = [i for i, r in enumerate(aa_records) if r["ss_type"] == "E"]
        for a_idx in range(len(e_indices)):
            i = e_indices[a_idx]
            si = aa_records[i]["sse_label"]
            for b_idx in range(a_idx + 1, len(e_indices)):
                j = e_indices[b_idx]
                if si != aa_records[j]["sse_label"] and D3_aa[i, j] < 5.6:
                    lines.append(
                        f'<line x1="{xs[i]:.2f}" y1="{ys[i]:.2f}" x2="{xs[j]:.2f}" y2="{ys[j]:.2f}" '
                        f'stroke="#94a3b8" stroke-width="0.22" stroke-dasharray="0.5,0.5" opacity="0.65"/>'
                    )

    # 2. Translucent SSE envelopes (ribbons/arrows) fitted over the 2D amino-acid positions
    sse_groups = OrderedDict()
    for i, rec in enumerate(aa_records):
        lbl = rec["sse_label"]
        if lbl:
            sse_groups.setdefault(lbl, []).append(i)

    for lbl, idxs in sse_groups.items():
        if len(idxs) < 2:
            continue
        col_hex = _sse_hex(aa_records[idxs[0]]["color"])
        ss_t = aa_records[idxs[0]]["ss_type"]
        p_start = np.array([xs[idxs[0]], ys[idxs[0]]])
        p_end = np.array([xs[idxs[-1]], ys[idxs[-1]]])
        if ss_t == "H":
            lines.append(
                f'<line x1="{p_start[0]:.2f}" y1="{p_start[1]:.2f}" x2="{p_end[0]:.2f}" y2="{p_end[1]:.2f}" '
                f'stroke="{col_hex}" stroke-width="2.2" stroke-linecap="round" opacity="0.28"/>'
            )
        elif ss_t == "E":
            vec = p_end - p_start
            vlen = float(np.linalg.norm(vec))
            if vlen > 0.8:
                u = vec / vlen
                n_vec = np.array([-u[1], u[0]])
                tip = p_end + u * 0.6
                base = p_end - u * 0.9
                p_left = base + n_vec * 0.95
                p_right = base - n_vec * 0.95
                lines.append(
                    f'<line x1="{p_start[0]:.2f}" y1="{p_start[1]:.2f}" x2="{base[0]:.2f}" y2="{base[1]:.2f}" '
                    f'stroke="{col_hex}" stroke-width="1.3" stroke-linecap="butt" opacity="0.38"/>'
                )
                lines.append(
                    f'<polygon points="{tip[0]:.2f},{tip[1]:.2f} {p_left[0]:.2f},{p_left[1]:.2f} {p_right[0]:.2f},{p_right[1]:.2f}" '
                    f'fill="{col_hex}" opacity="0.55"/>'
                )

    # 3. Continuous polypeptide C-alpha backbone trace
    for i in range(n - 1):
        r1, r2 = aa_records[i], aa_records[i + 1]
        if r1["chain_id"] == r2["chain_id"] and r1["dom_idx"] == r2["dom_idx"]:
            seq_gap = abs(r2["res_num"] - r1["res_num"])
            if seq_gap == 1 and D3_aa[i, i + 1] < 5.5:
                is_ss = (r1["ss_type"] != "L") and (r1["sse_label"] == r2["sse_label"])
                col_hex = _sse_hex(r1["color"]) if is_ss else "#94a3b8"
                sw = "0.52" if is_ss else "0.32"
                op = "0.90" if is_ss else "0.65"
                lines.append(
                    f'<line x1="{xs[i]:.2f}" y1="{ys[i]:.2f}" x2="{xs[i+1]:.2f}" y2="{ys[i+1]:.2f}" '
                    f'stroke="{col_hex}" stroke-width="{sw}" stroke-linecap="round" opacity="{op}"/>'
                )
            elif seq_gap <= 15:
                lines.append(
                    f'<line x1="{xs[i]:.2f}" y1="{ys[i]:.2f}" x2="{xs[i+1]:.2f}" y2="{ys[i+1]:.2f}" '
                    f'stroke="#cbd5e1" stroke-width="0.22" stroke-dasharray="0.6,0.6" opacity="0.55"/>'
                )

    # 4. Individual Amino Acid (C-alpha) positions
    r_ss = "0.52" if n <= 600 else "0.38"
    r_lp = "0.30" if n <= 600 else "0.22"
    for i, rec in enumerate(aa_records):
        ss_t = rec["ss_type"]
        if ss_t == "L":
            lines.append(
                f'<circle cx="{xs[i]:.2f}" cy="{ys[i]:.2f}" r="{r_lp}" fill="#94a3b8" opacity="0.72">'
                f'<title>{rec["res_name"]} {rec["res_num"]} (Chain {rec["chain_id"]}, Loop)</title></circle>'
            )
        else:
            col_hex = _sse_hex(rec["color"])
            stroke_col = "#0f172a" if ss_t == "E" else "#1e293b"
            lines.append(
                f'<circle cx="{xs[i]:.2f}" cy="{ys[i]:.2f}" r="{r_ss}" fill="{col_hex}" '
                f'stroke="{stroke_col}" stroke-width="0.10" opacity="0.95">'
                f'<title>{rec["res_name"]} {rec["res_num"]} ({rec["sse_label"]}, {ss_t})</title></circle>'
            )

    lines.append("</svg>")
    with open(svg_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def evaluate_induced_sse_metrics(ann, sse_layout, sse_to_indices, X2d_aa, h_type_map):
    """
    Evaluate 2DProts SSE-level metrics on the induced 2D SSE centers `sse_layout`:
    Rescales the induced 2D SSE centers so their weighted local RMS distance (<40 Å)
    matches `ann.get_distances()` (exact same normalization as `utils_dimred._rescale_to_target_distances`
    in `benchmark_dimred.py`), and then calls `benchmark_dimred.evaluate_layout_metrics`!
    """
    shown_sses = [
        s for s in sse_layout
        if s in ann.res_sizes and h_type_map.get(s, 3) in (1, 2)
    ]
    n = len(shown_sses)
    if n < 3:
        return {"local_err": 0.0, "stress_1": 0.0, "knn_pres": 100.0, "conflicts": 0, "n_sses": n}

    d3d_dict = ann.get_distances()
    D3 = np.zeros((n, n), dtype=float)
    D2_raw = np.zeros((n, n), dtype=float)
    weights = np.array([max(1.0, float(ann.res_sizes.get(s, 5))) for s in shown_sses])

    for i in range(n):
        si = shown_sses[i]
        pi = sse_layout[si]
        for j in range(i + 1, n):
            sj = shown_sses[j]
            pj = sse_layout[sj]
            D3[i, j] = D3[j, i] = float(d3d_dict[si][sj])
            D2_raw[i, j] = D2_raw[j, i] = math.hypot(pi[0] - pj[0], pi[1] - pj[1])

    mask = (D3 < 40.0) & np.triu(np.ones((n, n), dtype=bool), k=1)
    if not np.any(mask):
        mask = np.triu(np.ones((n, n), dtype=bool), k=1)
    W = (weights[:, None] + weights[None, :])[mask]
    rms_3d = float(np.sqrt(np.sum(W * (D3[mask] ** 2)) / max(1e-9, np.sum(W))))
    rms_2d = float(np.sqrt(np.sum(W * (D2_raw[mask] ** 2)) / max(1e-9, np.sum(W))))
    scale = (rms_3d / rms_2d) if (rms_2d > 1e-6 and rms_3d > 1e-6) else 1.0

    scaled_layout = {s: (sse_layout[s][0] * scale, sse_layout[s][1] * scale) for s in sse_layout}

    class _Holder:
        pass

    h = _Holder()
    h.layout = scaled_layout
    h.h_type = h_type_map
    return benchmark_dimred.evaluate_layout_metrics(ann, h)


def run_single_target_aa_method(target, method_id, ann, stats, ref_layout, aa_records, X3d, X3d_raw, D3_aa, backbone_bonds):
    """
    Run one (target, method_id) pair at the Amino-Acid level, align to the canonical reference
    frame `ref_layout`, save AA + induced SSE JSON and SVG, and compute both AA and SSE metrics.
    """
    tid = target["id"]
    json_path = os.path.join(JSON_DIR, f"{tid}_{method_id}.json")
    svg_path = os.path.join(SVG_DIR, f"{tid}_{method_id}.svg")

    t0 = time.perf_counter()
    X2d_aa = utils_dimred_aa.embed_amino_acids(
        D3_aa,
        X3d,
        X3d_raw=X3d_raw,
        backbone_bonds=backbone_bonds,
        method=method_id,
    )
    sse_layout, sse_angles, sse_to_indices = utils_dimred_aa.derive_sse_layout_from_aa(
        aa_records, X2d_aa, ann.vectors
    )
    elapsed_ms = (time.perf_counter() - t0) * 1000.0

    # Rigidly align the 2D amino-acid embedding (and induced SSE centers) onto `ref_layout`
    # so all 12 methods share the exact same global orientation as the SSE benchmark
    common = [k for k in sse_layout if k in ref_layout]
    if len(common) >= 2:
        A = np.array([[sse_layout[k][0], sse_layout[k][1]] for k in common], dtype=float).T
        B = np.array([[ref_layout[k][0], ref_layout[k][1]] for k in common], dtype=float).T
        try:
            R, t = superimpose2d.optimal_rotation_translation(A, B, allow_mirror=True)
            X2d_aa = (R @ X2d_aa.T + t).T
            sse_layout, sse_angles, sse_to_indices = utils_dimred_aa.derive_sse_layout_from_aa(
                aa_records, X2d_aa, ann.vectors
            )
        except Exception:
            pass

    # Compute Amino-Acid level metrics
    aa_metrics = evaluate_aa_metrics(D3_aa, X2d_aa, backbone_bonds)

    # Compute induced SSE-level metrics
    h_type_map = {s: stats.determine_types(s, ann.res_sizes) for s in ann.sses}
    sse_metrics = evaluate_induced_sse_metrics(ann, sse_layout, sse_to_indices, X2d_aa, h_type_map)

    # Save combined JSON (both SSE layout and full AA 2D trace for panel plotting)
    jdata = {
        "protein": tid,
        "method": method_id,
        "n_aa": len(aa_records),
        "layout": {k: [round(v[0], 3), round(v[1], 3)] for k, v in sse_layout.items()},
        "angles": {k: round(v, 2) for k, v in sse_angles.items()},
        "size": {k: ann.res_sizes.get(k, 5) for k in sse_layout},
        "h_type": h_type_map,
        "color": dict(ann.color),
        "aa_coords": [[round(float(p[0]), 3), round(float(p[1]), 3)] for p in X2d_aa],
        "aa_meta": [
            [r["res_num"], r["res_name"], r["chain_id"], r["dom_idx"], r["sse_label"] or "", r["ss_type"], r["color"]]
            for r in aa_records
        ],
    }
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(jdata, f)

    render_aa_svg(aa_records, X2d_aa, D3_aa, sse_layout, sse_angles, svg_path)

    return {
        "time_ms": round(elapsed_ms, 1),
        "n_aa": len(aa_records),
        "aa_local_err": aa_metrics["aa_local_err"],
        "aa_stress_1": aa_metrics["aa_stress_1"],
        "aa_knn_pres": aa_metrics["aa_knn_pres"],
        "backbone_err": aa_metrics["backbone_err"],
        # Induced SSE-level metrics for 1-to-1 comparison with benchmark_dimred.py
        "local_err": sse_metrics["local_err"],
        "stress_1": sse_metrics["stress_1"],
        "knn_pres": sse_metrics["knn_pres"],
        "conflicts": sse_metrics["conflicts"],
        "n_sses": sse_metrics["n_sses"],
    }


def main():
    only_targets = sys.argv[1:] if len(sys.argv) > 1 else None
    summary_path = os.path.join(OUT_DIR, "results_summary.json")
    prev_results = {}
    if os.path.exists(summary_path):
        try:
            with open(summary_path, "r", encoding="utf-8") as f:
                prev_results = json.load(f)
        except Exception:
            prev_results = {}

    print("Warming up Numba/UMAP JIT for Amino-Acid DimRed...", flush=True)
    _d_warm = np.eye(8, dtype=float)
    for i in range(8):
        for j in range(i + 1, 8):
            _d_warm[i, j] = _d_warm[j, i] = 3.8 * abs(i - j)
    _x_warm = np.column_stack([np.arange(8) * 3.8, np.zeros(8), np.zeros(8)])
    _ = utils_dimred_aa.embed_amino_acids(_d_warm, _x_warm, method="umap")

    all_results = dict(prev_results)

    for target in BENCHMARK_TARGETS:
        tid = target["id"]
        if only_targets and tid not in only_targets:
            continue
        if not only_targets and tid in all_results and len(all_results[tid].get("methods", {})) == len(METHODS):
            print(f"Skipping already-computed [{tid}] {target['name']}", flush=True)
            continue

        print(f"\n=== [AA-DimRed: {tid}] {target['name']} ({target['kind']}) ===", flush=True)
        ann, stats, ref_layout, ref_angles, chain_id, dom_id = benchmark_dimred.build_target_annotation(target)
        aa_records, X3d, X3d_raw, D3_aa, backbone_bonds = extract_target_amino_acids(target)
        n_clusters = len(ann.get_beta_connectivity_clusters())
        print(f"    Amino Acids (Cα): {len(aa_records)} | SSEs: {len(ann.sses)} | Clusters: {n_clusters}", flush=True)

        all_results[tid] = {
            "target": target,
            "n_aa": len(aa_records),
            "n_sses": len(ann.sses),
            "n_clusters": n_clusters,
            "methods": {},
        }
        for method_id, method_title, _ in METHODS:
            m = run_single_target_aa_method(
                target, method_id, ann, stats, ref_layout, aa_records, X3d, X3d_raw, D3_aa, backbone_bonds
            )
            all_results[tid]["methods"][method_id] = m
            print(
                f"    - {method_id:13s}: time={m['time_ms']:6.1f} ms | "
                f"AA_err={m['aa_local_err']:5.2f} Å | AA_str={m['aa_stress_1']:.3f} | AA_kNN={m['aa_knn_pres']:4.1f}% | "
                f"SSE_err={m['local_err']:5.2f} Å | SSE_kNN={m['knn_pres']:4.1f}%",
                flush=True,
            )

        with open(summary_path, "w", encoding="utf-8") as f:
            json.dump(all_results, f, indent=2, ensure_ascii=False)

    print(f"\nSaved Amino-Acid DimRed summary to {summary_path} ({len(all_results)} targets)", flush=True)


if __name__ == "__main__":
    main()
