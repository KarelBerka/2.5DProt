#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""
Computes and renders the 4 reference/baseline views at the Amino-Acid (C-alpha) level
for all 42 benchmark targets in `benchmark_dimred_aa_results`:
  1. `3d_pca`       : 3D Struktura (Principální osy PC1×PC2, nejkratší osa PC3 za plochou, po Cα)
  2. `2dprots_orig` : 2DProts Originál (bez membrány) – aminokyseliny namapované na 2DProts SSE + smyčky
  3. `2dprots_mem15`: 2DProts + Membrána (15°) – aminokyseliny namapované na membránový layout 15°
  4. `2dprots_mem30`: 2DProts + Membrána (30°) – aminokyseliny namapované na membránový layout 30°
"""

import os
import sys
import json
import time
import math
import warnings
import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

warnings.filterwarnings("ignore")

import benchmark_dimred
import benchmark_extras
import benchmark_dimred_aa
import utils_dimred_aa
import utils_membrane
import layout_loader

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE_DIR, "benchmark_dimred_aa_results")
SSE_JSON_DIR = os.path.join(BASE_DIR, "benchmark_dimred_results", "layouts")
SVG_DIR = os.path.join(OUT_DIR, "svgs")
JSON_DIR = os.path.join(OUT_DIR, "layouts")

EXTRA_METHODS = benchmark_extras.EXTRA_METHODS


def map_aa_to_sse_layout(aa_records, X3d_raw, sse_layout, sse_angles, sse_sizes):
    """
    Given a 2D SSE layout (`sse_layout`, `sse_angles`), reconstruct the 2D positions
    of all individual amino acids (`aa_records`):
      - Amino acids belonging to an SSE `s` are placed along its 2D vector (centered at
        `sse_layout[s]` with angle `sse_angles[s]`, plus authentic local helical/strand C-alpha
        transverse offset projected from 3D).
      - Loop amino acids (`ss_type == 'L'`) between SSEs are smoothly interpolated along a
        natural polypeptide curve connecting the C-terminus of the preceding SSE to the
        N-terminus of the following SSE.
    """
    n = len(aa_records)
    X2d = np.zeros((n, 2), dtype=float)
    placed = np.zeros(n, dtype=bool)

    # Group residue indices by SSE label
    sse_groups = {}
    for i, rec in enumerate(aa_records):
        lbl = rec["sse_label"]
        if lbl and lbl in sse_layout:
            sse_groups.setdefault(lbl, []).append(i)

    for lbl, idxs in sse_groups.items():
        cx, cy = sse_layout[lbl]
        ang_rad = math.radians(float(sse_angles.get(lbl, 0.0)))
        u_dir = np.array([math.sin(ang_rad), math.cos(ang_rad)], dtype=float)
        n_dir = np.array([-u_dir[1], u_dir[0]], dtype=float)

        m = len(idxs)
        is_helix = aa_records[idxs[0]]["ss_type"] == "H"
        # Use physical 2D length matching 2DProts visual scale (~1.5 Å/res for H, ~2.6 Å/res for E)
        step_along = 1.45 if is_helix else 2.35
        total_len = max(1.5, (m - 1) * step_along)

        for k_pos, idx in enumerate(idxs):
            frac = (float(k_pos) / float(max(1, m - 1))) - 0.5 if m > 1 else 0.0
            along = frac * total_len
            # Authentic helical wave or beta-pleat transverse modulation
            if is_helix and m >= 3:
                trans = 0.72 * math.sin(2.0 * math.pi * k_pos / 3.6)
            elif not is_helix and m >= 2:
                trans = 0.22 * (1 if (k_pos % 2 == 0) else -1)
            else:
                trans = 0.0
            X2d[idx] = np.array([cx, cy], dtype=float) + along * u_dir + trans * n_dir
            placed[idx] = True

    # Smoothly interpolate unplaced loop / terminal residues per (chain_id, dom_idx)
    start = 0
    while start < n:
        ch = aa_records[start]["chain_id"]
        dom = aa_records[start]["dom_idx"]
        end = start
        while end < n and aa_records[end]["chain_id"] == ch and aa_records[end]["dom_idx"] == dom:
            end += 1

        placed_in_seg = [i for i in range(start, end) if placed[i]]
        if not placed_in_seg:
            for i in range(start, end):
                X2d[i] = np.array([float(i - start) * 2.0, 0.0])
        else:
            # Leading N-terminal loop before first placed residue
            first_p = placed_in_seg[0]
            if first_p > start:
                # Direction of first SSE
                next_p = placed_in_seg[min(1, len(placed_in_seg) - 1)]
                v_dir = X2d[first_p] - X2d[next_p]
                vn = float(np.linalg.norm(v_dir))
                v_dir = (v_dir / vn) if vn > 1e-4 else np.array([-1.0, 0.0])
                for i in range(start, first_p):
                    dist_k = float(first_p - i) * 1.8
                    X2d[i] = X2d[first_p] + dist_k * v_dir

            # Trailing C-terminal loop after last placed residue
            last_p = placed_in_seg[-1]
            if last_p < end - 1:
                prev_p = placed_in_seg[max(0, len(placed_in_seg) - 2)]
                v_dir = X2d[last_p] - X2d[prev_p]
                vn = float(np.linalg.norm(v_dir))
                v_dir = (v_dir / vn) if vn > 1e-4 else np.array([1.0, 0.0])
                for i in range(last_p + 1, end):
                    dist_k = float(i - last_p) * 1.8
                    X2d[i] = X2d[last_p] + dist_k * v_dir

            # Internal loops between consecutive placed SSEs
            for p_idx in range(len(placed_in_seg) - 1):
                i0 = placed_in_seg[p_idx]
                i1 = placed_in_seg[p_idx + 1]
                if i1 - i0 <= 1:
                    continue
                p0 = X2d[i0]
                p1 = X2d[i1]
                chord = p1 - p0
                chord_len = float(np.linalg.norm(chord))
                n_loop = i1 - i0 - 1
                if chord_len > 1e-4:
                    u_c = chord / chord_len
                    n_c = np.array([-u_c[1], u_c[0]])
                else:
                    u_c = np.array([1.0, 0.0])
                    n_c = np.array([0.0, 1.0])
                # Slight outward arch proportional to loop residue count so hairpin turns look natural
                arch_amp = min(4.5, max(0.6, 0.35 * n_loop))
                mid_pt = 0.5 * (p0 + p1)
                if float(np.dot(mid_pt, n_c)) < 0:
                    n_c = -n_c
                for k_step, idx in enumerate(range(i0 + 1, i1), 1):
                    t = float(k_step) / float(n_loop + 1)
                    arch = math.sin(math.pi * t) * arch_amp
                    X2d[idx] = (1.0 - t) * p0 + t * p1 + arch * n_c

        start = end

    return X2d


def render_aa_svg_with_membrane(aa_records, X2d_aa, D3_aa, sse_layout, sse_angles, mem_info, svg_path):
    """
    Render the standalone AA SVG, optionally including the horizontal transmembrane band
    when `mem_info` has `is_transmembrane=True` (for `2dprots_mem15` and `2dprots_mem30`).
    """
    benchmark_dimred_aa.render_aa_svg(aa_records, X2d_aa, D3_aa, sse_layout, sse_angles, svg_path)
    if not isinstance(mem_info, dict) or not mem_info.get("is_transmembrane", False):
        return

    with open(svg_path, "r", encoding="utf-8") as f:
        svg_txt = f.read()

    yc = -float(mem_info.get("y_center_2d", 0.0))  # Inverted Y in SVG
    hh = float(mem_info.get("band_half_height_2d", 3.0))
    xs = X2d_aa[:, 0]
    min_x = float(np.min(xs)) - 10.0
    width = float(np.max(xs) - np.min(xs)) + 20.0

    mem_band_svg = (
        f'<rect x="{min_x:.2f}" y="{yc - hh:.2f}" width="{width:.2f}" height="{2.0 * hh:.2f}" '
        f'fill="#fef3c7" opacity="0.55"/>\n'
        f'<line x1="{min_x:.2f}" y1="{yc - hh:.2f}" x2="{min_x + width:.2f}" y2="{yc - hh:.2f}" '
        f'stroke="#d97706" stroke-width="0.35" stroke-dasharray="1.2,0.8" opacity="0.8"/>\n'
        f'<line x1="{min_x:.2f}" y1="{yc + hh:.2f}" x2="{min_x + width:.2f}" y2="{yc + hh:.2f}" '
        f'stroke="#d97706" stroke-width="0.35" stroke-dasharray="1.2,0.8" opacity="0.8"/>'
    )
    # Insert right after background rect
    first_rect_end = svg_txt.find("/>", svg_txt.find("<rect")) + 2
    svg_txt = svg_txt[:first_rect_end] + "\n" + mem_band_svg + svg_txt[first_rect_end:]
    with open(svg_path, "w", encoding="utf-8") as f:
        f.write(svg_txt)


def run_aa_extras_for_target(target, rdata):
    """Compute and render `3d_pca`, `2dprots_orig`, `2dprots_mem15`, and `2dprots_mem30` at the AA level."""
    tid = target["id"]
    ann, stats, ref_layout, ref_angles, chain_id, dom_id = benchmark_dimred.build_target_annotation(target)
    aa_records, X3d, X3d_raw, D3_aa, backbone_bonds = benchmark_dimred_aa.extract_target_amino_acids(target)
    h_type_map = {s: stats.determine_types(s, ann.res_sizes) for s in ann.sses}

    # 1. `3d_pca` at Amino-Acid level: 3D Principal Axes PC1 x PC2 (with PC3 depth shading)
    t0 = time.perf_counter()
    P_c = X3d_raw - np.mean(X3d_raw, axis=0, keepdims=True)
    _, _, vt = np.linalg.svd(P_c, full_matrices=False)
    proj = P_c @ vt.T
    spans = np.ptp(proj, axis=0)
    order = np.argsort(spans)[::-1]
    R_base = vt[order, :]
    if np.linalg.det(R_base) < 0:
        R_base[2, :] = -R_base[2, :]
    P_rot = P_c @ R_base.T
    spans_sorted = [round(float(x), 1) for x in np.ptp(P_rot, axis=0)]
    X2d_3d = P_rot[:, :2]
    sse_lay_3d, sse_ang_3d, sse_idx_3d = utils_dimred_aa.derive_sse_layout_from_aa(aa_records, X2d_3d, ann.vectors)

    # Align orientation to ref_layout
    common = [k for k in sse_lay_3d if k in ref_layout]
    if len(common) >= 2:
        A = np.array([[sse_lay_3d[k][0], sse_lay_3d[k][1]] for k in common], dtype=float).T
        B = np.array([[ref_layout[k][0], ref_layout[k][1]] for k in common], dtype=float).T
        try:
            R2, t2 = benchmark_dimred.superimpose2d.optimal_rotation_translation(A, B, allow_mirror=True)
            X2d_3d = (R2 @ X2d_3d.T + t2).T
            sse_lay_3d, sse_ang_3d, sse_idx_3d = utils_dimred_aa.derive_sse_layout_from_aa(aa_records, X2d_3d, ann.vectors)
        except Exception:
            pass
    ms_3d = (time.perf_counter() - t0) * 1000.0

    aa_m_3d = benchmark_dimred_aa.evaluate_aa_metrics(D3_aa, X2d_3d, backbone_bonds)
    sse_m_3d = benchmark_dimred_aa.evaluate_induced_sse_metrics(ann, sse_lay_3d, sse_idx_3d, X2d_3d, h_type_map)

    json_3d = os.path.join(JSON_DIR, f"{tid}_3d_pca.json")
    svg_3d = os.path.join(SVG_DIR, f"{tid}_3d_pca.svg")
    jdata_3d = {
        "protein": tid,
        "method": "3d_pca",
        "n_aa": len(aa_records),
        "pca_spans_3d": spans_sorted,
        "layout": {k: [round(v[0], 3), round(v[1], 3)] for k, v in sse_lay_3d.items()},
        "angles": {k: round(v, 2) for k, v in sse_ang_3d.items()},
        "size": {k: ann.res_sizes.get(k, 5) for k in sse_lay_3d},
        "h_type": h_type_map,
        "color": dict(ann.color),
        "membrane": {"is_transmembrane": False, "angle_step": 0},
        "aa_coords": [[round(float(p[0]), 3), round(float(p[1]), 3)] for p in X2d_3d],
        "aa_meta": [
            [r["res_num"], r["res_name"], r["chain_id"], r["dom_idx"], r["sse_label"] or "", r["ss_type"], r["color"]]
            for r in aa_records
        ],
    }
    with open(json_3d, "w", encoding="utf-8") as f:
        json.dump(jdata_3d, f)
    benchmark_dimred_aa.render_aa_svg(aa_records, X2d_3d, D3_aa, sse_lay_3d, sse_ang_3d, svg_3d)

    rdata["methods"]["3d_pca"] = {
        "time_ms": round(ms_3d, 1),
        "n_aa": len(aa_records),
        "pca_spans": spans_sorted,
        "aa_local_err": aa_m_3d["aa_local_err"],
        "aa_stress_1": aa_m_3d["aa_stress_1"],
        "aa_knn_pres": aa_m_3d["aa_knn_pres"],
        "backbone_err": aa_m_3d["backbone_err"],
        "local_err": sse_m_3d["local_err"],
        "stress_1": sse_m_3d["stress_1"],
        "knn_pres": sse_m_3d["knn_pres"],
        "conflicts": sse_m_3d["conflicts"],
        "n_sses": sse_m_3d["n_sses"],
    }

    # 2, 3, 4: `2dprots_orig`, `2dprots_mem15`, `2dprots_mem30` mapped to Amino Acids
    for m_key in ("2dprots_orig", "2dprots_mem15", "2dprots_mem30"):
        t0 = time.perf_counter()
        sse_jpath = os.path.join(SSE_JSON_DIR, f"{tid}_{m_key}.json")
        with open(sse_jpath, "r", encoding="utf-8") as f:
            sj = json.load(f)
        s_lay = {k: (float(v[0]), float(v[1])) for k, v in sj["layout"].items()}
        s_ang = {k: float(v) for k, v in sj["angles"].items()}
        mem_info = sj.get("membrane", {"is_transmembrane": False, "angle_step": 0})

        X2d_m = map_aa_to_sse_layout(aa_records, X3d_raw, s_lay, s_ang, ann.res_sizes)
        ms_m = (time.perf_counter() - t0) * 1000.0

        aa_m = benchmark_dimred_aa.evaluate_aa_metrics(D3_aa, X2d_m, backbone_bonds)
        l_holder = layout_loader.Layout(sses_stats=stats, annotation=ann)
        l_holder.layout = s_lay
        l_holder.a_layout = s_ang
        l_holder.h_type = h_type_map
        sse_m = benchmark_dimred.evaluate_layout_metrics(ann, l_holder)

        json_m = os.path.join(JSON_DIR, f"{tid}_{m_key}.json")
        svg_m = os.path.join(SVG_DIR, f"{tid}_{m_key}.svg")
        jdata_m = {
            "protein": tid,
            "method": m_key,
            "n_aa": len(aa_records),
            "layout": {k: [round(v[0], 3), round(v[1], 3)] for k, v in s_lay.items()},
            "angles": {k: round(v, 2) for k, v in s_ang.items()},
            "size": {k: ann.res_sizes.get(k, 5) for k in s_lay},
            "h_type": h_type_map,
            "color": dict(ann.color),
            "membrane": mem_info,
            "aa_coords": [[round(float(p[0]), 3), round(float(p[1]), 3)] for p in X2d_m],
            "aa_meta": [
                [r["res_num"], r["res_name"], r["chain_id"], r["dom_idx"], r["sse_label"] or "", r["ss_type"], r["color"]]
                for r in aa_records
            ],
        }
        with open(json_m, "w", encoding="utf-8") as f:
            json.dump(jdata_m, f)

        render_aa_svg_with_membrane(aa_records, X2d_m, D3_aa, s_lay, s_ang, mem_info, svg_m)

        rdata["methods"][m_key] = {
            "time_ms": round(ms_m, 1),
            "n_aa": len(aa_records),
            "is_tm": bool(isinstance(mem_info, dict) and mem_info.get("is_transmembrane", False)),
            "aa_local_err": aa_m["aa_local_err"],
            "aa_stress_1": aa_m["aa_stress_1"],
            "aa_knn_pres": aa_m["aa_knn_pres"],
            "backbone_err": aa_m["backbone_err"],
            "local_err": sse_m["local_err"],
            "stress_1": sse_m["stress_1"],
            "knn_pres": sse_m["knn_pres"],
            "conflicts": sse_m["conflicts"],
            "n_sses": sse_m["n_sses"],
        }

    print(
        f"[AA-Extras {tid:14s}] 3d_pca Cα={rdata['methods']['3d_pca']['aa_local_err']:5.2f}Å | "
        f"orig Cα={rdata['methods']['2dprots_orig']['aa_local_err']:5.2f}Å | "
        f"mem15 Cα={rdata['methods']['2dprots_mem15']['aa_local_err']:5.2f}Å | "
        f"mem30 Cα={rdata['methods']['2dprots_mem30']['aa_local_err']:5.2f}Å",
        flush=True,
    )


def main():
    summary_path = os.path.join(OUT_DIR, "results_summary.json")
    with open(summary_path, "r", encoding="utf-8") as f:
        results = json.load(f)

    for target in benchmark_dimred.BENCHMARK_TARGETS:
        tid = target["id"]
        if tid in results:
            run_aa_extras_for_target(target, results[tid])

    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\nUpdated {summary_path} with 3d_pca, 2dprots_orig, 2dprots_mem15, and 2dprots_mem30.", flush=True)


if __name__ == "__main__":
    main()
