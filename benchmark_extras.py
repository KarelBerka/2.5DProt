#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""
Computes and renders the 4 reference/baseline views for all benchmark targets:
1. `3d_pca`: 3D Aligned Protein (Principal Axes PC1 x PC2, shortest axis PC3 into depth behind plane)
2. `2dprots_orig`: Original 2DProts (without membrane alignment, continuous angles)
3. `2dprots_mem15`: 2DProts + Membrane alignment (15° angle quantization & barrel unrolling)
4. `2dprots_mem30`: 2DProts + Membrane alignment (30° angle quantization & barrel unrolling)
"""

import os
import sys
import json
import time
import math
import warnings
import numpy as np
import matplotlib
import gemmi

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

warnings.filterwarnings("ignore")

import benchmark_dimred
import layout_loader
import utils_membrane
import utils_matrix2layout
import utils_sses_types
import superimpose2d

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE_DIR, "benchmark_dimred_results")
SVG_DIR = os.path.join(OUT_DIR, "svgs")
JSON_DIR = os.path.join(OUT_DIR, "layouts")

EXTRA_METHODS = [
    ("3d_pca", "3D Struktura (Principální osy)", "3D průmět: PC1×PC2 v ploše, nejkratší osa PC3 za plochou"),
    ("2dprots_orig", "2DProts Originál (bez membrány)", "Původní 2DProts (konsenzus / doménová evoluce, spojité úhly)"),
    ("2dprots_mem15", "2DProts + Membrána (15°)", "Zarovnání na membránu + kvantizace úhlů SSE po 15°"),
    ("2dprots_mem30", "2DProts + Membrána (30°)", "Zarovnání na membránu + kvantizace úhlů SSE po 30°"),
]

# True transmembrane benchmark targets (GPCRs, SLC/ABC TM domains & complexes, outer-membrane beta-barrels, KcsA, AQP1, Na/K-ATPase, Proton pumps, Respiratory chain)
TM_TARGET_IDS = {
    "1f88A00",
    "1f88_dimer",
    "2rh1A02",
    "2rh1_complex",
    "1okcA00",
    "1pv6A01",
    "1pv6A02",
    "1pv6_complex",
    "2hydA01",
    "2hyd_dimer",
    "1bxwA00",
    "2omfA00",
    "2omf_trimer",
    "1bl8A00",
    "1bl8_tetramer",
    "1j4nA00",
    "1j4n_tetramer",
    "2zxeA02",
    "2zxe_complex",
    "1c3wA00",
    "2xnd_complex",
    "1zoy_complex",
    "1occA00",
    "1occ_complex",
}


def _color_to_hex(color_code):
    """Convert 2DProts rainbow/spectrum string or CSS color name to hex string #rrggbb."""
    if not color_code or not isinstance(color_code, str):
        return "#64748b"
    c = color_code.strip()
    if c.startswith("s") and c[1:].isdigit():
        val = float(c[1:]) / 900.0
        cmap = matplotlib.colormaps["turbo"]
        rgb = cmap(np.clip(val, 0.06, 0.94))[:3]
        return "#{:02x}{:02x}{:02x}".format(int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255))
    named = {
        "red": "#dc322f",
        "blue": "#268bd2",
        "green": "#2aa198",
        "grey": "#878e96",
        "gray": "#878e96",
        "orange": "#e67319",
        "purple": "#8c4dc4",
        "cyan": "#19b2cc",
        "magenta": "#cc40a6",
        "yellow": "#d9a61f",
        "brown": "#8c5933",
        "black": "#282c34",
    }
    if c.lower() in named:
        return named[c.lower()]
    try:
        rgb = matplotlib.colors.to_rgb(c)
        return "#{:02x}{:02x}{:02x}".format(int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255))
    except Exception:
        return "#3b82f6"


def _blend_hex_with_white(hex_col, alpha_depth):
    """Blend hex color with white/slate background according to depth factor alpha_depth in [0.35, 1.0]."""
    h = hex_col.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    bg_r, bg_g, bg_b = 248, 250, 252
    nr = int(r * alpha_depth + bg_r * (1.0 - alpha_depth))
    ng = int(g * alpha_depth + bg_g * (1.0 - alpha_depth))
    nb = int(b * alpha_depth + bg_b * (1.0 - alpha_depth))
    return f"#{nr:02x}{ng:02x}{nb:02x}"


def extract_target_ca_traces_in_ann_frame(target, ann):
    """
    Extract 3D CA backbone traces for each domain/subunit in the exact 3D coordinate
    frame of `ann.vectors`.
    Returns a list of (color_hex, np.ndarray of shape (M, 3)).
    """
    traces = []
    try:
        pdb_id = target["pdb"]
        use_au = target.get("use_au", True) if target["kind"] == "complex" else True
        model = benchmark_dimred.fetch_pdb_structure(pdb_id, use_au=use_au)
        poly_chains = benchmark_dimred.get_polymer_chains(model)
        ndom = len(target["domains"])

        for d_idx, (ch_id, dom_id, dom_code) in enumerate(target["domains"]):
            if ch_id not in poly_chains:
                continue
            ca_map = benchmark_dimred.extract_ca_map(poly_chains[ch_id])
            if not ca_map:
                continue
            dom_data = benchmark_dimred.fetch_domain_json(dom_code)
            hr = dom_data["helices_residues"]
            all_res = [r for pair in hr.values() for r in pair]
            if not all_res:
                continue
            res_offset = 0
            min_hr = min(all_res)
            min_ca = min(ca_map.keys())
            if min_hr <= 10 and min_ca > 15 and (min_ca - 1 + max(all_res)) in ca_map:
                res_offset = min_ca - 1
            r_min, r_max = max(1, min_hr + res_offset - 4), max(all_res) + res_offset + 4
            seq_nums = sorted(r for r in ca_map if r_min <= r <= r_max)
            if len(seq_nums) < 5:
                continue
            ca_pts = np.array([ca_map[r] for r in seq_nums], dtype=float)

            if target["kind"] == "domain":
                # Superpose crystal CA trace onto the domain's consensus 3dlayout frame
                l3d = dom_data["3dlayout"]
                src_pts, dst_pts = [], []
                for sse, pair in hr.items():
                    if sse not in l3d:
                        continue
                    r1, r2 = min(pair) + res_offset, max(pair) + res_offset
                    if r1 in ca_map and r2 in ca_map:
                        src_pts.extend([ca_map[r1].tolist(), ca_map[r2].tolist()])
                        dst_pts.extend([l3d[sse][0], l3d[sse][1]])
                if len(src_pts) >= 4:
                    sup = gemmi.superpose_positions(
                        [gemmi.Position(*p) for p in dst_pts],
                        [gemmi.Position(*p) for p in src_pts],
                    )
                    transformed = []
                    for pt in ca_pts:
                        gp = sup.transform.apply(gemmi.Position(float(pt[0]), float(pt[1]), float(pt[2])))
                        transformed.append([gp.x, gp.y, gp.z])
                    ca_pts = np.array(transformed, dtype=float)
                col_hex = "#64748b"
            else:
                dom_color = f"s{int(820.0 * d_idx / max(1, ndom))}"
                col_hex = _color_to_hex(dom_color)

            traces.append((col_hex, ca_pts))
    except Exception:
        pass
    return traces


def compute_3d_principal_axes_alignment(ann, ca_traces, ref_layout):
    """
    Align the 3D protein coordinates along its 3D Principal Axes (SVD) so that:
      - Axis 0 (x, horizontal): Principal axis with largest spatial extent (PC1)
      - Axis 1 (y, vertical):   Principal axis with second-largest spatial extent (PC2)
      - Axis 2 (z, depth):      Principal axis with SHORTEST spatial extent (PC3 behind plane)
    Preserves proper 3D right-handed geometry (det(R) = +1) while choosing one of the 4
    in-plane sign combinations (±PC1, ±PC2, ±PC3) that best matches `ref_layout`.
    """
    pts_list = []
    for _, tr in ca_traces:
        if len(tr) > 0:
            pts_list.append(tr)
    for sse, vec in ann.vectors.items():
        v0 = np.asarray(vec[0], dtype=float)
        v1 = np.asarray(vec[1], dtype=float)
        for t in np.linspace(0.0, 1.0, 5):
            pts_list.append((v0 * (1.0 - t) + v1 * t)[None, :])

    P = np.vstack(pts_list)
    center_3d = np.mean(P, axis=0)
    P_c = P - center_3d

    # SVD of centered 3D coordinates gives orthonormal principal axes
    _, _, vt = np.linalg.svd(P_c, full_matrices=False)
    # Sort the 3 principal axes by actual peak-to-peak 3D span so the shortest span is strictly axis 2 (z)
    proj = P_c @ vt.T
    spans = np.ptp(proj, axis=0)
    order = np.argsort(spans)[::-1]  # [largest_span_idx, middle_span_idx, shortest_span_idx]
    R_base = vt[order, :]
    if np.linalg.det(R_base) < 0:
        R_base[2, :] = -R_base[2, :]

    # Test the 4 right-handed sign flips (s1, s2, s1*s2) to best match ref_layout orientation
    common = [s for s in ann.sses if s in ref_layout and s in ann.vectors]
    best_R = R_base
    if len(common) >= 3:
        ref_arr = np.array([[ref_layout[s][0], ref_layout[s][1]] for s in common], dtype=float)
        ref_arr = ref_arr - np.mean(ref_arr, axis=0, keepdims=True)
        n_r = np.linalg.norm(ref_arr)
        if n_r > 1e-4:
            ref_arr = ref_arr / n_r
        mids_3d = np.array(
            [0.5 * (np.asarray(ann.vectors[s][0]) + np.asarray(ann.vectors[s][1])) - center_3d for s in common],
            dtype=float,
        )
        best_score = -1e9
        for s1 in (1.0, -1.0):
            for s2 in (1.0, -1.0):
                s3 = s1 * s2  # keeps det(R) == +1
                S_mat = np.diag([s1, s2, s3])
                R_cand = S_mat @ R_base
                xy = (mids_3d @ R_cand.T)[:, :2]
                xy = xy - np.mean(xy, axis=0, keepdims=True)
                n_xy = np.linalg.norm(xy)
                if n_xy > 1e-4:
                    xy = xy / n_xy
                score = float(np.sum(xy * ref_arr))
                if score > best_score:
                    best_score = score
                    best_R = R_cand

    P_rot = P_c @ best_R.T
    spans_sorted = np.ptp(P_rot, axis=0)  # [span_x, span_y, span_z_shortest]
    return center_3d, best_R, spans_sorted


def render_3d_principal_axes_svg(tid, target, ann, stats, ca_traces, center_3d, R_pca, spans_3d, svg_path, json_path, ligs_3d=None):
    """
    Render the 3D Principal-Axes view (`svgs/{tid}_3d_pca.svg`) with depth-sorted CA trace,
    3D-shaded helices, beta-strands, and 3D ligands, and save the projected 2D layout JSON (`layouts/{tid}_3d_pca.json`).
    """
    def tr3d(pt):
        return (np.asarray(pt, dtype=float) - center_3d) @ R_pca.T

    # Transform all SSE vectors and build 2D layout + angles for quantitative comparison
    sse_3d_rot = {}
    layout_2d = {}
    angles_2d = {}
    h_type = {}
    for s in ann.sses:
        ht_val = 1 if stats.is_prime(s) else (2 if stats.is_sec(s, ann.res_sizes) else 3)
        h_type[s] = ht_val
        v0 = tr3d(ann.vectors[s][0])
        v1 = tr3d(ann.vectors[s][1])
        mid = 0.5 * (v0 + v1)
        sse_3d_rot[s] = (v0, v1, mid)
        layout_2d[s] = [float(mid[0]), float(mid[1])]
        dx = float(v1[0] - v0[0])
        dy = float(v1[1] - v0[1])
        if math.hypot(dx, dy) > 1e-4:
            angles_2d[s] = math.degrees(math.atan2(dx, dy)) % 360.0
        else:
            angles_2d[s] = 0.0

    ligs_rot_3d = {}
    ligs_2d_json = {}
    if ligs_3d:
        for lk, lv in ligs_3d.items():
            c_rot = tr3d(lv["3d_layout"])
            ligs_rot_3d[lk] = (c_rot, lv)
            entry = dict(lv)
            entry["layout"] = [float(c_rot[0]), float(c_rot[1])]
            ligs_2d_json[lk] = entry

    # Save layout JSON for metric evaluation and matplotlib panel rendering
    jdata = {
        "layout": layout_2d,
        "angles": angles_2d,
        "helices_residues": ann.helices_residues,
        "size": ann.res_sizes,
        "protein": tid,
        "chain_id": ann.chain_id,
        "chain": "ALL" if target["kind"] == "complex" else target["domains"][0][0],
        "domain": "ALL" if target["kind"] == "complex" else target["domains"][0][1],
        "h_type": h_type,
        "beta_connectivity": ann.get_beta_connectivity_clusters(),
        "color": ann.color,
        "3dlayout": ann.layout_3d(),
        "membrane": {"is_transmembrane": False, "angle_step": 0},
        "pca_spans_3d": [round(float(x), 1) for x in spans_3d],
    }
    if ligs_2d_json:
        jdata["ligands"] = ligs_2d_json
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(jdata, f)

    # Collect primitives for depth-sorted SVG rendering (in SVG, +y is down, so y_svg = -y_3d)
    all_pts_xy = []
    all_z = []
    for s, (v0, v1, mid) in sse_3d_rot.items():
        if h_type.get(s, 3) in (1, 2):
            all_pts_xy.extend([[v0[0], v0[1]], [v1[0], v1[1]]])
            all_z.extend([v0[2], v1[2]])
    for lk, (c_rot, _) in ligs_rot_3d.items():
        all_pts_xy.append([float(c_rot[0]), float(c_rot[1])])
        all_z.append(float(c_rot[2]))
    rot_traces = []
    for col_hex, tr in ca_traces:
        tr_r = (tr - center_3d) @ R_pca.T
        rot_traces.append((col_hex, tr_r))
        all_pts_xy.extend(tr_r[:, :2].tolist())
        all_z.extend(tr_r[:, 2].tolist())

    if not all_pts_xy:
        all_pts_xy = [[-10, -10], [10, 10]]
        all_z = [-5.0, 5.0]

    xy_arr = np.array(all_pts_xy, dtype=float)
    z_min, z_max = float(min(all_z)), float(max(all_z))
    z_span = max(1.0, z_max - z_min)

    min_x, max_x = float(np.min(xy_arr[:, 0])), float(np.max(xy_arr[:, 0]))
    min_y, max_y = float(np.min(xy_arr[:, 1])), float(np.max(xy_arr[:, 1]))
    pad = max(6.0, 0.12 * max(max_x - min_x, max_y - min_y))
    vb_w = max(20.0, (max_x - min_x) + 2 * pad)
    vb_h = max(20.0, (max_y - min_y) + 2 * pad + 5.0)
    scale = 14.0  # pixels per Angstrom

    def to_svg_xy(x3, y3):
        sx = (x3 - (min_x - pad)) * scale
        sy = ((max_y + pad) - y3) * scale
        return sx, sy

    primitives = []

    # 1. Backbone CA trace segments (broken into chunks of 3 residues for smooth depth-sorting)
    for col_hex, tr_r in rot_traces:
        m = len(tr_r)
        step_c = 2 if m < 400 else 4
        for i in range(0, m - 1, step_c):
            chunk = tr_r[i : min(m, i + step_c + 1)]
            z_c = float(np.mean(chunk[:, 2]))
            depth_norm = (z_c - z_min) / z_span  # 0 = back, 1 = front
            alpha_d = 0.35 + 0.55 * depth_norm
            stroke_c = _blend_hex_with_white(col_hex, alpha_d * 0.75)
            lw = (1.3 + 1.2 * depth_norm)
            pts_str = " ".join(f"{to_svg_xy(p[0], p[1])[0]:.1f},{to_svg_xy(p[0], p[1])[1]:.1f}" for p in chunk)
            svg_el = (
                f'<polyline points="{pts_str}" fill="none" stroke="{stroke_c}" '
                f'stroke-width="{lw:.2f}" stroke-linecap="round" stroke-linejoin="round" opacity="{0.55 + 0.35 * depth_norm:.2f}"/>'
            )
            primitives.append((z_c - 0.5, svg_el))

    # 2. 3D SSE Cylinders (Helices) and 3D Arrows (Beta-strands)
    shown_sses = [s for s in ann.sses if h_type.get(s, 3) in (1, 2)]
    show_labels = len(shown_sses) <= 45
    for s in shown_sses:
        v0, v1, mid = sse_3d_rot[s]
        z_c = float(mid[2])
        depth_norm = (z_c - z_min) / z_span  # 0 = back, 1 = front
        base_hex = _color_to_hex(ann.color.get(s, "blue"))
        col_shaded = _blend_hex_with_white(base_hex, 0.50 + 0.50 * depth_norm)
        x0, y0 = to_svg_xy(v0[0], v0[1])
        x1, y1 = to_svg_xy(v1[0], v1[1])
        mx, my = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
        seg_len_px = math.hypot(x1 - x0, y1 - y0)

        is_helix = utils_sses_types.is_helic(s)
        if is_helix:
            w_cyl = (7.5 + 3.5 * depth_norm)
            if seg_len_px < 4.0:
                # Nearly end-on 3D helix -> draw 3D circular cylinder cross-section
                r_c = w_cyl * 0.68
                el = (
                    f'<g>'
                    f'<circle cx="{mx:.1f}" cy="{my:.1f}" r="{r_c + 1.2:.1f}" fill="#0f172a" opacity="{0.25 + 0.25 * depth_norm:.2f}"/>'
                    f'<circle cx="{mx:.1f}" cy="{my:.1f}" r="{r_c:.1f}" fill="{col_shaded}" stroke="#1e293b" stroke-width="1.2"/>'
                    f'<circle cx="{mx - 0.25 * r_c:.1f}" cy="{my - 0.25 * r_c:.1f}" r="{0.35 * r_c:.1f}" fill="#ffffff" opacity="0.45"/>'
                    f'</g>'
                )
            else:
                el = (
                    f'<g>'
                    f'<line x1="{x0:.1f}" y1="{y0:.1f}" x2="{x1:.1f}" y2="{y1:.1f}" '
                    f'stroke="#0f172a" stroke-width="{w_cyl + 2.4:.1f}" stroke-linecap="round" opacity="{0.22 + 0.25 * depth_norm:.2f}"/>'
                    f'<line x1="{x0:.1f}" y1="{y0:.1f}" x2="{x1:.1f}" y2="{y1:.1f}" '
                    f'stroke="{col_shaded}" stroke-width="{w_cyl:.1f}" stroke-linecap="round"/>'
                    f'<line x1="{x0:.1f}" y1="{y0:.1f}" x2="{x1:.1f}" y2="{y1:.1f}" '
                    f'stroke="#ffffff" stroke-width="{w_cyl * 0.28:.1f}" stroke-linecap="round" opacity="0.42"/>'
                    f'</g>'
                )
        else:
            # 3D Beta-strand arrow from v0 to v1
            w_str = (4.8 + 2.2 * depth_norm)
            if seg_len_px < 5.0:
                dx_u, dy_u = 0.0, -1.0
                x0, y0 = mx, my + 3.0
                x1, y1 = mx, my - 3.0
                seg_len_px = 6.0
            else:
                dx_u = (x1 - x0) / seg_len_px
                dy_u = (y1 - y0) / seg_len_px
            nx_u, ny_u = -dy_u, dx_u
            head_len = min(10.5, max(5.5, 0.35 * seg_len_px))
            head_w = w_str * 1.15
            half_w = w_str * 0.52
            bx = x1 - head_len * dx_u
            by = y1 - head_len * dy_u
            pts_poly = [
                (x0 + nx_u * half_w, y0 + ny_u * half_w),
                (bx + nx_u * half_w, by + ny_u * half_w),
                (bx + nx_u * head_w, by + ny_u * head_w),
                (x1, y1),
                (bx - nx_u * head_w, by - ny_u * head_w),
                (bx - nx_u * half_w, by - ny_u * half_w),
                (x0 - nx_u * half_w, y0 - ny_u * half_w),
            ]
            p_str = " ".join(f"{px:.1f},{py:.1f}" for px, py in pts_poly)
            el = (
                f'<polygon points="{p_str}" fill="{col_shaded}" '
                f'stroke="#0f172a" stroke-width="1.1" stroke-linejoin="round"/>'
            )

        if show_labels:
            el += (
                f'<text x="{mx:.1f}" y="{my + 3.2:.1f}" font-family="Verdana, sans-serif" '
                f'font-size="8.5" font-weight="bold" fill="#0f172a" stroke="#ffffff" '
                f'stroke-width="2.2" paint-order="stroke" text-anchor="middle">{s}</text>'
            )
        primitives.append((z_c, el))

    # 3. 3D Ligands & Cofactors (rendered as shaded 3D spheres with bold chemical codes)
    for lk, (c_rot, lv) in ligs_rot_3d.items():
        z_l = float(c_rot[2])
        depth_norm = (z_l - z_min) / z_span
        lx, ly = to_svg_xy(c_rot[0], c_rot[1])
        natoms = int(lv.get("atoms", 1))
        r_px = 5.2 if natoms < 10 else min(11.5, 4.2 + math.sqrt(natoms) * 0.95)
        f_col = lv.get("color", "#ef4444")
        s_col = lv.get("ocolor", "#7f1d1d")
        ecmp = lv.get("ecmp", lk)
        lig_el = (
            f'<g>'
            f'<circle cx="{lx:.1f}" cy="{ly:.1f}" r="{r_px + 2.0:.1f}" fill="#ffffff" opacity="0.88"/>'
            f'<circle cx="{lx:.1f}" cy="{ly:.1f}" r="{r_px:.1f}" fill="{f_col}" stroke="{s_col}" stroke-width="1.8"/>'
            f'<circle cx="{lx - 0.28 * r_px:.1f}" cy="{ly - 0.28 * r_px:.1f}" r="{0.34 * r_px:.1f}" fill="#ffffff" opacity="0.55"/>'
            f'<text x="{lx:.1f}" y="{ly - r_px - 3.0:.1f}" font-family="Verdana, sans-serif" '
            f'font-size="9.5" font-weight="bold" fill="{s_col}" stroke="#ffffff" '
            f'stroke-width="2.6" paint-order="stroke" text-anchor="middle">{ecmp}</text>'
            f'</g>'
        )
        primitives.append((z_l + 3.0, lig_el))

    # Sort all 3D primitives from back (z_min) to front (z_max)
    primitives.sort(key=lambda item: item[0])
    body_svg = "\n".join(item[1] for item in primitives)

    w_px = int(round(vb_w * scale))
    h_px = int(round(vb_h * scale))
    sx_ang, sy_ang, sz_ang = spans_3d[0], spans_3d[1], spans_3d[2]
    badge_text = f"3D Principální osy: PC1={sx_ang:.0f} Å × PC2={sy_ang:.0f} Å | hloubka PC3={sz_ang:.0f} Å"

    svg_doc = f"""<?xml version="1.0" encoding="utf-8" ?>
<svg xmlns="http://www.w3.org/2000/svg" width="{w_px}" height="{h_px}" viewBox="0 0 {w_px} {h_px}">
  <rect width="100%" height="100%" fill="#ffffff"/>
  {body_svg}
  <rect x="8" y="{h_px - 26}" width="{min(w_px - 16, 390)}" height="20" rx="5" fill="#f1f5f9" stroke="#cbd5e1" stroke-width="1"/>
  <text x="16" y="{h_px - 12}" font-family="Verdana, sans-serif" font-size="10.5" font-weight="bold" fill="#334155">{badge_text}</text>
</svg>
"""
    with open(svg_path, "w", encoding="utf-8") as f:
        f.write(svg_doc)


def compute_original_2dprots_layout(target, ann, stats, ref_layout):
    """
    Compute the authentic Original 2DProts layout (without membrane alignment):
      - For a single domain: uses the exact precomputed 2DProts domain layout & angles from 2dprots.ncbr.muni.cz.
      - For a multi-domain / multi-chain complex: takes each domain's precomputed 2DProts domain layout
        and assembles the complex using 2DProts' hierarchical `modify_protein_layout_via_evolution`.
    """
    if target["kind"] == "domain":
        _, _, dom_code = target["domains"][0]
        dom_data = benchmark_dimred.fetch_domain_json(dom_code)
        lay = {s: (float(v[0]), float(v[1])) for s, v in dom_data["layout"].items() if s in ann.sses}
        ang = {s: float(v) % 360.0 for s, v in dom_data["angles"].items() if s in ann.sses}
        return lay, ang

    # Complex: hierarchical domain assembly (`compute_whole_protein` logic)
    ndom = len(target["domains"])
    dom_keys = []
    domains_layout = {}
    domains_angles = {}
    sses_sizes_dom = {}
    layout_3d_dom = {}

    for d_idx, (ch_id, dom_id, dom_code) in enumerate(target["domains"]):
        dkey = f"D{d_idx + 1:02d}"
        dom_keys.append(dkey)
        dom_data = benchmark_dimred.fetch_domain_json(dom_code)
        prefix = f"{d_idx + 1:02d}" if ndom >= 10 else f"{d_idx + 1}"

        def rename_sse(s):
            if s.startswith("_"):
                return f"{s[:2]}{prefix}_{s[2:]}"
            return f"{s[0]}{prefix}_{s[1:]}"

        d_lay = {}
        d_ang = {}
        d_sz = {}
        d_3d = {}
        for s_orig, coord in dom_data["layout"].items():
            s_new = rename_sse(s_orig)
            if s_new in ann.sses:
                d_lay[s_new] = [float(coord[0]), float(coord[1])]
                d_ang[s_new] = float(dom_data["angles"].get(s_orig, 0.0))
                d_sz[s_new] = ann.res_sizes[s_new]
                d_3d[s_new] = ann.vectors[s_new]

        # Center each domain's 2D layout around (0, 0) before domain-level evolution
        if d_lay:
            cx = sum(v[0] for v in d_lay.values()) / len(d_lay)
            cy = sum(v[1] for v in d_lay.values()) / len(d_lay)
            for k in d_lay:
                d_lay[k] = (d_lay[k][0] - cx, d_lay[k][1] - cy)

        domains_layout[dkey] = d_lay
        domains_angles[dkey] = d_ang
        sses_sizes_dom[dkey] = d_sz
        layout_3d_dom[dkey] = d_3d

    # Build inter-domain 3D SSE-to-SSE distance matrix as expected by compute_protein_dimred_layout / count_protein_error
    dmat_all = ann.get_distances()
    s_dist_w = {"ALL": {dk1: {"ALL": {}} for dk1 in dom_keys}}
    for dk1 in dom_keys:
        for dk2 in dom_keys:
            if dk1 == dk2:
                continue
            s_dist_w["ALL"][dk1]["ALL"][dk2] = {
                s1: {s2: dmat_all[s1][s2] for s2 in layout_3d_dom[dk2]}
                for s1 in layout_3d_dom[dk1]
            }

    p_lay = {dk: (0.0, 0.0) for dk in dom_keys}
    p_ang = {dk: 0.0 for dk in dom_keys}

    p_lay_w = {"ALL": p_lay}
    p_ang_w = {"ALL": p_ang}
    d_lay_w = {"ALL": domains_layout}
    d_ang_w = {"ALL": domains_angles}
    sz_w = {"ALL": sses_sizes_dom}

    p_lay_w, p_ang_w, d_lay_w, d_ang_w = utils_matrix2layout.modify_protein_layout_via_evolution(
        p_lay_w, p_ang_w, d_lay_w, d_ang_w, s_dist_w, sz_w, method="sammon"
    )
    d_lay_w, d_ang_w = utils_matrix2layout.count_sse_protein_coords_based_on_cluster_coords(
        p_lay_w, p_ang_w, d_lay_w, d_ang_w, 1
    )

    out_lay = {}
    out_ang = {}
    for dk in dom_keys:
        for s, coord in d_lay_w["ALL"][dk].items():
            out_lay[s] = (float(coord[0]), float(coord[1]))
            out_ang[s] = float(d_ang_w["ALL"][dk][s]) % 360.0

    # Align to ref_layout via rigid 2D Procrustes
    common = [k for k in out_lay if k in ref_layout]
    if len(common) >= 2:
        A = np.array([[out_lay[k][0], out_lay[k][1]] for k in common], dtype=float).T
        B = np.array([[ref_layout[k][0], ref_layout[k][1]] for k in common], dtype=float).T
        R, t = superimpose2d.optimal_rotation_translation(A, B, allow_mirror=True)
        for k in list(out_lay.keys()):
            pt = R @ np.array([[out_lay[k][0]], [out_lay[k][1]]], dtype=float) + t
            out_lay[k] = (float(pt[0, 0]), float(pt[1, 0]))
        for k in list(out_ang.keys()):
            a_rad = math.radians(out_ang[k])
            v_rot = R @ np.array([math.sin(a_rad), math.cos(a_rad)], dtype=float)
            out_ang[k] = math.degrees(math.atan2(float(v_rot[0]), float(v_rot[1]))) % 360.0

    return out_lay, out_ang


def run_extras_for_target(target, rdata):
    """Compute and render `3d_pca`, `2dprots_orig`, `2dprots_mem15`, and `2dprots_mem30` for `target`."""
    import utils_benchmark_ligands
    tid = target["id"]
    ann, stats, ref_layout, ref_angles, chain_id, dom_id = benchmark_dimred.build_target_annotation(target)
    h_type = {s: (1 if stats.is_prime(s) else (2 if stats.is_sec(s, ann.res_sizes) else 3)) for s in ann.sses}
    all_c = ann.get_beta_connectivity_clusters()
    ligs_3d = benchmark_dimred.get_cached_target_ligands_3d(target, ann)

    # Ensure all 12 DR layout JSONs for this target have projected 2D ligands and updated SVGs
    for m_id, _, _ in benchmark_dimred.METHODS:
        m_json = os.path.join(JSON_DIR, f"{tid}_{m_id}.json")
        m_svg = os.path.join(SVG_DIR, f"{tid}_{m_id}.svg")
        if os.path.exists(m_json):
            with open(m_json, "r", encoding="utf-8") as f:
                mj = json.load(f)
            if ligs_3d and not mj.get("ligands"):
                mj_lay = {k: tuple(v) for k, v in mj.get("layout", {}).items()}
                mj_ang = {k: float(v) for k, v in mj.get("angles", {}).items()}
                mj["ligands"] = utils_benchmark_ligands.project_ligands_to_2d(ligs_3d, ann, mj_lay, mj_ang, mem_info=None)
                mj["membrane"] = {"is_transmembrane": False, "angle_step": 0}
                for _r in range(6):
                    try:
                        with open(m_json, "w", encoding="utf-8") as f:
                            json.dump(mj, f)
                        break
                    except OSError:
                        time.sleep(0.15)
                benchmark_dimred.render_svg_from_json(m_json, m_svg)

    # 1. 3D Principal Axes view (`3d_pca`)
    t0 = time.perf_counter()
    ca_traces = extract_target_ca_traces_in_ann_frame(target, ann)
    c3d, R_pca, spans_3d = compute_3d_principal_axes_alignment(ann, ca_traces, ref_layout)
    svg_3d = os.path.join(SVG_DIR, f"{tid}_3d_pca.svg")
    json_3d = os.path.join(JSON_DIR, f"{tid}_3d_pca.json")
    render_3d_principal_axes_svg(tid, target, ann, stats, ca_traces, c3d, R_pca, spans_3d, svg_3d, json_3d, ligs_3d=ligs_3d)
    ms_3d = (time.perf_counter() - t0) * 1000.0

    with open(json_3d, "r", encoding="utf-8") as f:
        j3d = json.load(f)
    dummy_l = layout_loader.Layout(sses_stats=stats, annotation=ann)
    dummy_l.layout = {k: tuple(v) for k, v in j3d["layout"].items()}
    dummy_l.a_layout = dict(j3d["angles"])
    dummy_l.h_type = h_type
    m_3d = benchmark_dimred.evaluate_layout_metrics(ann, dummy_l)
    m_3d["time_ms"] = round(ms_3d, 1)
    m_3d["pca_spans"] = [round(float(x), 1) for x in spans_3d]
    m_3d["n_ligands"] = len(ligs_3d)
    rdata["methods"]["3d_pca"] = m_3d

    # 2. Original 2DProts without membrane (`2dprots_orig`)
    t0 = time.perf_counter()
    orig_lay, orig_ang = compute_original_2dprots_layout(target, ann, stats, ref_layout)
    ms_orig = (time.perf_counter() - t0) * 1000.0

    l_orig = layout_loader.Layout(sses_stats=stats, annotation=ann)
    l_orig.layout = dict(orig_lay)
    l_orig.a_layout = dict(orig_ang)
    l_orig.h_type = h_type
    l_orig.membrane = {"is_transmembrane": False, "angle_step": 0}
    ligs_orig = utils_benchmark_ligands.project_ligands_to_2d(ligs_3d, ann, orig_lay, orig_ang, mem_info=None)
    l_orig.ligands_data = ligs_orig
    json_orig = os.path.join(JSON_DIR, f"{tid}_2dprots_orig.json")
    svg_orig = os.path.join(SVG_DIR, f"{tid}_2dprots_orig.svg")
    l_orig.save(chain_id, dom_id, json_orig)
    benchmark_dimred.render_svg_from_json(json_orig, svg_orig)
    m_orig = benchmark_dimred.evaluate_layout_metrics(ann, l_orig)
    m_orig["time_ms"] = round(ms_orig, 1)
    m_orig["n_ligands"] = len(ligs_orig)
    rdata["methods"]["2dprots_orig"] = m_orig

    # 3 & 4. 2DProts with Membrane alignment (15° and 30° variants)
    is_tm_target = tid in TM_TARGET_IDS
    for step_deg, m_key, b_mode in [
        (15.0, "2dprots_mem15", "unrolled_15"),
        (30.0, "2dprots_mem30", "unrolled_30"),
    ]:
        t0 = time.perf_counter()
        base_lay = {k: list(v) for k, v in orig_lay.items()}
        base_ang = dict(orig_ang)
        mem_lay, mem_ang, mem_info = utils_membrane.align_layout_to_membrane(
            base_lay,
            base_ang,
            ann.sizes,
            all_c,
            ann.vectors,
            helices_residues=ann.helices_residues,
            h_type=h_type,
            angle_step=step_deg,
            barrel_mode=b_mode,
            force_tm=is_tm_target,
        )
        ms_mem = (time.perf_counter() - t0) * 1000.0
        l_mem = layout_loader.Layout(sses_stats=stats, annotation=ann)
        l_mem.layout = {k: (float(v[0]), float(v[1])) for k, v in mem_lay.items()}
        l_mem.a_layout = {k: float(v) % 360.0 for k, v in mem_ang.items()}
        l_mem.h_type = h_type
        if mem_info:
            mem_info["angle_step"] = step_deg
            mem_info["is_transmembrane"] = True
            l_mem.membrane = mem_info
        else:
            l_mem.membrane = {"is_transmembrane": False, "angle_step": step_deg}

        ligs_mem = utils_benchmark_ligands.project_ligands_to_2d(
            ligs_3d, ann, l_mem.layout, l_mem.a_layout, mem_info=l_mem.membrane
        )
        l_mem.ligands_data = ligs_mem

        json_mem = os.path.join(JSON_DIR, f"{tid}_{m_key}.json")
        svg_mem = os.path.join(SVG_DIR, f"{tid}_{m_key}.svg")
        l_mem.save(chain_id, dom_id, json_mem)
        benchmark_dimred.render_svg_from_json(json_mem, svg_mem)
        m_mem = benchmark_dimred.evaluate_layout_metrics(ann, l_mem)
        m_mem["time_ms"] = round(ms_mem, 1)
        m_mem["is_tm"] = bool(mem_info)
        m_mem["n_ligands"] = len(ligs_mem)
        rdata["methods"][m_key] = m_mem

    lig_codes = sorted({v["ecmp"] for v in ligs_3d.values()}) if ligs_3d else []
    print(
        f"[{tid:14s}] 3d_pca={m_3d['local_err']:5.2f}Å ({spans_3d[0]:.0f}×{spans_3d[1]:.0f}×{spans_3d[2]:.0f}Å) | "
        f"orig={m_orig['local_err']:5.2f}Å | mem15={rdata['methods']['2dprots_mem15']['local_err']:5.2f}Å | "
        f"mem30={rdata['methods']['2dprots_mem30']['local_err']:5.2f}Å (TM={rdata['methods']['2dprots_mem30']['is_tm']}, ligs={len(ligs_3d)}:{','.join(lig_codes)})",
        flush=True,
    )


def main():
    summary_path = os.path.join(OUT_DIR, "results_summary.json")
    with open(summary_path, "r", encoding="utf-8") as f:
        results = json.load(f)

    only_targets = sys.argv[1:] if len(sys.argv) > 1 else None
    for target in benchmark_dimred.BENCHMARK_TARGETS:
        tid = target["id"]
        if only_targets and tid not in only_targets:
            continue
        if tid not in results:
            continue
        run_extras_for_target(target, results[tid])

    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\nUpdated {summary_path} with 3d_pca, 2dprots_orig, 2dprots_mem15, and 2dprots_mem30.", flush=True)


if __name__ == "__main__":
    main()
