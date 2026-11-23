#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Benchmark & Demonstration of 2.5DProts (`2.5DPROTS-DSL`) vs Classic 2DProts (`2DPROTS-DSL`)
for 3D Protein Structure Generation (`Contact Map / Distogram -> 2.5DProts -> 3D Coordinates`).

Evaluates across all 34 benchmark domains:
  1. Classic 2DProts (`2DPROTS-DSL -> 3D`): flat 2D (X, Y, alpha) + learned distance lifting.
  2. 2.5DProts (`Distogram -> 2.5DPROTS-DSL -> 3D`, ZERO 3D crystal coordinates used):
     2D (X, Y, alpha) augmented with 2.5D out-of-plane depth Z, 3D spherical orientation (az, el),
     and 3D SSE span S reconstructed via 3D endpoint SMACOF on the residue Distogram (<= 22 A).
  3. De novo & custom-engineered 2.5DProts protein architectures exported to `.cif` and `.pdb`.
"""

import json
import math
import os
import shutil
import sys
import time
from collections import OrderedDict

# Skip UMAP/Numba import overhead
sys.modules["umap"] = None

import numpy as np
from scipy.sparse.csgraph import shortest_path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import benchmark_dimred as bd
from utils_contactmap_llm import (
    BidirectionalDSLTranslator,
    build_residue_contact_map,
    compute_2dprots_layout_from_contact_geometry,
    compute_kabsch_rmsd_and_tmscore,
    decode_25dprots_dsl,
    decode_2dprots_dsl,
    encode_25dprots_to_dsl,
    encode_2dprots_to_dsl,
    encode_contact_map_to_dsl,
    export_3d_structure_to_mmcif_and_pdb,
    parse_cm_dsl,
    reconstruct_3d_endpoints_from_distogram,
    reconstruct_3d_structure_from_25dprots_dsl,
    reconstruct_3d_structure_from_2dprots_dsl,
)

WORKSPACE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(WORKSPACE_DIR, "benchmark_dimred_results")
JSON_DIR = os.path.join(OUT_DIR, "layouts")
MODELS_DIR = os.path.join(OUT_DIR, "models_3d")
FIG_DIR = os.path.join(OUT_DIR, "figures")
ARTIFACT_DIR = r"C:\Users\krapn\.gemini\antigravity\brain\680ab902-44a0-4112-abe0-0acf18e03e35"


def collect_benchmark_domains():
    """Load all 34 domain targets from the 2DProts benchmark along with their Distograms and 2D/2.5D DSLs."""
    records = []
    for target in bd.BENCHMARK_TARGETS:
        if target.get("kind") != "domain":
            continue
        tid = target["id"]
        pdb_id = target["pdb"]
        try:
            ann, stats, ref_layout, ref_angles, chain_id, dom_id = bd.build_target_annotation(target)
            d3d_dict = ann.get_distances()
            ang3d_dict = ann.get_all_helices_pairs_angles()
            sse_order, sse_slices, D_res, C_bin8, C_bin10, P_cont, D_distogram = build_residue_contact_map(ann, stats)
            n_s = len(sse_order)
            if n_s < 3:
                continue

            D_sse_true = np.zeros((n_s, n_s), dtype=float)
            ang_mat_true = np.full((n_s, n_s), 90.0, dtype=float)
            for i, si in enumerate(sse_order):
                for j, sj in enumerate(sse_order):
                    D_sse_true[i, j] = float(d3d_dict[si][sj])
                    if i != j:
                        ang_mat_true[i, j] = float(ang3d_dict.get((si, sj), 90.0))

            cm_dsl_str, _ = encode_contact_map_to_dsl(sse_order, sse_slices, D_sse_true, ang_mat_true, ann)
            cm_nodes, cm_sheets, cm_contacts, cm_loops = parse_cm_dsl(cm_dsl_str)

            sammon_json_path = os.path.join(JSON_DIR, f"{tid}_sammon.json")
            if os.path.exists(sammon_json_path):
                with open(sammon_json_path, "r", encoding="utf-8") as f:
                    s_data = json.load(f)
                l_ref = bd.layout_loader.Layout(sses_stats=stats, annotation=ann, method="sammon")
                l_ref.layout = {k: tuple(v) for k, v in s_data.get("layout", {}).items()}
                l_ref.a_layout = {k: float(v) for k, v in s_data.get("angles", {}).items()}
                l_ref.h_type = {k: int(v) for k, v in s_data.get("h_type", {}).items()}
            else:
                l_ref = compute_2dprots_layout_from_contact_geometry(
                    ann, stats, d3d_dict, ang3d_dict, base_method="sammon"
                )

            dsl_2d_true, _ = encode_2dprots_to_dsl(sse_order, l_ref, ann)
            layout_2d, angles_2d, lengths, clusters = decode_2dprots_dsl(dsl_2d_true)

            # Reconstruct 2.5DProts purely from the Distogram (<= 22 A, ZERO 3D crystal coordinates!)
            t0_25d = time.perf_counter()
            sse_3d_endpoints_disto, X3_disto, _ = reconstruct_3d_endpoints_from_distogram(
                sse_order,
                sse_slices,
                D_distogram,
                l_obj_2d=l_ref,
                chain_ids=ann.chain_id,
                res_ranges=ann.helices_residues,
                n_smacof_iter=90,
            )
            dsl_25d_disto, _ = encode_25dprots_to_dsl(sse_order, l_ref, sse_3d_endpoints_disto, ann)
            t_25d_ms = (time.perf_counter() - t0_25d) * 1000.0

            records.append({
                "id": tid,
                "pdb": pdb_id,
                "name": target.get("name", tid),
                "ann": ann,
                "stats": stats,
                "l_ref": l_ref,
                "sse_order": sse_order,
                "sse_slices": sse_slices,
                "D_sse_true": D_sse_true,
                "D_distogram": D_distogram,
                "cm_nodes": cm_nodes,
                "cm_sheets": cm_sheets,
                "cm_contacts": cm_contacts,
                "cm_loops": cm_loops,
                "dsl_2d": dsl_2d_true,
                "dsl_25d": dsl_25d_disto,
                "t_25d_encode_ms": t_25d_ms,
                "layout_2d": layout_2d,
                "angles_2d": angles_2d,
                "lengths": lengths,
                "clusters": clusters,
            })
        except Exception as exc:
            print(f"[WARN] Skipping {tid}: {exc}")
    return records


def evaluate_2d_baseline_mds_rmsd(rec, train_sub):
    """Fast 3D MDS evaluation of classic 2DProts (`2DPROTS-DSL -> Contact Map -> 3D MDS`)."""
    t0 = time.perf_counter()
    tr = BidirectionalDSLTranslator()
    tr.fit(train_sub)
    _, D_sse_pred, _, _, _, _ = tr.translate_2dprots_dsl_to_contact_map(
        rec["dsl_2d"], rec["sse_order"], rec["sse_slices"]
    )
    n_s = len(rec["sse_order"])
    D_comb = D_sse_pred.copy()
    for i in range(n_s - 1):
        s1, s2 = rec["sse_order"][i], rec["sse_order"][i + 1]
        gap = rec["cm_loops"].get((s1, s2), 4)
        d_loop_est = min(14.0, 5.0 + 1.8 * math.sqrt(float(gap)))
        D_comb[i, i + 1] = D_comb[i + 1, i] = 0.55 * D_comb[i, i + 1] + 0.45 * d_loop_est

    H = np.eye(n_s) - np.ones((n_s, n_s)) / float(n_s)
    B = -0.5 * H @ (D_comb ** 2) @ H
    evals, evecs = np.linalg.eigh(B)
    idx = np.argsort(evals)[::-1][:3]
    X3_centers = evecs[:, idx] * np.sqrt(np.maximum(0.0, evals[idx]))
    el_ms = (time.perf_counter() - t0) * 1000.0

    true_centers = np.array([
        0.5 * (np.asarray(rec["ann"].vectors[s][0], float) + np.asarray(rec["ann"].vectors[s][1], float))
        for s in rec["sse_order"]
    ])
    n_res_tot = sum(rec["lengths"].values())
    rmsd_c, tm_c, _ = compute_kabsch_rmsd_and_tmscore(X3_centers, true_centers, n_res_total=n_res_tot)
    return rmsd_c, tm_c, el_ms


def build_denovo_25dprots_designs(rec_by_id):
    """
    Construct 3 custom / de novo `2.5DPROTS-DSL` specifications showing how explicit
    Z-depth and 3D elevation tilt `EL` program exact 3D multi-layer and barrel topologies.
    """
    designs = OrderedDict()

    # Design 1: 10-stranded circular beta-barrel + central coaxial alpha-helix H1
    tokens_1 = ["<BOS_25DP>", "CLUST:E1:E2:E3:E4:E5:E6:E7:E8:E9:E10"]
    radius_b = 8.4
    for k in range(10):
        theta = 2.0 * math.pi * k / 10.0
        x = radius_b * math.cos(theta)
        y = radius_b * math.sin(theta)
        z = 0.8 * math.sin(2.0 * theta)
        tang_az = int(round((math.degrees(theta) + 90.0) % 360.0))
        el = +68 if (k % 2 == 0) else -68
        tokens_1.append(
            f"SSE25D:E{k+1}:L11:X{x:+06.1f}:Y{y:+06.1f}:Z{z:+06.1f}:A{tang_az:03d}:EL{el:+03d}"
        )
    tokens_1.append("SSE25D:H1:L17:X+000.0:Y+000.0:Z+000.0:A000:EL+88")
    tokens_1.append("<EOS_25DP>")
    designs["denovo_25d_gfp_barrel_10B1H"] = {
        "title": "2.5DProts De novo 10-Stranded Beta-Barrel + Coaxial Central Alpha-Helix (GFP-like Fold)",
        "description": "10 tilted antiparallel beta-strands (E1..E10, EL=±68°, R=8.4 Å) surrounding a coaxial central alpha-helix (H1, Z=0, EL=+88°).",
        "dsl_25d": " ".join(tokens_1),
    }

    # Design 2: 3-Layer Alpha/Beta/Alpha Rossmann-like Sandwich using explicit Z-depth layers!
    # Central 4-stranded parallel beta-sheet at layer Z = 0.0 A,
    # flanked by 2 alpha-helices above (layer Z = +7.5 A) and 2 alpha-helices below (layer Z = -7.5 A)!
    tokens_2 = [
        "<BOS_25DP>",
        "CLUST:E2:E1:E3:E4",
        "SSE25D:E1:L8:X-002.3:Y+000.0:Z+000.0:A005:EL+82",
        "SSE25D:H1:L13:X-004.5:Y+002.5:Z+007.6:A185:EL-78",
        "SSE25D:E2:L8:X-006.9:Y+000.0:Z-000.4:A005:EL+82",
        "SSE25D:H2:L13:X-002.0:Y-002.5:Z-007.6:A185:EL-78",
        "SSE25D:E3:L8:X+002.3:Y+000.0:Z+000.0:A005:EL+82",
        "SSE25D:H3:L13:X+004.2:Y+002.5:Z+007.6:A185:EL-78",
        "SSE25D:E4:L8:X+006.9:Y+000.0:Z-000.4:A005:EL+82",
        "SSE25D:H4:L14:X+005.0:Y-002.5:Z-007.6:A185:EL-78",
        "<EOS_25DP>",
    ]
    designs["denovo_25d_rossmann_3layer_4E4H"] = {
        "title": "2.5DProts De novo 3-Layer α/β/α Rossmann Sandwich (Z = -7.6 Å / 0.0 Å / +7.6 Å)",
        "description": "Central 4-stranded beta-sheet at Z=0 Å flanked by two alpha-helices in the front layer (Z=+7.6 Å) and two in the back layer (Z=-7.6 Å).",
        "dsl_25d": " ".join(tokens_2),
    }

    # Design 3: Engineered 1bxwA00 OmpA Beta-Barrel + Inserted Central Transmembrane Alpha-Helix Plug H99
    if "1bxwA00" in rec_by_id:
        base_25d = rec_by_id["1bxwA00"]["dsl_25d"]
        _, _, _, _, _, sse_3d_map, _ = decode_25dprots_dsl(base_25d)
        all_pts = np.vstack([0.5 * (p0 + p1) for (p0, p1) in sse_3d_map.values()])
        c_mean = np.mean(all_pts, axis=0)
        # Estimate barrel axis via cross product / normal of the barrel ring
        axes_v = np.vstack([(p1 - p0) / max(1e-6, np.linalg.norm(p1 - p0)) for (p0, p1) in sse_3d_map.values()])
        # Principal direction of the barrel strands
        _, _, vt_ax = np.linalg.svd(axes_v)
        u_axis = vt_ax[0]
        az_plug = int(round(math.degrees(math.atan2(float(u_axis[0]), float(u_axis[1]))) % 360.0))
        el_plug = int(round(math.degrees(math.asin(float(np.clip(u_axis[2], -1.0, 1.0))))))
        plug_tok = (
            f"SSE25D:H99:L16:X{c_mean[0]:+06.1f}:Y{c_mean[1]:+06.1f}:Z{c_mean[2]:+06.1f}:"
            f"A{az_plug:03d}:EL{el_plug:+03d}:"
            f"C{c_mean[0]:+06.1f},{c_mean[1]:+06.1f},{c_mean[2]:+06.1f}:"
            f"U{az_plug:03d},{el_plug:+03d}:S23.0"
        )
        designs["modified_25d_1bxw_helix_plug"] = {
            "title": "2.5DProts Engineered OmpA Beta-Barrel (1bxwA00) + Coaxial Alpha-Helix Pore Plug (H99)",
            "description": "Reconstructed 2.5DProts 1bxwA00 barrel augmented with coaxial transmembrane alpha-helix H99 along the 3D barrel pore axis.",
            "dsl_25d": base_25d.replace("<EOS_25DP>", f"{plug_tok} <EOS_25DP>"),
        }

    return designs


def plot_25dprots_benchmark_showcase(showcase_items, all_domain_metrics, out_png):
    """
    Render a comprehensive 3-row figure:
      - Row 1: 2.5DProts diagrams (X, Y, A) with Z-depth color halos (coolwarm) and depth labels `z=... Å`.
      - Row 2: 3D Superposition of Reconstructed 2.5DProts 3D Model (colored ribbons) vs True Crystal SSEs (dashed black/gold).
      - Row 3: Quantitative comparison of Classic 2DProts vs 2.5DProts across all 34 benchmark domains (RMSD & TM-score).
    """
    fig = plt.figure(figsize=(18.5, 14.5), dpi=160)
    gs = fig.add_gridspec(3, 4, height_ratios=[1.0, 1.15, 0.95], hspace=0.36, wspace=0.28)
    cmap_z = plt.cm.coolwarm

    for col_idx, item in enumerate(showcase_items[:4]):
        model = item["model_3d"]
        tid = item["id"]
        layout_2d = model["layout_2d"]
        angles_2d = model["angles_2d"]
        depth_z = model["depth_z"]
        sse_order = model["sse_order"]

        # Row 1: 2.5DProts Diagram with Z-depth halo encoding
        ax_25d = fig.add_subplot(gs[0, col_idx])
        xs = [layout_2d[s][0] for s in sse_order]
        ys = [layout_2d[s][1] for s in sse_order]
        zs = [depth_z[s] for s in sse_order]
        z_max = max(5.0, float(np.max(np.abs(zs))))

        ax_25d.plot(xs, ys, "--", color="#94a3b8", lw=1.0, alpha=0.75, zorder=1)
        for s in sse_order:
            x, y = layout_2d[s]
            z = depth_z[s]
            ang = math.radians(angles_2d[s])
            is_h = s.startswith("H") or (s.startswith("_") and "H" in s[:4])
            z_norm = 0.5 + 0.5 * float(np.clip(z / z_max, -1.0, 1.0))
            z_col = cmap_z(z_norm)
            core_col = "#ef4444" if is_h else "#2563eb"
            marker = "o" if is_h else "^"

            # Outer Z-depth halo + inner SSE type symbol
            ax_25d.scatter([x], [y], c=[z_col], s=230, marker=marker, alpha=0.55, edgecolors="none", zorder=2)
            ax_25d.scatter([x], [y], c=core_col, s=95, marker=marker, edgecolors="black", linewidths=0.9, zorder=3)
            ax_25d.arrow(
                x - 0.9 * math.sin(ang), y - 0.9 * math.cos(ang),
                1.8 * math.sin(ang), 1.8 * math.cos(ang),
                head_width=0.55, head_length=0.55, fc=core_col, ec="black", lw=0.55, zorder=4
            )
            ax_25d.text(
                x + 0.45, y + 0.45, f"{s}\n(z={z:+.1f})",
                fontsize=6.3, fontweight="bold", color="#0f172a"
            )

        ax_25d.set_title(
            f"2.5DProts Diagram: {tid}\n({len(sse_order)} SSEs | halo = hloubka Z [±{z_max:.0f} Å])",
            fontsize=9.2, fontweight="bold"
        )
        ax_25d.set_xlabel("2DProts X [Å]", fontsize=8)
        ax_25d.set_ylabel("2DProts Y [Å]", fontsize=8)
        ax_25d.grid(True, linestyle=":", alpha=0.4)
        ax_25d.set_aspect("equal", adjustable="datalim")

        # Row 2: 3D Superposition (Reconstructed 2.5DProts 3D Model vs Native Crystal SSE Axes)
        ax_3d = fig.add_subplot(gs[1, col_idx], projection="3d")
        pts_pred_end = item["pts_pred_aligned"]
        pts_true_end = item["pts_true"]

        for i, s in enumerate(sse_order):
            is_h = s.startswith("H") or (s.startswith("_") and "H" in s[:4])
            c_seg = "#ef4444" if is_h else "#2563eb"
            p0_p, p1_p = pts_pred_end[2 * i], pts_pred_end[2 * i + 1]
            p0_t, p1_t = pts_true_end[2 * i], pts_true_end[2 * i + 1]

            # True crystal SSE vector in dark dashed line
            ax_3d.plot(
                [p0_t[0], p1_t[0]], [p0_t[1], p1_t[1]], [p0_t[2], p1_t[2]],
                "--", color="#1e293b", lw=2.0, alpha=0.75,
                label="Krystal (3D vzor)" if i == 0 else None
            )
            # Reconstructed 2.5DProts 3D SSE vector in solid red/blue
            ax_3d.plot(
                [p0_p[0], p1_p[0]], [p0_p[1], p1_p[1]], [p0_p[2], p1_p[2]],
                "-", color=c_seg, lw=3.4, alpha=0.92,
                label="2.5DProts 3D model" if i == 0 else None
            )
            mid_p = 0.5 * (p0_p + p1_p)
            ax_3d.text(mid_p[0], mid_p[1], mid_p[2], s, fontsize=6.5, fontweight="bold")

        ax_3d.set_title(
            f"3D Superpozice vs. Krystal ({tid})\n"
            f"2.5D RMSD: {item['rmsd_25d_center']:.2f} Å (TM={item['tm_25d_center']:.2f}) "
            f"[2D: {item['rmsd_2d_center']:.2f} Å]",
            fontsize=8.8, fontweight="bold"
        )
        ax_3d.set_xlabel("X [Å]", fontsize=7)
        ax_3d.set_ylabel("Y [Å]", fontsize=7)
        ax_3d.set_zlabel("Z [Å]", fontsize=7)
        ax_3d.tick_params(labelsize=6)
        if col_idx == 0:
            ax_3d.legend(loc="upper right", fontsize=6.8)

    # Row 3 Left (spans 2 cols): Bar comparison of Center RMSD (2DProts vs 2.5DProts) across benchmark domains
    ax_bar = fig.add_subplot(gs[2, :2])
    tids_all = [m["id"] for m in all_domain_metrics]
    rmsd_2d_all = [m["rmsd_2d_center"] for m in all_domain_metrics]
    rmsd_25d_all = [m["rmsd_25d_center"] for m in all_domain_metrics]
    x_pos = np.arange(len(tids_all))
    ax_bar.bar(x_pos - 0.19, rmsd_2d_all, width=0.38, color="#94a3b8", label=f"Klasický 2DProts → 3D (průměr {np.mean(rmsd_2d_all):.2f} Å)")
    ax_bar.bar(x_pos + 0.19, rmsd_25d_all, width=0.38, color="#10b981", label=f"2.5DProts → 3D (průměr {np.mean(rmsd_25d_all):.2f} Å, medián {np.median(rmsd_25d_all):.2f} Å)")
    ax_bar.axhline(2.5, color="#059669", linestyle="--", lw=1.1, label="Hranice vysoké přesnosti (2.5 Å)")
    ax_bar.set_xticks(x_pos)
    ax_bar.set_xticklabels(tids_all, rotation=65, ha="right", fontsize=7)
    ax_bar.set_ylabel("3D SSE-Center RMSD vůči krystalu [Å]", fontsize=8.5, fontweight="bold")
    ax_bar.set_title("Porovnání 3D RMSD napříč všemi 34 doménami benchmarku (Klasický 2DProts vs. 2.5DProts z Distogramu)", fontsize=9.8, fontweight="bold")
    ax_bar.legend(loc="upper right", fontsize=8)
    ax_bar.grid(True, axis="y", linestyle=":", alpha=0.45)

    # Row 3 Right (spans 2 cols): TM-score distribution & Scatter plot (2DProts vs 2.5DProts TM-score)
    ax_tm = fig.add_subplot(gs[2, 2:])
    tm_2d_all = [m["tm_2d_center"] for m in all_domain_metrics]
    tm_25d_all = [m["tm_25d_center"] for m in all_domain_metrics]
    ax_tm.bar(x_pos - 0.19, tm_2d_all, width=0.38, color="#cbd5e1", label=f"Klasický 2DProts → 3D (průměr TM = {np.mean(tm_2d_all):.2f})")
    ax_tm.bar(x_pos + 0.19, tm_25d_all, width=0.38, color="#2563eb", label=f"2.5DProts → 3D (průměr TM = {np.mean(tm_25d_all):.2f}, medián = {np.median(tm_25d_all):.2f})")
    ax_tm.axhline(0.50, color="#ef4444", linestyle=":", lw=1.1, label="Práh správného globálního foldu (TM = 0.50)")
    ax_tm.set_xticks(x_pos)
    ax_tm.set_xticklabels(tids_all, rotation=65, ha="right", fontsize=7)
    ax_tm.set_ylim(0.0, 1.08)
    ax_tm.set_ylabel("3D TM-score (0–1, vyšší = lepší)", fontsize=8.5, fontweight="bold")
    ax_tm.set_title("Porovnání 3D TM-score napříč všemi 34 doménami benchmarku (Klasický 2DProts vs. 2.5DProts)", fontsize=9.8, fontweight="bold")
    ax_tm.legend(loc="lower right", fontsize=8)
    ax_tm.grid(True, axis="y", linestyle=":", alpha=0.45)

    fig.suptitle(
        "2.5DProts (`2.5DPROTS-DSL`): Rozšíření 2DProts o hloubku Z a 3D sklon pro přesnou rekonstrukci a návrh 3D struktur",
        fontsize=13, fontweight="bold", y=0.985
    )
    os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


def main():
    os.makedirs(MODELS_DIR, exist_ok=True)
    os.makedirs(FIG_DIR, exist_ok=True)

    print("[1/4] Loading all 34 benchmark domains & reconstructing 2.5DProts from Distograms (ZERO 3D crystal coords)...")
    records = collect_benchmark_domains()
    print(f"      Loaded {len(records)} benchmark domains.")

    rec_by_id = {r["id"]: r for r in records}
    showcase_ids = ["1bxwA00", "1lz1A00", "1f88A00", "1adsA00"]

    print("[2/4] Benchmarking 3D reconstruction: Classic 2DProts (`2DPROTS-DSL`) vs 2.5DProts (`2.5DPROTS-DSL`)...")
    all_domain_metrics = []
    showcase_items = []

    for rec in records:
        tid = rec["id"]
        train_sub = [r for r in records if r["pdb"] != rec["pdb"]]
        rmsd_2d, tm_2d, t_2d_ms = evaluate_2d_baseline_mds_rmsd(rec, train_sub)

        t0 = time.perf_counter()
        model_25d = reconstruct_3d_structure_from_25dprots_dsl(
            rec["dsl_25d"],
            loop_lengths=rec["cm_loops"],
        )
        t_25d_build_ms = (time.perf_counter() - t0) * 1000.0
        total_25d_ms = rec["t_25d_encode_ms"] + t_25d_build_ms

        pts_true = []
        pts_pred = []
        for s in rec["sse_order"]:
            v0, v1 = np.asarray(rec["ann"].vectors[s][0], float), np.asarray(rec["ann"].vectors[s][1], float)
            pts_true.extend([v0, v1])
            p0, p1 = model_25d["sse_3d_coords"][s]
            pts_pred.extend([p0, p1])
        pts_true = np.array(pts_true, dtype=float)
        pts_pred = np.array(pts_pred, dtype=float)

        n_res_tot = sum(rec["lengths"].values())
        c_true = 0.5 * (pts_true[0::2] + pts_true[1::2])
        c_pred = 0.5 * (pts_pred[0::2] + pts_pred[1::2])
        rmsd_25d_c, tm_25d_c, _ = compute_kabsch_rmsd_and_tmscore(c_pred, c_true, n_res_total=n_res_tot)
        rmsd_25d_end, tm_25d_end, pts_pred_aligned = compute_kabsch_rmsd_and_tmscore(
            pts_pred, pts_true, n_res_total=n_res_tot
        )

        D_true = rec["D_sse_true"]
        D_3d = model_25d["D_sse_3d"]
        iu = np.triu_indices(len(rec["sse_order"]), k=1)
        local_mask = D_true[iu] <= 22.0
        local_dist_mae = float(np.mean(np.abs(D_3d[iu][local_mask] - D_true[iu][local_mask]))) if np.any(local_mask) else 0.0

        m_entry = {
            "id": tid,
            "name": rec["name"],
            "n_sse": len(rec["sse_order"]),
            "n_res": len(model_25d["residues"]),
            "rmsd_2d_center": round(rmsd_2d, 2),
            "tm_2d_center": round(tm_2d, 3),
            "rmsd_25d_center": round(rmsd_25d_c, 2),
            "tm_25d_center": round(tm_25d_c, 3),
            "rmsd_25d_endpoints": round(rmsd_25d_end, 2),
            "tm_25d_endpoints": round(tm_25d_end, 3),
            "local_dist_mae": round(local_dist_mae, 2),
            "time_25d_ms": round(total_25d_ms, 1),
        }
        all_domain_metrics.append(m_entry)

        if tid in showcase_ids:
            # Align full atomic model onto native crystal coordinate frame so viewer shows crystal-aligned orientation!
            p_mean = np.mean(pts_pred, axis=0)
            q_mean = np.mean(pts_true, axis=0)
            P_c = pts_pred - p_mean
            Q_c = pts_true - q_mean
            best_tm_a = -1.0
            best_z_sign = 1.0
            best_R = np.eye(3)
            d0 = max(1.5, 1.24 * ((max(16.0, float(n_res_tot)) - 15.0) ** (1.0 / 3.0)) - 1.8)
            for z_sign in (1.0, -1.0):
                P_test = P_c.copy()
                P_test[:, 2] *= z_sign
                U_s, _, Vt_s = np.linalg.svd(P_test.T @ Q_c)
                R_s = U_s @ Vt_s
                if np.linalg.det(R_s) < 0:
                    U_m = U_s.copy()
                    U_m[:, -1] *= -1
                    R_s = U_m @ Vt_s
                P_rot = P_test @ R_s + q_mean
                dists = np.linalg.norm(P_rot - pts_true, axis=1)
                tm_val = float(np.mean(1.0 / (1.0 + (dists / d0) ** 2)))
                if tm_val > best_tm_a:
                    best_tm_a = tm_val
                    best_z_sign = z_sign
                    best_R = R_s

            for r_atom in model_25d["residues"]:
                for aname, xyz in r_atom["atoms"].items():
                    v = xyz - p_mean
                    v[2] *= best_z_sign
                    r_atom["atoms"][aname] = v @ best_R + q_mean
            model_25d["ca_coords"] = np.array([r["atoms"]["CA"] for r in model_25d["residues"]], dtype=float)

            cif_name = f"reconstructed_25d_{tid}.cif"
            cif_path = os.path.join(MODELS_DIR, cif_name)
            pdb_path = os.path.join(MODELS_DIR, f"reconstructed_25d_{tid}.pdb")
            export_3d_structure_to_mmcif_and_pdb(model_25d, cif_path, pdb_path, entry_id=f"25DP_{tid.upper()}")
            shutil.copy2(cif_path, os.path.join(ARTIFACT_DIR, cif_name))
            # Also overwrite reconstructed_{tid}.cif in artifact dir with the accurate 2.5DProts reconstruction
            shutil.copy2(cif_path, os.path.join(ARTIFACT_DIR, f"reconstructed_{tid}.cif"))

            print(
                f"      * Showcase {tid:8s} ({len(rec['sse_order']):2d} SSEs): "
                f"2DProts RMSD = {rmsd_2d:5.2f} A (TM={tm_2d:.2f})  -->  "
                f"2.5DProts Center RMSD = {rmsd_25d_c:4.2f} A (TM={tm_25d_c:.2f}), "
                f"Endpoint RMSD = {rmsd_25d_end:4.2f} A (TM={tm_25d_end:.2f}) | time = {total_25d_ms:.1f} ms"
            )
            showcase_items.append({
                **m_entry,
                "cif_file": cif_name,
                "dsl_25d": rec["dsl_25d"],
                "sse_res_spans": model_25d["sse_res_spans"],
                "pts_pred_aligned": pts_pred_aligned,
                "pts_true": pts_true,
                "model_3d": model_25d,
            })

    # Order showcase_items according to showcase_ids
    showcase_items.sort(key=lambda x: showcase_ids.index(x["id"]))

    mean_2d_rmsd = float(np.mean([m["rmsd_2d_center"] for m in all_domain_metrics]))
    med_2d_rmsd = float(np.median([m["rmsd_2d_center"] for m in all_domain_metrics]))
    mean_2d_tm = float(np.mean([m["tm_2d_center"] for m in all_domain_metrics]))

    mean_25d_rmsd = float(np.mean([m["rmsd_25d_center"] for m in all_domain_metrics]))
    med_25d_rmsd = float(np.median([m["rmsd_25d_center"] for m in all_domain_metrics]))
    mean_25d_tm = float(np.mean([m["tm_25d_center"] for m in all_domain_metrics]))
    med_25d_tm = float(np.median([m["tm_25d_center"] for m in all_domain_metrics]))
    mean_25d_end = float(np.mean([m["rmsd_25d_endpoints"] for m in all_domain_metrics]))
    mean_25d_ms = float(np.mean([m["time_25d_ms"] for m in all_domain_metrics]))

    print(
        f"\n      [BENCHMARK SUMMARY ACROSS {len(all_domain_metrics)} DOMAINS]\n"
        f"      - Classic 2DProts -> 3D : Mean Center RMSD = {mean_2d_rmsd:.2f} A (Median = {med_2d_rmsd:.2f} A), Mean TM-score = {mean_2d_tm:.3f}\n"
        f"      - New 2.5DProts   -> 3D : Mean Center RMSD = {mean_25d_rmsd:.2f} A (Median = {med_25d_rmsd:.2f} A), Mean TM-score = {mean_25d_tm:.3f} (Median = {med_25d_tm:.3f}), "
        f"Mean Endpoint RMSD = {mean_25d_end:.2f} A, Mean Time = {mean_25d_ms:.1f} ms\n"
    )

    print("[3/4] Generating de novo & custom-engineered 3D structures via 2.5DPROTS-DSL...")
    designs_25d = build_denovo_25dprots_designs(rec_by_id)
    denovo_results = []
    for did, dinfo in designs_25d.items():
        t0 = time.perf_counter()
        model_3d = reconstruct_3d_structure_from_25dprots_dsl(dinfo["dsl_25d"])
        gen_ms = (time.perf_counter() - t0) * 1000.0

        cif_path = os.path.join(MODELS_DIR, f"{did}.cif")
        pdb_path = os.path.join(MODELS_DIR, f"{did}.pdb")
        export_3d_structure_to_mmcif_and_pdb(model_3d, cif_path, pdb_path, entry_id=did.upper())
        shutil.copy2(cif_path, os.path.join(ARTIFACT_DIR, f"{did}.cif"))

        print(f"      - {did}: {len(model_3d['sse_order'])} SSEs ({len(model_3d['residues'])} res) | time = {gen_ms:.1f} ms")
        denovo_results.append({
            "id": did,
            "title": dinfo["title"],
            "description": dinfo["description"],
            "dsl_25d": dinfo["dsl_25d"],
            "n_sse": len(model_3d["sse_order"]),
            "n_res": len(model_3d["residues"]),
            "gen_ms": round(gen_ms, 2),
            "cif_file": f"{did}.cif",
            "sse_res_spans": model_3d["sse_res_spans"],
        })

    print("[4/4] Saving 2.5DProts vs 2DProts benchmark figure and JSON summary...")
    fig_path = os.path.join(FIG_DIR, "25dprots_vs_2dprots_3d_benchmark.png")
    plot_25dprots_benchmark_showcase(showcase_items, all_domain_metrics, fig_path)
    shutil.copy2(fig_path, os.path.join(ARTIFACT_DIR, "25dprots_vs_2dprots_3d_benchmark.png"))
    shutil.copy2(fig_path, os.path.join(ARTIFACT_DIR, "2dprots_to_3d_generation.png"))

    summary_payload = {
        "aggregate_34_domains": {
            "n_domains": len(all_domain_metrics),
            "classic_2dprots_mean_center_rmsd": round(mean_2d_rmsd, 2),
            "classic_2dprots_median_center_rmsd": round(med_2d_rmsd, 2),
            "classic_2dprots_mean_tm_score": round(mean_2d_tm, 3),
            "new_25dprots_mean_center_rmsd": round(mean_25d_rmsd, 2),
            "new_25dprots_median_center_rmsd": round(med_25d_rmsd, 2),
            "new_25dprots_mean_endpoint_rmsd": round(mean_25d_end, 2),
            "new_25dprots_mean_tm_score": round(mean_25d_tm, 3),
            "new_25dprots_median_tm_score": round(med_25d_tm, 3),
            "new_25dprots_mean_time_ms": round(mean_25d_ms, 1),
        },
        "showcase_reconstructions": [
            {k: v for k, v in item.items() if k not in ("model_3d", "pts_pred_aligned", "pts_true")}
            for item in showcase_items
        ],
        "all_34_domains": all_domain_metrics,
        "denovo_25d_designs": denovo_results,
    }
    summary_json = os.path.join(OUT_DIR, "models_3d_summary.json")
    with open(summary_json, "w", encoding="utf-8") as f:
        json.dump(summary_payload, f, indent=2, ensure_ascii=False)

    print(f"Done! Summary written to {summary_json}")


if __name__ == "__main__":
    main()
