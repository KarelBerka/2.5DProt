#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""
Generates for Amino-Acid (C-alpha) Level Dimensionality Reduction (`benchmark_dimred_aa.py`):
1. Interactive HTML comparison gallery (`benchmark_dimred_aa_results/index.html`) with SSE vs AA toggle
2. Quantitative summary charts (`benchmark_dimred_aa_results/figures/benchmark_aa_summary_charts.png`)
3. Direct SSE-DimRed vs AA-DimRed comparison chart (`benchmark_dimred_aa_results/figures/comparison_sse_vs_aa_charts.png`)
4. Visual 2D Amino-Acid layout comparison panels (`benchmark_dimred_aa_results/figures/panel_aa_*.png`)
5. Side-by-side SSE vs Amino-Acid comparison panels (`benchmark_dimred_aa_results/figures/sidebyside_*.png`)
"""

import os
import json
import math
from collections import OrderedDict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import benchmark_dimred
import benchmark_extras
import benchmark_dimred_aa

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE_DIR, "benchmark_dimred_aa_results")
SSE_OUT_DIR = os.path.join(BASE_DIR, "benchmark_dimred_results")
FIG_DIR = os.path.join(OUT_DIR, "figures")
JSON_DIR = os.path.join(OUT_DIR, "layouts")
SSE_JSON_DIR = os.path.join(SSE_OUT_DIR, "layouts")
os.makedirs(FIG_DIR, exist_ok=True)

METHODS = benchmark_extras.EXTRA_METHODS + benchmark_dimred_aa.METHODS


def _sse_rgb(color_code):
    """Convert 2DProts rainbow/spectrum string or CSS color to matplotlib RGB."""
    if not color_code or not isinstance(color_code, str):
        return (0.55, 0.58, 0.62)
    c = color_code.strip()
    if c.startswith("s") and c[1:].isdigit():
        val = float(c[1:]) / 900.0
        cmap = matplotlib.colormaps["turbo"]
        return cmap(np.clip(val, 0.05, 0.95))[:3]
    try:
        return matplotlib.colors.to_rgb(c)
    except Exception:
        return (0.35, 0.55, 0.80)


def compute_aggregate_stats(results):
    """Compute overall, domain-only, and complex-only aggregate statistics per method."""
    agg = {}
    for m_id, m_title, m_desc in METHODS:
        times_all = []
        aa_lerr_all, aa_lerr_dom, aa_lerr_cmp = [], [], []
        aa_str_all, aa_str_dom, aa_str_cmp = [], [], []
        aa_knn_all, aa_knn_dom, aa_knn_cmp = [], [], []
        bb_err_all = []

        # Induced SSE-level metrics
        sse_lerr_all, sse_lerr_dom, sse_lerr_cmp = [], [], []
        sse_str_all, sse_knn_all, sse_conf_all = [], [], []

        wins_aa_lerr = 0
        wins_aa_str = 0
        wins_aa_knn = 0

        for tid, rdata in results.items():
            m_dict = rdata.get("methods", {})
            if m_id not in m_dict:
                continue
            m = m_dict[m_id]
            kind = rdata["target"]["kind"]

            times_all.append(m["time_ms"])
            aa_lerr_all.append(m["aa_local_err"])
            aa_str_all.append(m["aa_stress_1"])
            aa_knn_all.append(m["aa_knn_pres"])
            bb_err_all.append(m["backbone_err"])

            sse_lerr_all.append(m["local_err"])
            sse_str_all.append(m["stress_1"])
            sse_knn_all.append(m["knn_pres"])
            sse_conf_all.append(m["conflicts"])

            if kind == "domain":
                aa_lerr_dom.append(m["aa_local_err"])
                aa_str_dom.append(m["aa_stress_1"])
                aa_knn_dom.append(m["aa_knn_pres"])
                sse_lerr_dom.append(m["local_err"])
            else:
                aa_lerr_cmp.append(m["aa_local_err"])
                aa_str_cmp.append(m["aa_stress_1"])
                aa_knn_cmp.append(m["aa_knn_pres"])
                sse_lerr_cmp.append(m["local_err"])

            best_l = min(v["aa_local_err"] for v in m_dict.values())
            best_s = min(v["aa_stress_1"] for v in m_dict.values())
            best_k = max(v["aa_knn_pres"] for v in m_dict.values())
            if abs(m["aa_local_err"] - best_l) < 1e-6:
                wins_aa_lerr += 1
            if abs(m["aa_stress_1"] - best_s) < 1e-6:
                wins_aa_str += 1
            if abs(m["aa_knn_pres"] - best_k) < 1e-6:
                wins_aa_knn += 1

        agg[m_id] = {
            "id": m_id,
            "title": m_title,
            "desc": m_desc,
            "mean_time_ms": float(np.mean(times_all)) if times_all else 0.0,
            "mean_aa_lerr": float(np.mean(aa_lerr_all)) if aa_lerr_all else 0.0,
            "mean_aa_lerr_dom": float(np.mean(aa_lerr_dom)) if aa_lerr_dom else 0.0,
            "mean_aa_lerr_cmp": float(np.mean(aa_lerr_cmp)) if aa_lerr_cmp else 0.0,
            "mean_aa_stress": float(np.mean(aa_str_all)) if aa_str_all else 0.0,
            "mean_aa_knn": float(np.mean(aa_knn_all)) if aa_knn_all else 0.0,
            "mean_bb_err": float(np.mean(bb_err_all)) if bb_err_all else 0.0,
            "mean_sse_lerr": float(np.mean(sse_lerr_all)) if sse_lerr_all else 0.0,
            "mean_sse_lerr_dom": float(np.mean(sse_lerr_dom)) if sse_lerr_dom else 0.0,
            "mean_sse_lerr_cmp": float(np.mean(sse_lerr_cmp)) if sse_lerr_cmp else 0.0,
            "mean_sse_stress": float(np.mean(sse_str_all)) if sse_str_all else 0.0,
            "mean_sse_knn": float(np.mean(sse_knn_all)) if sse_knn_all else 0.0,
            "mean_sse_conflicts": float(np.mean(sse_conf_all)) if sse_conf_all else 0.0,
            "wins_aa_lerr": wins_aa_lerr,
            "wins_aa_str": wins_aa_str,
            "wins_aa_knn": wins_aa_knn,
        }

    evol_time = max(1.0, agg.get("evolution", {}).get("mean_time_ms", 1.0))
    for m_id in agg:
        agg[m_id]["speedup"] = evol_time / max(0.5, agg[m_id]["mean_time_ms"])
    return agg


def plot_aa_summary_charts(agg, out_path):
    """Generate a 4-panel quantitative benchmark chart for Amino-Acid level DimRed."""
    m_ids = [m[0] for m in METHODS if m[0] in agg]
    labels = [
        agg[m]["title"] if m in ("3d_pca", "2dprots_orig", "2dprots_mem15", "2dprots_mem30")
        else agg[m]["title"].split("(")[0].strip()
        for m in m_ids
    ]

    lerr_dom = [agg[m]["mean_aa_lerr_dom"] for m in m_ids]
    lerr_cmp = [agg[m]["mean_aa_lerr_cmp"] for m in m_ids]
    stress_all = [agg[m]["mean_aa_stress"] for m in m_ids]
    knn_all = [agg[m]["mean_aa_knn"] for m in m_ids]
    times_ms = [max(0.5, agg[m]["mean_time_ms"]) for m in m_ids]

    fig, axes = plt.subplots(2, 2, figsize=(16.5, 12.0), dpi=160)
    fig.patch.set_facecolor("#ffffff")
    y_pos = np.arange(len(m_ids))

    # 1. Local Amino-Acid Distance Error (<40 Å) - Domains vs Complexes
    ax = axes[0, 0]
    h = 0.38
    ax.barh(y_pos - h / 2, lerr_dom, height=h, color="#0284c7", label="Jednotlivé domény (Cα)")
    ax.barh(y_pos + h / 2, lerr_cmp, height=h, color="#ea580c", label="Celé komplexy (multi-3D Cα)")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9.5, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlabel("Průměrná lokální chyba vzdáleností Cα <40 Å [Å] (nižší = lepší)", fontsize=10)
    ax.set_title("A. Lokální přesnost na úrovni aminokyselin (<40 Å): Domény vs. Komplexy", fontsize=11.2, fontweight="bold")
    ax.grid(axis="x", linestyle="--", alpha=0.4)
    ax.legend(loc="lower right", fontsize=9)
    for i, (vd, vc) in enumerate(zip(lerr_dom, lerr_cmp)):
        ax.text(vd + 0.08, i - h / 2, f"{vd:.2f} Å", va="center", fontsize=8, color="#0c4a6e")
        ax.text(vc + 0.08, i + h / 2, f"{vc:.2f} Å", va="center", fontsize=8, color="#7c2d12")
    ax.set_xlim(0, max(max(lerr_dom), max(lerr_cmp)) * 1.22)

    # 2. Global Metric Stress-1 on C-alpha
    ax = axes[0, 1]
    colors_s = ["#10b981" if s == min(stress_all) else "#3b82f6" for s in stress_all]
    ax.barh(y_pos, stress_all, color=colors_s, height=0.62)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9.5, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlabel("Kruskalův globální Stress-1 nad maticí Cα (nižší = lepší)", fontsize=10)
    ax.set_title("B. Globální metrické zkreslení tvaru po aminokyselinách (Stress-1)", fontsize=11.2, fontweight="bold")
    ax.grid(axis="x", linestyle="--", alpha=0.4)
    for i, v in enumerate(stress_all):
        ax.text(v + 0.005, i, f"{v:.3f}", va="center", fontsize=8.5)
    ax.set_xlim(0, max(stress_all) * 1.18)

    # 3. 3D Contact k-NN Preservation (%) on Amino Acids
    ax = axes[1, 0]
    colors_k = ["#059669" if k == max(knn_all) else "#6366f1" for k in knn_all]
    ax.barh(y_pos, knn_all, color=colors_k, height=0.62)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9.5, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlabel("Zachování prostorových 3D kontaktů Cα v 2D [%] (vyšší = lepší)", fontsize=10)
    ax.set_title("C. Zachování terciárních kontaktů aminokyselin (k-NN, |i−j|>2)", fontsize=11.2, fontweight="bold")
    ax.grid(axis="x", linestyle="--", alpha=0.4)
    for i, v in enumerate(knn_all):
        ax.text(v + 0.6, i, f"{v:.1f}%", va="center", fontsize=8.5)
    ax.set_xlim(0, max(knn_all) * 1.16)

    # 4. Execution Time (ms) & Backbone Error
    ax = axes[1, 1]
    colors_t = ["#dc2626" if m == "evolution" else "#0ea5e9" for m in m_ids]
    ax.barh(y_pos, times_ms, color=colors_t, height=0.62)
    ax.set_xscale("log")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9.5, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlabel("Průměrný čas výpočtu [ms, log měřítko] (nižší = rychlejší)", fontsize=10)
    ax.set_title("D. Rychlost výpočtu po aminokyselinách (průměr ~430 Cα / cíl)", fontsize=11.2, fontweight="bold")
    ax.grid(axis="x", linestyle="--", alpha=0.4)
    for i, (v, m) in enumerate(zip(times_ms, m_ids)):
        bb = agg[m]["mean_bb_err"]
        ax.text(v * 1.12, i, f"{v:.0f} ms (páteř Δ={bb:.2f}Å)", va="center", fontsize=8.2)
    ax.set_xlim(1, max(times_ms) * 8.0)

    plt.suptitle(
        "Redukce dimenzionality po jednotlivých aminokyselinách (Cα): Srovnání 12 algoritmů nad 38 cíli",
        fontsize=13.5,
        fontweight="bold",
        y=0.99,
    )
    plt.tight_layout(rect=[0, 0, 1, 0.965])
    plt.savefig(out_path, dpi=160)
    plt.close(fig)


def plot_sse_vs_aa_comparison_chart(agg_aa, agg_sse, out_path):
    """
    Generate a 3-panel direct comparison chart contrasting SSE-level DimRed vs
    Amino-Acid-level DimRed on the exact same 2DProts SSE evaluation metrics!
    """
    m_ids = [m[0] for m in METHODS]
    labels = [m[1].split("(")[0].replace("-AA", "").strip() for m in METHODS]
    y_pos = np.arange(len(m_ids))
    h = 0.38

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 6.8), dpi=160)
    fig.patch.set_facecolor("#ffffff")

    # Panel 1: Induced 2DProts SSE Local Error (<40 Å)
    ax = axes[0]
    sse_l = [agg_sse[m]["mean_lerr"] for m in m_ids]
    aa_l = [agg_aa[m]["mean_sse_lerr"] for m in m_ids]
    ax.barh(y_pos - h / 2, sse_l, height=h, color="#64748b", label="Po prvcích SSE (Benchmark DimRed)")
    ax.barh(y_pos + h / 2, aa_l, height=h, color="#16a34a", label="Po aminokyselinách Cα (AA-DimRed)")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9.5, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlabel("Lokální chyba 2DProts <40 Å [Å] (nižší = lepší)", fontsize=9.5)
    ax.set_title("1. Lokální přesnost SSE (<40 Å)\nSSE-DimRed vs. AA-DimRed", fontsize=11, fontweight="bold")
    ax.grid(axis="x", linestyle="--", alpha=0.4)
    ax.legend(loc="lower right", fontsize=8.5)
    for i, (vs, va) in enumerate(zip(sse_l, aa_l)):
        ax.text(vs + 0.08, i - h / 2, f"{vs:.2f}", va="center", fontsize=7.8, color="#334155")
        ax.text(va + 0.08, i + h / 2, f"{va:.2f}", va="center", fontsize=7.8, color="#14532d", fontweight="bold")
    ax.set_xlim(0, max(max(sse_l), max(aa_l)) * 1.22)

    # Panel 2: Global Metric Stress-1
    ax = axes[1]
    sse_s = [agg_sse[m]["mean_stress"] for m in m_ids]
    aa_s = [agg_aa[m]["mean_sse_stress"] for m in m_ids]
    ax.barh(y_pos - h / 2, sse_s, height=h, color="#64748b", label="Po prvcích SSE")
    ax.barh(y_pos + h / 2, aa_s, height=h, color="#2563eb", label="Po aminokyselinách Cα")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9.5, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlabel("Globální Kruskalův Stress-1 (nižší = lepší)", fontsize=9.5)
    ax.set_title("2. Globální Stress-1 (celkový tvar)\nSSE-DimRed vs. AA-DimRed", fontsize=11, fontweight="bold")
    ax.grid(axis="x", linestyle="--", alpha=0.4)
    ax.legend(loc="lower right", fontsize=8.5)
    for i, (vs, va) in enumerate(zip(sse_s, aa_s)):
        ax.text(vs + 0.005, i - h / 2, f"{vs:.3f}", va="center", fontsize=7.8, color="#334155")
        ax.text(va + 0.005, i + h / 2, f"{va:.3f}", va="center", fontsize=7.8, color="#1e3a8a", fontweight="bold")
    ax.set_xlim(0, max(max(sse_s), max(aa_s)) * 1.22)

    # Panel 3: 4-NN Topology Preservation (%)
    ax = axes[2]
    sse_k = [agg_sse[m]["mean_knn"] for m in m_ids]
    aa_k = [agg_aa[m]["mean_sse_knn"] for m in m_ids]
    ax.barh(y_pos - h / 2, sse_k, height=h, color="#64748b", label="Po prvcích SSE")
    ax.barh(y_pos + h / 2, aa_k, height=h, color="#7c3aed", label="Po aminokyselinách Cα")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9.5, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlabel("Zachování 4-NN sousedství SSE [%] (vyšší = lepší)", fontsize=9.5)
    ax.set_title("3. Zachování lokální 3D topologie (4-NN)\nSSE-DimRed vs. AA-DimRed", fontsize=11, fontweight="bold")
    ax.grid(axis="x", linestyle="--", alpha=0.4)
    ax.legend(loc="lower right", fontsize=8.5)
    for i, (vs, va) in enumerate(zip(sse_k, aa_k)):
        ax.text(vs + 0.7, i - h / 2, f"{vs:.1f}%", va="center", fontsize=7.8, color="#334155")
        ax.text(va + 0.7, i + h / 2, f"{va:.1f}%", va="center", fontsize=7.8, color="#4c1d95", fontweight="bold")
    ax.set_xlim(0, max(max(sse_k), max(aa_k)) * 1.20)

    plt.suptitle(
        "Přímé srovnání: Redukce dimenze po prvcích sekundární struktury (SSE) vs. po jednotlivých aminokyselinách (Cα)",
        fontsize=13,
        fontweight="bold",
        y=0.99,
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    plt.savefig(out_path, dpi=160)
    plt.close(fig)


def _draw_aa_on_axis(ax, jdata, m_id=None):
    """Draw the 2D Amino-Acid trace, SSE ribbons/arrows, and individual C-alpha beads on a matplotlib axis."""
    aa_coords = np.asarray(jdata.get("aa_coords", []), dtype=float)
    aa_meta = jdata.get("aa_meta", [])
    mem = jdata.get("membrane", {})
    n = len(aa_coords)
    if n == 0:
        return

    if isinstance(mem, dict) and mem.get("is_transmembrane", False) and m_id in ("2dprots_mem15", "2dprots_mem30"):
        yc = float(mem.get("y_center_2d", 0.0))
        hh = float(mem.get("band_half_height_2d", 3.0))
        ax.axhspan(yc - hh, yc + hh, color="#fde68a", alpha=0.38, zorder=0)
        ax.axhline(yc + hh, color="#d97706", linestyle="--", linewidth=0.9, alpha=0.7, zorder=1)
        ax.axhline(yc - hh, color="#d97706", linestyle="--", linewidth=0.9, alpha=0.7, zorder=1)

    xs = aa_coords[:, 0]
    ys = aa_coords[:, 1]

    # 1. Draw continuous polypeptide backbone segments per (chain_id, dom_idx)
    seg_x, seg_y = [xs[0]], [ys[0]]
    for i in range(1, n):
        m_prev = aa_meta[i - 1]
        m_curr = aa_meta[i]
        same_chain = (m_prev[2] == m_curr[2]) and (m_prev[3] == m_curr[3])
        seq_gap = abs(int(m_curr[0]) - int(m_prev[0]))
        if same_chain and seq_gap <= 15:
            seg_x.append(xs[i])
            seg_y.append(ys[i])
        else:
            if len(seg_x) >= 2:
                ax.plot(seg_x, seg_y, color="#94a3b8", linewidth=0.9, alpha=0.55, zorder=1)
            seg_x, seg_y = [xs[i]], [ys[i]]
    if len(seg_x) >= 2:
        ax.plot(seg_x, seg_y, color="#94a3b8", linewidth=0.9, alpha=0.55, zorder=1)

    # 2. Group amino acids by SSE to draw translucent SSE ribbons / directional arrows + C-alpha trace
    sse_groups = OrderedDict()
    for i, meta in enumerate(aa_meta):
        s_lbl = meta[4]
        if s_lbl:
            sse_groups.setdefault(s_lbl, []).append(i)

    for s_lbl, idxs in sse_groups.items():
        if len(idxs) < 2:
            continue
        ss_type = aa_meta[idxs[0]][5]
        rgb = _sse_rgb(aa_meta[idxs[0]][6])
        gx = xs[idxs]
        gy = ys[idxs]
        # Connect consecutive C-alpha atoms inside the SSE
        ax.plot(gx, gy, color=rgb, linewidth=1.8 if ss_type == "H" else 1.5, alpha=0.92, zorder=3)
        if ss_type == "H":
            # Translucent helix axis cylinder
            ax.plot([gx[0], gx[-1]], [gy[0], gy[-1]], color=rgb, linewidth=5.0, solid_capstyle="round", alpha=0.25, zorder=2)
        elif ss_type == "E":
            # Beta-strand directional arrow at C-terminal residue
            ax.annotate(
                "",
                xy=(gx[-1], gy[-1]),
                xytext=(gx[0], gy[0]),
                arrowprops=dict(arrowstyle="-|>,head_length=0.42,head_width=0.32", lw=2.0, color=rgb, alpha=0.55),
                zorder=2,
            )

    # 3. Scatter individual Amino Acid (C-alpha) positions
    loop_mask = np.array([m[5] == "L" for m in aa_meta], dtype=bool)
    ss_mask = ~loop_mask
    pt_size_ss = 10.0 if n <= 350 else (5.5 if n <= 900 else 2.5)
    pt_size_lp = 4.5 if n <= 350 else (2.5 if n <= 900 else 1.2)

    if np.any(loop_mask):
        ax.scatter(xs[loop_mask], ys[loop_mask], s=pt_size_lp, color="#94a3b8", alpha=0.65, edgecolors="none", zorder=3)
    if np.any(ss_mask):
        ss_colors = [_sse_rgb(aa_meta[i][6]) for i in np.where(ss_mask)[0]]
        ax.scatter(
            xs[ss_mask],
            ys[ss_mask],
            s=pt_size_ss,
            c=ss_colors,
            edgecolors="#0f172a",
            linewidths=0.22 if n <= 600 else 0.1,
            alpha=0.95,
            zorder=4,
        )

    cx = 0.5 * (float(np.min(xs)) + float(np.max(xs)))
    cy = 0.5 * (float(np.min(ys)) + float(np.max(ys)))
    span = max(float(np.max(xs) - np.min(xs)), float(np.max(ys) - np.min(ys)), 10.0) * 0.58
    ax.set_xlim(cx - span, cx + span)
    ax.set_ylim(cy - span, cy + span)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_edgecolor("#cbd5e1")


def render_aa_comparison_panel_png(tid, title, results_entry, out_path):
    """
    Render a clean 4x4 (16-panel) matplotlib multi-panel figure showing:
    Row 1: 4 Reference / 2DProts views (`3d_pca`, `2dprots_orig`, `2dprots_mem15`, `2dprots_mem30`) mapped to Cα
    Rows 2-4: All 12 Amino-Acid dimensionality reduction algorithms (`evolution` .. `umap`)
    """
    fig, axes = plt.subplots(4, 4, figsize=(16.5, 15.0), dpi=140)
    fig.patch.set_facecolor("#ffffff")

    for idx, (m_id, m_title, _) in enumerate(METHODS):
        ax = axes[idx // 4, idx % 4]
        is_ref_row = idx < 4
        ax.set_facecolor("#f1f5f9" if is_ref_row else "#f8fafc")
        jpath = os.path.join(JSON_DIR, f"{tid}_{m_id}.json")
        if not os.path.exists(jpath):
            ax.set_visible(False)
            continue
        with open(jpath, "r", encoding="utf-8") as f:
            jdata = json.load(f)

        _draw_aa_on_axis(ax, jdata, m_id=m_id)
        for spine in ax.spines.values():
            spine.set_edgecolor("#6366f1" if is_ref_row else "#cbd5e1")
            spine.set_linewidth(1.4 if is_ref_row else 1.0)

        m_info = results_entry.get("methods", {}).get(m_id, {})
        short_name = m_title if is_ref_row else m_title.split("(")[0].strip()
        if m_id == "3d_pca" and "pca_spans_3d" in jdata:
            sp3 = jdata["pca_spans_3d"]
            sub_line = f"PC1×2×3: {sp3[0]:.0f}×{sp3[1]:.0f}×{sp3[2]:.0f}Å | Cα err={m_info.get('aa_local_err', 0):.2f}Å"
        elif m_id in ("2dprots_mem15", "2dprots_mem30"):
            tm_tag = "TM pás" if m_info.get("is_tm") else "Glob."
            sub_line = f"[{tm_tag}] Cα err={m_info.get('aa_local_err', 0):.2f}Å | SSE={m_info.get('local_err', 0):.2f}Å | kNN={m_info.get('aa_knn_pres', 0):.0f}%"
        else:
            sub_line = (
                f"Cα err={m_info.get('aa_local_err', 0):.2f}Å | SSE={m_info.get('local_err', 0):.2f}Å | "
                f"kNN={m_info.get('aa_knn_pres', 0):.0f}% | {m_info.get('time_ms', 0):.0f}ms"
            )
        ax.set_title(
            f"{short_name}\n{sub_line}",
            fontsize=8.1,
            fontweight="bold",
            pad=3.5,
        )

    n_aa = results_entry.get("n_aa", 0)
    plt.suptitle(f"{title} — {n_aa} aminokyselin (Cα) [16 variant]", fontsize=13.5, fontweight="bold", y=0.992)
    plt.tight_layout(rect=[0, 0, 1, 0.968])
    plt.savefig(out_path, dpi=140)
    plt.close(fig)


def render_sidebyside_sse_vs_aa_png(tid, title, sse_res_entry, aa_res_entry, out_path):
    """
    Render a 2x4 side-by-side visual comparison figure for key representative methods
    (PCA, Metric MDS, Sammon wMDS, Kamada-Kawai / UMAP):
      Top row: Traditional SSE-level DimRed (`benchmark_dimred.py`)
      Bottom row: Amino-Acid (C-alpha) level DimRed (`benchmark_dimred_aa.py`)
    """
    sel_methods = [
        ("pca", "PCA (3D→2D)"),
        ("mds", "Metric MDS (SMACOF)"),
        ("sammon", "Sammon / Local wMDS"),
        ("kamada_kawai", "Kamada-Kawai (Pružinový)"),
    ]
    fig, axes = plt.subplots(2, 4, figsize=(16, 8.2), dpi=145)
    fig.patch.set_facecolor("#ffffff")

    for col_i, (m_id, m_label) in enumerate(sel_methods):
        # Top row: SSE-level layout
        ax_top = axes[0, col_i]
        ax_top.set_facecolor("#f8fafc")
        sse_jpath = os.path.join(SSE_JSON_DIR, f"{tid}_{m_id}.json")
        if os.path.exists(sse_jpath):
            with open(sse_jpath, "r", encoding="utf-8") as f:
                sj = json.load(f)
            lay = sj.get("layout", {})
            ang = sj.get("angles", {})
            sz = sj.get("size", {})
            ht = sj.get("h_type", {})
            col = sj.get("color", {})
            shown = [s for s in lay if ht.get(s, 3) in (1, 2)]
            if len(shown) >= 2:
                ax_top.plot([lay[s][0] for s in shown], [lay[s][1] for s in shown], color="#94a3b8", lw=0.8, alpha=0.45, zorder=1)
            all_x, all_y = [], []
            for s in shown:
                x, y = lay[s]
                a_rad = math.radians(float(ang.get(s, 0.0)))
                length = max(1.2, float(sz.get(s, 5)) * 0.22)
                dx, dy = length * math.sin(a_rad), length * math.cos(a_rad)
                x0, y0, x1, y1 = x - dx, y - dy, x + dx, y + dy
                all_x.extend([x0, x1])
                all_y.extend([y0, y1])
                rgb = _sse_rgb(col.get(s, "blue"))
                if s.startswith("H") or s.startswith("_H"):
                    ax_top.plot([x0, x1], [y0, y1], color=rgb, lw=4.8, solid_capstyle="round", zorder=3)
                else:
                    ax_top.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>,head_length=0.45,head_width=0.35", lw=2.4, color=rgb), zorder=4)
            if all_x and all_y:
                cx, cy = 0.5 * (min(all_x) + max(all_x)), 0.5 * (min(all_y) + max(all_y))
                span = max(max(all_x) - min(all_x), max(all_y) - min(all_y), 10.0) * 0.60
                ax_top.set_xlim(cx - span, cx + span)
                ax_top.set_ylim(cy - span, cy + span)
            ax_top.set_aspect("equal")
            ax_top.set_xticks([])
            ax_top.set_yticks([])
            for spine in ax_top.spines.values():
                spine.set_edgecolor("#cbd5e1")
            sm = sse_res_entry.get("methods", {}).get(m_id, {})
            ax_top.set_title(
                f"[SSE] {m_label}\nSSE err={sm.get('local_err', 0):.2f}Å | Str={sm.get('stress_1', 0):.2f} | 4-NN={sm.get('knn_pres', 0):.0f}%",
                fontsize=8.8,
                fontweight="bold",
                color="#334155",
            )

        # Bottom row: Amino-Acid level layout
        ax_bot = axes[1, col_i]
        ax_bot.set_facecolor("#f0fdf4" if col_i in (1, 2) else "#f8fafc")
        aa_jpath = os.path.join(JSON_DIR, f"{tid}_{m_id}.json")
        if os.path.exists(aa_jpath):
            with open(aa_jpath, "r", encoding="utf-8") as f:
                aj = json.load(f)
            _draw_aa_on_axis(ax_bot, aj)
            am = aa_res_entry.get("methods", {}).get(m_id, {})
            ax_bot.set_title(
                f"[Aminokyseliny Cα] {m_label}\nSSE err={am.get('local_err', 0):.2f}Å | Str={am.get('stress_1', 0):.2f} | 4-NN={am.get('knn_pres', 0):.0f}%",
                fontsize=8.8,
                fontweight="bold",
                color="#166534",
            )

    plt.suptitle(
        f"Srovnání reprezentace po SSE (nahoře) vs. po jednotlivých aminokyselinách Cα (dole): {title}",
        fontsize=12.5,
        fontweight="bold",
        y=0.99,
    )
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(out_path, dpi=145)
    plt.close(fig)


def generate_html_gallery(results, agg, agg_sse, html_path):
    """Generate the interactive HTML gallery (`benchmark_dimred_aa_results/index.html`) with AA vs SSE switcher."""
    categories = OrderedDict()
    for t in benchmark_dimred.BENCHMARK_TARGETS:
        tid = t["id"]
        if tid not in results:
            continue
        categories.setdefault(t["category"], []).append((tid, results[tid]))

    best_aa_lerr = min(a["mean_aa_lerr"] for a in agg.values())
    best_sse_lerr = min(a["mean_sse_lerr"] for a in agg.values())
    best_aa_str = min(a["mean_aa_stress"] for a in agg.values())
    best_aa_knn = max(a["mean_aa_knn"] for a in agg.values())
    best_time = min(a["mean_time_ms"] for a in agg.values())

    agg_rows_html = []
    for m_id, m_title, m_desc in METHODS:
        a = agg[m_id]
        s_ref = agg_sse.get(m_id, {})
        sse_orig_lerr = s_ref.get("mean_lerr", 0.0)
        delta_lerr = a["mean_sse_lerr"] - sse_orig_lerr
        delta_str = f"{delta_lerr:+.2f} Å"

        c_aal = "best-cell" if abs(a["mean_aa_lerr"] - best_aa_lerr) < 0.03 else ""
        c_ssl = "best-cell" if abs(a["mean_sse_lerr"] - best_sse_lerr) < 0.03 else ""
        c_str = "best-cell" if abs(a["mean_aa_stress"] - best_aa_str) < 0.005 else ""
        c_knn = "best-cell" if abs(a["mean_aa_knn"] - best_aa_knn) < 0.4 else ""
        c_tim = "best-cell" if abs(a["mean_time_ms"] - best_time) < 5.0 else ""

        agg_rows_html.append(f"""
        <tr>
            <td><strong>{m_title}</strong><div class="sub">{m_desc}</div></td>
            <td class="{c_tim}">{a['mean_time_ms']:.1f} ms</td>
            <td class="{c_aal}"><strong>{a['mean_aa_lerr']:.2f} Å</strong></td>
            <td class="{c_str}">{a['mean_aa_stress']:.3f}</td>
            <td class="{c_knn}">{a['mean_aa_knn']:.1f}%</td>
            <td>{a['mean_bb_err']:.2f} Å</td>
            <td class="{c_ssl}"><strong>{a['mean_sse_lerr']:.2f} Å</strong> <span class="badge">{delta_str} vs SSE</span></td>
            <td>{a['mean_sse_knn']:.1f}%</td>
            <td>{a['wins_aa_lerr']} / {a['wins_aa_str']} / {a['wins_aa_knn']}</td>
        </tr>
        """)

    sections_html = []
    for cat_name, items in categories.items():
        cards_html = []
        for tid, rdata in items:
            t = rdata["target"]
            kind_badge = (
                '<span class="kind-badge complex">Komplex (Multi-3D)</span>'
                if t["kind"] == "complex"
                else '<span class="kind-badge domain">Doména</span>'
            )
            m_dict = rdata["methods"]
            min_lerr = min(v["aa_local_err"] for v in m_dict.values()) if m_dict else 0
            min_stress = min(v["aa_stress_1"] for v in m_dict.values()) if m_dict else 0
            max_knn = max(v["aa_knn_pres"] for v in m_dict.values()) if m_dict else 0

            grid_items = []
            for m_id, m_title, _ in METHODS:
                if m_id not in m_dict:
                    continue
                m = m_dict[m_id]
                svg_aa_rel = f"svgs/{tid}_{m_id}.svg"
                svg_sse_rel = f"../benchmark_dimred_results/svgs/{tid}_{m_id}.svg"
                is_best_l = abs(m["aa_local_err"] - min_lerr) < 1e-5
                border_cls = "method-card best-card" if is_best_l else "method-card"
                short_title = m_title.split("(")[0].strip()
                grid_items.append(f"""
                <div class="{border_cls}" data-method="{m_id}">
                    <div class="method-header">
                        <span class="method-name">{short_title}</span>
                        <span class="method-time">{m['time_ms']:.0f} ms</span>
                    </div>
                    <a href="{svg_aa_rel}" target="_blank" class="svg-link" data-aa="{svg_aa_rel}" data-sse="{svg_sse_rel}">
                        <img src="{svg_aa_rel}" class="preview-img" data-aa="{svg_aa_rel}" data-sse="{svg_sse_rel}" alt="{tid} {m_id}" loading="lazy" />
                    </a>
                    <div class="method-metrics">
                        <span class="{'m-win' if is_best_l else ''}" title="Lokální chyba Cα <40 Å">Cα: {m['aa_local_err']:.2f}Å</span>
                        <span title="Odvozená lokální chyba SSE <40 Å">SSE: {m['local_err']:.2f}Å</span>
                        <span class="{'m-win' if abs(m['aa_knn_pres'] - max_knn) < 1e-5 else ''}" title="Zachování Cα k-NN">kNN: {m['aa_knn_pres']:.0f}%</span>
                    </div>
                </div>
                """)

            cards_html.append(f"""
            <div class="target-block" id="target-{tid}" data-kind="{t['kind']}">
                <div class="target-header">
                    <div>
                        <h3>{t['name']} <code class="tid-code">{tid}</code> {kind_badge}</h3>
                        <p class="target-desc">{t['desc']} — <strong>{rdata['n_aa']} aminokyselin (Cα)</strong> v {rdata['n_sses']} SSE ({rdata['n_clusters']} klastrů)</p>
                    </div>
                </div>
                <div class="methods-grid">
                    {''.join(grid_items)}
                </div>
            </div>
            """)

        sections_html.append(f"""
        <section class="category-section">
            <h2 class="category-title">{cat_name}</h2>
            {''.join(cards_html)}
        </section>
        """)

    html = f"""<!DOCTYPE html>
<html lang="cs">
<head>
    <meta charset="UTF-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0" />
    <title>2DProts – Redukce dimenze po jednotlivých aminokyselinách (Cα) – 12 algoritmů</title>
    <style>
        :root {{
            --bg: #f8fafc;
            --card: #ffffff;
            --text: #0f172a;
            --muted: #475569;
            --border: #e2e8f0;
            --accent: #0284c7;
            --win: #dcfce7;
            --win-text: #166534;
        }}
        * {{ box-sizing: border-box; }}
        body {{
            margin: 0;
            padding: 24px 32px;
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            background: var(--bg);
            color: var(--text);
            line-height: 1.5;
        }}
        header {{
            background: linear-gradient(135deg, #0f172a 0%, #0369a1 100%);
            color: #fff;
            padding: 28px 32px;
            border-radius: 14px;
            margin-bottom: 28px;
            box-shadow: 0 10px 25px rgba(15, 23, 42, 0.12);
        }}
        header h1 {{ margin: 0 0 8px 0; font-size: 1.75rem; }}
        header p {{ margin: 0; color: #e0f2fe; font-size: 1rem; }}
        .summary-card {{
            background: var(--card);
            border: 1px solid var(--border);
            border-radius: 12px;
            padding: 20px 24px;
            margin-bottom: 28px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.03);
            overflow-x: auto;
        }}
        table.summary-table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 0.9rem;
        }}
        table.summary-table th, table.summary-table td {{
            padding: 10px 12px;
            border-bottom: 1px solid var(--border);
            text-align: left;
        }}
        table.summary-table th {{
            background: #f1f5f9;
            font-weight: 700;
            color: #1e293b;
        }}
        .sub {{ font-size: 0.78rem; color: var(--muted); }}
        .best-cell {{ background: var(--win); color: var(--win-text); font-weight: 700; }}
        .badge {{
            display: inline-block;
            background: #dcfce7;
            color: #166534;
            font-size: 0.75rem;
            padding: 1px 6px;
            border-radius: 99px;
            margin-left: 4px;
        }}
        .filter-bar {{
            position: sticky;
            top: 0;
            z-index: 50;
            background: rgba(255,255,255,0.95);
            backdrop-filter: blur(8px);
            border: 1px solid var(--border);
            padding: 12px 20px;
            border-radius: 10px;
            margin-bottom: 24px;
            display: flex;
            gap: 12px;
            align-items: center;
            flex-wrap: wrap;
        }}
        .filter-btn {{
            border: 1px solid #cbd5e1;
            background: #fff;
            padding: 6px 14px;
            border-radius: 8px;
            cursor: pointer;
            font-weight: 600;
            font-size: 0.85rem;
        }}
        .filter-btn.active {{
            background: var(--accent);
            color: #fff;
            border-color: var(--accent);
        }}
        .category-section {{ margin-bottom: 40px; }}
        .category-title {{
            font-size: 1.35rem;
            padding-bottom: 8px;
            border-bottom: 2px solid #cbd5e1;
            margin-bottom: 18px;
            color: #0f172a;
        }}
        .target-block {{
            background: var(--card);
            border: 1px solid var(--border);
            border-radius: 12px;
            padding: 18px 20px;
            margin-bottom: 22px;
            box-shadow: 0 2px 6px rgba(0,0,0,0.03);
        }}
        .target-header h3 {{ margin: 0 0 4px 0; font-size: 1.12rem; }}
        .target-desc {{ margin: 0 0 14px 0; color: var(--muted); font-size: 0.9rem; }}
        .tid-code {{
            background: #f1f5f9;
            padding: 2px 7px;
            border-radius: 5px;
            font-size: 0.88rem;
            color: #334155;
        }}
        .kind-badge {{
            display: inline-block;
            padding: 2px 9px;
            border-radius: 99px;
            font-size: 0.76rem;
            font-weight: 700;
            margin-left: 6px;
        }}
        .kind-badge.domain {{ background: #dbeafe; color: #1d4ed8; }}
        .kind-badge.complex {{ background: #fef3c7; color: #b45309; }}
        .methods-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(215px, 1fr));
            gap: 12px;
        }}
        .method-card {{
            border: 1px solid var(--border);
            border-radius: 9px;
            padding: 8px;
            background: #fff;
            display: flex;
            flex-direction: column;
            justify-content: space-between;
            transition: transform 0.12s ease, box-shadow 0.12s ease;
        }}
        .method-card:hover {{
            transform: translateY(-2px);
            box-shadow: 0 6px 16px rgba(0,0,0,0.08);
        }}
        .method-card.best-card {{
            border: 2px solid #22c55e;
            background: #f0fdf4;
        }}
        .method-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            font-size: 0.8rem;
            margin-bottom: 6px;
        }}
        .method-name {{ font-weight: 700; color: #1e293b; }}
        .method-time {{ color: #64748b; font-size: 0.74rem; }}
        .method-card img {{
            width: 100%;
            height: 175px;
            object-fit: contain;
            background: #ffffff;
            border-radius: 6px;
            border: 1px solid #f1f5f9;
        }}
        .method-metrics {{
            display: flex;
            justify-content: space-between;
            font-size: 0.74rem;
            color: #475569;
            margin-top: 6px;
            padding-top: 5px;
            border-top: 1px solid #f1f5f9;
        }}
        .m-win {{ color: #15803d; font-weight: 700; }}
    </style>
</head>
<body>
    <header>
        <h1>2DProts: 2D reprezentace po jednotlivých aminokyselinách (Cα) – 12 algoritmů DimRed</h1>
        <p>Komplexní benchmark nad 38 cíli (cytochromy P450, enzymy EC 1–7, GPCR, SLC a ABC transportéry, β-barely, 14-3-3, hemoglobin, KcsA, akvaporin, F1-ATPáza a GroEL) počítaný přímo z 3D pozic aminokyselin</p>
    </header>

    <div class="summary-card">
        <h2 style="margin-top:0;">Celkové kvantitativní srovnání algoritmů na úrovni aminokyselin (průměr přes {len(results)} struktur)</h2>
        <table class="summary-table">
            <thead>
                <tr>
                    <th>Algoritmus</th>
                    <th>Průměrný čas</th>
                    <th>Lokální chyba Cα &lt;40 Å</th>
                    <th>Globální Stress-1 (Cα)</th>
                    <th>3D kontakty k-NN (Cα)</th>
                    <th>Chyba vazby páteře</th>
                    <th>Odvozená chyba SSE &lt;40 Å</th>
                    <th>Zachování 4-NN (SSE)</th>
                    <th>Počet 1. míst (Err / Str / kNN)</th>
                </tr>
            </thead>
            <tbody>
                {''.join(agg_rows_html)}
            </tbody>
        </table>
    </div>

    <div class="filter-bar">
        <strong>Filtr struktur:</strong>
        <button class="filter-btn kind-btn active" onclick="filterKind('all', this)">Vše ({len(results)})</button>
        <button class="filter-btn kind-btn" onclick="filterKind('domain', this)">Pouze jednotlivé domény</button>
        <button class="filter-btn kind-btn" onclick="filterKind('complex', this)">Pouze celé komplexy (Multi-3D)</button>
        <span style="margin-left:18px;"><strong>Režim zobrazení SVG:</strong></span>
        <button class="filter-btn mode-btn active" onclick="switchMode('aa', this)">Po aminokyselinách (Cα + smyčky)</button>
        <button class="filter-btn mode-btn" onclick="switchMode('sse', this)">Původní po SSE (pro porovnání)</button>
    </div>

    {''.join(sections_html)}

    <script>
        function filterKind(kind, btn) {{
            document.querySelectorAll('.kind-btn').forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
            document.querySelectorAll('.target-block').forEach(el => {{
                if (kind === 'all' || el.getAttribute('data-kind') === kind) {{
                    el.style.display = '';
                }} else {{
                    el.style.display = 'none';
                }}
            }});
        }}
        function switchMode(mode, btn) {{
            document.querySelectorAll('.mode-btn').forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
            document.querySelectorAll('.preview-img').forEach(img => {{
                img.src = img.getAttribute('data-' + mode);
            }});
            document.querySelectorAll('.svg-link').forEach(a => {{
                a.href = a.getAttribute('data-' + mode);
            }});
        }}
    </script>
</body>
</html>
"""
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html)


def main():
    summary_path = os.path.join(OUT_DIR, "results_summary.json")
    with open(summary_path, "r", encoding="utf-8") as f:
        results = json.load(f)

    sse_summary_path = os.path.join(SSE_OUT_DIR, "results_summary.json")
    with open(sse_summary_path, "r", encoding="utf-8") as f:
        sse_results = json.load(f)

    sse_agg_path = os.path.join(SSE_OUT_DIR, "aggregate_stats.json")
    with open(sse_agg_path, "r", encoding="utf-8") as f:
        agg_sse = json.load(f)

    agg = compute_aggregate_stats(results)
    with open(os.path.join(OUT_DIR, "aggregate_stats.json"), "w", encoding="utf-8") as f:
        json.dump(agg, f, indent=2, ensure_ascii=False)

    chart_path = os.path.join(FIG_DIR, "benchmark_aa_summary_charts.png")
    plot_aa_summary_charts(agg, chart_path)
    print(f"Saved AA summary chart to {chart_path}")

    comp_chart_path = os.path.join(FIG_DIR, "comparison_sse_vs_aa_charts.png")
    plot_sse_vs_aa_comparison_chart(agg, agg_sse, comp_chart_path)
    print(f"Saved SSE vs AA comparison chart to {comp_chart_path}")

    featured = [
        ("1tqnA00", "Cytochrom P450 3A4 (1tqnA00 – globulární enzym, 37 SSE)"),
        ("1adsA00", "EC 1: Aldózareduktáza – kanonický (β/α)8 TIM barel (1adsA00)"),
        ("1f88_dimer", "GPCR: Bovinní rodopsin – dimerní 7TM komplex (1f88_dimer)"),
        ("1pv6_complex", "SLC / MFS: Laktózová permeáza LacY – celý 12TM komplex N+C (1pv6_complex)"),
        ("2hyd_dimer", "ABC transportér Sav1866 – celý homodimerní komplex (2hyd_dimer, EC 7)"),
        ("1emaA00", "Beta-barel: Zelený fluorescenční protein GFP – 11-vláknový barel (1emaA00)"),
        ("2omf_trimer", "Beta-barel: Porin OmpF – celý biologický trimer 3×16 β-vláken (2omf_trimer)"),
        ("1qjb_dimer", "14-3-3ζ protein – homodimerní komplex ve tvaru W (1qjb_dimer)"),
        ("2hhb_tetramer", "Hemoglobin – celý kvartérní tetramer α2β2 (2hhb_tetramer)"),
        ("1bl8_tetramer", "Draslíkový kanál KcsA – celý biologický tetramer (1bl8_tetramer)"),
        ("1j4n_tetramer", "Akvaporin AQP1 – celý biologický tetramer 4× kanál (1j4n_tetramer)"),
        ("1bmf_complex", "Extra velký komplex: F1-ATPáza α3β3γ motor – 19 domén (1bmf_complex)"),
        ("1aon_ring", "Extra velký komplex: Chaperonin GroEL – heptamerní kruh 7 podjednotek (1aon_ring)"),
    ]
    for tid, title in featured:
        if tid in results:
            p_out = os.path.join(FIG_DIR, f"panel_aa_{tid}.png")
            render_aa_comparison_panel_png(tid, title, results[tid], p_out)
            print(f"Saved AA panel {p_out}")
            if tid in sse_results:
                sbs_out = os.path.join(FIG_DIR, f"sidebyside_{tid}.png")
                render_sidebyside_sse_vs_aa_png(tid, title, sse_results[tid], results[tid], sbs_out)
                print(f"Saved side-by-side panel {sbs_out}")

    html_path = os.path.join(OUT_DIR, "index.html")
    generate_html_gallery(results, agg, agg_sse, html_path)
    print(f"Saved interactive AA HTML gallery to {html_path}")


if __name__ == "__main__":
    main()
