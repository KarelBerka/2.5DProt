#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""
Generates:
1. Interactive HTML comparison gallery (`benchmark_dimred_results/index.html`) with all 16 views:
   - 4 Reference / 2DProts views (`3d_pca`, `2dprots_orig`, `2dprots_mem15`, `2dprots_mem30`)
   - 12 Dimensionality Reduction algorithms (`evolution`, `pca`, `pcoa`, `mds`, `nmmds`, `sammon`,
     `isomap`, `spectral`, `lle`, `kamada_kawai`, `tsne`, `umap`)
2. Quantitative summary charts (`benchmark_dimred_results/figures/benchmark_summary_charts.png`)
3. 4x4 (16-panel) visual comparison figures (`benchmark_dimred_results/figures/panel_*.png`)
"""

import os
import sys
import json
import math
import shutil
from collections import OrderedDict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import benchmark_dimred
import benchmark_extras

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE_DIR, "benchmark_dimred_results")
FIG_DIR = os.path.join(OUT_DIR, "figures")
JSON_DIR = os.path.join(OUT_DIR, "layouts")
os.makedirs(FIG_DIR, exist_ok=True)

ALL_METHODS = benchmark_extras.EXTRA_METHODS + benchmark_dimred.METHODS


def compute_aggregate_stats(results):
    """Compute overall, domain-only, and complex-only aggregate statistics across all 16 views."""
    agg = OrderedDict()
    for m_id, m_title, m_desc in ALL_METHODS:
        times_all, lerr_all, stress_all, knn_all, conf_all = [], [], [], [], []
        times_dom, lerr_dom, stress_dom, knn_dom = [], [], [], []
        times_cmp, lerr_cmp, stress_cmp, knn_cmp = [], [], [], []
        wins_lerr = 0
        wins_stress = 0
        wins_knn = 0

        for tid, rdata in results.items():
            m_dict = rdata.get("methods", {})
            if m_id not in m_dict:
                continue
            m = m_dict[m_id]
            kind = rdata["target"]["kind"]
            times_all.append(m["time_ms"])
            lerr_all.append(m["local_err"])
            stress_all.append(m["stress_1"])
            knn_all.append(m["knn_pres"])
            conf_all.append(m["conflicts"])

            if kind == "domain":
                times_dom.append(m["time_ms"])
                lerr_dom.append(m["local_err"])
                stress_dom.append(m["stress_1"])
                knn_dom.append(m["knn_pres"])
            else:
                times_cmp.append(m["time_ms"])
                lerr_cmp.append(m["local_err"])
                stress_cmp.append(m["stress_1"])
                knn_cmp.append(m["knn_pres"])

            best_l = min(v["local_err"] for v in m_dict.values())
            best_s = min(v["stress_1"] for v in m_dict.values())
            best_k = max(v["knn_pres"] for v in m_dict.values())
            if abs(m["local_err"] - best_l) < 1e-6:
                wins_lerr += 1
            if abs(m["stress_1"] - best_s) < 1e-6:
                wins_stress += 1
            if abs(m["knn_pres"] - best_k) < 1e-6:
                wins_knn += 1

        agg[m_id] = {
            "id": m_id,
            "title": m_title,
            "desc": m_desc,
            "mean_time_ms": float(np.mean(times_all)) if times_all else 0.0,
            "median_time_ms": float(np.median(times_all)) if times_all else 0.0,
            "mean_lerr": float(np.mean(lerr_all)) if lerr_all else 0.0,
            "mean_lerr_dom": float(np.mean(lerr_dom)) if lerr_dom else 0.0,
            "mean_lerr_cmp": float(np.mean(lerr_cmp)) if lerr_cmp else 0.0,
            "mean_stress": float(np.mean(stress_all)) if stress_all else 0.0,
            "mean_stress_dom": float(np.mean(stress_dom)) if stress_dom else 0.0,
            "mean_stress_cmp": float(np.mean(stress_cmp)) if stress_cmp else 0.0,
            "mean_knn": float(np.mean(knn_all)) if knn_all else 0.0,
            "mean_knn_dom": float(np.mean(knn_dom)) if knn_dom else 0.0,
            "mean_knn_cmp": float(np.mean(knn_cmp)) if knn_cmp else 0.0,
            "mean_conflicts": float(np.mean(conf_all)) if conf_all else 0.0,
            "wins_lerr": wins_lerr,
            "wins_stress": wins_stress,
            "wins_knn": wins_knn,
        }

    evol_time = max(1.0, agg.get("evolution", {}).get("mean_time_ms", 1.0))
    for m_id in agg:
        agg[m_id]["speedup"] = evol_time / max(0.5, agg[m_id]["mean_time_ms"])
    return agg


def plot_summary_charts(agg, out_path):
    """Generate a 4-panel quantitative benchmark chart comparing all 16 views/algorithms."""
    m_ids = [m[0] for m in ALL_METHODS if m[0] in agg]
    labels = [agg[m]["title"].split("(")[0].strip() if m not in ("2dprots_mem15", "2dprots_mem30", "3d_pca", "2dprots_orig") else agg[m]["title"] for m in m_ids]

    lerr_dom = [agg[m]["mean_lerr_dom"] for m in m_ids]
    lerr_cmp = [agg[m]["mean_lerr_cmp"] for m in m_ids]
    stress_all = [agg[m]["mean_stress"] for m in m_ids]
    knn_all = [agg[m]["mean_knn"] for m in m_ids]
    times_ms = [max(0.5, agg[m]["mean_time_ms"]) for m in m_ids]

    fig, axes = plt.subplots(2, 2, figsize=(16.5, 12.0), dpi=160)
    fig.patch.set_facecolor("#ffffff")
    y_pos = np.arange(len(m_ids))

    # 1. Local Weighted Error (<40 Å) - Domains vs Complexes
    ax = axes[0, 0]
    h = 0.38
    ax.barh(y_pos - h / 2, lerr_dom, height=h, color="#2563eb", label="Jednotlivé domény")
    ax.barh(y_pos + h / 2, lerr_cmp, height=h, color="#f59e0b", label="Celé komplexy (multi-3D)")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8.8, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlabel("Průměrná lokální chyba vzdáleností <40 Å [Å] (nižší = lepší)", fontsize=10)
    ax.set_title("A. Lokální přesnost (<40 Å): Domény vs. Komplexy", fontsize=11.5, fontweight="bold")
    ax.grid(axis="x", linestyle="--", alpha=0.4)
    ax.legend(loc="lower right", fontsize=9)
    for i, (vd, vc) in enumerate(zip(lerr_dom, lerr_cmp)):
        ax.text(vd + 0.08, i - h / 2, f"{vd:.2f} Å", va="center", fontsize=7.5, color="#1e3a8a")
        ax.text(vc + 0.08, i + h / 2, f"{vc:.2f} Å", va="center", fontsize=7.5, color="#92400e")
    ax.set_xlim(0, max(max(lerr_dom), max(lerr_cmp)) * 1.22)

    # 2. Global Metric Stress-1
    ax = axes[0, 1]
    colors_s = ["#10b981" if s == min(stress_all) else ("#8b5cf6" if m.startswith("2dprots") or m == "3d_pca" else "#3b82f6") for s, m in zip(stress_all, m_ids)]
    ax.barh(y_pos, stress_all, color=colors_s, height=0.64)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8.8, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlabel("Kruskalův globální Stress-1 (nižší = lepší)", fontsize=10)
    ax.set_title("B. Globální metrické zkreslení tvaru (Stress-1)", fontsize=11.5, fontweight="bold")
    ax.grid(axis="x", linestyle="--", alpha=0.4)
    for i, v in enumerate(stress_all):
        ax.text(v + 0.005, i, f"{v:.3f}", va="center", fontsize=8.0)
    ax.set_xlim(0, max(stress_all) * 1.18)

    # 3. k-NN Neighborhood Preservation (%)
    ax = axes[1, 0]
    colors_k = ["#059669" if k == max(knn_all) else ("#8b5cf6" if m.startswith("2dprots") or m == "3d_pca" else "#6366f1") for k, m in zip(knn_all, m_ids)]
    ax.barh(y_pos, knn_all, color=colors_k, height=0.64)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8.8, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlabel("Zachování 3D k-NN sousedství v 2D [%] (vyšší = lepší)", fontsize=10)
    ax.set_title("C. Zachování lokální topologie a kontaktů (k-NN, k=4)", fontsize=11.5, fontweight="bold")
    ax.grid(axis="x", linestyle="--", alpha=0.4)
    for i, v in enumerate(knn_all):
        ax.text(v + 0.6, i, f"{v:.1f}%", va="center", fontsize=8.0)
    ax.set_xlim(0, max(knn_all) * 1.16)

    # 4. Execution Time (ms, log scale)
    ax = axes[1, 1]
    colors_t = ["#dc2626" if m == "evolution" else ("#8b5cf6" if m.startswith("2dprots") or m == "3d_pca" else "#0ea5e9") for m in m_ids]
    ax.barh(y_pos, times_ms, color=colors_t, height=0.64)
    ax.set_xscale("log")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8.8, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlabel("Průměrný čas výpočtu [ms, log měřítko] (nižší = rychlejší)", fontsize=10)
    ax.set_title("D. Výpočetní rychlost metod (průměr na protein)", fontsize=11.5, fontweight="bold")
    ax.grid(axis="x", linestyle="--", alpha=0.4)
    for i, (v, m) in enumerate(zip(times_ms, m_ids)):
        sp = agg[m]["speedup"]
        txt = f"{v:.0f} ms ({sp:.0f}×)" if m != "evolution" else f"{v:.0f} ms (1× ref)"
        ax.text(v * 1.12, i, txt, va="center", fontsize=8.0)
    ax.set_xlim(0.3, max(times_ms) * 7.0)

    plt.suptitle(
        "Srovnání 3D Principálních os, 2DProts (Originál / 15° / 30° s membránou) a 12 DR algoritmů",
        fontsize=13.5,
        fontweight="bold",
        y=0.99,
    )
    plt.tight_layout(rect=[0, 0, 1, 0.965])
    plt.savefig(out_path, dpi=160)
    plt.close(fig)


def _sse_rgb(color_code):
    """Convert 2DProts rainbow/spectrum string or CSS color to matplotlib RGB."""
    if not color_code or not isinstance(color_code, str):
        return (0.55, 0.58, 0.62)
    c = color_code.strip()
    if c.startswith("s") and c[1:].isdigit():
        val = float(c[1:]) / 900.0
        cmap = matplotlib.colormaps["turbo"]
        return cmap(np.clip(val, 0.05, 0.95))[:3]
    named = {
        "red": (0.86, 0.20, 0.18),
        "blue": (0.15, 0.55, 0.82),
        "green": (0.16, 0.63, 0.60),
        "grey": (0.55, 0.58, 0.62),
        "gray": (0.55, 0.58, 0.62),
        "orange": (0.90, 0.45, 0.10),
        "purple": (0.55, 0.30, 0.75),
        "cyan": (0.10, 0.70, 0.80),
        "magenta": (0.80, 0.25, 0.65),
        "yellow": (0.85, 0.65, 0.12),
        "brown": (0.55, 0.35, 0.20),
    }
    if c.lower() in named:
        return named[c.lower()]
    try:
        return matplotlib.colors.to_rgb(c)
    except Exception:
        return (0.35, 0.55, 0.80)


def render_comparison_panel_png(tid, title, results_entry, out_path):
    """
    Render a 4x4 (16-panel) matplotlib figure showing:
    Row 1: 3D Principal Axes (`3d_pca`), `2dprots_orig`, `2dprots_mem15`, `2dprots_mem30`
    Rows 2-4: All 12 dimensionality reduction algorithms (`evolution` .. `umap`)
    """
    fig, axes = plt.subplots(4, 4, figsize=(16.5, 15.0), dpi=140)
    fig.patch.set_facecolor("#ffffff")

    for idx, (m_id, m_title, _) in enumerate(ALL_METHODS):
        ax = axes[idx // 4, idx % 4]
        is_ref_row = idx < 4
        ax.set_facecolor("#f1f5f9" if is_ref_row else "#f8fafc")
        jpath = os.path.join(JSON_DIR, f"{tid}_{m_id}.json")
        if not os.path.exists(jpath):
            ax.set_visible(False)
            continue
        with open(jpath, "r", encoding="utf-8") as f:
            jdata = json.load(f)

        lay = jdata.get("layout", {})
        ang = jdata.get("angles", {})
        sz = jdata.get("size", {})
        ht = jdata.get("h_type", {})
        col = jdata.get("color", {})
        mem = jdata.get("membrane", {})

        shown = [s for s in lay if ht.get(s, 3) in (1, 2)]

        # If membrane band is active (in 2dprots_mem15 / 2dprots_mem30), draw horizontal membrane slab
        if isinstance(mem, dict) and mem.get("is_transmembrane", False) and m_id in ("2dprots_mem15", "2dprots_mem30"):
            yc = float(mem.get("y_center_2d", 0.0))
            hh = float(mem.get("band_half_height_2d", 3.0))
            ax.axhspan(yc - hh, yc + hh, color="#fde68a", alpha=0.38, zorder=0)
            ax.axhline(yc + hh, color="#d97706", linestyle="--", linewidth=0.9, alpha=0.7, zorder=1)
            ax.axhline(yc - hh, color="#d97706", linestyle="--", linewidth=0.9, alpha=0.7, zorder=1)

        if len(shown) >= 2:
            xs = [lay[s][0] for s in shown]
            ys = [lay[s][1] for s in shown]
            ax.plot(xs, ys, color="#94a3b8", linewidth=0.85, alpha=0.48, zorder=1)

        all_x, all_y = [], []
        for s in shown:
            x, y = lay[s]
            a_deg = float(ang.get(s, 0.0))
            a_rad = math.radians(a_deg)
            length = max(1.2, float(sz.get(s, 5)) * 0.22)
            dx = length * math.sin(a_rad)
            dy = length * math.cos(a_rad)
            x0, y0 = x - dx, y - dy
            x1, y1 = x + dx, y + dy
            all_x.extend([x0, x1])
            all_y.extend([y0, y1])
            rgb = _sse_rgb(col.get(s, "blue"))
            is_helix = s.startswith("H") or s.startswith("_H")
            if is_helix:
                ax.plot([x0, x1], [y0, y1], color=rgb, linewidth=4.6, solid_capstyle="round", zorder=3)
                ax.plot([x0, x1], [y0, y1], color="#1e293b", linewidth=5.6, solid_capstyle="round", zorder=2, alpha=0.25)
            else:
                ax.annotate(
                    "",
                    xy=(x1, y1),
                    xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>,head_length=0.45,head_width=0.35", lw=2.3, color=rgb),
                    zorder=4,
                )

        # Draw bound ligands & cofactors (if present in layout JSON)
        ligs = jdata.get("ligands", {})
        if isinstance(ligs, dict) and ligs:
            import matplotlib.patheffects as patheffects
            import utils_benchmark_ligands
            for lk, lv in ligs.items():
                if not isinstance(lv, dict) or not lv.get("layout"):
                    continue
                lx, ly = float(lv["layout"][0]), float(lv["layout"][1])
                all_x.append(lx)
                all_y.append(ly)
                ecmp = str(lv.get("ecmp") or lk)
                natoms = int(lv.get("atoms", 1))
                f_col, s_col = utils_benchmark_ligands.get_ligand_colors(ecmp)
                if lv.get("color"):
                    f_col = lv["color"]
                if lv.get("ocolor"):
                    s_col = lv["ocolor"]
                ms = 5.6 if natoms < 10 else min(9.5, 4.6 + math.sqrt(natoms) * 0.72)
                ax.plot(lx, ly, marker="o", markersize=ms + 2.2, color="#ffffff", alpha=0.90, zorder=5)
                ax.plot(
                    lx,
                    ly,
                    marker="o",
                    markersize=ms,
                    markerfacecolor=f_col,
                    markeredgecolor=s_col,
                    markeredgewidth=1.2,
                    zorder=6,
                )
                txt_el = ax.annotate(
                    ecmp,
                    xy=(lx, ly),
                    xytext=(0, ms * 0.55 + 1.2),
                    textcoords="offset points",
                    fontsize=5.8,
                    fontweight="bold",
                    color="#0f172a",
                    ha="center",
                    va="bottom",
                    zorder=7,
                )
                txt_el.set_path_effects([patheffects.withStroke(linewidth=2.0, foreground="#ffffff")])

        if all_x and all_y:
            cx = 0.5 * (min(all_x) + max(all_x))
            cy = 0.5 * (min(all_y) + max(all_y))
            span = max(max(all_x) - min(all_x), max(all_y) - min(all_y), 10.0) * 0.62
            ax.set_xlim(cx - span, cx + span)
            ax.set_ylim(cy - span, cy + span)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_edgecolor("#6366f1" if is_ref_row else "#cbd5e1")
            spine.set_linewidth(1.4 if is_ref_row else 1.0)

        m_info = results_entry.get("methods", {}).get(m_id, {})
        short_name = m_title if is_ref_row else m_title.split("(")[0].strip()
        if m_id == "3d_pca" and "pca_spans_3d" in jdata:
            sp3 = jdata["pca_spans_3d"]
            sub_line = f"PC1×2×3: {sp3[0]:.0f}×{sp3[1]:.0f}×{sp3[2]:.0f}Å | err={m_info.get('local_err', 0):.2f}Å"
        elif m_id in ("2dprots_mem15", "2dprots_mem30"):
            tm_tag = "TM pás" if m_info.get("is_tm") else "Glob."
            sub_line = f"[{tm_tag}] err={m_info.get('local_err', 0):.2f}Å | str={m_info.get('stress_1', 0):.2f} | kNN={m_info.get('knn_pres', 0):.0f}%"
        else:
            sub_line = (
                f"err={m_info.get('local_err', 0):.2f}Å | str={m_info.get('stress_1', 0):.2f} | "
                f"kNN={m_info.get('knn_pres', 0):.0f}% | {m_info.get('time_ms', 0):.0f}ms"
            )
        ax.set_title(f"{short_name}\n{sub_line}", fontsize=8.2, fontweight="bold", pad=3.5)

    plt.suptitle(title, fontsize=13.5, fontweight="bold", y=0.992)
    plt.tight_layout(rect=[0, 0, 1, 0.968])
    plt.savefig(out_path, dpi=140)
    plt.close(fig)


def generate_html_gallery(results, agg, html_path):
    """Generate the self-contained interactive HTML gallery (`benchmark_dimred_results/index.html`)."""
    categories = OrderedDict()
    for t in benchmark_dimred.BENCHMARK_TARGETS:
        tid = t["id"]
        if tid not in results:
            continue
        cat = t["category"]
        categories.setdefault(cat, []).append((tid, results[tid]))

    agg_rows_html = []
    best_lerr = min(a["mean_lerr"] for a in agg.values())
    best_lerr_dom = min(a["mean_lerr_dom"] for a in agg.values())
    best_lerr_cmp = min(a["mean_lerr_cmp"] for a in agg.values())
    best_stress = min(a["mean_stress"] for a in agg.values())
    best_knn = max(a["mean_knn"] for a in agg.values())
    best_time = min(a["mean_time_ms"] for a in agg.values())

    for m_id, m_title, m_desc in ALL_METHODS:
        if m_id not in agg:
            continue
        a = agg[m_id]
        c_lerr = "best-cell" if abs(a["mean_lerr"] - best_lerr) < 0.02 else ""
        c_ldom = "best-cell" if abs(a["mean_lerr_dom"] - best_lerr_dom) < 0.02 else ""
        c_lcmp = "best-cell" if abs(a["mean_lerr_cmp"] - best_lerr_cmp) < 0.02 else ""
        c_str = "best-cell" if abs(a["mean_stress"] - best_stress) < 0.005 else ""
        c_knn = "best-cell" if abs(a["mean_knn"] - best_knn) < 0.3 else ""
        c_tim = "best-cell" if abs(a["mean_time_ms"] - best_time) < 5.0 else ""
        row_cls = "ref-row" if m_id in ("3d_pca", "2dprots_orig", "2dprots_mem15", "2dprots_mem30") else ""

        agg_rows_html.append(f"""
        <tr class="{row_cls}">
            <td><strong>{m_title}</strong><div class="sub">{m_desc}</div></td>
            <td class="{c_tim}">{a['mean_time_ms']:.1f} ms <span class="badge">{a['speedup']:.1f}×</span></td>
            <td class="{c_lerr}"><strong>{a['mean_lerr']:.2f} Å</strong></td>
            <td class="{c_ldom}">{a['mean_lerr_dom']:.2f} Å</td>
            <td class="{c_lcmp}">{a['mean_lerr_cmp']:.2f} Å</td>
            <td class="{c_str}">{a['mean_stress']:.3f}</td>
            <td class="{c_knn}">{a['mean_knn']:.1f}%</td>
            <td>{a['mean_conflicts']:.2f}</td>
            <td>{a['wins_lerr']} / {a['wins_stress']} / {a['wins_knn']}</td>
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
            min_lerr = min(v["local_err"] for v in m_dict.values()) if m_dict else 0
            min_stress = min(v["stress_1"] for v in m_dict.values()) if m_dict else 0
            max_knn = max(v["knn_pres"] for v in m_dict.values()) if m_dict else 0

            grid_items = []
            for m_id, m_title, _ in ALL_METHODS:
                if m_id not in m_dict:
                    continue
                m = m_dict[m_id]
                svg_rel = f"svgs/{tid}_{m_id}.svg"
                is_best_l = abs(m["local_err"] - min_lerr) < 1e-5
                is_ref = m_id in ("3d_pca", "2dprots_orig", "2dprots_mem15", "2dprots_mem30")
                border_cls = "method-card best-card" if is_best_l else ("method-card ref-card" if is_ref else "method-card")
                short_title = m_title if is_ref else m_title.split("(")[0].strip()
                extra_badge = ""
                if m_id == "3d_pca" and "pca_spans" in m:
                    sp = m["pca_spans"]
                    extra_badge = f'<span class="method-time">{sp[0]:.0f}×{sp[1]:.0f}×{sp[2]:.0f}Å</span>'
                elif m_id in ("2dprots_mem15", "2dprots_mem30"):
                    extra_badge = f'<span class="method-time">{"TM" if m.get("is_tm") else "Glob."}</span>'
                else:
                    extra_badge = f'<span class="method-time">{m["time_ms"]:.0f} ms</span>'

                grid_items.append(f"""
                <div class="{border_cls}" data-method="{m_id}">
                    <div class="method-header">
                        <span class="method-name">{short_title}</span>
                        {extra_badge}
                    </div>
                    <a href="{svg_rel}" target="_blank" title="Klikněte pro otevření plného SVG ({tid} - {m_title})">
                        <img src="{svg_rel}" alt="{tid} {m_id}" loading="lazy" />
                    </a>
                    <div class="method-metrics">
                        <span class="{'m-win' if is_best_l else ''}" title="Lokální chyba <40 Å">Err: {m['local_err']:.2f}Å</span>
                        <span class="{'m-win' if abs(m['stress_1'] - min_stress) < 1e-5 else ''}" title="Kruskalův Stress-1">Str: {m['stress_1']:.2f}</span>
                        <span class="{'m-win' if abs(m['knn_pres'] - max_knn) < 1e-5 else ''}" title="Zachování 4-NN">kNN: {m['knn_pres']:.0f}%</span>
                    </div>
                </div>
                """)

            cards_html.append(f"""
            <div class="target-block" id="target-{tid}" data-kind="{t['kind']}">
                <div class="target-header">
                    <div>
                        <h3>{t['name']} <code class="tid-code">{tid}</code> {kind_badge}</h3>
                        <p class="target-desc">{t['desc']} — <strong>{rdata['n_sses']} SSE</strong> ({rdata['n_clusters']} klastrů)</p>
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
    <title>2DProts – 3D Principální osy, 2DProts (Orig / 15° / 30°) a 12 DR algoritmů</title>
    <style>
        :root {{
            --bg: #f8fafc;
            --card: #ffffff;
            --text: #0f172a;
            --muted: #475569;
            --border: #e2e8f0;
            --accent: #2563eb;
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
            background: linear-gradient(135deg, #0f172a 0%, #1e3a8a 100%);
            color: #fff;
            padding: 28px 32px;
            border-radius: 14px;
            margin-bottom: 28px;
            box-shadow: 0 10px 25px rgba(15, 23, 42, 0.12);
        }}
        header h1 {{ margin: 0 0 8px 0; font-size: 1.7rem; }}
        header p {{ margin: 0; color: #cbd5e1; font-size: 0.98rem; }}
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
            font-size: 0.88rem;
        }}
        table.summary-table th, table.summary-table td {{
            padding: 9px 12px;
            border-bottom: 1px solid var(--border);
            text-align: left;
        }}
        table.summary-table th {{
            background: #f1f5f9;
            font-weight: 700;
            color: #1e293b;
        }}
        tr.ref-row {{ background: #f5f3ff; }}
        .sub {{ font-size: 0.76rem; color: var(--muted); }}
        .best-cell {{ background: var(--win); color: var(--win-text); font-weight: 700; }}
        .badge {{
            display: inline-block;
            background: #e0f2fe;
            color: #0369a1;
            font-size: 0.74rem;
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
            grid-template-columns: repeat(auto-fill, minmax(210px, 1fr));
            gap: 11px;
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
        .method-card.ref-card {{
            border: 1.5px solid #8b5cf6;
            background: #faf5ff;
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
            font-size: 0.78rem;
            margin-bottom: 6px;
        }}
        .method-name {{ font-weight: 700; color: #1e293b; }}
        .method-time {{ color: #64748b; font-size: 0.72rem; }}
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
            font-size: 0.73rem;
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
        <h1>2DProts: 3D Principální osy, 2DProts (Orig / 15° / 30° s membránou) a 12 algoritmů redukce dimenze</h1>
        <p>Komplexní benchmark nad {len(results)} strukturami (cytochromy P450, enzymy EC 1–7, GPCR, SLC a ABC transportéry, β-barely, 14-3-3, hemoglobin, KcsA, akvaporin, F1-ATPáza a GroEL – jako jednotlivé domény i celé komplexy)</p>
    </header>

    <div class="summary-card">
        <h2 style="margin-top:0;">Celkové kvantitativní srovnání všech 16 pohledů a algoritmů (průměr přes {len(results)} struktur)</h2>
        <table class="summary-table">
            <thead>
                <tr>
                    <th>Pohled / Algoritmus</th>
                    <th>Průměrný čas (zrychlení)</th>
                    <th>Lokální chyba &lt;40 Å (Vše)</th>
                    <th>Lokální chyba (Domény)</th>
                    <th>Lokální chyba (Komplexy)</th>
                    <th>Globální Stress-1</th>
                    <th>Zachování 4-NN</th>
                    <th>Kolize SSE</th>
                    <th>Počet 1. míst (Err / Str / kNN)</th>
                </tr>
            </thead>
            <tbody>
                {''.join(agg_rows_html)}
            </tbody>
        </table>
    </div>

    <div class="filter-bar">
        <strong>Filtr zobrazení:</strong>
        <button class="filter-btn active" onclick="filterKind('all', this)">Vše ({len(results)})</button>
        <button class="filter-btn" onclick="filterKind('domain', this)">Pouze jednotlivé domény</button>
        <button class="filter-btn" onclick="filterKind('complex', this)">Pouze celé komplexy (Multi-3D)</button>
    </div>

    {''.join(sections_html)}

    <script>
        function filterKind(kind, btn) {{
            document.querySelectorAll('.filter-btn').forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
            document.querySelectorAll('.target-block').forEach(el => {{
                if (kind === 'all' || el.getAttribute('data-kind') === kind) {{
                    el.style.display = '';
                }} else {{
                    el.style.display = 'none';
                }}
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

    agg = compute_aggregate_stats(results)
    with open(os.path.join(OUT_DIR, "aggregate_stats.json"), "w", encoding="utf-8") as f:
        json.dump(agg, f, indent=2, ensure_ascii=False)

    chart_path = os.path.join(FIG_DIR, "benchmark_summary_charts.png")
    plot_summary_charts(agg, chart_path)
    print(f"Saved summary chart to {chart_path}")

    featured = [
        ("1tqnA00", "Cytochrom P450 3A4 (1tqnA00 – globulární enzym, 37 SSE)"),
        ("1adsA00", "EC 1: Aldózareduktáza – kanonický (β/α)8 TIM barel (1adsA00)"),
        ("1f88_dimer", "GPCR: Bovinní rodopsin – dimerní 7TM komplex (1f88_dimer)"),
        ("1pv6_complex", "SLC / MFS: Laktózová permeáza LacY – celý 12TM komplex N+C (1pv6_complex)"),
        ("2hyd_dimer", "ABC transportér Sav1866 – celý homodimerní komplex (2hyd_dimer, EC 7)"),
        ("1emaA00", "Beta-barel (solubilní): Zelený fluorescenční protein GFP – 11-vláknový barel (1emaA00)"),
        ("1rbpA00", "Beta-barel (solubilní lipokalin): Retinol-vázající protein RBP4 s retinolem RTL (1rbpA00)"),
        ("1cbsA00", "Beta-barel (solubilní iLBP): CRABP-II s kyselinou all-trans-retinovou REA (1cbsA00)"),
        ("1stpA00", "Beta-barel (solubilní): Streptavidin – monomer s biotinem BTN (1stpA00)"),
        ("1stp_tetramer", "Beta-barel (solubilní komplex): Streptavidin – homotetramer 4×8 vláken + 4×BTN (1stp_tetramer)"),
        ("2omf_trimer", "Beta-barel (TM komplex): Porin OmpF – celý biologický trimer 3×16 β-vláken (2omf_trimer)"),
        ("1qjb_dimer", "14-3-3ζ protein – homodimerní komplex ve tvaru W (1qjb_dimer)"),
        ("2hhb_tetramer", "Hemoglobin – celý kvartérní tetramer α2β2 (2hhb_tetramer)"),
        ("1bl8_tetramer", "Draslíkový kanál KcsA – celý biologický tetramer (1bl8_tetramer)"),
        ("1j4n_tetramer", "Akvaporin AQP1 – celý biologický tetramer 4× kanál (1j4n_tetramer)"),
        ("1bmf_complex", "Extra velký komplex: F1-ATPáza α3β3γ motor – 19 domén (1bmf_complex)"),
        ("1aon_ring", "Extra velký komplex: Chaperonin GroEL – heptamerní kruh 7 podjednotek (1aon_ring)"),
        ("2zxeA02", "Na+/K+-ATPáza: 10TM transmembránová doména α-podjednotky s K+ a CLR (2zxeA02)"),
        ("2zxe_complex", "Na+/K+-ATPáza: Kompletní αβ-heterokomplex – 6 domén s K+, MgF4, Mg2+, CLR (2zxe_complex)"),
        ("1c3wA00", "Protonová pumpa: Bakteriorhodopsin – 7TM pumpa s retinalem RET a lipidy (1c3wA00)"),
        ("2xnd_complex", "Kompletní protonová pumpa: F1c8-ATP syntáza – α3β3γδε + TM rotor c8, 31 domén (2xnd_complex)"),
        ("1zoy_complex", "Dýchací řetězec – Komplex II: Sukcinátdehydrogenáza s FAD, 3×FeS, UQ1, HEM (1zoy_complex)"),
        ("1occA00", "Dýchací řetězec – Komplex IV: Cytochrom c oxidáza – podjednotka I s 2×HEA, CuB, Mg (1occA00)"),
        ("1occ_complex", "Dýchací řetězec – Komplex IV: Cytochrom c oxidáza – celý 13-podjednotkový komplex (1occ_complex)"),
    ]
    only_panels = set(sys.argv[1:]) if len(sys.argv) > 1 else None
    for tid, title in featured:
        if only_panels and tid not in only_panels:
            continue
        if tid in results:
            p_out = os.path.join(FIG_DIR, f"panel_{tid}.png")
            render_comparison_panel_png(tid, title, results[tid], p_out)
            print(f"Saved panel {p_out}")

    html_path = os.path.join(OUT_DIR, "index.html")
    generate_html_gallery(results, agg, html_path)
    print(f"Saved interactive HTML gallery to {html_path}")


if __name__ == "__main__":
    main()
