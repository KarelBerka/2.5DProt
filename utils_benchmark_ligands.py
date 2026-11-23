#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""
Ligand extraction (from 3D mmCIF crystal structures and 2DProts domain JSONs)
and 3D-to-2D binding-pocket projection for all 16 benchmark views.
"""

import os
import math
import gzip
import json
import numpy as np
import gemmi
import utils_membrane
import utils_sses_types

# Repetitive crystallization detergents, buffer ions, or bulk lipids capped at max 4 copies per target
BULK_DETERGENT_OR_BUFFER_CODES = {
    "LI1", "C8E", "BNG", "PC1", "LDM", "SO4", "GOL", "BU1", "12P", "ACM", "CAU", "PLM",
}

# Distinct, high-contrast chemical/functional color palette for ligands & cofactors
LIGAND_PALETTE = {
    # Hemes & porphyrins
    "HEM": ("#ef4444", "#991b1b"),
    "HEA": ("#dc2626", "#7f1d1d"),
    # Flavins & quinones & retinoids (retinal, retinol, retinoic acid) & biotin
    "FAD": ("#f59e0b", "#92400e"),
    "FMN": ("#f59e0b", "#92400e"),
    "UQ1": ("#d97706", "#78350f"),
    "RET": ("#f97316", "#9a3412"),
    "RTL": ("#f97316", "#9a3412"),
    "REA": ("#fb923c", "#9a3412"),
    "BTN": ("#ec4899", "#831843"),
    "BLV": ("#06b6d4", "#155e75"),
    # Iron-sulfur clusters
    "FES": ("#ea580c", "#7c2d12"),
    "SF4": ("#c2410c", "#7c2d12"),
    "F3S": ("#b45309", "#78350f"),
    # Nucleotides & nicotinamides & phosphate analogs
    "ATP": ("#3b82f6", "#1d4ed8"),
    "ADP": ("#60a5fa", "#1e40af"),
    "ANP": ("#2563eb", "#1e3a8a"),
    "NAD": ("#0ea5e9", "#0369a1"),
    "NAP": ("#0284c7", "#075985"),
    "MF4": ("#a855f7", "#6b21a8"),
    "PO4": ("#c084fc", "#6b21a8"),
    # Catalytic & channel ions
    "K":   ("#8b5cf6", "#4c1d95"),
    "NA":  ("#a78bfa", "#5b21b6"),
    "MG":  ("#10b981", "#065f46"),
    "MN":  ("#14b8a6", "#115e59"),
    "CU":  ("#06b6d4", "#155e75"),
    "ZN":  ("#64748b", "#334155"),
    "HG":  ("#94a3b8", "#475569"),
    # Inhibitors / substrates
    "CXT": ("#ec4899", "#9d174d"),
    "BHF": ("#f43f5e", "#9f1239"),
    "OXM": ("#fb7185", "#9f1239"),
    "SCN": ("#f472b6", "#831843"),
    # Membrane lipids, sterols & glycans
    "CLR": ("#eab308", "#854d0e"),
    "CDL": ("#ca8a04", "#713f12"),
    "EPH": ("#d97706", "#78350f"),
    "SQU": ("#facc15", "#854d0e"),
    "NAG": ("#22c55e", "#166534"),
    "NDG": ("#4ade80", "#166534"),
}


def get_ligand_colors(ecmp):
    """Return (fill_hex, stroke_hex) for a chemical component code `ecmp`."""
    code = (ecmp or "").upper().strip()
    if code in LIGAND_PALETTE:
        return LIGAND_PALETTE[code]
    if code in ("LI1", "PC1", "LDM", "PLM", "C8E", "BNG", "12P", "BU1"):
        return ("#fde047", "#854d0e")
    return ("#f472b6", "#9d174d")


def _index_to_letters(idx):
    """Convert 0-based index to pure uppercase letter key ('B', 'C', ..., 'Z', 'AA', 'AB', ...)."""
    idx = idx + 1  # Start at 'B' like standard 2DProts ligand keys
    chars = []
    while True:
        chars.append(chr(ord("A") + (idx % 26)))
        idx = idx // 26 - 1
        if idx < 0:
            break
    return "".join(reversed(chars))


def _min_dist_point_to_sse_3d(pt_3d, vec_pair):
    """Return (min_dist_3d, t_frac, closest_pt_3d) from `pt_3d` to 3D segment `vec_pair`."""
    v0 = np.asarray(vec_pair[0], dtype=float)
    v1 = np.asarray(vec_pair[1], dtype=float)
    dv = v1 - v0
    l2 = float(np.dot(dv, dv))
    if l2 < 1e-8:
        return float(np.linalg.norm(pt_3d - v0)), 0.5, v0
    t = float(np.clip(np.dot(pt_3d - v0, dv) / l2, 0.0, 1.0))
    q = v0 + t * dv
    return float(np.linalg.norm(pt_3d - q)), t, q


def extract_target_ligands_3d(target, ann, cache_dir):
    """
    Extract all biologically relevant ligands for `target` in the exact 3D coordinate
    frame of `ann.vectors`.
    Returns an OrderedDict mapping pure-letter ligand keys ('B', 'C', ...) to ligand dicts.
    """
    import benchmark_dimred

    pdb_id = target["pdb"]
    use_au = target.get("use_au", True) if target["kind"] == "complex" else True
    model = benchmark_dimred.fetch_pdb_structure(pdb_id, use_au=use_au)

    # Read chemical component descriptions from cached mmCIF or domain JSONs
    comp_desc = {}
    cif_path = os.path.join(cache_dir, f"pdb_{pdb_id}.cif.gz")
    if os.path.exists(cif_path):
        try:
            with gzip.open(cif_path, "rb") as f:
                cif_text = f.read().decode("utf-8", errors="replace")
            doc = gemmi.cif.read_string(cif_text)
            block = doc.sole_block()
            for cid, cnm in zip(block.find_values("_chem_comp.id"), block.find_values("_chem_comp.name")):
                comp_desc[cid.strip().strip("'\"")] = cnm.strip().strip("'\"")
        except Exception:
            pass

    for _, _, dom_code in target["domains"]:
        try:
            ddata = benchmark_dimred.fetch_domain_json(dom_code)
            for lval in ddata.get("ligands", {}).values():
                if lval.get("ecmp") and lval.get("edsc"):
                    comp_desc[lval["ecmp"]] = lval["edsc"]
        except Exception:
            pass

    # For single-domain targets, compute the Kabsch transform from crystal chain onto `ann.vectors`
    sup_transform = None
    if target["kind"] == "domain":
        ch_id, _, dom_code = target["domains"][0]
        poly_chains = benchmark_dimred.get_polymer_chains(model)
        if ch_id in poly_chains:
            ca_map = benchmark_dimred.extract_ca_map(poly_chains[ch_id])
            ddata = benchmark_dimred.fetch_domain_json(dom_code)
            hr = ddata.get("helices_residues", {})
            l3d = ddata.get("3dlayout", {})
            if ca_map and hr:
                all_res = [r for pair in hr.values() for r in pair]
                res_offset = 0
                min_hr = min(all_res)
                min_ca = min(ca_map.keys())
                if min_hr <= 10 and min_ca > 15 and (min_ca - 1 + max(all_res)) in ca_map:
                    res_offset = min_ca - 1
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
                    sup_transform = sup.transform

    # Collect candidate non-polymer, non-water residues from the crystal structure
    raw_candidates = []
    target_chains = {ch_id for ch_id, _, _ in target["domains"]}
    for ch in model:
        for res in ch:
            if res.entity_type == gemmi.EntityType.Polymer or res.entity_type == gemmi.EntityType.Water:
                continue
            ecmp = res.name.strip()
            if ecmp in ("HOH", "DOD", "WAT"):
                continue
            atom_pts = []
            for a in res:
                pos = a.pos
                if sup_transform is not None:
                    pos = sup_transform.apply(pos)
                atom_pts.append([float(pos.x), float(pos.y), float(pos.z)])
            if not atom_pts:
                continue
            pts_arr = np.array(atom_pts, dtype=float)
            c3d = pts_arr.mean(axis=0)

            # Minimum distance from any ligand atom (or centroid) to any SSE in `ann.vectors`
            d_min = 1e9
            for s, vec_pair in ann.vectors.items():
                d_c, _, _ = _min_dist_point_to_sse_3d(c3d, vec_pair)
                if d_c < d_min:
                    d_min = d_c
            if d_min <= 14.5:
                # Prefer ligands on target chains or very close to SSEs
                if target["kind"] == "domain" and ch.name not in target_chains and d_min > 9.5:
                    continue
                raw_candidates.append({
                    "ecmp": ecmp,
                    "edsc": comp_desc.get(ecmp, ecmp),
                    "atoms": len(atom_pts),
                    "e_type": "non-polymer",
                    "ept": None,
                    "auth_chain": ch.name,
                    "res": res.seqid.num,
                    "center_3d": [round(float(c3d[0]), 3), round(float(c3d[1]), 3), round(float(c3d[2]), 3)],
                    "d_min": float(d_min),
                })

    # Sort by distance to SSEs so the most tightly coordinated ligands come first
    raw_candidates.sort(key=lambda x: (x["ecmp"] in BULK_DETERGENT_OR_BUFFER_CODES, x["d_min"]))

    filtered = []
    bulk_counts = {}
    for cand in raw_candidates:
        code = cand["ecmp"]
        if code in BULK_DETERGENT_OR_BUFFER_CODES:
            if bulk_counts.get(code, 0) >= 3:
                continue
            bulk_counts[code] = bulk_counts.get(code, 0) + 1
        # Deduplicate symmetry-overlapping copies at the exact same 3D center (< 1.2 A)
        if any(
            math.dist(cand["center_3d"], prev["center_3d"]) < 1.2 and cand["ecmp"] == prev["ecmp"]
            for prev in filtered
        ):
            continue
        filtered.append(cand)

    out = {}
    for idx, item in enumerate(filtered):
        lkey = _index_to_letters(idx)
        fill_col, stroke_col = get_ligand_colors(item["ecmp"])
        out[lkey] = {
            "layout": [0.0, 0.0],
            "atoms": item["atoms"],
            "e_type": item["e_type"],
            "ept": item["ept"],
            "edsc": item["edsc"],
            "ecmp": item["ecmp"],
            "auth_chain": item["auth_chain"],
            "res": item["res"],
            "3d_layout": item["center_3d"],
            "color": fill_col,
            "ocolor": stroke_col,
        }
    return out


def project_ligands_to_2d(ligs_3d, ann, layout_2d, angles_2d, mem_info=None):
    """
    Map 3D ligands (`ligs_3d`) into the 2D coordinate space of (`layout_2d`, `angles_2d`).
    Preserves local binding-site topology relative to coordinating SSE segments in 2D
    and resolves ligand-SSE and ligand-ligand overlaps.
    """
    if not ligs_3d or not layout_2d:
        return {}

    valid_sses = [s for s in layout_2d if s in ann.vectors]
    if not valid_sses:
        return {}

    # Precompute 2D start/end segments for all SSEs in `layout_2d`
    sse_2d_segs = {}
    for s in valid_sses:
        cx, cy = float(layout_2d[s][0]), float(layout_2d[s][1])
        ang_rad = math.radians(float(angles_2d.get(s, 0.0)))
        # Half-length in 2D canvas units
        hl = float(ann.res_sizes.get(s, 8)) * 0.145
        v0_2d = np.array([cx - hl * math.sin(ang_rad), cy - hl * math.cos(ang_rad)], dtype=float)
        v1_2d = np.array([cx + hl * math.sin(ang_rad), cy + hl * math.cos(ang_rad)], dtype=float)
        sse_2d_segs[s] = (v0_2d, v1_2d, np.array([cx, cy], dtype=float))

    # Global 2D centroid of all SSEs
    all_centers_2d = np.array([sse_2d_segs[s][2] for s in valid_sses], dtype=float)
    global_center_2d = all_centers_2d.mean(axis=0)

    is_tm = bool(mem_info and mem_info.get("is_transmembrane", False))
    band_half = float(mem_info.get("band_half_height_2d", 3.0)) if is_tm else 3.0
    tm_set = set(mem_info.get("tm_helices", [])) if is_tm else set()

    out_ligs = {}
    anchors_2d = {}

    for lkey, ldata in ligs_3d.items():
        c3d = np.asarray(ldata["3d_layout"], dtype=float)
        ecmp = (ldata.get("ecmp") or "").upper().strip()

        # Find closest point on every SSE segment in 3D and corresponding fractional point in 2D
        sse_hits = []
        for s in valid_sses:
            d3d, t_frac, q3d = _min_dist_point_to_sse_3d(c3d, ann.vectors[s])
            v0_2d, v1_2d, _ = sse_2d_segs[s]
            q2d = (1.0 - t_frac) * v0_2d + t_frac * v1_2d
            sse_hits.append((d3d, s, t_frac, q3d, q2d))

        sse_hits.sort(key=lambda x: x[0])
        d_min = sse_hits[0][0]
        s_near = sse_hits[0][1]
        q2d_near = sse_hits[0][4]

        # Coordinating SSEs (up to 5 within local neighborhood)
        cutoff = max(d_min * 1.65, d_min + 6.5)
        coord_hits = [h for h in sse_hits[:5] if h[0] <= cutoff]
        if not coord_hits:
            coord_hits = [sse_hits[0]]

        weights = []
        for d3d, s, t_frac, q3d, q2d in coord_hits:
            w = math.exp(-((d3d - d_min) ** 2) / (2.0 * (3.2 ** 2))) / ((d3d + 1.2) ** 2)
            weights.append((w, d3d, s, q3d, q2d))
        w_sum = sum(w for w, _, _, _, _ in weights) or 1.0

        anchor_2d = sum(w * q2d for w, _, _, _, q2d in weights) / w_sum
        anchor_3d = sum(w * q3d for w, _, _, q3d, _ in weights) / w_sum

        # Check if enclosed inside a multi-SSE binding pocket
        r_3d_from_anchor = float(np.linalg.norm(c3d - anchor_3d))
        if len(coord_hits) >= 2 and r_3d_from_anchor < max(4.0, d_min * 0.85):
            # Inside binding pocket: place at weighted pocket center with slight bias toward nearest SSE point
            pos_2d = 0.78 * anchor_2d + 0.22 * q2d_near
        else:
            # Peripheral ligand: offset outward from nearest SSE point away from local centroid
            dir_vec = q2d_near - anchor_2d
            if np.linalg.norm(dir_vec) < 0.25:
                dir_vec = q2d_near - global_center_2d
            norm_d = float(np.linalg.norm(dir_vec))
            u_out = dir_vec / norm_d if norm_d > 1e-4 else np.array([1.0, 0.0])
            r_off = min(1.65, max(0.75, 0.16 * d_min))
            pos_2d = 0.65 * q2d_near + 0.35 * anchor_2d + u_out * r_off

        if is_tm:
            f_tm = sum(w for w, _, s, _, _ in weights if s in tm_set) / w_sum
            if ecmp in utils_membrane.PORE_ION_CODES and f_tm >= 0.45:
                pos_2d[1] = float(np.clip(pos_2d[1], -band_half * 0.75, band_half * 0.75))
            elif ecmp in utils_membrane.MEMBRANE_LIPID_CODES or ecmp in BULK_DETERGENT_OR_BUFFER_CODES:
                if f_tm >= 0.45:
                    pos_2d[1] = float(np.clip(pos_2d[1], -band_half + 0.35, band_half - 0.35))

        new_entry = dict(ldata)
        new_entry["layout"] = [float(pos_2d[0]), float(pos_2d[1])]
        out_ligs[lkey] = new_entry
        anchors_2d[lkey] = np.array([float(pos_2d[0]), float(pos_2d[1])], dtype=float)

    # Clear SSE segments using utils_membrane.deoverlap_ligands_from_sses
    utils_membrane.deoverlap_ligands_from_sses(
        out_ligs,
        {s: list(layout_2d[s]) for s in valid_sses},
        angles_2d,
        ann.res_sizes,
        mem_info if is_tm else None,
    )

    # Resolve pairwise ligand-ligand overlaps while keeping ligands tethered to their binding site
    lkeys = list(out_ligs.keys())
    if len(lkeys) >= 2:
        for _ in range(18):
            for i in range(len(lkeys)):
                k1 = lkeys[i]
                p1 = out_ligs[k1]["layout"]
                r1 = 0.32 if out_ligs[k1].get("atoms", 1) < 10 else min(0.68, math.sqrt(out_ligs[k1]["atoms"]) / 10.0 + 0.12)
                for j in range(i + 1, len(lkeys)):
                    k2 = lkeys[j]
                    p2 = out_ligs[k2]["layout"]
                    r2 = 0.32 if out_ligs[k2].get("atoms", 1) < 10 else min(0.68, math.sqrt(out_ligs[k2]["atoms"]) / 10.0 + 0.12)
                    min_sep = r1 + r2 + 0.65
                    dx = p1[0] - p2[0]
                    dy = p1[1] - p2[1]
                    dist = math.hypot(dx, dy)
                    if dist < min_sep:
                        if dist > 1e-4:
                            ux, uy = dx / dist, dy / dist
                        else:
                            ang = (i * 1.1 + j * 2.3)
                            ux, uy = math.cos(ang), math.sin(ang)
                        push = 0.52 * (min_sep - dist)
                        p1[0] += ux * push
                        p1[1] += uy * push
                        p2[0] -= ux * push
                        p2[1] -= uy * push
            # Gentle tether back toward pocket anchor
            for k in lkeys:
                p = out_ligs[k]["layout"]
                anc = anchors_2d[k]
                ddx = p[0] - float(anc[0])
                ddy = p[1] - float(anc[1])
                d_anc = math.hypot(ddx, ddy)
                if d_anc > 2.2:
                    p[0] = float(anc[0]) + (ddx / d_anc) * 2.2
                    p[1] = float(anc[1]) + (ddy / d_anc) * 2.2

    return out_ligs
