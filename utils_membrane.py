#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""
Transmembrane helix bundle detection from 3D SSE vectors, 15-degree angle
quantization, and 2D layout membrane alignment for 2DProts.
"""

import math
import numpy as np
import utils_sses_types


ANGLE_STEP = 30.0


def snap_angle(angle, step=None):
    """Quantize an angle in degrees to the nearest multiple of `step` (mod 360)."""
    if step is None:
        step = ANGLE_STEP
    if step <= 0:
        return float(angle % 360.0)
    return float((round(float(angle) / step) * step) % 360.0)


def snap_angles_dict(angles_dict, beta_connectivity=None, step=None):
    """
    Quantize all SSE angles in `angles_dict` to multiples of `step` degrees,
    preserving exact parallel (0 deg) / antiparallel (180 deg) relationships
    within each beta-sheet cluster.
    """
    if step is None:
        step = ANGLE_STEP
    if not angles_dict:
        return angles_dict

    snapped = {}
    handled = set()

    if beta_connectivity:
        for cluster in beta_connectivity:
            sses = [s for s in cluster.split(":") if s in angles_dict]
            if len(sses) > 1:
                ref_sse = sses[0]
                ref_orig = angles_dict[ref_sse] % 360.0
                ref_snap = snap_angle(ref_orig, step)
                snapped[ref_sse] = ref_snap
                handled.add(ref_sse)
                for s in sses[1:]:
                    diff = (angles_dict[s] - ref_orig) % 360.0
                    if 90.0 < diff < 270.0:
                        snapped[s] = (ref_snap + 180.0) % 360.0
                    else:
                        snapped[s] = ref_snap
                    handled.add(s)

    for sse, ang in angles_dict.items():
        if sse not in handled:
            snapped[sse] = snap_angle(ang, step)

    return snapped


def _detect_transmembrane_beta_barrel(vectors_3d, helices_residues=None, h_type=None, angle_step=None, beta_connectivity=None):
    """
    Detect transmembrane beta-barrels (e.g. 8-stranded OmpA/OmpX, 16-stranded Porin, 18-stranded OmpF)
    from 3D SSE start/end vectors and beta-sheet connectivity, while excluding soluble globular
    beta-barrels (such as GFP, lipocalins, and TIM barrels).
    """
    step = angle_step if angle_step is not None else ANGLE_STEP
    if not vectors_3d:
        return None

    sheets = {}
    med_helices_count = 0
    for sse, vec_pair in vectors_3d.items():
        if not vec_pair or len(vec_pair) < 2:
            continue
        if h_type is not None and sse in h_type and h_type[sse] not in (1, 2):
            continue
        s = np.array(vec_pair[0], dtype=float)
        e = np.array(vec_pair[1], dtype=float)
        diff = e - s
        L = float(np.linalg.norm(diff))
        if L < 1e-3:
            continue
        if utils_sses_types.is_sheet(sse):
            sheets[sse] = {
                "start": s,
                "end": e,
                "mid": 0.5 * (s + e),
                "vec": diff,
                "unit": diff / L,
                "length": L,
            }
        elif utils_sses_types.is_helic(sse) and L >= 15.0:
            med_helices_count += 1

    # Find transmembrane beta-barrel cluster(s) from beta_connectivity, or fallback to long sheets
    barrel_cl = []
    if beta_connectivity:
        multi_barrels = []
        for cl_str in beta_connectivity:
            parts = [s for s in cl_str.split(":") if s in sheets]
            if len(parts) >= 8:
                multi_barrels.extend(parts)
            if len(parts) > len(barrel_cl):
                barrel_cl = parts
        if len(multi_barrels) > len(barrel_cl):
            barrel_cl = multi_barrels
    else:
        barrel_cl = [s for s, v in sheets.items() if v["length"] >= 18.0]

    if len(barrel_cl) < 8:
        return None

    # Exclude alpha/beta globular folds (e.g. TIM barrels with many helices)
    if med_helices_count > 2 and not (len(barrel_cl) >= 14 and med_helices_count <= max(4, len(barrel_cl) // 4)):
        return None

    lengths = sorted([sheets[s]["length"] for s in barrel_cl], reverse=True)
    max_L = lengths[0]
    top6_mean_L = float(np.mean(lengths[: min(6, len(lengths))]))
    long_strands_count = sum(1 for L in lengths if L >= 22.0)

    # Transmembrane beta-barrels have very long membrane-spanning strands (max >= 42 A, top-6 mean >= 36 A),
    # whereas globular beta-barrels (GFP 1ema max=37.2 A, Lipocalin 1bbp max=33.5 A, TIM barrel max=22.3 A) are shorter.
    if max_L < 41.5 or top6_mean_L < 35.5 or long_strands_count < 6:
        return None

    def sse_sort_key(sse_name):
        if helices_residues and sse_name in helices_residues:
            hr = helices_residues[sse_name]
            try:
                return min(int(str(hr[0]).split("_")[-1]), int(str(hr[1]).split("_")[-1]))
            except ValueError:
                pass
        digits = "".join(c for c in sse_name if c.isdigit())
        return int(digits) if digits else 0

    seq_strands = sorted(barrel_cl, key=sse_sort_key)

    # Reference strand: first full-length barrel strand (L >= 20 A).
    # In outer-membrane beta-barrels, Strand 1 starts in the periplasm (bottom, -y) and goes UP to extracellular (+y).
    # Since +n_mem points from TOP (+y) to BOTTOM (-y) in 2DProts, Strand 1 has u . n_mem < 0 -> ref_u = -u_strand1.
    ref_s = next((s for s in seq_strands if sheets[s]["length"] >= 20.0), seq_strands[0])
    ref_u = -sheets[ref_s]["unit"]

    weighted_u = np.zeros(3, dtype=float)
    for s in seq_strands:
        u = sheets[s]["unit"]
        w = sheets[s]["length"]
        weighted_u += w * (u if np.dot(u, ref_u) >= 0.0 else -u)

    norm_val = float(np.linalg.norm(weighted_u))
    if norm_val < 1e-6:
        return None
    n_mem = weighted_u / norm_val

    long_seq_strands = [s for s in seq_strands if sheets[s]["length"] >= 20.0]
    mids = np.array([sheets[s]["mid"] for s in (long_seq_strands or seq_strands)], dtype=float)
    c_mem = mids.mean(axis=0)
    z_vals = np.dot(mids - c_mem, -n_mem)
    if float(np.std(z_vals)) > 8.5:
        return None

    abs_dots = [abs(float(np.dot(sheets[s]["unit"], n_mem))) for s in seq_strands if sheets[s]["length"] >= 20.0]
    mean_abs_dot = float(np.mean(abs_dots)) if abs_dots else 0.0
    if mean_abs_dot < 0.70:
        return None

    mean_tilt_deg = math.degrees(math.acos(min(1.0, max(0.0, mean_abs_dot))))
    max_tm_tilt = 30.0 if step >= 30.0 else 45.0
    snapped_tilt = min(max_tm_tilt, snap_angle(mean_tilt_deg, step))

    arbitrary = np.array([0.0, 1.0, 0.0]) if abs(n_mem[1]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e_lat1 = np.cross(n_mem, arbitrary)
    e_lat1 = e_lat1 / np.linalg.norm(e_lat1)
    e_lat2 = np.cross(n_mem, e_lat1)
    e_lat2 = e_lat2 / np.linalg.norm(e_lat2)

    sse_z_3d = {}
    sse_x_3d = {}
    sse_y_depth_3d = {}
    strand_cos_n = {}
    tm_tilts_deg = {}
    tm_target_angles = {}

    for sse, vec_pair in vectors_3d.items():
        if not vec_pair or len(vec_pair) < 2:
            continue
        s = np.array(vec_pair[0], dtype=float)
        e = np.array(vec_pair[1], dtype=float)
        mid = 0.5 * (s + e)
        vec = e - s
        L = float(np.linalg.norm(vec))
        if L < 1e-3:
            continue
        u = vec / L
        rel_mid = mid - c_mem
        sse_z_3d[sse] = float(np.dot(rel_mid, -n_mem))
        sse_x_3d[sse] = float(np.dot(rel_mid, e_lat1))
        sse_y_depth_3d[sse] = float(np.dot(rel_mid, e_lat2))
        cos_n = float(np.dot(u, n_mem))
        strand_cos_n[sse] = cos_n
        if sse in barrel_cl:
            tm_tilts_deg[sse] = round(math.degrees(math.acos(min(1.0, abs(cos_n)))), 1)
            if cos_n >= 0.0:
                tm_target_angles[sse] = (0.0 + snapped_tilt) % 360.0
            else:
                tm_target_angles[sse] = (180.0 + snapped_tilt) % 360.0

    # Exclude short N/C-terminal pre/post-strands that lie completely outside the membrane slab (e.g. E6 in 2omfA00)
    while len(seq_strands) > 8 and sheets[seq_strands[0]]["length"] < 15.0 and abs(sse_z_3d.get(seq_strands[0], 0.0)) > 14.0:
        seq_strands.pop(0)
    while len(seq_strands) > 8 and sheets[seq_strands[-1]]["length"] < 15.0 and abs(sse_z_3d.get(seq_strands[-1], 0.0)) > 14.0:
        seq_strands.pop(-1)

    return {
        "is_transmembrane": True,
        "is_beta_barrel": True,
        "tm_helices": seq_strands,
        "ref_strand": ref_s,
        "barrel_tilt_deg": float(snapped_tilt),
        "mean_3d_tilt_deg": round(mean_tilt_deg, 1),
        "strand_cos_n": strand_cos_n,
        "tm_tilts_deg": tm_tilts_deg,
        "tm_target_angles": tm_target_angles,
        "surface_target_angles": {},
        "normal_3d": [float(x) for x in n_mem],
        "center_3d": [float(x) for x in c_mem],
        "e_lat1_3d": [float(x) for x in e_lat1],
        "e_lat2_3d": [float(x) for x in e_lat2],
        "sse_z_3d": sse_z_3d,
        "sse_x_3d": sse_x_3d,
        "sse_y_depth_3d": sse_y_depth_3d,
        "thickness_3d": 28.0,
        "band_half_height_2d": 2.8,
    }


def detect_transmembrane_from_3d(vectors_3d, helices_residues=None, h_type=None, angle_step=None, beta_connectivity=None, force_tm=False):
    """
    Automatically detect transmembrane (TM) helix bundles or transmembrane beta-barrels
    and the 3D membrane plane solely from 3D SSE start/end vectors (`3dlayout` / `Annotation.vectors`).

    Returns a dictionary with membrane metadata if a TM bundle or TM beta-barrel is detected,
    or None if the protein/domain is globular.
    """
    step = angle_step if angle_step is not None else ANGLE_STEP
    if not vectors_3d:
        return None

    helices = {}
    sheets_count = 0
    med_helices_count = 0

    for sse, vec_pair in vectors_3d.items():
        if not vec_pair or len(vec_pair) < 2:
            continue
        if h_type is not None and sse in h_type and h_type[sse] not in (1, 2):
            continue
        s = np.array(vec_pair[0], dtype=float)
        e = np.array(vec_pair[1], dtype=float)
        diff = e - s
        length = float(np.linalg.norm(diff))
        if length < 1e-3:
            continue

        if utils_sses_types.is_sheet(sse):
            if length >= 3.0:
                sheets_count += 1
        elif utils_sses_types.is_helic(sse):
            if length >= 15.0:
                med_helices_count += 1
            helices[sse] = {
                "start": s,
                "end": e,
                "mid": 0.5 * (s + e),
                "vec": diff,
                "unit": diff / length,
                "length": length,
            }

    # Candidate TM helices must span at least ~25 A in 3D (~17 residues), or >=15.5 A when force_tm=True on pure helical channels
    min_cand_len = 15.5 if (force_tm and sheets_count <= 4) else 25.0
    candidates = {k: v for k, v in helices.items() if v["length"] >= min_cand_len}
    if len(candidates) < 2:
        return _detect_transmembrane_beta_barrel(
            vectors_3d, helices_residues=helices_residues, h_type=h_type, angle_step=step, beta_connectivity=beta_connectivity
        )

    # Find largest mutually aligned bundle of long helices lying in a common slab
    cos_max_angle = math.cos(math.radians(40.0 if force_tm else 38.0))
    max_z_off = 18.0 if force_tm else 14.0
    best_cluster = []

    for seed_k, seed_v in candidates.items():
        u0 = seed_v["unit"]
        m0 = seed_v["mid"]
        aligned = []
        for k, v in candidates.items():
            cos_a = abs(float(np.dot(v["unit"], u0)))
            if cos_a >= cos_max_angle:
                # Check midpoint projection along u0 so helices are side-by-side in a slab
                z_offset = abs(float(np.dot(v["mid"] - m0, u0)))
                if z_offset <= max_z_off:
                    aligned.append(k)
        if len(aligned) > len(best_cluster):
            best_cluster = aligned

    # Also check funnel/cone-tilted TM barrels or wide oligomeric TM complexes (e.g., 1okc, 1bl8_tetramer, 1j4n_tetramer)
    if (len(best_cluster) < 4 or (force_tm and sheets_count <= 4)) and len(candidates) >= 4 and (sheets_count <= 4 or force_tm):
        cand_keys = list(candidates.keys())
        mids_all = np.array([candidates[k]["mid"] for k in cand_keys])
        us_all = np.array([candidates[k]["unit"] for k in cand_keys])
        mids_c = mids_all - mids_all.mean(axis=0)
        _, _, vh_plane = np.linalg.svd(mids_c, full_matrices=False)
        n_plane = vh_plane[2]
        dots_plane = np.abs(us_all @ n_plane)
        z_offsets_plane = np.abs(mids_c @ n_plane)
        z_tol = 14.0 if force_tm else 8.0
        cone_aligned = [
            cand_keys[idx]
            for idx in range(len(cand_keys))
            if dots_plane[idx] >= (0.50 if force_tm else 0.56) and z_offsets_plane[idx] <= z_tol
        ]
        if len(cone_aligned) >= 4 and float(np.mean([dots_plane[cand_keys.index(k)] for k in cone_aligned])) >= (0.62 if force_tm else 0.70):
            if len(cone_aligned) > len(best_cluster):
                best_cluster = cone_aligned

    n_tm = len(best_cluster)
    if n_tm < 2:
        return _detect_transmembrane_beta_barrel(
            vectors_3d, helices_residues=helices_residues, h_type=h_type, angle_step=step, beta_connectivity=beta_connectivity
        )

    avg_tm_length = sum(candidates[k]["length"] for k in best_cluster) / n_tm
    tm_fraction = n_tm / max(1, med_helices_count)

    # Distinguish true transmembrane bundles from globular alpha/beta domains
    is_valid_tm = bool(force_tm and n_tm >= 2)
    if n_tm >= 4 and tm_fraction >= 0.45 and avg_tm_length >= 28.0:
        is_valid_tm = True
    elif n_tm >= 2 and tm_fraction >= 0.50 and avg_tm_length >= 32.0 and sheets_count <= n_tm:
        is_valid_tm = True

    if not is_valid_tm:
        return _detect_transmembrane_beta_barrel(
            vectors_3d, helices_residues=helices_residues, h_type=h_type, angle_step=step, beta_connectivity=beta_connectivity
        )

    # Order TM helices by sequence order if residue numbers are available
    def sse_sort_key(sse_name):
        if helices_residues and sse_name in helices_residues:
            hr = helices_residues[sse_name]
            try:
                return min(int(str(hr[0]).split("_")[-1]), int(str(hr[1]).split("_")[-1]))
            except ValueError:
                pass
        digits = "".join(c for c in sse_name if c.isdigit())
        return int(digits) if digits else 0

    best_cluster = sorted(best_cluster, key=sse_sort_key)

    # Compute length-weighted membrane normal vector n_mem
    tm_mids = np.array([candidates[k]["mid"] for k in best_cluster])
    c_mem = tm_mids.mean(axis=0)
    ref_u = candidates[best_cluster[0]]["unit"]
    if len(best_cluster) >= 4:
        _, _, vh_tm = np.linalg.svd(tm_mids - c_mem, full_matrices=False)
        n_slab = vh_tm[2]
        if np.dot(n_slab, ref_u) < 0:
            n_slab = -n_slab
        ref_u = n_slab

    weighted_sum = np.zeros(3, dtype=float)
    for k in best_cluster:
        u = candidates[k]["unit"]
        w = candidates[k]["length"]
        if np.dot(u, ref_u) < 0:
            weighted_sum -= w * u
        else:
            weighted_sum += w * u

    norm_val = np.linalg.norm(weighted_sum)
    if norm_val < 1e-6:
        return None
    # n_mem points in the direction of the first TM helix (start -> end, i.e. top -> bottom in 2DProts)
    n_mem = weighted_sum / norm_val

    # Verify planar slab tightness
    z_tm = np.dot(tm_mids - c_mem, -n_mem)
    if float(np.std(z_tm)) > (14.0 if force_tm else 8.5):
        return None

    # Find principal lateral axis e_lat1 in the membrane plane (perpendicular to n_mem)
    tm_mids_centered = tm_mids - c_mem
    tm_mids_in_plane = tm_mids_centered - np.outer(np.dot(tm_mids_centered, n_mem), n_mem)
    k_half = max(1, len(best_cluster) // 2)
    seq_vec_3d = tm_mids_in_plane[len(best_cluster) - k_half :].mean(axis=0) - tm_mids_in_plane[:k_half].mean(axis=0)
    if len(best_cluster) >= 2 and np.linalg.norm(tm_mids_in_plane) > 1e-3:
        _, _, vh = np.linalg.svd(tm_mids_in_plane, full_matrices=False)
        e_lat1 = vh[0]
        e_lat1 = e_lat1 - np.dot(e_lat1, n_mem) * n_mem
        e_lat1 = e_lat1 / np.linalg.norm(e_lat1)
        if np.linalg.norm(seq_vec_3d) > 1.0:
            if np.dot(seq_vec_3d, e_lat1) < 0:
                e_lat1 = -e_lat1
        elif np.dot(candidates[best_cluster[-1]]["mid"] - candidates[best_cluster[0]]["mid"], e_lat1) < 0:
            e_lat1 = -e_lat1
    else:
        arbitrary = np.array([0.0, 1.0, 0.0]) if abs(n_mem[1]) < 0.9 else np.array([1.0, 0.0, 0.0])
        e_lat1 = np.cross(n_mem, arbitrary)
        e_lat1 = e_lat1 / np.linalg.norm(e_lat1)
    e_lat2 = np.cross(n_mem, e_lat1)
    e_lat2 = e_lat2 / np.linalg.norm(e_lat2)

    # Compute 3D height z_3d (along -n_mem, where +z_3d is up/+y in 2D) and target 2D angles
    sse_z_3d = {}
    sse_x_3d = {}
    sse_y_depth_3d = {}
    tm_target_angles = {}
    surface_target_angles = {}
    tm_tilts_deg = {}

    for sse, vec_pair in vectors_3d.items():
        if not vec_pair or len(vec_pair) < 2:
            continue
        s = np.array(vec_pair[0], dtype=float)
        e = np.array(vec_pair[1], dtype=float)
        mid = 0.5 * (s + e)
        vec = e - s
        L = float(np.linalg.norm(vec))
        if L < 1e-3:
            continue
        u = vec / L
        rel_mid = mid - c_mem
        sse_z_3d[sse] = float(np.dot(rel_mid, -n_mem))
        sse_x_3d[sse] = float(np.dot(rel_mid, e_lat1))
        sse_y_depth_3d[sse] = float(np.dot(rel_mid, e_lat2))

        if sse in best_cluster:
            # True 3D tilt angle from membrane normal (0..45 deg)
            cos_n = float(np.dot(u, n_mem))
            tilt_deg = math.degrees(math.acos(min(1.0, abs(cos_n))))
            tm_tilts_deg[sse] = round(tilt_deg, 1)

            # Determine tilt direction sign in lateral plane
            lat_proj = float(np.dot(u, e_lat1))
            if abs(lat_proj) < 0.08:
                lat_proj = float(np.dot(u, e_lat2))
            tilt_sign = -1.0 if lat_proj >= 0 else 1.0

            # In 2DProts, angle=0 means start(N) at +y (top) and end(C) at -y (bottom).
            # Since n_mem points from top(+y) to bottom(-y), cos_n > 0 means base angle 0 deg,
            # while cos_n < 0 means base angle 180 deg.
            # Cap TM 2D tilt so TM helices always clearly traverse the membrane (0 or +/-30 deg in 30-deg mode).
            max_tm_tilt = 30.0 if step >= 30.0 else 45.0
            snapped_tilt = min(max_tm_tilt, snap_angle(tilt_deg, step))
            if cos_n >= 0:
                tm_target_angles[sse] = (0.0 + tilt_sign * snapped_tilt) % 360.0
            else:
                tm_target_angles[sse] = (180.0 - tilt_sign * snapped_tilt) % 360.0
        elif utils_sses_types.is_helic(sse) and L >= 10.0:
            # Check for amphipathic surface helices lying flat on the membrane (like GPCR Helix 8)
            cos_n = float(np.dot(u, n_mem))
            if abs(cos_n) <= 0.30 and 12.0 <= abs(sse_z_3d[sse]) <= 28.0:
                lat_proj = float(np.dot(u, e_lat1))
                base_horiz = 270.0 if lat_proj >= 0 else 90.0
                surface_target_angles[sse] = snap_angle(base_horiz, step)

    # Estimate membrane hydrophobic thickness in 3D (default ~30 A -> 6.0 canvas units)
    thickness_3d = min(34.0, max(26.0, avg_tm_length * 0.75))

    return {
        "is_transmembrane": True,
        "tm_helices": best_cluster,
        "tm_tilts_deg": tm_tilts_deg,
        "tm_target_angles": tm_target_angles,
        "surface_target_angles": surface_target_angles,
        "normal_3d": [float(x) for x in n_mem],
        "center_3d": [float(x) for x in c_mem],
        "e_lat1_3d": [float(x) for x in e_lat1],
        "e_lat2_3d": [float(x) for x in e_lat2],
        "sse_z_3d": sse_z_3d,
        "sse_x_3d": sse_x_3d,
        "sse_y_depth_3d": sse_y_depth_3d,
        "thickness_3d": float(thickness_3d),
        "band_half_height_2d": float(thickness_3d * 0.1),  # matches Image.CIRCLE_MULTIPLIER = 0.1
    }


def _align_beta_barrel_to_membrane(layout, angles, sizes, beta_connectivity, vectors_3d,
                                   mem_info, helices_residues=None, h_type=None,
                                   ligands=None, channels=None, angle_step=None,
                                   barrel_mode="unrolled_30"):
    """
    Aligns a transmembrane beta-barrel to the horizontal membrane band:
      - In `barrel_mode="unrolled"` (default), unrolls the beta-barrel strands sequentially
        left-to-right across the membrane band at their quantized 3D barrel shear angle
        (e.g. 30/210 deg in 30-deg mode, or 0/180 deg when barrel_mode="unrolled_0").
      - In `barrel_mode="arch"` (or `"arch_0"`), preserves the 2D cylindrical front/back wall
        projection inside the horizontal membrane band.
      - Groups split collinear strand segments (e.g. beta-bulge split strands) onto the same
        tilted strand track with guaranteed non-overlapping vertical separation.
      - Positions extramembrane loops, turns, and hairpins directly above/below their
        sequence-flanking barrel strand endpoints.
    """
    step = angle_step if angle_step is not None else ANGLE_STEP
    seq_strands = [s for s in mem_info["tm_helices"] if s in layout]
    barrel_set = set(seq_strands)

    def seq_key(sse_name):
        if helices_residues and sse_name in helices_residues:
            hr = helices_residues[sse_name]
            try:
                return min(int(str(hr[0]).split("_")[-1]), int(str(hr[1]).split("_")[-1]))
            except ValueError:
                pass
        digits = "".join(c for c in sse_name if c.isdigit())
        return int(digits) if digits else 0

    seq_strands = sorted(seq_strands, key=seq_key)
    ref_s = mem_info.get("ref_strand", seq_strands[0])
    if ref_s not in angles:
        ref_s = seq_strands[0]

    # Rotate 2D layout so ref_s points at 180 deg (N-entry at bottom -y, C-exit at top +y)
    ref_ang_rad = math.radians((angles[ref_s] - 180.0) % 360.0)
    ca = math.cos(ref_ang_rad)
    sa = math.sin(ref_ang_rad)
    pts_b = np.array([layout[s] for s in seq_strands], dtype=float)
    cb2d = pts_b.mean(axis=0)

    def rot_to_vert(pt):
        dx = pt[0] - cb2d[0]
        dy = pt[1] - cb2d[1]
        return [ca * dx - sa * dy, sa * dx + ca * dy]

    rot_l = {s: rot_to_vert(pt) for s, pt in layout.items()}
    if rot_l[seq_strands[-1]][0] < rot_l[seq_strands[0]][0]:
        for s in rot_l:
            rot_l[s][0] = -rot_l[s][0]

    if barrel_mode.endswith("_0"):
        tilt = 0.0
    elif barrel_mode.endswith("_30"):
        tilt = 30.0
    elif barrel_mode.endswith("_15"):
        tilt = 15.0
    else:
        tilt = float(mem_info.get("barrel_tilt_deg", 30.0 if step >= 30.0 else step))

    tilt_rad = math.radians(tilt)
    strand_cos_n = mem_info.get("strand_cos_n", {})
    sse_z_3d = mem_info.get("sse_z_3d", {})
    band_half = float(mem_info.get("band_half_height_2d", 2.8))

    new_angles = {}
    for s, ang in angles.items():
        if s in barrel_set:
            if strand_cos_n.get(s, 1.0) >= 0.0:
                new_angles[s] = (0.0 + tilt) % 360.0
            else:
                new_angles[s] = (180.0 + tilt) % 360.0
        else:
            if abs(strand_cos_n.get(s, 0.0)) >= 0.55 and utils_sses_types.is_sheet(s):
                base_a = 0.0 if strand_cos_n[s] >= 0.0 else 180.0
                new_angles[s] = (base_a + tilt) % 360.0
            else:
                rot_a = (ang - math.degrees(ref_ang_rad)) % 360.0
                new_angles[s] = snap_angle(rot_a, step)

    new_layout = {}

    if not barrel_mode.startswith("arch"):
        # Unrolled barrel mode: group split collinear strand pieces onto the same column track
        columns = []
        for s in seq_strands:
            sign_s = 1 if strand_cos_n.get(s, 1.0) >= 0.0 else -1
            if columns:
                prev_s = columns[-1][-1]
                sign_prev = 1 if strand_cos_n.get(prev_s, 1.0) >= 0.0 else -1
                L_s = float(sizes.get(s, 20.0)) if sizes else 20.0
                L_prev = float(sizes.get(prev_s, 20.0)) if sizes else 20.0
                if sign_s == sign_prev and min(L_s, L_prev) < 15.0:
                    columns[-1].append(s)
                    continue
            columns.append([s])

        def _b_dom_prefix(name):
            s = name.lstrip("_")
            if s and s[0] in ("H", "E") and "_" in s:
                return s[1 : s.index("_")]
            return ""

        spacing = (1.68 if tilt >= 25.0 else 1.62) if len(columns) > 22 else (2.55 if tilt >= 25.0 else 2.45)
        n_cols = len(columns)
        strand_col_x = {}
        col_xs = []
        curr_x = 0.0
        for col_idx, col_sses in enumerate(columns):
            if col_idx > 0:
                prev_dom = _b_dom_prefix(columns[col_idx - 1][0])
                curr_dom = _b_dom_prefix(col_sses[0])
                curr_x += spacing + (2.0 if prev_dom != curr_dom else 0.0)
            col_xs.append(curr_x)
        mean_col_x = sum(col_xs) / max(1, len(col_xs))
        col_xs = [cx - mean_col_x for cx in col_xs]

        for col_idx, col_sses in enumerate(columns):
            base_x = col_xs[col_idx]
            if len(col_sses) == 1:
                s = col_sses[0]
                s_offset = 0.14 * sse_z_3d.get(s, 0.0)
                nx = base_x + s_offset * math.sin(tilt_rad)
                ny = s_offset * math.cos(tilt_rad)
                new_layout[s] = [float(nx), float(ny)]
                strand_col_x[s] = base_x
            else:
                # Multiple collinear pieces in the same strand track (e.g. E12 + E19 in 2omfA00):
                # sort from bottom (-z_3d) to top (+z_3d) and enforce non-overlapping 1D intervals along track
                col_sorted = sorted(col_sses, key=lambda x: sse_z_3d.get(x, 0.0))
                s_offsets = [0.14 * sse_z_3d.get(x, 0.0) for x in col_sorted]
                hls = [0.1 * float(sizes.get(x, 10.0)) if sizes else 1.0 for x in col_sorted]
                for i_c in range(1, len(col_sorted)):
                    min_sep = hls[i_c - 1] + hls[i_c] + 0.75
                    if s_offsets[i_c] - s_offsets[i_c - 1] < min_sep:
                        if hls[i_c] < hls[i_c - 1]:
                            s_offsets[i_c] = s_offsets[i_c - 1] + min_sep
                        else:
                            s_offsets[i_c - 1] = s_offsets[i_c] - min_sep
                for s, s_off in zip(col_sorted, s_offsets):
                    nx = base_x + s_off * math.sin(tilt_rad)
                    ny = s_off * math.cos(tilt_rad)
                    new_layout[s] = [float(nx), float(ny)]
                    strand_col_x[s] = base_x

        # Compute true vertical extent of barrel strands so extramembrane SSEs never collide with strand tips
        barrel_top_y = float(band_half)
        barrel_bot_y = -float(band_half)
        for b in barrel_set:
            if b in new_layout:
                hl_b = 0.10 * float(sizes.get(b, 12.0)) if sizes else 1.2
                ang_b = math.radians(new_angles.get(b, 0.0))
                vy_b = abs(math.cos(ang_b)) * hl_b
                barrel_top_y = max(barrel_top_y, new_layout[b][1] + vy_b)
                barrel_bot_y = min(barrel_bot_y, new_layout[b][1] - vy_b)

        # Position non-barrel SSE clusters anchored to sequence-flanking barrel strands of the same subunit
        sse_to_cl = {}
        if beta_connectivity:
            for cl_str in beta_connectivity:
                parts = [s for s in cl_str.split(":") if s in layout and s not in barrel_set]
                for s in parts:
                    sse_to_cl[s] = parts

        processed = set(barrel_set)
        non_barrel_clusters = []
        cluster_side = {}
        for s in layout:
            if s in processed:
                continue
            cl = sse_to_cl.get(s, [s])
            cl = [x for x in cl if x not in processed]
            if not cl:
                continue

            cl_sorted = sorted(cl, key=seq_key)
            cl_dom = _b_dom_prefix(cl_sorted[0])
            dom_strands = [b for b in seq_strands if _b_dom_prefix(b) == cl_dom] or seq_strands
            cl_seq = sum(seq_key(x) for x in cl_sorted) / len(cl_sorted)
            before_b = [b for b in dom_strands if seq_key(b) <= cl_seq]
            after_b = [b for b in dom_strands if seq_key(b) >= cl_seq]
            b1 = before_b[-1] if before_b else dom_strands[0]
            b2 = after_b[0] if after_b else dom_strands[-1]
            base_anchor_x = 0.5 * (strand_col_x.get(b1, new_layout[b1][0]) + strand_col_x.get(b2, new_layout[b2][0]))

            mean_z = sum(sse_z_3d.get(x, 0.0) for x in cl_sorted) / len(cl_sorted)
            side_sign = 1 if mean_z >= 0.0 else -1
            if side_sign > 0:
                target_cy = max(barrel_top_y + 2.5, barrel_top_y + 1.5 + 0.55 * max(0.0, abs(mean_z) - 6.0))
            else:
                target_cy = min(barrel_bot_y - 2.5, barrel_bot_y - 1.5 - 0.55 * max(0.0, abs(mean_z) - 6.0))

            anchor_x = base_anchor_x + np.clip(target_cy, -5.0, 5.0) * math.tan(tilt_rad)
            n_c = len(cl_sorted)
            for idx_c, x in enumerate(cl_sorted):
                dx = (idx_c - 0.5 * (n_c - 1)) * 2.15
                dy = 0.40 * (sse_z_3d.get(x, mean_z) - mean_z)
                new_layout[x] = [float(anchor_x + dx), float(target_cy + dy)]
                processed.add(x)
            non_barrel_clusters.append(cl_sorted)
            cluster_side[tuple(cl_sorted)] = side_sign

        if non_barrel_clusters:
            _resolve_extramembrane_overlaps(
                new_layout,
                new_angles,
                sizes,
                fixed_sses=barrel_set,
                movable_clusters=non_barrel_clusters,
                cluster_side=cluster_side,
                top_bound_y=barrel_top_y + 2.2,
                bot_bound_y=barrel_bot_y - 2.2,
                min_sep=3.1,
            )
    else:
        # Cylindrical arch projection mode
        ordered_by_x = sorted(seq_strands, key=lambda s: (rot_l[s][0], rot_l[s][1]))
        min_dx = 2.05 if len(seq_strands) <= 10 else 1.80
        xs = [rot_l[s][0] * 1.35 for s in ordered_by_x]
        for i in range(1, len(xs)):
            if xs[i] - xs[i - 1] < min_dx:
                shift = min_dx - (xs[i] - xs[i - 1])
                for j in range(i, len(xs)):
                    xs[j] += shift
        mean_x = sum(xs) / len(xs)
        xs = [x - mean_x for x in xs]

        ys_raw = [rot_l[s][1] for s in ordered_by_x]
        max_abs_y = max(1.0, max(abs(y) for y in ys_raw))
        arch_amp = 1.15 if len(seq_strands) <= 10 else 1.55
        for s, base_x in zip(ordered_by_x, xs):
            z3d = sse_z_3d.get(s, 0.0)
            s_offset = (rot_l[s][1] / max_abs_y) * arch_amp + 0.08 * z3d
            nx = base_x + s_offset * math.sin(tilt_rad)
            ny = s_offset * math.cos(tilt_rad)
            new_layout[s] = [float(nx), float(ny)]

        for s in layout:
            if s in new_layout:
                continue
            z3d = sse_z_3d.get(s, 0.0)
            nx = rot_l[s][0] * 1.25
            ny = 0.45 * z3d
            if abs(z3d) > 8.0:
                ny = max(ny, band_half + 2.5) if z3d > 0.0 else min(ny, -(band_half + 2.5))
            new_layout[s] = [float(nx), float(ny)]

    if ligands:
        localize_ligands_2d(
            ligands,
            old_layout=rot_l,
            new_layout=new_layout,
            new_angles=new_angles,
            sizes=sizes,
            mem_info=mem_info,
            pre_transform=rot_to_vert,
        )

    if channels:
        localize_channels_2d(
            channels,
            layout_2d=new_layout,
            angles_2d=new_angles,
            sizes=sizes,
            vectors_3d=vectors_3d,
            h_type=h_type,
            ligands=ligands,
            mem_info=mem_info,
            old_layout=rot_l,
            pre_transform=rot_to_vert,
        )

    mem_info["y_center_2d"] = 0.0
    mem_info["barrel_mode"] = barrel_mode
    mem_info["angle_step"] = float(step)
    return new_layout, new_angles, mem_info


def _sse_half_len_2d(sse, sizes):
    sz = float(sizes.get(sse, 12.0)) if sizes else 12.0
    return (0.145 if utils_sses_types.is_helic(sse) else 0.10) * sz


def _seg_dist_and_vec_2d(p1, q1, p2, q2):
    """Returns (min_dist, midpoint_delta_vector_from_2_to_1) between 2D segments p1-q1 and p2-q2."""
    u = q1 - p1
    v = q2 - p2
    w = p1 - p2
    a = float(np.dot(u, u))
    b = float(np.dot(u, v))
    c = float(np.dot(v, v))
    d = float(np.dot(u, w))
    e = float(np.dot(v, w))
    D = a * c - b * b
    sc, sN, sD = D, D, D
    tc, tN, tD = D, D, D
    if D < 1e-8:
        sN = 0.0
        sD = 1.0
        tN = e
        tD = c
    else:
        sN = b * e - c * d
        tN = a * e - b * d
        if sN < 0.0:
            sN = 0.0
            tN = e
            tD = c
        elif sN > sD:
            sN = sD
            tN = e + b
            tD = c
    if tN < 0.0:
        tN = 0.0
        if -d < 0.0:
            sN = 0.0
        elif -d > a:
            sN = sD
        else:
            sN = -d
            sD = a
    elif tN > tD:
        tN = tD
        if (-d + b) < 0.0:
            sN = 0.0
        elif (-d + b) > a:
            sN = sD
        else:
            sN = -d + b
            sD = a
    sc = 0.0 if abs(sN) < 1e-8 else sN / sD
    tc = 0.0 if abs(tN) < 1e-8 else tN / tD
    dP = w + (sc * u) - (tc * v)
    dist = float(np.linalg.norm(dP))
    return dist, dP


def _resolve_extramembrane_overlaps(layout, angles, sizes, fixed_sses, movable_clusters,
                                    cluster_side, top_bound_y, bot_bound_y, n_iters=65, min_sep=3.25):
    """
    Separates overlapping extramembrane SSE segments by shifting movable clusters (rigid beta-sheets
    or single extramembrane helices) in the open extramembrane space while strictly keeping upper
    clusters above `top_bound_y` and lower clusters below `bot_bound_y`.
    """
    if not movable_clusters:
        return

    sse_to_cid = {}
    for cid, cl in enumerate(movable_clusters):
        for s in cl:
            sse_to_cid[s] = cid

    all_sses = [s for s in layout]
    seg_vecs = {}
    for s in all_sses:
        hl = _sse_half_len_2d(s, sizes)
        ang_r = math.radians(angles.get(s, 0.0))
        seg_vecs[s] = np.array([math.sin(ang_r) * hl, math.cos(ang_r) * hl], dtype=float)

    def enforce_vertical_bounds():
        for cl in movable_clusters:
            side = cluster_side.get(tuple(cl), 1)
            if side > 0:
                min_tip_y = min(layout[s][1] - abs(seg_vecs[s][1]) for s in cl)
                if min_tip_y < top_bound_y:
                    dy = top_bound_y - min_tip_y
                    for s in cl:
                        layout[s][1] += dy
            else:
                max_tip_y = max(layout[s][1] + abs(seg_vecs[s][1]) for s in cl)
                if max_tip_y > bot_bound_y:
                    dy = bot_bound_y - max_tip_y
                    for s in cl:
                        layout[s][1] += dy

    enforce_vertical_bounds()

    for it in range(n_iters):
        disp = {cid: np.zeros(2, dtype=float) for cid in range(len(movable_clusters))}
        max_overlap = 0.0
        n = len(all_sses)
        for i in range(n):
            s1 = all_sses[i]
            cid1 = sse_to_cid.get(s1, None)
            c1 = np.array(layout[s1], dtype=float)
            v1 = seg_vecs[s1]
            p1, q1 = c1 - v1, c1 + v1
            r1 = float(np.linalg.norm(v1))
            for j in range(i + 1, n):
                s2 = all_sses[j]
                cid2 = sse_to_cid.get(s2, None)
                if cid1 is None and cid2 is None:
                    continue
                if cid1 is not None and cid1 == cid2:
                    continue
                c2 = np.array(layout[s2], dtype=float)
                v2 = seg_vecs[s2]
                r2 = float(np.linalg.norm(v2))
                # Quick bounding-circle reject
                dc = c1 - c2
                if abs(dc[0]) > r1 + r2 + min_sep or abs(dc[1]) > r1 + r2 + min_sep:
                    continue
                p2, q2 = c2 - v2, c2 + v2
                dist, dP = _seg_dist_and_vec_2d(p1, q1, p2, q2)
                if dist < min_sep:
                    ov = min_sep - dist
                    if ov > max_overlap:
                        max_overlap = ov
                    if dist > 1e-4:
                        dir_vec = dP / dist
                    else:
                        if np.linalg.norm(dc) > 1e-4:
                            dir_vec = dc / float(np.linalg.norm(dc))
                        else:
                            dir_vec = np.array([1.0 if (i + j) % 2 == 0 else -1.0, 0.3], dtype=float)
                            dir_vec /= float(np.linalg.norm(dir_vec))
                    # Encourage horizontal fanning so extramembrane domains use lateral space without inflating total height
                    dir_vec[0] *= 1.75
                    dir_vec[1] *= 0.65
                    norm_dv = float(np.linalg.norm(dir_vec))
                    if norm_dv > 1e-6:
                        dir_vec /= norm_dv
                    step_mag = 0.48 * ov
                    if cid1 is not None and cid2 is not None:
                        disp[cid1] += step_mag * dir_vec
                        disp[cid2] -= step_mag * dir_vec
                    elif cid1 is not None:
                        disp[cid1] += 1.8 * step_mag * dir_vec
                    elif cid2 is not None:
                        disp[cid2] -= 1.8 * step_mag * dir_vec

        if max_overlap < 0.05:
            break
        for cid, cl in enumerate(movable_clusters):
            dv = disp[cid]
            mag = float(np.linalg.norm(dv))
            if mag > 1e-5:
                if mag > 2.5:
                    dv = dv * (2.5 / mag)
                for s in cl:
                    layout[s][0] = float(layout[s][0] + dv[0])
                    layout[s][1] = float(layout[s][1] + dv[1])
        enforce_vertical_bounds()


def align_layout_to_membrane(layout, angles, sizes, beta_connectivity, vectors_3d,
                             helices_residues=None, h_type=None, ligands=None, channels=None,
                             angle_step=None, barrel_mode="unrolled_30", force_tm=None):
    """
    If a transmembrane helix bundle or transmembrane beta-barrel is detected in `vectors_3d`:
      1. Rotates/reflects the 2D layout so the lateral spread of TM SSEs is horizontal (x-axis)
         and extramembrane domains align with the membrane normal (y-axis).
      2. Aligns the vertical (y-axis) positions of SSEs with the membrane normal so all TM SSEs
         pass through a horizontal membrane band centered at `y_mem = 0.0`, while extracellular
         and intracellular SSEs sit in the open space above and below the TM helix tips at 1:1 scale.
      3. Sets TM SSE 2D angles to vertical (0/180 deg) plus their quantized 3D tilt (in `angle_step` steps).
      4. Smoothly transforms ligands and channels to preserve their relative neighborhood to SSEs.

    Always quantizes all SSE angles to multiples of `angle_step` degrees.
    Returns (new_layout, new_angles, membrane_info).
    """
    step = angle_step if angle_step is not None else ANGLE_STEP
    if force_tm is False:
        mem_info = None
    else:
        mem_info = detect_transmembrane_from_3d(
            vectors_3d,
            helices_residues=helices_residues,
            h_type=h_type,
            angle_step=step,
            beta_connectivity=beta_connectivity,
            force_tm=bool(force_tm),
        )
    if not mem_info or not layout:
        snapped_angles = snap_angles_dict(angles, beta_connectivity=beta_connectivity, step=step)
        if ligands and layout:
            localize_ligands_2d(ligands, old_layout=layout, new_layout=layout, new_angles=snapped_angles, sizes=sizes, mem_info=None)
        if channels and layout:
            localize_channels_2d(
                channels,
                layout_2d=layout,
                angles_2d=snapped_angles,
                sizes=sizes,
                vectors_3d=vectors_3d,
                h_type=h_type,
                ligands=ligands,
                mem_info=None,
            )
        return layout, snapped_angles, None

    if mem_info.get("is_beta_barrel"):
        return _align_beta_barrel_to_membrane(
            layout,
            angles,
            sizes,
            beta_connectivity,
            vectors_3d,
            mem_info,
            helices_residues=helices_residues,
            h_type=h_type,
            ligands=ligands,
            channels=channels,
            angle_step=step,
            barrel_mode=barrel_mode,
        )

    tm_helices = [h for h in mem_info["tm_helices"] if h in layout]
    if len(tm_helices) < 2:
        snapped_angles = snap_angles_dict(angles, beta_connectivity=beta_connectivity, step=step)
        if ligands and layout:
            localize_ligands_2d(ligands, old_layout=layout, new_layout=layout, new_angles=snapped_angles, sizes=sizes, mem_info=None)
        if channels and layout:
            localize_channels_2d(
                channels,
                layout_2d=layout,
                angles_2d=snapped_angles,
                sizes=sizes,
                vectors_3d=vectors_3d,
                h_type=h_type,
                ligands=ligands,
                mem_info=None,
            )
        return layout, snapped_angles, None

    sse_z_3d = mem_info["sse_z_3d"]
    sse_x_3d = mem_info["sse_x_3d"]
    band_half = mem_info["band_half_height_2d"]  # ~3.0 canvas units
    tm_set = set(tm_helices)
    non_tm_sses = [s for s in layout if s not in tm_set]
    ext_sses = [s for s in non_tm_sses if abs(sse_z_3d.get(s, 0.0)) > 8.0]

    # 1. Center 2D layout on TM bundle centroid and find the 2D orthogonal transformation R2d in O(2)
    #    that aligns TM lateral spread with +x and extramembrane domains with +y (z_3d > 0) / -y (z_3d < 0)
    tm_pts = np.array([layout[h] for h in tm_helices], dtype=float)
    tm_center_2d = tm_pts.mean(axis=0)
    tm_centered = tm_pts - tm_center_2d

    if len(ext_sses) >= 2:
        P_list = []
        Q_list = []
        for s, coord in layout.items():
            w = 2.0 if s in tm_set else (2.5 if abs(sse_z_3d.get(s, 0.0)) > 8.0 else 1.0)
            P_list.append(w * (np.array(coord, dtype=float) - tm_center_2d))
            Q_list.append(w * np.array([sse_x_3d.get(s, 0.0), sse_z_3d.get(s, 0.0)], dtype=float))
        P_mat = np.array(P_list, dtype=float)
        Q_mat = np.array(Q_list, dtype=float)
        U_o, _, Vt_o = np.linalg.svd(P_mat.T @ Q_mat)
        R2d = U_o @ Vt_o  # 2x2 orthogonal matrix mapping row vector p -> p @ R2d
    else:
        k_half_2d = max(1, len(tm_helices) // 2)
        seq_vec_2d = tm_centered[len(tm_helices) - k_half_2d :].mean(axis=0) - tm_centered[:k_half_2d].mean(axis=0)
        if np.linalg.norm(tm_centered) > 1e-3:
            _, _, vh2d = np.linalg.svd(tm_centered, full_matrices=False)
            ax2d = vh2d[0]
            if np.linalg.norm(seq_vec_2d) > 1.0:
                if np.dot(seq_vec_2d, ax2d) < 0:
                    ax2d = -ax2d
            elif np.dot(np.array(layout[tm_helices[-1]]) - np.array(layout[tm_helices[0]]), ax2d) < 0:
                ax2d = -ax2d
        else:
            ax2d = np.array([1.0, 0.0])
        rot_angle_rad = -math.atan2(ax2d[1], ax2d[0])
        cos_r = math.cos(rot_angle_rad)
        sin_r = math.sin(rot_angle_rad)
        R2d = np.array([[cos_r, sin_r], [-sin_r, cos_r]], dtype=float)
        if ext_sses:
            test_y = sum(((np.array(layout[s]) - tm_center_2d) @ R2d)[1] * sse_z_3d.get(s, 0.0) for s in ext_sses)
            if test_y < 0.0:
                R2d[:, 1] = -R2d[:, 1]

    def rot2d(pt):
        v = (np.array(pt, dtype=float) - tm_center_2d) @ R2d
        return [float(v[0]), float(v[1])]

    def rot_ang_2d(ang):
        rad = math.radians(ang)
        v = np.array([math.sin(rad), math.cos(rad)], dtype=float) @ R2d
        return math.degrees(math.atan2(float(v[0]), float(v[1]))) % 360.0

    rotated_layout = {sse: rot2d(coord) for sse, coord in layout.items()}
    rotated_angles = {sse: rot_ang_2d(ang) for sse, ang in angles.items()}

    # 2. Ensure TM helices have clean horizontal separation along x while aligning their y to the membrane band
    new_layout = {}
    sse_to_cluster = {}
    if beta_connectivity:
        for cl in beta_connectivity:
            parts = [s for s in cl.split(":") if s in rotated_layout]
            for s in parts:
                sse_to_cluster[s] = parts

    def _sse_dom_prefix(name):
        s = name.lstrip("_")
        if s and s[0] in ("H", "E") and "_" in s:
            return s[1 : s.index("_")]
        return ""

    # Check if TM helices belong to multiple distinct subunits/domains (e.g. homo-tetramers 1j4n_tetramer, 1bl8_tetramer, 1f88_dimer)
    tm_by_dom = {}
    for h in tm_helices:
        d_pref = _sse_dom_prefix(h)
        tm_by_dom.setdefault(d_pref, []).append(h)

    all_doms = {_sse_dom_prefix(s) for s in rotated_layout}
    pure_tm_oligomer = (
        len(tm_by_dom) >= 2
        and all(len(v) >= 2 for v in tm_by_dom.values())
        and all_doms == set(tm_by_dom.keys())
        and max((abs(sse_z_3d.get(s, 0.0)) for s in rotated_layout), default=0.0) < 36.0
    )

    # Position TM helices first
    for h in tm_helices:
        rx, ry = rotated_layout[h]
        x3d = sse_x_3d.get(h, rx)
        z3d = sse_z_3d.get(h, 0.0)
        nx = 0.60 * rx * 1.25 + 0.40 * x3d * 1.15
        ny = float(np.clip(0.14 * z3d + 0.08 * ry, -1.85, 1.85))
        new_layout[h] = [float(nx), float(ny)]

    tm_targ = mem_info.get("tm_target_angles", {})
    if pure_tm_oligomer:
        # Order subunits around their 2D arrangement so each monomer's TM bundle stays contiguous
        dom_order = sorted(
            tm_by_dom.keys(),
            key=lambda d: (
                round(sum(rotated_layout[h][1] for h in tm_by_dom[d]) / len(tm_by_dom[d]) / 12.0),
                sum(rotated_layout[h][0] for h in tm_by_dom[d]) / len(tm_by_dom[d]),
            ),
        )
        tm_sorted_by_x = []
        for d in dom_order:
            d_cx = sum(rotated_layout[h][0] for h in tm_by_dom[d]) / len(tm_by_dom[d])
            d_sorted = sorted(tm_by_dom[d], key=lambda h: (rotated_layout[h][0] - d_cx))
            tm_sorted_by_x.extend(d_sorted)
            for h in d_sorted:
                new_layout[h][0] = float(rotated_layout[h][0] - d_cx)
    else:
        tm_sorted_by_x = sorted(tm_helices, key=lambda h: new_layout[h][0])

    # Enforce a minimum horizontal separation between TM helices (accounting for 30-deg tilts at top/bottom)
    for idx in range(1, len(tm_sorted_by_x)):
        h_prev = tm_sorted_by_x[idx - 1]
        h_curr = tm_sorted_by_x[idx]
        ang_prev = math.radians(tm_targ.get(h_prev, 0.0))
        ang_curr = math.radians(tm_targ.get(h_curr, 0.0))
        tan_prev = math.tan(ang_prev) if abs(math.cos(ang_prev)) > 0.2 else 0.0
        tan_curr = math.tan(ang_curr) if abs(math.cos(ang_curr)) > 0.2 else 0.0
        same_dom = _sse_dom_prefix(h_prev) == _sse_dom_prefix(h_curr)
        base_gap = 2.65 if same_dom else (3.40 if pure_tm_oligomer else 2.75)
        min_dx = max(base_gap, 1.75 + 2.15 * abs(tan_curr - tan_prev))
        if new_layout[h_curr][0] - new_layout[h_prev][0] < min_dx:
            shift = min_dx - (new_layout[h_curr][0] - new_layout[h_prev][0])
            for j in range(idx, len(tm_sorted_by_x)):
                new_layout[tm_sorted_by_x[j]][0] += shift
    # Re-center TM helices horizontally around x = 0
    mean_tm_x = sum(new_layout[h][0] for h in tm_helices) / len(tm_helices)
    for h in tm_helices:
        new_layout[h][0] -= mean_tm_x

    dom_old_tm_center = {}
    dom_new_tm_cx = {}
    if pure_tm_oligomer:
        for d, d_tms in tm_by_dom.items():
            dom_old_tm_center[d] = (
                sum(rotated_layout[h][0] for h in d_tms) / len(d_tms),
                sum(rotated_layout[h][1] for h in d_tms) / len(d_tms),
            )
            dom_new_tm_cx[d] = sum(new_layout[h][0] for h in d_tms) / len(d_tms)

    # Compute actual top and bottom vertical reach of the TM helix bars so extramembrane domains
    # sit completely clear of the TM helix ends
    tm_top_y = float(band_half)
    tm_bot_y = -float(band_half)
    for h in tm_helices:
        hl_h = _sse_half_len_2d(h, sizes)
        ang_h = math.radians(tm_targ.get(h, 0.0))
        vy_h = abs(math.cos(ang_h)) * hl_h
        tm_top_y = max(tm_top_y, new_layout[h][1] + vy_h)
        tm_bot_y = min(tm_bot_y, new_layout[h][1] - vy_h)

    # 3. Set TM helix angles (preserving 3D tilt in quantized steps) & snap all other SSE angles
    new_angles = dict(rotated_angles)
    for h, target_ang in mem_info["tm_target_angles"].items():
        if h in new_angles:
            new_angles[h] = target_ang
    for h, surf_ang in mem_info["surface_target_angles"].items():
        if h in new_angles:
            new_angles[h] = surf_ang

    new_angles = snap_angles_dict(new_angles, beta_connectivity=beta_connectivity, step=step)

    # 4. Position non-TM SSE clusters in the open extramembrane space above (+y) and below (-y)
    processed_non_tm = set(tm_helices)
    upper_clusters = []
    lower_clusters = []
    movable_clusters = []
    cluster_side = {}

    for sse in rotated_layout:
        if sse in processed_non_tm:
            continue
        cl = sse_to_cluster.get(sse, [sse])
        cl = [s for s in cl if s not in processed_non_tm]
        if not cl:
            continue
        for s in cl:
            processed_non_tm.add(s)

        z_vals = [sse_z_3d.get(s, 0.0) for s in cl]
        mean_z = sum(z_vals) / len(z_vals)
        orig_cx = sum(rotated_layout[s][0] for s in cl) / len(cl)
        orig_cy = sum(rotated_layout[s][1] for s in cl) / len(cl)
        x3d_c = sum(sse_x_3d.get(s, orig_cx) for s in cl) / len(cl)

        d_pref = _sse_dom_prefix(cl[0])
        if pure_tm_oligomer and d_pref in dom_new_tm_cx:
            old_dx = orig_cx - dom_old_tm_center[d_pref][0]
            old_dy = orig_cy - dom_old_tm_center[d_pref][1]
            target_cx = dom_new_tm_cx[d_pref] + 1.10 * old_dx
            side_sign = 1 if mean_z >= 0.0 else -1
            if side_sign > 0:
                raw_cy = 0.45 * max(0.0, mean_z) + 0.35 * abs(old_dy)
                upper_clusters.append((cl, target_cx, raw_cy, orig_cx, orig_cy))
            else:
                raw_cy = 0.45 * min(0.0, mean_z) - 0.35 * abs(old_dy)
                lower_clusters.append((cl, target_cx, raw_cy, orig_cx, orig_cy))
        else:
            # Preserve 90% of 2DProts horizontal unfolding expanded by 1.38x so extramembrane
            # domains spread generously across the open horizontal space outside the membrane
            # while keeping natural 1:1 vertical height so the square viewport stays zoomed in!
            target_cx = 1.38 * (0.90 * orig_cx + 0.10 * x3d_c)
            side_sign = 1 if (mean_z > 2.0 or (abs(mean_z) <= 2.0 and orig_cy >= 0.0)) else -1
            if side_sign > 0:
                raw_cy = 0.45 * max(0.0, mean_z) + 0.55 * orig_cy
                upper_clusters.append((cl, target_cx, raw_cy, orig_cx, orig_cy))
            else:
                raw_cy = 0.45 * min(0.0, mean_z) + 0.55 * orig_cy
                lower_clusters.append((cl, target_cx, raw_cy, orig_cx, orig_cy))
        movable_clusters.append(cl)
        cluster_side[tuple(cl)] = side_sign

    # Place upper extramembrane clusters and shift them as a group so the lowest SSE tip clears tm_top_y + 2.2
    if upper_clusters:
        for cl, t_cx, r_cy, o_cx, o_cy in upper_clusters:
            for s in cl:
                dx = (rotated_layout[s][0] - o_cx) * 1.25
                dy = (rotated_layout[s][1] - o_cy) * 1.00
                new_layout[s] = [float(t_cx + dx), float(r_cy + dy)]
        min_upper_tip = min(
            new_layout[s][1] - abs(math.cos(math.radians(new_angles.get(s, 0.0)))) * _sse_half_len_2d(s, sizes)
            for cl, _, _, _, _ in upper_clusters
            for s in cl
        )
        shift_up = (tm_top_y + 2.2) - min_upper_tip
        for cl, _, _, _, _ in upper_clusters:
            for s in cl:
                new_layout[s][1] = float(new_layout[s][1] + shift_up)

    # Place lower extramembrane clusters and shift them as a group so the highest SSE tip clears tm_bot_y - 2.2
    if lower_clusters:
        for cl, t_cx, r_cy, o_cx, o_cy in lower_clusters:
            for s in cl:
                dx = (rotated_layout[s][0] - o_cx) * 1.25
                dy = (rotated_layout[s][1] - o_cy) * 1.00
                new_layout[s] = [float(t_cx + dx), float(r_cy + dy)]
        max_lower_tip = max(
            new_layout[s][1] + abs(math.cos(math.radians(new_angles.get(s, 0.0)))) * _sse_half_len_2d(s, sizes)
            for cl, _, _, _, _ in lower_clusters
            for s in cl
        )
        shift_dn = (tm_bot_y - 2.2) - max_lower_tip
        for cl, _, _, _, _ in lower_clusters:
            for s in cl:
                new_layout[s][1] = float(new_layout[s][1] + shift_dn)

    # Resolve any remaining segment-segment overlaps among extramembrane SSEs while keeping TM clearance
    if movable_clusters:
        _resolve_extramembrane_overlaps(
            new_layout,
            new_angles,
            sizes,
            fixed_sses=tm_set,
            movable_clusters=movable_clusters,
            cluster_side=cluster_side,
            top_bound_y=tm_top_y + 2.0,
            bot_bound_y=tm_bot_y - 2.0,
            min_sep=2.95,
        )

    # 5. Localize ligand and channel 2D coordinates accurately relative to SSEs and membrane
    localize_ligands_2d(
        ligands,
        old_layout=rotated_layout,
        new_layout=new_layout,
        new_angles=new_angles,
        sizes=sizes,
        mem_info=mem_info,
        pre_transform=rot2d,
    )

    if channels:
        localize_channels_2d(
            channels,
            layout_2d=new_layout,
            angles_2d=new_angles,
            sizes=sizes,
            vectors_3d=vectors_3d,
            h_type=h_type,
            ligands=ligands,
            mem_info=mem_info,
            old_layout=rotated_layout,
            pre_transform=rot2d,
        )

    mem_info["y_center_2d"] = 0.0
    return new_layout, new_angles, mem_info


MEMBRANE_LIPID_CODES = {
    "LI1", "SQU", "PC1", "CDL", "PAL", "LDM", "PLM", "CLR", "OLA", "LFA",
    "BNG", "C8E", "LDA", "PEF", "POP", "POG", "LHG", "DGD", "MGD", "BOG",
    "SDS", "LMT", "HTO", "HTG", "MYR", "STE", "SPC", "PEK", "PGT", "cardiolipin",
}

MEMBRANE_LIPID_KEYWORDS = (
    "LIPID", "PHOSPHATIDYL", "CARDIOLIPIN", "SQUALENE", "PALMIT",
    "CHOLESTEROL", "DIACYLGLYCEROL", "PHYTANYL", "OCTYL", "DODECYL",
    "GLUCOSIDE", "MALTOSIDE", "SPHINGO", "CERAMIDE",
)

PORE_ION_CODES = {"K", "NA", "CA", "CL", "BR", "CS", "RB", "MG"}

GLYCAN_CODES = {"NAG", "MAN", "BMA", "GAL", "FUC", "SIA", "NDG", "GLC", "XYP"}


def _is_membrane_lipid(ecmp, edsc):
    if ecmp in MEMBRANE_LIPID_CODES:
        return True
    for kw in MEMBRANE_LIPID_KEYWORDS:
        if kw in edsc:
            return True
    return False


def localize_ligands_2d(ligands, old_layout, new_layout, new_angles=None, sizes=None,
                        mem_info=None, pre_transform=None):
    """
    Re-localizes 2D ligand coordinates so that:
      - Unscaled 3D Angstrom offsets (8-18 units in raw layout.json) are compressed
        to true 2D canvas scale (~1.2-2.4 canvas units from coordinating SSEs).
      - Internal pocket cofactors/inhibitors (e.g. RET, HEM, CXT) sit directly inside
        their coordinating helix bundle pocket.
      - Membrane-associated lipids (e.g. LI1, SQU, PC1, CDL, PAL) sit inside the
        horizontal membrane band (upper/lower leaflet) alongside their coordinating TM helices.
      - Channel/pore ions (e.g. K+ in KcsA) align along the central transmembrane pore
        between the pore-lining helices.
      - Pairwise ligand-ligand and ligand-SSE center collisions are smoothly resolved.
    """
    if not ligands or not old_layout or not new_layout:
        return

    valid_sses = [s for s in new_layout if s in old_layout]
    if not valid_sses:
        return

    tm_set = set(mem_info.get("tm_helices", [])) if mem_info else set()
    band_half = float(mem_info.get("band_half_height_2d", 3.0)) if mem_info else 3.0
    sse_z_3d = mem_info.get("sse_z_3d", {}) if mem_info else {}

    # Identify intramembrane SSEs (TM helices + pore helices inside the membrane band)
    intramem_sses = set(tm_set)
    if mem_info:
        for s in valid_sses:
            if abs(new_layout[s][1]) <= band_half + 0.5 and abs(sse_z_3d.get(s, 99.0)) <= 14.0:
                intramem_sses.add(s)

    lig_keys = [k for k in ligands if "layout" in ligands[k] and ligands[k]["layout"] is not None]
    if not lig_keys:
        return

    old_lig_pts = {}
    for k in lig_keys:
        pt = ligands[k]["layout"]
        old_lig_pts[k] = list(pre_transform(pt)) if pre_transform else [float(pt[0]), float(pt[1])]

    target_pts = {}
    lig_role = {}      # "pocket", "pore_ion", "lipid", "glycan", "surface"
    pore_ion_keys = []

    for k in lig_keys:
        l_obj = ligands[k]
        ecmp = (l_obj.get("ecmp") or "").upper().strip()
        edsc = (l_obj.get("edsc") or "").upper().strip()
        opt = old_lig_pts[k]

        # Distances to all SSEs in old layout
        dists = []
        for s in valid_sses:
            sp = old_layout[s]
            d = math.hypot(opt[0] - sp[0], opt[1] - sp[1])
            dists.append((d, s))
        dists.sort(key=lambda t: t[0])

        d_min, s_near = dists[0]
        # Coordinating SSE neighborhood (up to 5 closest SSEs within local cutoff)
        cutoff = max(d_min * 1.55, d_min + 4.2)
        coord_sses = [(d, s) for d, s in dists[:5] if d <= cutoff]
        if not coord_sses:
            coord_sses = [dists[0]]

        # Exponentially decayed inverse-distance weights
        weights = []
        for d, s in coord_sses:
            w = math.exp(-((d - d_min) ** 2) / (2.0 * (2.6 ** 2))) / ((d + 1.0) ** 2)
            weights.append((w, s, d))
        w_sum = sum(w for w, _, _ in weights) or 1.0

        anchor_old_x = sum(w * old_layout[s][0] for w, s, _ in weights) / w_sum
        anchor_old_y = sum(w * old_layout[s][1] for w, s, _ in weights) / w_sum
        anchor_new_x = sum(w * new_layout[s][0] for w, s, _ in weights) / w_sum
        anchor_new_y = sum(w * new_layout[s][1] for w, s, _ in weights) / w_sum

        dx_anchor = opt[0] - anchor_old_x
        dy_anchor = opt[1] - anchor_old_y
        r_anchor = math.hypot(dx_anchor, dy_anchor)

        dx_near = opt[0] - old_layout[s_near][0]
        dy_near = opt[1] - old_layout[s_near][1]
        r_near = math.hypot(dx_near, dy_near)

        # Fraction of coordinating weight on intramembrane SSEs
        f_intramem = sum(w for w, s, _ in weights if s in intramem_sses) / w_sum
        f_tm = sum(w for w, s, _ in weights if s in tm_set) / w_sum

        # Check if ligand sits inside a multi-SSE binding pocket
        is_enclosed_pocket = (
            len(coord_sses) >= 2
            and r_anchor < max(2.2, d_min * 0.82)
            and d_min < 7.0
            and ecmp not in MEMBRANE_LIPID_CODES
            and ecmp not in GLYCAN_CODES
        )

        if mem_info is not None:
            # 1. Channel / pore ions coordinated by intramembrane helices
            if ecmp in PORE_ION_CODES and f_intramem >= 0.6:
                lig_role[k] = "pore_ion"
                pore_ion_keys.append(k)
                target_pts[k] = [anchor_new_x, anchor_new_y]
                continue

            # 2. Internal transmembrane binding pocket cofactor / inhibitor (e.g. RET, CXT)
            if is_enclosed_pocket and f_intramem >= 0.5:
                lig_role[k] = "pocket"
                px = anchor_new_x + 0.35 * dx_anchor
                py = anchor_new_y + 0.30 * dy_anchor
                py = max(-band_half + 0.65, min(band_half - 0.65, py))
                target_pts[k] = [float(px), float(py)]
                continue

            # 3. Membrane lipids / detergents (or peripheral TM-bound hydrophobic ligands)
            if _is_membrane_lipid(ecmp, edsc) or (f_tm >= 0.65 and ecmp not in GLYCAN_CODES and l_obj.get("atoms", 0) >= 8):
                lig_role[k] = "lipid"
                # Anchor horizontally to the closest 1-2 TM helices
                tm_coord = [(w, s) for w, s, _ in weights if s in intramem_sses][:2]
                if not tm_coord:
                    tm_coord = [(1.0, s_near)]
                tw_sum = sum(w for w, _ in tm_coord) or 1.0
                tm_new_x = sum(w * new_layout[s][0] for w, s in tm_coord) / tw_sum
                tm_new_y = sum(w * new_layout[s][1] for w, s in tm_coord) / tw_sum
                tm_old_x = sum(w * old_layout[s][0] for w, s in tm_coord) / tw_sum
                tm_old_y = sum(w * old_layout[s][1] for w, s in tm_coord) / tw_sum

                ldx = opt[0] - tm_old_x
                ldy = opt[1] - tm_old_y
                sign_x = 1.0 if (ldx + 0.25 * ldy) >= 0.0 else -1.0
                lat_mag = min(1.85, max(0.95, 0.28 * abs(ldx) + 0.75))
                px = tm_new_x + sign_x * lat_mag

                # Map old radial y-offset into upper or lower leaflet inside the membrane band
                leaflet_y = tm_new_y + (band_half - 0.65) * math.tanh(ldy / 5.5)
                py = max(-band_half + 0.35, min(band_half - 0.35, leaflet_y))
                target_pts[k] = [float(px), float(py)]
                continue

        # 4. Globular internal pocket ligand (e.g. HEM in Cytochrome P450)
        if is_enclosed_pocket:
            lig_role[k] = "pocket"
            px = anchor_new_x + 0.45 * dx_anchor
            py = anchor_new_y + 0.45 * dy_anchor
            target_pts[k] = [float(px), float(py)]
            continue

        # 5. Peripheral / extramembrane / glycan / surface ion ligand:
        # Compress inflated 3D Angstrom distance to clean 2D canvas scale around nearest SSE(s)
        lig_role[k] = "glycan" if ecmp in GLYCAN_CODES else "surface"
        r_scaled = 1.15 + 1.25 * math.tanh(max(0.0, r_near - 1.2) / 6.0)
        if r_near > 1e-4:
            ux = dx_near / r_near
            uy = dy_near / r_near
        else:
            ux, uy = 1.0, 0.0

        # Blend 70% nearest SSE + 30% weighted neighborhood anchor
        base_x = 0.70 * new_layout[s_near][0] + 0.30 * anchor_new_x
        base_y = 0.70 * new_layout[s_near][1] + 0.30 * anchor_new_y
        px = base_x + ux * r_scaled
        py = base_y + uy * r_scaled

        if mem_info is not None and s_near not in intramem_sses:
            # Keep extramembrane ligands on the same side of the membrane as their coordinating SSE
            if new_layout[s_near][1] > band_half:
                py = max(band_half + 0.8, py)
            elif new_layout[s_near][1] < -band_half:
                py = min(-band_half - 0.8, py)

        target_pts[k] = [float(px), float(py)]

    # Special handling for transmembrane pore ions (e.g. K+ selectivity filter in 1bl8A00):
    # Place them in a clean vertical single-file along the central channel axis between pore & TM helices
    if pore_ion_keys and mem_info is not None:
        intramem_list = sorted(intramem_sses, key=lambda s: new_layout[s][0])
        if len(intramem_list) >= 2:
            # Place pore axis in the open channel gap between the inner TM helix and pore helix
            x_right = new_layout[intramem_list[-1]][0]
            x_mid = new_layout[intramem_list[-2]][0]
            pore_x = 0.52 * x_mid + 0.48 * x_right
        else:
            pore_x = 0.0
        n_ions = len(pore_ion_keys)
        for idx, k in enumerate(sorted(pore_ion_keys)):
            if n_ions == 1:
                py = 0.4
            else:
                # Stagger vertically with >=1.75 canvas units spacing so K labels never touch circles
                frac = idx / float(n_ions - 1)
                py = -1.55 + frac * 3.50
            target_pts[k] = [float(pore_x), float(py)]

    # Precompute 2D line segments for all SSEs so ligands never sit on top of helix bars or labels
    sse_segments = {}
    for s in valid_sses:
        sx, sy = new_layout[s]
        ang_rad = math.radians(new_angles.get(s, 0.0)) if new_angles else 0.0
        hl = float(sizes.get(s, 10.0)) * 0.1 if sizes else 1.2
        dx_s = math.sin(ang_rad) * hl
        dy_s = math.cos(ang_rad) * hl
        sse_segments[s] = (sx - dx_s, sy - dy_s, sx + dx_s, sy + dy_s, sx, sy)

    # Approximate visual circle radius in canvas units for each ligand
    lig_radii = {}
    for k in lig_keys:
        atoms = max(1, int(ligands[k].get("atoms", 5)))
        lig_radii[k] = min(0.75, max(0.28, math.sqrt(atoms) * 0.095))

    for k in lig_keys:
        ligands[k]["layout"] = [float(target_pts[k][0]), float(target_pts[k][1])]
        ligands[k]["_role"] = lig_role.get(k, "surface")

    deoverlap_ligands_from_sses(ligands, new_layout, new_angles, sizes, mem_info)


def deoverlap_ligands_from_sses(ligands, layout, angles, sizes, mem_info=None):
    """
    Ensures ligand circles and their upper text labels do not collide with SSE bars,
    SSE center labels, or other ligands in the final resolved 2D layout.
    """
    if not ligands or not layout:
        return

    lig_keys = [k for k in ligands if "layout" in ligands[k] and ligands[k]["layout"] is not None]
    if not lig_keys:
        return

    band_half = float(mem_info.get("band_half_height_2d", 3.0)) if mem_info else 3.0

    sse_segments = {}
    for s, (sx, sy) in layout.items():
        ang_rad = math.radians(angles.get(s, 0.0)) if angles else 0.0
        hl = float(sizes.get(s, 10.0)) * 0.1 if sizes else 1.2
        dx_s = math.sin(ang_rad) * hl
        dy_s = math.cos(ang_rad) * hl
        sse_segments[s] = (sx - dx_s, sy - dy_s, sx + dx_s, sy + dy_s, sx, sy)

    lig_radii = {}
    for k in lig_keys:
        atoms = max(1, int(ligands[k].get("atoms", 5)))
        lig_radii[k] = min(0.75, max(0.30, math.sqrt(atoms) * 0.10 + 0.10))

    target_pts = {k: [float(ligands[k]["layout"][0]), float(ligands[k]["layout"][1])] for k in lig_keys}
    sorted_keys = sorted(lig_keys)
    n_lig = len(sorted_keys)

    def seg_closest(px, py, ax, ay, bx, by, sx, sy):
        abx = bx - ax
        aby = by - ay
        ab2 = abx * abx + aby * aby
        if ab2 > 1e-6:
            t = max(0.0, min(1.0, ((px - ax) * abx + (py - ay) * aby) / ab2))
            return ax + t * abx, ay + t * aby
        return sx, sy

    for _ in range(65):
        # 1. Push ligands (both circle center and upper text label) clear of SSE segments and SSE labels
        for k in sorted_keys:
            role = ligands[k].get("_role", "surface")
            r_k = lig_radii[k]
            min_seg_d = (0.90 + r_k) if role in ("pocket", "pore_ion") else (1.18 + r_k)
            min_center_d = 1.35 + r_k
            px, py = target_pts[k]

            for s, (ax, ay, bx, by, sx, sy) in sse_segments.items():
                qx, qy = seg_closest(px, py, ax, ay, bx, by, sx, sy)
                ddx = px - qx
                ddy = py - qy
                dist = math.hypot(ddx, ddy)
                if dist < min_seg_d:
                    if dist > 1e-4:
                        push = (min_seg_d - dist) * 0.65
                        if role == "pore_ion":
                            px += (1.0 if ddx >= 0 else -1.0) * push
                        elif abs(by - ay) > abs(bx - ax):
                            # Primarily vertical/tilted helix: push horizontally off the helix bar
                            sign_x = 1.0 if ddx >= 0 else -1.0
                            px += sign_x * push * 0.85 + (ddx / dist) * push * 0.25
                            py += (ddy / dist) * push * 0.25
                        else:
                            px += (ddx / dist) * push
                            py += (ddy / dist) * push
                    else:
                        px += min_seg_d * 0.7

                # Also keep the text label right above the ligand circle clear of the SSE bar
                lbl_y = py + r_k + 0.45
                lqx, lqy = seg_closest(px, lbl_y, ax, ay, bx, by, sx, sy)
                ldx = px - lqx
                ldy = lbl_y - lqy
                ldist = math.hypot(ldx, ldy)
                if ldist < 0.95 and role != "pore_ion":
                    if ldist > 1e-4:
                        lpush = (0.95 - ldist) * 0.50
                        px += (ldx / ldist) * lpush
                        py += (ldy / ldist) * lpush

                # Extra clearance around SSE midpoint where the 'H..'/'E..' label is drawn
                cdx = px - sx
                cdy = py - sy
                cdist = math.hypot(cdx, cdy)
                if cdist < min_center_d and role != "pore_ion":
                    if cdist > 1e-4:
                        cpush = (min_center_d - cdist) * 0.55
                        px += (cdx / cdist) * cpush
                        py += (cdy / cdist) * cpush

                # Also keep the ligand's upper text label (px, py + r_k + 0.32) clear of the SSE center label (sx, sy - 0.14)
                lbl_cdx = px - sx
                lbl_cdy = (py + r_k + 0.32) - (sy - 0.14)
                lbl_cdist = math.hypot(lbl_cdx, lbl_cdy)
                if lbl_cdist < 1.38 and role != "pore_ion":
                    if lbl_cdist > 1e-4:
                        lcpush = (1.38 - lbl_cdist) * 0.55
                        px += (lbl_cdx / lbl_cdist) * lcpush * 0.65
                        py += (lbl_cdy / lbl_cdist) * lcpush * 0.85

            target_pts[k] = [px, py]

        # 2. Pairwise repulsion between nearby ligands (accounting for circle radii + upper text label)
        for i in range(n_lig):
            k1 = sorted_keys[i]
            r1 = ligands[k1].get("_role", "surface")
            for j in range(i + 1, n_lig):
                k2 = sorted_keys[j]
                r2 = ligands[k2].get("_role", "surface")
                p1 = target_pts[k1]
                p2 = target_pts[k2]
                ddx = p1[0] - p2[0]
                ddy = p1[1] - p2[1]
                dist = math.hypot(ddx, ddy)
                r_sum = lig_radii[k1] + lig_radii[k2]
                if r1 == "pore_ion" and r2 == "pore_ion":
                    min_d = 1.75
                else:
                    vert_factor = 0.45 if abs(ddx) < 1.45 else 0.0
                    min_d = r_sum + 1.32 + vert_factor

                if dist < min_d:
                    if dist > 1e-4:
                        ux = ddx / dist
                        uy = ddy / dist
                    else:
                        angle = (i * 1.7 + j * 0.9) % (2.0 * math.pi)
                        ux = math.cos(angle)
                        uy = math.sin(angle)
                    push = (min_d - dist) * 0.50
                    if r1 == "pore_ion" and r2 == "pore_ion":
                        p1[1] += (1.0 if ddy >= 0 else -1.0) * push
                        p2[1] -= (1.0 if ddy >= 0 else -1.0) * push
                    else:
                        p1[0] += ux * push
                        p1[1] += uy * push
                        p2[0] -= ux * push
                        p2[1] -= uy * push

        # 3. Clamp membrane lipids and pocket ligands inside the membrane band
        if mem_info is not None:
            is_bb = bool(mem_info.get("is_beta_barrel"))
            for k in sorted_keys:
                role = ligands[k].get("_role", "surface")
                if role == "lipid":
                    if is_bb and -1.18 < target_pts[k][1] < 0.68:
                        # Keep beta-barrel detergent/lipid labels out of the y~0 midline where E.. labels sit
                        target_pts[k][1] = 0.72 if target_pts[k][1] >= -0.25 else -1.22
                    target_pts[k][1] = max(-band_half + 0.32, min(band_half - 0.32, target_pts[k][1]))
                elif role == "pocket":
                    target_pts[k][1] = max(-band_half + 0.55, min(band_half - 0.55, target_pts[k][1]))
                elif role == "pore_ion":
                    target_pts[k][1] = max(-band_half + 0.15, min(band_half - 0.15, target_pts[k][1]))

    for k in sorted_keys:
        ligands[k]["layout"] = [float(target_pts[k][0]), float(target_pts[k][1])]



def superimpose_membrane_layouts(layout_all_lists):
    """
    Superimpose multiple membrane-aligned layouts across a protein family while
    strictly preserving the horizontal membrane plane (y = 0) and 15-degree quantized angles.
    Only tests membrane-preserving reflections R = diag(sx, sy) with sx, sy in {-1, +1}
    and horizontal translation tx (ty = 0).
    """
    n_obj = len(layout_all_lists)
    if n_obj == 0:
        return []
    if n_obj == 1:
        return [(np.eye(2), np.zeros((2, 1)))]

    # Pick consensus reference as the layout with the most points
    i_ref = max(range(n_obj), key=lambda i: len(layout_all_lists[i]))
    ref_dict = dict(layout_all_lists[i_ref])

    rot_trans = [(np.eye(2), np.zeros((2, 1))) for _ in range(n_obj)]
    candidate_signs = [(1.0, 1.0), (-1.0, 1.0), (1.0, -1.0), (-1.0, -1.0)]

    for _ in range(3):
        accum_pts = {}
        accum_cnt = {}
        for i in range(n_obj):
            obj = layout_all_lists[i]
            common = [k for k in obj if k in ref_dict]
            if len(common) < 2:
                R = np.eye(2)
                t = np.zeros((2, 1))
            else:
                A = np.array([obj[k] for k in common], dtype=float)
                B = np.array([ref_dict[k] for k in common], dtype=float)
                best_err = float("inf")
                best_R = np.eye(2)
                best_tx = 0.0
                for sx, sy in candidate_signs:
                    tx = float(np.mean(B[:, 0] - sx * A[:, 0]))
                    dx = sx * A[:, 0] + tx - B[:, 0]
                    dy = sy * A[:, 1] - B[:, 1]
                    err = float(np.mean(dx * dx + dy * dy))
                    if err < best_err:
                        best_err = err
                        best_R = np.array([[sx, 0.0], [0.0, sy]], dtype=float)
                        best_tx = tx
                R = best_R
                t = np.array([[best_tx], [0.0]], dtype=float)

            rot_trans[i] = (R, t)
            for k, pt in obj.items():
                tx_pt = R[0, 0] * pt[0] + t[0, 0]
                ty_pt = R[1, 1] * pt[1]
                if k not in accum_pts:
                    accum_pts[k] = [tx_pt, ty_pt]
                    accum_cnt[k] = 1
                else:
                    accum_pts[k][0] += tx_pt
                    accum_pts[k][1] += ty_pt
                    accum_cnt[k] += 1

        ref_dict = {
            k: [accum_pts[k][0] / accum_cnt[k], accum_pts[k][1] / accum_cnt[k]]
            for k in accum_pts
        }

    return rot_trans


TUNNEL_PALETTE = [
    "#00838f",  # Teal / Cyan (Solvent / Water / Pore)
    "#1565c0",  # Royal Blue (Channel 2a / Ion channel)
    "#7b1fa2",  # Violet / Purple (Channel 2b / Substrate access)
    "#ef6c00",  # Amber / Orange (Channel 2e)
    "#ad1457",  # Rose / Magenta (Channel 2c)
    "#2e7d32",  # Emerald Green (Channel 5)
    "#4527a0",  # Deep Indigo
    "#00695c",  # Dark Teal
]

COFACTOR_CODES = {
    "HEM", "HEC", "HEA", "HEB", "HEO", "RET", "FAD", "FMN",
    "NAP", "NDP", "NAD", "NAI", "SAM", "SAH", "PLP", "TPP", "COB",
}


def format_tunnel_label(raw_name, idx=0, is_tm_pore=False, has_pore_ions=False):
    """Format a compact, human-readable label for a protein tunnel/channel."""
    if not raw_name:
        if is_tm_pore:
            return "Ion channel" if has_pore_ions else "TM pore"
        return f"Ch. {idx + 1}"
    s = str(raw_name).strip()
    low = s.lower()
    if "solvent" in low:
        return "Solvent ch."
    if "water" in low:
        return "Water ch."
    if "substrate" in low:
        return "Substrate ch."
    if "ion" in low and "channel" in low:
        return "Ion channel"
    if "pore" in low or low == "path":
        return "Ion channel" if has_pore_ions else "TM pore"
    if low.startswith("channel "):
        suffix = s[8:].strip()
        return f"Ch. {suffix}"
    # Check if purely numeric (legacy unannotated MOLE length ID like '23.56' or Caver ID '1')
    if s.replace(".", "").isdigit():
        if is_tm_pore:
            return "Ion channel" if has_pore_ions else "TM pore"
        return f"Ch. {idx + 1}"
    if len(s) > 16:
        return s[:14] + ".."
    return s


def pick_tunnel_color(label, idx=0):
    """Assign a consistent, high-contrast biological color to each tunnel."""
    low = (label or "").lower()
    if "solvent" in low or "pore" in low:
        return "#00838f"
    if "water" in low:
        return "#0097a7"
    if "2a" in low or "ion" in low:
        return "#1565c0"
    if "2b" in low or "substrate" in low:
        return "#7b1fa2"
    if "2e" in low:
        return "#ef6c00"
    if "2c" in low:
        return "#ad1457"
    if "ch. 5" in low or low.endswith(" 5"):
        return "#2e7d32"
    return TUNNEL_PALETTE[idx % len(TUNNEL_PALETTE)]


def _split_channel_key(key):
    """Split a channel point key like '23.561A0004' into (prefix, index)."""
    d = -1
    for i, c in enumerate(key):
        if c.isalpha():
            d = i
            break
    if d == -1:
        return key, 0
    prefix = key[:d]
    suffix = key[d + 1 :]
    try:
        idx = int(suffix)
    except ValueError:
        idx = 0
    return prefix, idx


def localize_channels_2d(channels, layout_2d, angles_2d=None, sizes=None,
                         vectors_3d=None, h_type=None, ligands=None,
                         mem_info=None, old_layout=None, pre_transform=None):
    """
    Deduplicates, names, and localizes 3D tunnel profiles (`channels`) into the 2D layout:
      1. Splits any concatenated tunnel jumps (> 14 A in 3D) and filters redundant duplicate
         tunnels across ChannelsDB calculation methods (prioritizing curated MOLE channels).
      2. Maps 3D tunnel coordinates onto 2D SSE segments using inverse-distance segment
         interpolation + damped 3D-to-2D Jacobian projection.
      3. For transmembrane pores (`mem_info` present and vertical span >= 15 A), projects
         vertical height directly from the 3D membrane normal (y = 0.18 * z_3d) and centers
         lateral x along the pore ion axis or central inter-helix gap.
      4. For active-site cofactor tunnels (e.g. P450 tunnels starting at HEM/HEC), anchors
         the tunnel origin at the cofactor's 2D center, fans out closely overlapping exit
         trajectories, and clears SSE bars at the tunnel mouth.
      5. Smooths each tunnel trajectory so it renders as a clean, kink-free biological pathway.
    """
    if not channels or not layout_2d:
        return

    # 1. Group channel points into ordered tunnel chains
    raw_groups = {}
    for k, val in list(channels.items()):
        if not isinstance(val, dict) or "layout" not in val:
            continue
        fn = val.get("filename", "")
        prefix, pt_idx = _split_channel_key(k)
        raw_groups.setdefault((fn, prefix), []).append((pt_idx, k, val))

    if not raw_groups:
        return

    candidate_tunnels = []
    for (fn, prefix), items in raw_groups.items():
        items.sort(key=lambda x: x[0])
        # Truncate if there is a discontinuous 3D jump (> 14 A) from a legacy ID collision
        clean_items = [items[0]]
        for it in items[1:]:
            p_prev = clean_items[-1][2].get("3d_layout")
            p_curr = it[2].get("3d_layout")
            if p_prev is not None and p_curr is not None:
                dist3 = math.dist(p_prev, p_curr)
                if dist3 > 14.5:
                    break
            clean_items.append(it)

        if len(clean_items) < 2:
            continue

        pts_3d = [
            np.array(it[2]["3d_layout"], dtype=float)
            for it in clean_items
            if it[2].get("3d_layout") is not None and len(it[2]["3d_layout"]) == 3
        ]
        has_3d = len(pts_3d) == len(clean_items)
        len_3d = 0.0
        if has_3d:
            pts_arr = np.array(pts_3d, dtype=float)
            len_3d = float(np.sum(np.linalg.norm(np.diff(pts_arr, axis=0), axis=1)))
        else:
            pts_arr = None

        first_val = clean_items[0][2]
        cid = str(first_val.get("id", prefix))
        tname = str(first_val.get("text_names", cid)).strip()
        is_descriptive = bool(
            tname
            and not tname.replace(".", "").isdigit()
            and not tname.lower().startswith("tunnel")
            and tname.lower() not in ("path", "none")
        )
        is_mole = ("." in cid) or is_descriptive
        mean_rad = float(np.mean([it[2].get("radius", 1.5) for it in clean_items]))

        candidate_tunnels.append({
            "fn": fn,
            "prefix": prefix,
            "items": clean_items,
            "has_3d": has_3d,
            "pts_3d": pts_arr,
            "len_3d": len_3d,
            "cid": cid,
            "raw_name": tname,
            "is_descriptive": is_descriptive,
            "is_mole": is_mole,
            "mean_rad": mean_rad,
        })

    if not candidate_tunnels:
        return

    # 2. Deduplicate tunnels per file (prioritizing descriptive / MOLE tunnels over Caver fragments)
    by_fn = {}
    for ct in candidate_tunnels:
        by_fn.setdefault(ct["fn"], []).append(ct)

    kept_tunnels = []
    for fn, t_list in by_fn.items():
        has_mole = any(t["is_mole"] and t["len_3d"] >= 8.0 for t in t_list)
        # Sort priority: descriptive > MOLE > Caver, then by length/radius
        t_list.sort(
            key=lambda t: (
                0 if t["is_descriptive"] else (1 if t["is_mole"] else 2),
                -round(t["mean_rad"], 1),
                t["len_3d"] if t["len_3d"] > 0 else 25.0,
            )
        )
        fn_kept = []
        for t in t_list:
            if has_mole and not t["is_mole"]:
                continue
            if t["has_3d"] and t["len_3d"] < 7.5 and len(t_list) > 1:
                continue
            if t["has_3d"] and t["len_3d"] > 65.0 and not t["is_descriptive"] and len(fn_kept) >= 1:
                continue

            is_dup = False
            if t["has_3d"]:
                p_start = t["pts_3d"][0]
                p_end = t["pts_3d"][-1]
                for kt in fn_kept:
                    if not kt["has_3d"]:
                        continue
                    if t["is_descriptive"] and t["raw_name"] == kt["raw_name"]:
                        is_dup = True
                        break
                    k_start = kt["pts_3d"][0]
                    k_end = kt["pts_3d"][-1]
                    if np.linalg.norm(p_end - k_end) < 6.2 and np.linalg.norm(p_start - k_start) < 6.8:
                        is_dup = True
                        break
                    # Check if exit is close and trajectory mostly overlaps
                    if np.linalg.norm(p_end - k_end) < 7.5:
                        dists = [float(np.min(np.linalg.norm(kt["pts_3d"] - p, axis=1))) for p in t["pts_3d"]]
                        if np.percentile(dists, 80) < 5.0:
                            is_dup = True
                            break
            if not is_dup and len(fn_kept) < 6:
                fn_kept.append(t)

        if not fn_kept and t_list:
            fn_kept.append(t_list[0])
        kept_tunnels.extend(fn_kept)

    # Rebuild `channels` dictionary in-place with only the kept tunnels
    kept_keys_order = []
    kept_map = {}
    for t in kept_tunnels:
        for _, k, val in t["items"]:
            kept_keys_order.append(k)
            kept_map[k] = val

    for k in list(channels.keys()):
        if k not in kept_map:
            del channels[k]

    # 3. Build 3D and 2D SSE segment geometry for 3D->2D mapping and 2D obstacle avoidance
    import heapq
    import utils_sses_types

    segs = {}
    vis_segs_2d = []
    pts3d_list, pts2d_list = [], []
    for s, coord2d in layout_2d.items():
        if h_type and sizes:
            if utils_sses_types.sse_should_be_shown(s, h_type, sizes) != 0 and len(layout_2d) > 6:
                continue
        ang_deg = angles_2d.get(s, 0.0) if angles_2d else 0.0
        th = math.radians(ang_deg)
        sz = sizes.get(s, 8) if sizes else 8
        hl = 0.1 * float(sz)
        cx, cy = float(coord2d[0]), float(coord2d[1])
        # Match Image.secstr: angle_start is at 180+angle (-sin, -cos), angle_target is at angle (+sin, +cos)
        a2 = np.array([cx - hl * math.sin(th), cy - hl * math.cos(th)], dtype=float)
        b2 = np.array([cx + hl * math.sin(th), cy + hl * math.cos(th)], dtype=float)
        c2_pt = np.array([cx, cy], dtype=float)
        vis_segs_2d.append((s, a2, b2, c2_pt, hl))
        if vectors_3d and s in vectors_3d and vectors_3d[s] and len(vectors_3d[s]) >= 2:
            a3 = np.array(vectors_3d[s][0], dtype=float)
            b3 = np.array(vectors_3d[s][1], dtype=float)
            segs[s] = (a3, b3, a2, b2, c2_pt, hl)
            pts3d_list.extend([a3, b3, 0.5 * (a3 + b3)])
            pts2d_list.extend([a2, b2, c2_pt])

    has_3d_map = len(segs) >= 2
    if has_3d_map:
        pts3d_arr = np.array(pts3d_list, dtype=float)
        pts2d_arr = np.array(pts2d_list, dtype=float)
        c3 = pts3d_arr.mean(axis=0)
        c2 = pts2d_arr.mean(axis=0)
        X = pts3d_arr - c3
        Y = pts2d_arr - c2
        J_T, _, _, _ = np.linalg.lstsq(X.T @ X + 40.0 * np.eye(3), X.T @ Y, rcond=None)
        J = J_T.T
        U, S, Vt = np.linalg.svd(J, full_matrices=False)
        S = np.clip(S, 0.22, 0.36)
        J = U @ np.diag(S) @ Vt
    else:
        J = np.array([[0.30, 0.0, 0.0], [0.0, 0.30, 0.0]], dtype=float)

    def _pt_to_seg(p, a, b):
        u = b - a
        L2 = float(np.dot(u, u))
        if L2 < 1e-9:
            return float(np.linalg.norm(p - a)), 0.0, a
        t_p = float(np.clip(np.dot(p - a, u) / L2, 0.0, 1.0))
        q = a + t_p * u
        return float(np.linalg.norm(p - q)), t_p, q

    def _segs_intersect(p0, p1, q0, q1):
        def _ccw(A, B, C):
            return (C[1] - A[1]) * (B[0] - A[0]) > (B[1] - A[1]) * (C[0] - A[0])
        return (_ccw(p0, q0, q1) != _ccw(p1, q0, q1)) and (_ccw(p0, p1, q0) != _ccw(p0, p1, q1))

    def _seg_to_seg_dist(p0, p1, q0, q1):
        if _segs_intersect(p0, p1, q0, q1):
            return 0.0
        d1, _, _ = _pt_to_seg(p0, q0, q1)
        d2, _, _ = _pt_to_seg(p1, q0, q1)
        d3, _, _ = _pt_to_seg(q0, p0, p1)
        d4, _, _ = _pt_to_seg(q1, p0, p1)
        return min(d1, d2, d3, d4)

    def _convex_hull_2d(pts):
        u_pts = sorted(set((round(float(p[0]), 4), round(float(p[1]), 4)) for p in pts))
        if len(u_pts) <= 2:
            return [np.array(p, dtype=float) for p in u_pts]
        def _cross(o, a, b):
            return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
        lower = []
        for p in u_pts:
            while len(lower) >= 2 and _cross(lower[-2], lower[-1], p) <= 0:
                lower.pop()
            lower.append(p)
        upper = []
        for p in reversed(u_pts):
            while len(upper) >= 2 and _cross(upper[-2], upper[-1], p) <= 0:
                upper.pop()
            upper.append(p)
        return [np.array(p, dtype=float) for p in (lower[:-1] + upper[:-1])]

    def _ray_hull_distance(origin, u_dir, hull):
        best_t = 0.0
        nh = len(hull)
        for i in range(nh):
            a = hull[i]
            b = hull[(i + 1) % nh]
            v = b - a
            det = -u_dir[0] * v[1] + u_dir[1] * v[0]
            if abs(det) < 1e-8:
                continue
            rhs = a - origin
            t_r = (-rhs[0] * v[1] + rhs[1] * v[0]) / det
            s_r = (u_dir[0] * rhs[1] - u_dir[1] * rhs[0]) / det
            if t_r > 0 and -1e-4 <= s_r <= 1.0001 and t_r > best_t:
                best_t = float(t_r)
        return best_t

    def map_3d_point_to_2d(p3):
        items = []
        for s, (a3, b3, a2, b2, _, _) in segs.items():
            u3 = b3 - a3
            t_p = float(np.clip(np.dot(p3 - a3, u3) / (np.dot(u3, u3) + 1e-9), 0.0, 1.0))
            q3 = a3 + t_p * u3
            q2 = a2 + t_p * (b2 - a2)
            dist = float(np.linalg.norm(p3 - q3))
            est2 = q2 + 0.65 * (J @ (p3 - q3))
            items.append((dist, est2, s))
        items.sort(key=lambda x: x[0])
        top = items[: min(4, len(items))]
        ws = np.array([1.0 / ((it[0] + 1.8) ** 3.4) for it in top], dtype=float)
        ws /= ws.sum()
        res = sum(w * it[1] for w, it in zip(ws, top))
        return res, items[0][0]

    # Pre-sample all visible 2D SSE segments for fast boundary/envelope queries
    all_sse_sample_pts = []
    for _, a2, b2, _, _ in vis_segs_2d:
        for alpha in np.linspace(0.0, 1.0, 9):
            all_sse_sample_pts.append(a2 + alpha * (b2 - a2))
    hull_2d = _convex_hull_2d(all_sse_sample_pts) if all_sse_sample_pts else []

    # Identify active-site cofactor ligand (e.g. HEM, HEC, RET) and pore ions (e.g. K+) if present
    cofactor_2d = None
    pore_ion_pts = []
    if ligands:
        best_cof_atoms = -1
        for lk, lval in ligands.items():
            if not isinstance(lval, dict) or "layout" not in lval or lval["layout"] is None:
                continue
            ecmp = str(lval.get("ecmp") or "").upper().strip()
            edsc = str(lval.get("edsc") or "").upper().strip()
            role = lval.get("_role", "")
            atoms = int(lval.get("atoms", 1) or 1)
            lpt = np.array(lval["layout"], dtype=float)
            if role == "pore_ion" or ecmp in PORE_ION_CODES:
                pore_ion_pts.append(lpt)
            elif ecmp in COFACTOR_CODES or "HEME" in edsc or "RETINAL" in edsc or "FLAVIN" in edsc:
                if atoms > best_cof_atoms:
                    best_cof_atoms = atoms
                    cofactor_2d = lpt
            elif role == "pocket" and atoms >= 15 and not _is_membrane_lipid(ecmp, edsc):
                if atoms > best_cof_atoms:
                    best_cof_atoms = atoms
                    cofactor_2d = lpt

    # 4. Initial classification and natural direction extraction
    hub_tunnels = []
    for t_idx, t in enumerate(kept_tunnels):
        n_pts = len(t["items"])
        if t["has_3d"] and has_3d_map:
            pts_3d = t["pts_3d"]
            mapped = [map_3d_point_to_2d(p3) for p3 in pts_3d]
            raw_2d = np.array([m[0] for m in mapped], dtype=float)
            start_min_d3 = mapped[0][1]

            is_tm_pore = False
            if mem_info and "normal_3d" in mem_info and "center_3d" in mem_info:
                c_mem = np.array(mem_info["center_3d"], dtype=float)
                n_mem = np.array(mem_info["normal_3d"], dtype=float)
                z_prof = np.dot(pts_3d - c_mem, -n_mem)
                z_span = float(np.max(z_prof) - np.min(z_prof))
                z_min_abs = float(np.min(np.abs(z_prof)))
                if z_span >= 15.0 and z_min_abs <= 12.0:
                    is_tm_pore = True

            t["is_tm_pore"] = is_tm_pore
            label = format_tunnel_label(
                t["raw_name"], idx=t_idx, is_tm_pore=is_tm_pore, has_pore_ions=bool(pore_ion_pts)
            )
            color = pick_tunnel_color(label, idx=t_idx)
            t["label"] = label
            t["color"] = color

            if is_tm_pore:
                c_mem = np.array(mem_info["center_3d"], dtype=float)
                n_mem = np.array(mem_info["normal_3d"], dtype=float)
                y_center_2d = float(mem_info.get("y_center_2d", 0.0))
                band_hh = float(mem_info.get("band_half_height_2d", 3.2))
                z_prof = np.dot(pts_3d - c_mem, -n_mem)

                mean_seg_x = float(np.mean(raw_2d[:, 0]))
                if pore_ion_pts:
                    ion_arr = np.array(pore_ion_pts, dtype=float)
                    x_anchor = float(np.mean(ion_arr[:, 0]))
                else:
                    tm_hs = [h for h in mem_info.get("tm_helices", []) if h in layout_2d]
                    tm_xs = sorted(float(layout_2d[h][0]) for h in tm_hs)
                    if len(tm_xs) >= 2:
                        gaps = [0.5 * (tm_xs[i] + tm_xs[i + 1]) for i in range(len(tm_xs) - 1) if tm_xs[i + 1] - tm_xs[i] >= 3.0]
                        x_anchor = min(gaps, key=lambda g: abs(g - mean_seg_x)) if gaps else float(np.mean(tm_xs))
                    else:
                        x_anchor = mean_seg_x

                # Determine full vertical extent of the structure near the pore so both ends exit outside
                y_min_struct = y_center_2d - band_hh
                y_max_struct = y_center_2d + band_hh
                for _, a2, b2, c2_pt, _ in vis_segs_2d:
                    if abs(c2_pt[0] - x_anchor) <= 7.5 or abs(a2[0] - x_anchor) <= 6.0 or abs(b2[0] - x_anchor) <= 6.0:
                        y_min_struct = min(y_min_struct, float(a2[1]), float(b2[1]))
                        y_max_struct = max(y_max_struct, float(a2[1]), float(b2[1]))

                y_bot_ext = y_min_struct - 1.35
                y_top_ext = y_max_struct + 1.35
                if float(z_prof[-1]) >= float(z_prof[0]):
                    y_prof = np.linspace(y_bot_ext, y_top_ext, n_pts)
                else:
                    y_prof = np.linspace(y_top_ext, y_bot_ext, n_pts)

                # Route x_prof at each y_k through the open horizontal corridor between visible SSE bars
                def _open_x_corridor_at_y(yk, ref_x):
                    blocked = []
                    for _, a2, b2, _, _ in vis_segs_2d:
                        ya, yb = float(a2[1]), float(b2[1])
                        xa, xb = float(a2[0]), float(b2[0])
                        ymin_s, ymax_s = min(ya, yb), max(ya, yb)
                        if yk < ymin_s - 0.85 or yk > ymax_s + 0.85:
                            continue
                        # Find sub-segment of [a2, b2] within |y - yk| <= 1.05
                        if abs(yb - ya) < 0.8:
                            x_lo, x_hi = min(xa, xb), max(xa, xb)
                        else:
                            t0 = np.clip(((yk - 1.05) - ya) / (yb - ya), 0.0, 1.0)
                            t1 = np.clip(((yk + 1.05) - ya) / (yb - ya), 0.0, 1.0)
                            x0 = xa + t0 * (xb - xa)
                            x1 = xa + t1 * (xb - xa)
                            x_lo, x_hi = min(x0, x1), max(x0, x1)
                        blocked.append((float(x_lo) - 1.05, float(x_hi) + 1.05))
                    if not blocked:
                        return ref_x
                    blocked.sort()
                    merged = [list(blocked[0])]
                    for b_lo, b_hi in blocked[1:]:
                        if b_lo <= merged[-1][1] + 0.15:
                            merged[-1][1] = max(merged[-1][1], b_hi)
                        else:
                            merged.append([b_lo, b_hi])
                    # Build open intervals between [-30, +30]
                    open_ints = []
                    prev_hi = -35.0
                    for b_lo, b_hi in merged:
                        if b_lo - prev_hi >= 0.6:
                            open_ints.append((prev_hi, b_lo))
                        prev_hi = b_hi
                    open_ints.append((prev_hi, 35.0))
                    # Pick open interval closest to ref_x (with preference for wider intervals near x_anchor)
                    best_x = ref_x
                    best_score = 1e9
                    for o_lo, o_hi in open_ints:
                        w_gap = o_hi - o_lo
                        if o_lo <= ref_x <= o_hi:
                            cand_x = 0.5 * (o_lo + o_hi) if w_gap < 3.6 else np.clip(ref_x, o_lo + 0.6, o_hi - 0.6)
                        else:
                            cand_x = 0.5 * (o_lo + o_hi) if w_gap < 3.6 else (o_hi - 0.6 if ref_x > o_hi else o_lo + 0.6)
                        score = abs(cand_x - ref_x) + 0.45 * abs(cand_x - x_anchor)
                        if score < best_score:
                            best_score = score
                            best_x = float(cand_x)
                    return best_x

                x_prof = np.zeros(n_pts, dtype=float)
                mid_k = int(np.argmin(np.abs(y_prof - y_center_2d)))
                x_prof[mid_k] = _open_x_corridor_at_y(y_prof[mid_k], x_anchor)
                for k in range(mid_k + 1, n_pts):
                    x_prof[k] = _open_x_corridor_at_y(y_prof[k], 0.7 * x_prof[k - 1] + 0.3 * x_anchor)
                for k in range(mid_k - 1, -1, -1):
                    x_prof[k] = _open_x_corridor_at_y(y_prof[k], 0.7 * x_prof[k + 1] + 0.3 * x_anchor)

                # Smooth x_prof while respecting SSE clearance and pore ion alignment
                traj_2d = np.column_stack([x_prof, y_prof])
                for _ in range(24):
                    new_x = traj_2d[:, 0].copy()
                    for k in range(1, n_pts - 1):
                        new_x[k] = 0.25 * traj_2d[k - 1, 0] + 0.50 * traj_2d[k, 0] + 0.25 * traj_2d[k + 1, 0]
                    traj_2d[:, 0] = new_x
                    for k in range(n_pts):
                        pk = traj_2d[k]
                        if pore_ion_pts:
                            for ion_p in pore_ion_pts:
                                if abs(pk[1] - ion_p[1]) < 1.8:
                                    pk[0] = 0.65 * pk[0] + 0.35 * float(ion_p[0])
                        for _, a2, b2, _, _ in vis_segs_2d:
                            d_s, _, q_s = _pt_to_seg(pk, a2, b2)
                            if d_s < 1.32:
                                dx_s = pk[0] - q_s[0]
                                sign_x = 1.0 if dx_s >= 0 else -1.0
                                pk[0] += sign_x * (1.32 - d_s) * 0.75
                        traj_2d[k] = pk

                t["traj_2d"] = traj_2d
            else:
                # Globular / Active-site tunnel
                if cofactor_2d is not None and start_min_d3 <= 9.5:
                    origin_2d = cofactor_2d.copy()
                    t["is_hub_tunnel"] = True
                else:
                    origin_2d = raw_2d[0].copy()
                    t["is_hub_tunnel"] = False

                shift0 = origin_2d - raw_2d[0]
                guide_2d = np.zeros_like(raw_2d)
                p3_0 = pts_3d[0]
                for k in range(n_pts):
                    tau = float(k) / max(1, n_pts - 1)
                    seg_anchored = raw_2d[k] + ((1.0 - tau) ** 1.15) * shift0
                    direct_3d = origin_2d + 1.25 * (J @ (pts_3d[k] - p3_0))
                    blend = 0.35 * tau
                    guide_2d[k] = (1.0 - blend) * seg_anchored + blend * direct_3d

                # Natural exit direction from outer half of 3D-mapped trajectory + 3D exit vector
                outer_vecs = guide_2d[max(1, n_pts // 2) :] - origin_2d
                mean_outer = np.mean(outer_vecs, axis=0)
                end_v = guide_2d[-1] - origin_2d
                jac_v = J @ (pts_3d[-1] - p3_0)
                dir_v = np.zeros(2, dtype=float)
                if float(np.linalg.norm(end_v)) > 1e-3:
                    dir_v += 0.45 * (end_v / np.linalg.norm(end_v))
                if float(np.linalg.norm(mean_outer)) > 1e-3:
                    dir_v += 0.35 * (mean_outer / np.linalg.norm(mean_outer))
                if float(np.linalg.norm(jac_v)) > 1e-3:
                    dir_v += 0.20 * (jac_v / np.linalg.norm(jac_v))
                nat_angle = math.atan2(dir_v[1], dir_v[0]) if float(np.linalg.norm(dir_v)) > 1e-4 else 0.0

                t["origin_2d"] = origin_2d
                t["guide_2d"] = guide_2d
                t["nat_angle"] = nat_angle
                t["traj_2d"] = guide_2d.copy()
                if t.get("is_hub_tunnel"):
                    hub_tunnels.append(t)
        else:
            pts_2d = []
            for _, _, val in t["items"]:
                pt = val["layout"]
                if pre_transform is not None and not val.get("_localized_2d"):
                    pt = pre_transform(pt)
                pts_2d.append(np.array(pt, dtype=float))
            traj_2d = np.array(pts_2d, dtype=float)
            label = format_tunnel_label(t["raw_name"], idx=t_idx, is_tm_pore=False)
            color = pick_tunnel_color(label, idx=t_idx)
            t["label"] = label
            t["color"] = color
            t["is_tm_pore"] = False
            t["traj_2d"] = traj_2d

    # 5. De-crowd hub tunnel exit angles and route each non-TM tunnel through inter-SSE corridors to outside the structure
    if len(hub_tunnels) >= 2:
        angles = [ht["nat_angle"] for ht in hub_tunnels]
        min_sep = math.radians(34.0)
        for _ in range(35):
            for i in range(len(hub_tunnels)):
                for j in range(i + 1, len(hub_tunnels)):
                    diff = (angles[i] - angles[j] + math.pi) % (2.0 * math.pi) - math.pi
                    if abs(diff) < min_sep:
                        push = 0.5 * (min_sep - abs(diff))
                        sign = 1.0 if diff >= 0 else -1.0
                        angles[i] += sign * push
                        angles[j] -= sign * push
        for ht, ang in zip(hub_tunnels, angles):
            ht["target_angle"] = ang
    elif len(hub_tunnels) == 1:
        hub_tunnels[0]["target_angle"] = hub_tunnels[0]["nat_angle"]

    def _compute_exterior_exit_pt(o2, theta):
        u_dir = np.array([math.cos(theta), math.sin(theta)], dtype=float)
        n_dir = np.array([-u_dir[1], u_dir[0]], dtype=float)
        r_cone = 0.0
        sector_pars = []
        for q in all_sse_sample_pts:
            rel = q - o2
            d_par = float(np.dot(rel, u_dir))
            if d_par <= 0:
                continue
            d_perp = abs(float(np.dot(rel, n_dir)))
            ang_q = math.atan2(rel[1], rel[0])
            da = abs((ang_q - theta + math.pi) % (2.0 * math.pi) - math.pi)
            if da <= math.radians(24.0) or d_perp <= 3.8:
                r_cone = max(r_cone, d_par)
            if da <= math.radians(42.0):
                sector_pars.append(d_par)
        r_sector_wide = float(np.percentile(sector_pars, 85)) if sector_pars else 10.5
        r_ext = max(r_cone, 0.85 * r_sector_wide, 9.5) + 2.0
        p_ext = o2 + r_ext * u_dir
        # Ensure p_ext is at least 1.85 units clear of every SSE segment
        for _ in range(10):
            for _, a2, b2, _, _ in vis_segs_2d:
                d_s, _, q_s = _pt_to_seg(p_ext, a2, b2)
                if d_s < 1.85:
                    diff_v = p_ext - q_s
                    dn = float(np.linalg.norm(diff_v))
                    p_ext += (diff_v / dn if dn > 1e-4 else u_dir) * (1.85 - d_s)
        return p_ext

    # Pre-generate candidate corridor graph nodes in the open spaces between SSE bars
    corridor_nodes_base = []
    for s, a2, b2, c2_pt, hl in vis_segs_2d:
        u_s = b2 - a2
        un = float(np.linalg.norm(u_s))
        t_hat = u_s / un if un > 1e-5 else np.array([0.0, 1.0])
        n_hat = np.array([-t_hat[1], t_hat[0]])
        for d_tip in (1.45, 2.35):
            corridor_nodes_base.append(a2 - d_tip * t_hat)
            corridor_nodes_base.append(b2 + d_tip * t_hat)
        for sgn in (-1.0, 1.0):
            corridor_nodes_base.append(a2 - 1.05 * t_hat + sgn * 1.55 * n_hat)
            corridor_nodes_base.append(b2 + 1.05 * t_hat + sgn * 1.55 * n_hat)
            corridor_nodes_base.append(c2_pt + sgn * 1.75 * n_hat)

    # Add midpoints of gaps between nearby SSE endpoints/centers
    sse_key_pts = []
    for _, a2, b2, c2_pt, _ in vis_segs_2d:
        sse_key_pts.extend([a2, b2, c2_pt])
    for i_p in range(len(sse_key_pts)):
        for j_p in range(i_p + 1, len(sse_key_pts)):
            dist_ij = float(np.linalg.norm(sse_key_pts[i_p] - sse_key_pts[j_p]))
            if 2.2 <= dist_ij <= 9.5:
                corridor_nodes_base.append(0.5 * (sse_key_pts[i_p] + sse_key_pts[j_p]))

    # Keep only candidate nodes that lie in open space (>= 1.25 units from all visible SSE bars)
    valid_base_nodes = []
    for cand in corridor_nodes_base:
        if all(_pt_to_seg(cand, a2, b2)[0] >= 1.25 for _, a2, b2, _, _ in vis_segs_2d):
            valid_base_nodes.append(cand)

    routed_tunnel_pts = []
    for t in kept_tunnels:
        if t.get("is_tm_pore") or not (t.get("has_3d") and has_3d_map):
            continue
        o2 = t["origin_2d"]
        base_theta = t.get("target_angle", t.get("nat_angle", 0.0))

        # Refine exit angle within +-14 deg to maximize corridor clearance while keeping separation from other hub tunnels
        best_theta = base_theta
        best_ang_score = 1e9
        other_thetas = [ht.get("target_angle", ht.get("nat_angle", 0.0)) for ht in hub_tunnels if ht is not t]
        for d_deg in np.linspace(-14.0, 14.0, 15):
            cand_th = base_theta + math.radians(float(d_deg))
            if any(abs((cand_th - oth + math.pi) % (2.0 * math.pi) - math.pi) < math.radians(27.0) for oth in other_thetas):
                continue
            cand_ext = _compute_exterior_exit_pt(o2, cand_th)
            # Score ray from o2 to cand_ext by how many SSE bar mid-sections it crosses
            ray_pen = 0.15 * abs(d_deg)
            for _, a2, b2, _, _ in vis_segs_2d:
                if _pt_to_seg(o2, a2, b2)[0] < 1.2:
                    continue
                d_ray = _seg_to_seg_dist(o2, cand_ext, a2, b2)
                if d_ray < 0.2:
                    ray_pen += 12.0
                elif d_ray < 1.5:
                    ray_pen += 4.5 * (1.5 - d_ray)
            if ray_pen < best_ang_score:
                best_ang_score = ray_pen
                best_theta = cand_th

        t["target_angle"] = best_theta
        p_exit = _compute_exterior_exit_pt(o2, best_theta)
        u_exit = p_exit - o2
        R_total = float(np.linalg.norm(u_exit))
        u_hat = u_exit / max(1e-5, R_total)
        n_hat = np.array([-u_hat[1], u_hat[0]])

        # Add radial guide nodes pushed outside SSE bars
        nodes = [o2, p_exit]
        for frac in np.linspace(0.14, 0.88, 8):
            for lat in (0.0, -2.0, 2.0, -4.0, 4.0, -6.0, 6.0):
                gp = o2 + (frac * R_total) * u_hat + lat * math.sin(math.pi * frac) * n_hat
                if all(_pt_to_seg(gp, a2, b2)[0] >= 1.25 for _, a2, b2, _, _ in vis_segs_2d):
                    nodes.append(gp)
        for bn in valid_base_nodes:
            rel = bn - o2
            d_par = float(np.dot(rel, u_hat))
            d_perp = abs(float(np.dot(rel, n_hat)))
            if -1.5 <= d_par <= R_total + 2.0 and d_perp <= max(10.5, 0.58 * R_total):
                nodes.append(bn)

        # Identify pocket SSEs immediately adjacent to o2 so the first step out of the pocket isn't deadlocked
        pocket_sses = {s for s, a2, b2, _, _ in vis_segs_2d if _pt_to_seg(o2, a2, b2)[0] < 2.25}

        n_nodes = len(nodes)
        dist_arr = [1e18] * n_nodes
        parent = [-1] * n_nodes
        dist_arr[0] = 0.0
        pq = [(0.0, 0)]

        while pq:
            c_d, u_idx = heapq.heappop(pq)
            if c_d > dist_arr[u_idx] + 1e-6:
                continue
            if u_idx == 1:
                break
            pu = nodes[u_idx]
            ru = float(np.dot(pu - o2, u_hat))
            for v_idx in range(1, n_nodes):
                if v_idx == u_idx:
                    continue
                pv = nodes[v_idx]
                step_len = float(np.linalg.norm(pv - pu))
                if step_len > 12.5 or step_len < 0.2:
                    continue
                rv = float(np.dot(pv - o2, u_hat))
                # Prevent backward zig-zags toward origin
                if rv < ru - 0.8:
                    continue

                # Check clearance of segment [pu, pv] from all visible SSE bars
                seg_pen = 0.0
                p_start_chk = pu + 0.18 * (pv - pu) if u_idx == 0 else pu
                for s, a2, b2, _, _ in vis_segs_2d:
                    d_seg = _seg_to_seg_dist(p_start_chk, pv, a2, b2)
                    if d_seg < 0.55:
                        if u_idx == 0 and s in pocket_sses and _pt_to_seg(pv, a2, b2)[0] >= 1.25:
                            seg_pen += 85.0 * (0.65 - d_seg)
                        else:
                            seg_pen += 25000.0 * (0.75 - d_seg)
                    elif d_seg < 2.25:
                        seg_pen += 38.0 * ((2.25 - d_seg) ** 2)

                # Lateral deviation from tunnel's target corridor ray + preference for multi-hop corridor steps
                mid_uv = 0.5 * (pu + pv)
                lat_dist = abs(float(np.dot(mid_uv - o2, n_hat)))
                dev_pen = 0.24 * lat_dist * step_len + 0.05 * (step_len ** 2)

                # Repulsion from already routed tunnels so branches use distinct inter-SSE gaps
                tun_pen = 0.0
                for prev_traj in routed_tunnel_pts:
                    d_prev = float(np.min(np.linalg.norm(prev_traj - mid_uv, axis=1)))
                    if d_prev < 2.4:
                        tun_pen += 18.0 * ((2.4 - d_prev) ** 2)

                nd = c_d + step_len + seg_pen + dev_pen + tun_pen
                if nd < dist_arr[v_idx]:
                    dist_arr[v_idx] = nd
                    parent[v_idx] = u_idx
                    heapq.heappush(pq, (nd, v_idx))

        # Reconstruct polygonal path
        poly = []
        curr = 1
        while curr != -1:
            poly.append(nodes[curr])
            curr = parent[curr]
        poly.reverse()
        if len(poly) < 2:
            poly = [o2, p_exit]

        # Resample polygonal path into n_pts evenly spaced points
        n_pts = len(t["items"])
        poly_arr = np.array(poly, dtype=float)
        seg_lens = np.linalg.norm(np.diff(poly_arr, axis=0), axis=1)
        cum_lens = np.concatenate([[0.0], np.cumsum(seg_lens)])
        tot_len = max(1e-5, float(cum_lens[-1]))
        target_s = np.linspace(0.0, tot_len, n_pts)
        traj = np.zeros((n_pts, 2), dtype=float)
        traj[:, 0] = np.interp(target_s, cum_lens, poly_arr[:, 0])
        traj[:, 1] = np.interp(target_s, cum_lens, poly_arr[:, 1])

        # Corridor-centered elastic relaxation: balances gentle curve smoothing with soft repulsion
        # from nearby SSE bars so the tunnel glides smoothly down the center of each inter-SSE corridor
        # instead of pulling taut against SSE corners.
        for it_smooth in range(28):
            new_traj = traj.copy()
            for k in range(1, n_pts - 1):
                # 5-point neighborhood blend for graceful, wide-radius curvature
                k_m2 = max(0, k - 2)
                k_p2 = min(n_pts - 1, k + 2)
                sm_pt = (
                    0.12 * traj[k_m2]
                    + 0.24 * traj[k - 1]
                    + 0.28 * traj[k]
                    + 0.24 * traj[k + 1]
                    + 0.12 * traj[k_p2]
                )
                # Soft corridor-centering push away from any SSE bar within 2.05 units
                rep_vec = np.zeros(2, dtype=float)
                for s, a2, b2, _, _ in vis_segs_2d:
                    if k <= 2 and s in pocket_sses:
                        continue
                    d_s, _, q_s = _pt_to_seg(sm_pt, a2, b2)
                    if d_s < 2.05:
                        diff_v = sm_pt - q_s
                        dn = float(np.linalg.norm(diff_v))
                        if dn > 1e-4:
                            rep_vec += (diff_v / dn) * ((2.05 - d_s) * 0.45)
                cand = sm_pt + rep_vec
                # Hard safety floor at 1.38 units from non-pocket SSEs
                for s, a2, b2, _, _ in vis_segs_2d:
                    if k <= 2 and s in pocket_sses:
                        continue
                    d_s, _, q_s = _pt_to_seg(cand, a2, b2)
                    if d_s < 1.38:
                        diff_v = cand - q_s
                        dn = float(np.linalg.norm(diff_v))
                        if dn > 1e-4:
                            cand = q_s + (diff_v / dn) * 1.38
                # Forbid smoothing from jumping across an SSE bar
                crosses = False
                for s, a2, b2, _, _ in vis_segs_2d:
                    if k <= 2 and s in pocket_sses:
                        continue
                    if _segs_intersect(traj[k], cand, a2, b2):
                        crosses = True
                        break
                    d_old = min(_seg_to_seg_dist(new_traj[k - 1], traj[k], a2, b2), _seg_to_seg_dist(traj[k], traj[k + 1], a2, b2))
                    d_new = min(_seg_to_seg_dist(new_traj[k - 1], cand, a2, b2), _seg_to_seg_dist(cand, traj[k + 1], a2, b2))
                    if d_new < 0.85 and d_new < d_old - 1e-4:
                        crosses = True
                        break
                if not crosses:
                    new_traj[k] = cand
            traj = new_traj

        t["traj_2d"] = traj
        routed_tunnel_pts.append(traj[max(1, n_pts // 5) :])

    # 6. Write updated coordinates, label, and color back to `channels`
    for t in kept_tunnels:
        traj = t["traj_2d"]
        for idx_p, (_, k, val) in enumerate(t["items"]):
            val["layout"] = [float(traj[idx_p, 0]), float(traj[idx_p, 1])]
            val["text_names"] = t["label"]
            val["color"] = t["color"]
            val["_localized_2d"] = True


def compute_ghost_sses(layout_2d, angles_2d, sizes, vectors_3d, h_type=None,
                       beta_connectivity=None, mem_info=None, helices_residues=None,
                       ligands=None):
    """
    Identifies strong local 3D SSE interactions that could not fit together in the 2D layout
    and generates non-overlapping 'ghost' SSE placements (`duch sekundární struktury`):
      1. `barrel_seam`: Periodic closure of an unrolled transmembrane beta-barrel, placing
         a ghost of the C-terminal barrel strand next to the N-terminal barrel strand on the left
         and a ghost of the N-terminal barrel strand next to the C-terminal strand on the right.
      2. `barrel_loop`: Inward-folded pore-constricting loop SSEs in beta-barrels (e.g. Loop L3
         helix H31 in 16-stranded porins 2porA00 and 2omfA00) that fold into the barrel lumen in 3D
         to contact central barrel wall strands (E42/E52/E54), shown as a docked ghost with a
         subtle folding arc from the loop's sequence position.
      3. `local_contact`: Split beta-sheets or strong non-sequential 3D SSE contacts (d3D <= 5.2 A)
         that were separated in 2D (d2D >= 7.5 canvas units), such as split beta-hairpins E44:E107
         and E97:E108 in Cytochrome P450 1tqnA00 or the extracellular beta-plug in GPCR 1f88A00.
    """
    if not layout_2d or not angles_2d or not sizes or not vectors_3d:
        return []

    vis_sses = [
        s for s in layout_2d
        if s in vectors_3d and s in angles_2d and s in sizes
        and (not h_type or utils_sses_types.sse_should_be_shown(s, h_type, sizes) == 0)
    ]
    if len(vis_sses) < 2:
        return []

    def _seg2d(coord, ang_deg, sz):
        hl = max(0.45, float(sz) * 0.10)
        rad = math.radians(float(ang_deg))
        v = np.array([math.sin(rad) * hl, math.cos(rad) * hl], dtype=float)
        c = np.array(coord, dtype=float)
        return c - v, c + v

    def _seg_dist_2d(a0, a1, b0, b1):
        ts = np.linspace(0.0, 1.0, 9)
        pa = a0[None, :] + ts[:, None] * (a1 - a0)[None, :]
        pb = b0[None, :] + ts[:, None] * (b1 - b0)[None, :]
        return float(np.min(np.linalg.norm(pa[:, None, :] - pb[None, :, :], axis=-1)))

    def _seg_dist_3d(s1, s2):
        v1, v2 = vectors_3d.get(s1), vectors_3d.get(s2)
        if not v1 or not v2 or len(v1) < 2 or len(v2) < 2:
            return 999.0
        a0, a1 = np.array(v1[0], dtype=float), np.array(v1[1], dtype=float)
        b0, b1 = np.array(v2[0], dtype=float), np.array(v2[1], dtype=float)
        ts = np.linspace(0.0, 1.0, 15)
        pa = a0[None, :] + ts[:, None] * (a1 - a0)[None, :]
        pb = b0[None, :] + ts[:, None] * (b1 - b0)[None, :]
        return float(np.min(np.linalg.norm(pa[:, None, :] - pb[None, :, :], axis=-1)))

    def _res_range(s):
        if helices_residues and s in helices_residues:
            hr = helices_residues[s]
            try:
                r0 = int(str(hr[0]).split("_")[-1])
                r1 = int(str(hr[1]).split("_")[-1])
                return min(r0, r1), max(r0, r1)
            except ValueError:
                pass
        digits = "".join(c for c in s if c.isdigit())
        val = int(digits) if digits else 0
        return val, val

    def _seq_gap(s1, s2):
        a0, a1 = _res_range(s1)
        b0, b1 = _res_range(s2)
        if a1 < b0:
            return b0 - a1
        if b1 < a0:
            return a0 - b1
        return 0

    existing_segs = {s: _seg2d(layout_2d[s], angles_2d[s], sizes[s]) for s in vis_sses}
    seq_ordered_vis = sorted(vis_sses, key=lambda x: _res_range(x)[0])
    connector_segs = []
    for k_vis in range(len(seq_ordered_vis) - 1):
        s_a, s_b = seq_ordered_vis[k_vis], seq_ordered_vis[k_vis + 1]
        _, a_end = existing_segs[s_a]
        b_start, _ = existing_segs[s_b]
        connector_segs.append((s_a, s_b, a_end, b_start))
    placed_ghosts = []

    def _ghost_clearance(g_coord, g_ang, g_sz, ignore_sses=()):
        ga0, ga1 = _seg2d(g_coord, g_ang, g_sz)
        min_d = 999.0
        for s, (sa0, sa1) in existing_segs.items():
            if s in ignore_sses:
                continue
            d = _seg_dist_2d(ga0, ga1, sa0, sa1)
            if d < min_d:
                min_d = d
        for s_a, s_b, ca0, ca1 in connector_segs:
            # Trim 20% from each end of the connector segment so endpoints at SSE tips don't false-trigger
            c_trim0 = ca0 + 0.20 * (ca1 - ca0)
            c_trim1 = ca0 + 0.80 * (ca1 - ca0)
            d_conn = _seg_dist_2d(ga0, ga1, c_trim0, c_trim1) + 0.55
            if d_conn < min_d:
                min_d = d_conn
        for pg in placed_ghosts:
            pa0, pa1 = _seg2d(pg["coord"], pg["angle"], pg["size"])
            d = _seg_dist_2d(ga0, ga1, pa0, pa1)
            if d < min_d:
                min_d = d
        if ligands:
            for lg in ligands.values():
                lpt = lg.get("layout")
                if lpt:
                    lp = np.array(lpt, dtype=float)
                    d_lg = _seg_dist_2d(ga0, ga1, lp, lp) - 0.35
                    if d_lg < min_d:
                        min_d = d_lg
        return min_d

    is_bb = isinstance(mem_info, dict) and bool(mem_info.get("is_beta_barrel"))
    barrel_mode = mem_info.get("barrel_mode", "unrolled_30") if is_bb else ""
    barrel_set = set(mem_info.get("tm_helices", [])) if is_bb else set()

    # 1. Unrolled transmembrane beta-barrel seam closure & inward-folded pore loops
    if is_bb and not barrel_mode.startswith("arch"):
        vis_barrel = sorted(
            [s for s in barrel_set if s in vis_sses],
            key=lambda x: layout_2d[x][0],
        )
        full_barrel = [s for s in vis_barrel if float(sizes.get(s, 0.0)) >= 18.0]
        if len(full_barrel) >= 6:
            s_left = full_barrel[0]
            s_right = full_barrel[-1]
            # Verify 3D barrel closure proximity
            if _seg_dist_3d(s_left, s_right) <= 7.5:
                xs_full = [layout_2d[s][0] for s in full_barrel]
                spacing = float(np.median(np.diff(xs_full))) if len(xs_full) > 1 else 2.55
                spacing = max(2.25, min(2.85, spacing))

                # Right seam: ghost of s_left placed one column to the right of s_right
                right_most_x = max(layout_2d[s][0] for s in vis_barrel)
                gx_r = right_most_x + spacing
                gy_r = float(layout_2d[s_left][1])
                placed_ghosts.append({
                    "source_sse": s_left,
                    "partner_sse": s_right,
                    "coord": [float(gx_r), float(gy_r)],
                    "angle": float(angles_2d[s_left]),
                    "size": float(sizes[s_left]),
                    "kind": "barrel_seam",
                    "show_arc": False,
                })

                # Left seam: ghost of s_right placed one column to the left of the leftmost barrel strand
                left_most_s = vis_barrel[0]
                gx_l = layout_2d[left_most_s][0] - spacing
                gy_l = float(layout_2d[s_right][1])
                placed_ghosts.append({
                    "source_sse": s_right,
                    "partner_sse": s_left,
                    "coord": [float(gx_l), float(gy_l)],
                    "angle": float(angles_2d[s_right]),
                    "size": float(sizes[s_right]),
                    "kind": "barrel_seam",
                    "show_arc": False,
                })

        # 2. Inward-folded pore-constricting loop SSEs (e.g. H31 in 2porA00 and 2omfA00)
        sse_z_3d = mem_info.get("sse_z_3d", {})
        band_half = float(mem_info.get("band_half_height_2d", 2.8))
        barrel_top_y = max(
            [band_half] + [
                layout_2d[b][1] + abs(math.cos(math.radians(angles_2d[b]))) * 0.10 * float(sizes[b])
                for b in vis_barrel
            ]
        )
        for s_loop in vis_sses:
            if s_loop in barrel_set:
                continue
            z_loop = abs(float(sse_z_3d.get(s_loop, 99.0)))
            if z_loop > 9.5:
                continue
            # Find closest barrel wall strands in 3D that are far in 2D
            close_bar = []
            for b in full_barrel:
                d3 = _seg_dist_3d(s_loop, b)
                d2 = float(np.linalg.norm(np.array(layout_2d[s_loop]) - np.array(layout_2d[b])))
                if d3 <= 5.8 and d2 >= 7.5 and _seq_gap(s_loop, b) >= 15:
                    close_bar.append((d3, b))
            if close_bar:
                close_bar.sort(key=lambda t: t[0])
                placed_loop = False
                for _, b_cand in close_bar:
                    if placed_loop:
                        break
                    b_a0, b_a1 = existing_segs[b_cand]
                    b_top = b_a0 if b_a0[1] >= b_a1[1] else b_a1
                    for dx_try in [0.0, -0.85, 0.85]:
                        if placed_loop:
                            break
                        for dy_try in [1.25, 1.65, 2.05]:
                            c_pt = [float(b_top[0] + dx_try), float(max(b_top[1] + dy_try, barrel_top_y + 0.95))]
                            if _ghost_clearance(c_pt, angles_2d[s_loop], sizes[s_loop], ignore_sses=(b_cand,)) >= 1.25:
                                placed_ghosts.append({
                                    "source_sse": s_loop,
                                    "partner_sse": b_cand,
                                    "coord": c_pt,
                                    "angle": float(angles_2d[s_loop]),
                                    "size": float(sizes[s_loop]),
                                    "kind": "barrel_loop",
                                    "show_arc": True,
                                    "origin_coord": [float(layout_2d[s_loop][0]), float(layout_2d[s_loop][1])],
                                })
                                placed_loop = True
                                break

    # 3. General split beta-sheets and strong non-sequential 3D local contacts separated in 2D
    sse_to_cl_id = {}
    if beta_connectivity:
        for idx_cl, cl_str in enumerate(beta_connectivity):
            for s in cl_str.split(":"):
                sse_to_cl_id[s] = idx_cl

    candidate_pairs = []
    for i in range(len(vis_sses)):
        for j in range(i + 1, len(vis_sses)):
            s1, s2 = vis_sses[i], vis_sses[j]
            if s1 in barrel_set and s2 in barrel_set:
                continue
            if sse_to_cl_id.get(s1, -1) == sse_to_cl_id.get(s2, -2):
                continue
            if _seq_gap(s1, s2) < 8:
                continue
            d2 = float(np.linalg.norm(np.array(layout_2d[s1]) - np.array(layout_2d[s2])))
            if d2 < 7.5:
                continue
            d3 = _seg_dist_3d(s1, s2)
            both_sheets = utils_sses_types.is_sheet(s1) and utils_sses_types.is_sheet(s2)
            one_sheet = utils_sses_types.is_sheet(s1) or utils_sses_types.is_sheet(s2)
            if both_sheets and d3 <= 5.2:
                candidate_pairs.append((d3 - 1.5, d3, s1, s2, True))
            elif one_sheet and not is_bb and d3 <= 4.65 and d2 >= 9.0:
                candidate_pairs.append((d3, d3, s1, s2, False))

    candidate_pairs.sort(key=lambda t: t[0])
    used_cl_pairs = set()
    local_count = 0

    def _try_place_lateral_ghost(src_s, part_s):
        p_c = np.array(layout_2d[part_s], dtype=float)
        p_ang = float(angles_2d[part_s])
        p_rad = math.radians(p_ang)
        u_par = np.array([math.sin(p_rad), math.cos(p_rad)], dtype=float)
        u_perp = np.array([math.cos(p_rad), -math.sin(p_rad)], dtype=float)
        g_sz = min(float(sizes[src_s]), 18.0)
        # If both are sheets, align ghost parallel/antiparallel to partner sheet
        if utils_sses_types.is_sheet(src_s) and utils_sses_types.is_sheet(part_s):
            v1 = np.array(vectors_3d[src_s][1], dtype=float) - np.array(vectors_3d[src_s][0], dtype=float)
            v2 = np.array(vectors_3d[part_s][1], dtype=float) - np.array(vectors_3d[part_s][0], dtype=float)
            dot3 = float(np.dot(v1, v2))
            g_ang = p_ang if dot3 >= 0.0 else (p_ang + 180.0) % 360.0
        else:
            g_ang = float(angles_2d[src_s])

        best_cand = None
        best_score = -999.0
        p_hl = max(0.5, float(sizes[part_s]) * 0.10)
        for d_perp in [2.05, 2.40, 2.75]:
            for sign_p in [1.0, -1.0]:
                for frac_par in [0.0, 0.25, -0.25]:
                    c_pt = p_c + sign_p * d_perp * u_perp + (frac_par * p_hl) * u_par
                    cl_val = _ghost_clearance(c_pt, g_ang, g_sz, ignore_sses=(part_s,))
                    score = cl_val - 1.5 * abs(frac_par) - 0.35 * (d_perp - 2.05)
                    if cl_val >= 1.25 and score > best_score:
                        best_score = score
                        best_cand = ([float(c_pt[0]), float(c_pt[1])], float(g_ang), float(g_sz))
        return best_cand

    for _, d3, s1, s2, both_sheets in candidate_pairs:
        if local_count >= 3:
            break
        cl1 = sse_to_cl_id.get(s1, s1)
        cl2 = sse_to_cl_id.get(s2, s2)
        pair_key = tuple(sorted([str(cl1), str(cl2)]))
        if pair_key in used_cl_pairs:
            continue

        if both_sheets:
            # Place mutual ghosts on both split beta-sheet parts if space permits
            cand1 = _try_place_lateral_ghost(s2, s1)
            if cand1:
                placed_ghosts.append({
                    "source_sse": s2,
                    "partner_sse": s1,
                    "coord": cand1[0],
                    "angle": cand1[1],
                    "size": cand1[2],
                    "kind": "local_contact",
                    "show_arc": False,
                })
            cand2 = _try_place_lateral_ghost(s1, s2)
            if cand2:
                placed_ghosts.append({
                    "source_sse": s1,
                    "partner_sse": s2,
                    "coord": cand2[0],
                    "angle": cand2[1],
                    "size": cand2[2],
                    "kind": "local_contact",
                    "show_arc": False,
                })
            if cand1 or cand2:
                used_cl_pairs.add(pair_key)
                local_count += 1
        else:
            # Place a ghost of the smaller SSE next to the larger partner SSE
            if float(sizes[s1]) <= float(sizes[s2]):
                src_s, part_s = s1, s2
            else:
                src_s, part_s = s2, s1
            cand = _try_place_lateral_ghost(src_s, part_s)
            if cand:
                placed_ghosts.append({
                    "source_sse": src_s,
                    "partner_sse": part_s,
                    "coord": cand[0],
                    "angle": cand[1],
                    "size": cand[2],
                    "kind": "local_contact",
                    "show_arc": False,
                })
                used_cl_pairs.add(pair_key)
                local_count += 1

    return placed_ghosts




