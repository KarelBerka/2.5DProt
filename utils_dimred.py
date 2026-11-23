#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""
Pluggable Dimensionality Reduction Core for 2DProts.

Provides multiple dimensionality reduction algorithms for projecting
multi-3D protein secondary structure elements (SSEs), beta-sheet clusters,
and multi-chain/multi-domain complexes into 2D layouts:

Supported algorithms (`method`):
  - "evolution"    : Original 2DProts stochastic hill-climbing / evolutionary optimizer
                     (with fast NumPy acceleration when C binary is unavailable)
  - "pca"          : Principal Component Analysis on 3D coordinates (fallback to PCoA)
  - "pcoa" / "cmds": Classical Multidimensional Scaling (Principal Coordinates Analysis)
  - "mds"          : Metric Multidimensional Scaling (SMACOF)
  - "nmmds"        : Non-Metric Multidimensional Scaling
  - "sammon"       : Locally-Weighted MDS / Sammon Mapping tailored to the 2DProts
                     local-contact objective (< 40 A, weighted by SSE size, L-BFGS-B)
  - "isomap"       : Isometric Feature Mapping (geodesic k-NN graph distances + MDS)
  - "spectral"     : Laplacian Eigenmaps / Spectral Embedding on RBF contact graph
  - "lle"          : Locally Linear Embedding
  - "kamada_kawai" : Force-directed spring layout on the 3D distance graph
  - "tsne"         : t-Distributed Stochastic Neighbor Embedding (precomputed distances)
  - "umap"         : Uniform Manifold Approximation and Projection (precomputed distances)

Any algorithm can also be combined with local evolutionary/gradient refinement
by appending "+refine" or "+evolution" (e.g. "pca+refine", "umap+refine").
"""

import os
import math
import logging
from collections import OrderedDict
import numpy as np
from scipy.optimize import minimize
from scipy.sparse.csgraph import shortest_path
from sklearn.decomposition import PCA
from sklearn.manifold import MDS, TSNE, Isomap, SpectralEmbedding, LocallyLinearEmbedding
import networkx as nx

import superimpose2d

try:
    import umap
    HAS_UMAP = True
except ImportError:
    HAS_UMAP = False


# Global default dimensionality reduction method (can be overridden via env var TWODPROTS_DIMRED or function args)
DEFAULT_METHOD = os.environ.get("TWODPROTS_DIMRED", "evolution").strip().lower()

AVAILABLE_METHODS = [
    "evolution",
    "pca",
    "pcoa",
    "mds",
    "nmmds",
    "sammon",
    "isomap",
    "spectral",
    "lle",
    "kamada_kawai",
    "tsne",
    "umap",
]


def get_active_method(flags=None, method=None):
    """Resolve the active dimensionality reduction method from explicit arg, flags, or environment."""
    if method:
        return str(method).strip().lower()
    if isinstance(flags, dict) and "dimred" in flags and flags["dimred"]:
        return str(flags["dimred"]).strip().lower()
    return os.environ.get("TWODPROTS_DIMRED", DEFAULT_METHOD).strip().lower()


def _build_cluster_distance_matrix(clusters, sses_distance, sses_sizes, sses_positions=None, c_layout=None):
    """
    Build a symmetric NxN 3D distance matrix D and optional Nx3 3D midpoint array X3d
    for a list of clusters (where each cluster is 'SSE1' or 'SSE1:SSE2:...').
    Uses parallel-axis variance correction when `c_layout` is provided for multi-strand clusters.
    """
    n = len(clusters)
    D = np.zeros((n, n), dtype=float)
    weights = np.ones(n, dtype=float)
    cluster_sses = [c.split(":") for c in clusters]

    # Internal squared radius R_i^2 of each cluster in 2D (0 for single-SSE clusters)
    rad_sq = np.zeros(n, dtype=float)
    if c_layout:
        for i, cl_name in enumerate(clusters):
            sses_i = cluster_sses[i]
            if len(sses_i) > 1 and cl_name in c_layout:
                cl_dict = c_layout[cl_name]
                pts = [cl_dict[s] for s in sses_i if s in cl_dict]
                if len(pts) > 1:
                    arr = np.asarray(pts, dtype=float)
                    arr = arr - np.mean(arr, axis=0, keepdims=True)
                    rad_sq[i] = float(np.mean(np.sum(arr ** 2, axis=1)))

    for i, sses_i in enumerate(cluster_sses):
        w_i = sum(max(1.0, float(sses_sizes.get(s, 5.0))) for s in sses_i)
        weights[i] = w_i

    for i in range(n):
        sses_i = cluster_sses[i]
        for j in range(i + 1, n):
            sses_j = cluster_sses[j]
            num = 0.0
            den = 0.0
            for s1 in sses_i:
                w1 = max(1.0, float(sses_sizes.get(s1, 5.0)))
                d_row = sses_distance.get(s1, {})
                for s2 in sses_j:
                    if s2 in d_row:
                        w2 = max(1.0, float(sses_sizes.get(s2, 5.0)))
                        w = w1 * w2
                        num += w * (float(d_row[s2]) ** 2)
                        den += w
            dist_sq = (num / den) if den > 0.0 else 225.0
            if len(sses_i) > 1 or len(sses_j) > 1:
                dist_sq = max(2.25, dist_sq - rad_sq[i] - rad_sq[j])
            dist_ij = max(math.sqrt(dist_sq), 1.5)
            D[i, j] = dist_ij
            D[j, i] = dist_ij

    X3d = None
    if sses_positions:
        coords_3d = []
        valid_3d = True
        for sses_i in cluster_sses:
            pts = []
            for s in sses_i:
                if s in sses_positions and len(sses_positions[s]) >= 2:
                    v0 = np.asarray(sses_positions[s][0], dtype=float)
                    v1 = np.asarray(sses_positions[s][1], dtype=float)
                    pts.append(0.5 * (v0 + v1))
            if pts:
                coords_3d.append(np.mean(pts, axis=0))
            else:
                valid_3d = False
                break
        if valid_3d and len(coords_3d) == n:
            X3d = np.vstack(coords_3d)

    return D, weights, X3d


def _classical_mds(D):
    """
    Classical Multidimensional Scaling (PCoA / Torgerson) from an NxN distance matrix D.
    Deterministic and closed-form via eigendecomposition of the double-centered matrix.
    """
    n = D.shape[0]
    if n <= 1:
        return np.zeros((n, 2), dtype=float)
    if n == 2:
        d = max(2.0, float(D[0, 1]))
        return np.array([[-0.5 * d, 0.0], [0.5 * d, 0.0]], dtype=float)

    D2 = D ** 2
    H = np.eye(n) - np.ones((n, n)) / float(n)
    B = -0.5 * (H @ D2 @ H)
    B = 0.5 * (B + B.T)

    eigvals, eigvecs = np.linalg.eigh(B)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    top_vals = np.maximum(eigvals[:2], 1e-6)
    X2d = eigvecs[:, :2] * np.sqrt(top_vals)[None, :]
    return X2d


def _rescale_to_target_distances(X2d, D, weights=None):
    """
    Rescale 2D coordinates X2d so their weighted RMS local pairwise distance
    matches the target 3D distance matrix D, and orient canonically onto Classical MDS.
    """
    n = X2d.shape[0]
    if n <= 1:
        return np.zeros((n, 2), dtype=float)

    X2d = X2d - np.mean(X2d, axis=0, keepdims=True)
    diff = X2d[:, None, :] - X2d[None, :, :]
    d2d = np.sqrt(np.sum(diff ** 2, axis=-1))

    mask = np.triu(np.ones((n, n), dtype=bool), k=1)
    local_mask = mask & (D < 40.0)
    if not np.any(local_mask):
        local_mask = mask

    if weights is not None:
        W = (weights[:, None] + weights[None, :])[local_mask]
    else:
        W = np.ones(np.sum(local_mask), dtype=float)

    rms_3d = float(np.sqrt(np.sum(W * (D[local_mask] ** 2)) / max(1e-9, np.sum(W))))
    rms_2d = float(np.sqrt(np.sum(W * (d2d[local_mask] ** 2)) / max(1e-9, np.sum(W))))

    if rms_2d > 1e-6 and rms_3d > 1e-6:
        X2d = X2d * (rms_3d / rms_2d)

    # Align orientation with Classical MDS so all algorithms share a consistent reference frame
    X_ref = _classical_mds(D)
    if n >= 2 and np.linalg.norm(X_ref) > 1e-6 and np.linalg.norm(X2d) > 1e-6:
        try:
            R, t = superimpose2d.optimal_rotation_translation(X2d.T, X_ref.T, allow_mirror=True)
            X2d = (R @ X2d.T + t).T
        except Exception:
            pass

    return X2d


def _sammon_2dprots_optimize(D, weights, X_init, fixed_indices=None, fixed_coords=None,
                             start_indices=None, start_coords=None, max_iter=250):
    """
    Deterministic smooth optimization (L-BFGS-B) of the 2DProts weighted local distance stress
    (Sammon / Weighted MDS with < 40 A cutoff and optional fixed/template cluster constraints).
    """
    n = D.shape[0]
    if n <= 1:
        return X_init.copy()

    fixed_set = set(fixed_indices or [])
    movable_indices = [i for i in range(n) if i not in fixed_set]
    if not movable_indices:
        return X_init.copy()

    X_full = X_init.copy()
    if fixed_indices and fixed_coords is not None:
        for idx, coord in zip(fixed_indices, fixed_coords):
            X_full[idx] = coord

    # Precompute pair weights matching 2DProts philosophy:
    # Emphasize distances < 40 A, weighted by (size_i + size_j) / max(D_ij, 4.0)
    W_Mat = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            d3 = D[i, j]
            if d3 < 45.0:
                decay = max(0.05, 1.0 - (d3 / 50.0) ** 2)
                w = (weights[i] + weights[j]) * decay / max(d3, 4.0)
            else:
                w = 0.02 * (weights[i] + weights[j]) / max(d3, 10.0)
            W_Mat[i, j] = w
            W_Mat[j, i] = w

    has_start = bool(start_indices) and start_coords is not None
    if has_start:
        start_idx_arr = np.array(start_indices, dtype=int)
        start_pts_arr = np.asarray(start_coords, dtype=float)
        start_w_arr = np.array([weights[i] * n / 20.0 for i in start_indices], dtype=float)

    mov_arr = np.array(movable_indices, dtype=int)

    def obj_and_grad(x_flat):
        X = X_full.copy()
        X[mov_arr] = x_flat.reshape(-1, 2)

        diff = X[:, None, :] - X[None, :, :]
        dist2 = np.sum(diff ** 2, axis=-1)
        dist = np.sqrt(dist2 + 1e-4)

        # Residual between 2D distance and 3D distance
        res = dist - D
        # Extra repulsion if 2D distance is smaller than 2.2 units to prevent overlaps
        repulse_mask = (dist < 2.2) & (np.eye(n) == 0)
        loss = 0.5 * np.sum( np.triu(W_Mat * (res ** 2), k=1) )

        # Gradient w.r.t. X
        coeff = W_Mat * res / dist
        np.fill_diagonal(coeff, 0.0)

        if np.any(repulse_mask):
            pen = 2.2 - dist
            loss += 5.0 * 0.5 * np.sum(np.triu(repulse_mask * (pen ** 2), k=1))
            coeff = coeff - 10.0 * repulse_mask * pen / dist

        grad = np.sum(coeff[:, :, None] * diff, axis=1)

        if has_start:
            s_diff = X[start_idx_arr] - start_pts_arr
            s_dist = np.sqrt(np.sum(s_diff ** 2, axis=-1) + 1e-4)
            loss += np.sum(start_w_arr * s_dist)
            grad[start_idx_arr] += (start_w_arr / s_dist)[:, None] * s_diff

        return float(loss), grad[mov_arr].ravel()

    x0 = X_full[mov_arr].ravel()
    res = minimize(obj_and_grad, x0, method="L-BFGS-B", jac=True, options={"maxiter": max_iter, "ftol": 1e-6})
    X_full[mov_arr] = res.x.reshape(-1, 2)
    return X_full


def embed_distance_matrix(D, weights=None, X3d=None, method="mds", random_state=42):
    """
    Core embedding dispatcher: takes an NxN distance matrix D (and optional Nx3 coordinates X3d)
    and returns an Nx2 array of 2D coordinates in Angstrom scale.
    """
    n = D.shape[0]
    if n == 0:
        return np.zeros((0, 2), dtype=float)
    if n == 1:
        return np.zeros((1, 2), dtype=float)
    if n == 2:
        d = max(2.2, float(D[0, 1]))
        return np.array([[-0.5 * d, 0.0], [0.5 * d, 0.0]], dtype=float)

    base_method = method.split("+")[0].strip().lower()
    if weights is None:
        weights = np.ones(n, dtype=float)

    # Always compute Classical MDS as a reliable fallback / initialization
    X_pcoa = _classical_mds(D)

    try:
        if base_method == "pca":
            if X3d is not None and X3d.shape[0] == n:
                X_centered = X3d - np.mean(X3d, axis=0, keepdims=True)
                pca = PCA(n_components=2, random_state=random_state)
                X2d = pca.fit_transform(X_centered)
                return X2d
            else:
                return X_pcoa

        elif base_method in ("pcoa", "cmds", "classical_mds"):
            return X_pcoa

        elif base_method == "mds":
            mds = MDS(
                n_components=2,
                metric=True,
                dissimilarity="precomputed",
                n_init=1,
                max_iter=300,
                random_state=random_state,
            )
            X2d = mds.fit_transform(D, init=X_pcoa)
            return _rescale_to_target_distances(X2d, D, weights)

        elif base_method == "nmmds":
            mds = MDS(
                n_components=2,
                metric=False,
                dissimilarity="precomputed",
                n_init=1,
                max_iter=300,
                random_state=random_state,
            )
            X2d = mds.fit_transform(D, init=X_pcoa)
            return _rescale_to_target_distances(X2d, D, weights)

        elif base_method in ("sammon", "wmds"):
            X2d = _sammon_2dprots_optimize(D, weights, X_pcoa, max_iter=300)
            return X2d

        elif base_method == "isomap":
            k = max(2, min(6, n - 1))
            # Build k-NN graph from D and compute shortest path geodesic distances
            G_dist = np.full((n, n), np.inf, dtype=float)
            for i in range(n):
                nbrs = np.argsort(D[i])[1 : k + 1]
                for j in nbrs:
                    G_dist[i, j] = D[i, j]
                    G_dist[j, i] = D[j, i]
            geo_D = shortest_path(G_dist, directed=False)
            # Replace any disconnected components with direct distances * 1.5
            inf_mask = ~np.isfinite(geo_D)
            if np.any(inf_mask):
                geo_D[inf_mask] = D[inf_mask] * 1.5
            X2d = _classical_mds(geo_D)
            return _rescale_to_target_distances(X2d, D, weights)

        elif base_method == "spectral":
            if n <= 3:
                return X_pcoa
            med_d = float(np.median(D[D > 0])) if np.any(D > 0) else 15.0
            sigma = max(4.0, med_d * 0.85)
            affinity = np.exp(-(D ** 2) / (2.0 * (sigma ** 2)))
            np.fill_diagonal(affinity, 1.0)
            spec = SpectralEmbedding(
                n_components=2,
                affinity="precomputed",
                random_state=random_state,
            )
            X2d = spec.fit_transform(affinity)
            return _rescale_to_target_distances(X2d, D, weights)

        elif base_method == "lle":
            if n <= 4:
                return X_pcoa
            X_in = X3d if (X3d is not None and X3d.shape[0] == n) else _classical_mds_nd(D, min(3, n - 1))
            k = max(2, min(6, n - 2))
            lle = LocallyLinearEmbedding(
                n_neighbors=k,
                n_components=2,
                method="standard",
                reg=1e-2,
                random_state=random_state,
            )
            X2d = lle.fit_transform(X_in)
            return _rescale_to_target_distances(X2d, D, weights)

        elif base_method in ("kamada_kawai", "kk"):
            G = nx.Graph()
            for i in range(n):
                G.add_node(i)
            for i in range(n):
                for j in range(i + 1, n):
                    G.add_edge(i, j, weight=max(1.0, float(D[i, j])))
            init_pos = {i: X_pcoa[i] for i in range(n)}
            pos = nx.kamada_kawai_layout(G, dist={i: {j: float(D[i, j]) for j in range(n)} for i in range(n)}, pos=init_pos)
            X2d = np.array([pos[i] for i in range(n)], dtype=float)
            return _rescale_to_target_distances(X2d, D, weights)

        elif base_method == "tsne":
            if n <= 3:
                return X_pcoa
            perp = max(1.5, min(8.0, float(n - 1) / 2.2))
            tsne = TSNE(
                n_components=2,
                perplexity=perp,
                metric="precomputed",
                init=X_pcoa,
                learning_rate="auto",
                max_iter=500,
                random_state=random_state,
            )
            X2d = tsne.fit_transform(D)
            return _rescale_to_target_distances(X2d, D, weights)

        elif base_method == "umap":
            if not HAS_UMAP or n <= 4:
                return _sammon_2dprots_optimize(D, weights, X_pcoa, max_iter=150)
            n_neighbors = max(2, min(8, n - 1))
            reducer = umap.UMAP(
                n_components=2,
                n_neighbors=n_neighbors,
                min_dist=0.35,
                metric="precomputed",
                init=X_pcoa,
                random_state=random_state,
                n_epochs=250,
            )
            X2d = reducer.fit_transform(D)
            return _rescale_to_target_distances(X2d, D, weights)

        else:
            logging.warning(f"Unknown dimensionality reduction method '{method}', falling back to 'mds'")
            return embed_distance_matrix(D, weights=weights, X3d=X3d, method="mds", random_state=random_state)

    except Exception as exc:
        logging.warning(f"Dimensionality reduction '{method}' encountered an issue ({exc}); using PCoA/Sammon fallback.")
        return _sammon_2dprots_optimize(D, weights, X_pcoa, max_iter=150)


def _classical_mds_nd(D, ndim=3):
    """Helper for N-dimensional Classical MDS coordinates when 3D coords are absent."""
    n = D.shape[0]
    ndim = min(ndim, n)
    D2 = D ** 2
    H = np.eye(n) - np.ones((n, n)) / float(n)
    B = -0.5 * (H @ D2 @ H)
    B = 0.5 * (B + B.T)
    eigvals, eigvecs = np.linalg.eigh(B)
    idx = np.argsort(eigvals)[::-1]
    eigvals = np.maximum(eigvals[idx[:ndim]], 1e-6)
    return eigvecs[:, idx[:ndim]] * np.sqrt(eigvals)[None, :]


def align_embedding_to_references(clusters, X2d, layout, fixed_clusters, start_layout):
    """
    Align the newly computed 2D embedding `X2d` (for `clusters`) onto:
      1. `fixed_clusters` (if non-empty, e.g. second pass in 2-step prime/sec optimization), or
      2. `start_layout` (if provided, e.g. protein family consensus template).
    """
    cluster_to_idx = {c: i for i, c in enumerate(clusters)}

    # 1. If fixed_clusters exist, anchor strictly to fixed_clusters
    fixed_list = [c for c in fixed_clusters if c in cluster_to_idx and c in layout]
    if len(fixed_list) >= 2:
        idx_f = [cluster_to_idx[c] for c in fixed_list]
        A = X2d[idx_f].T  # 2 x k
        B = np.array([layout[c] for c in fixed_list], dtype=float).T  # 2 x k
        try:
            R, t = superimpose2d.optimal_rotation_translation(A, B, allow_mirror=True)
            X2d = (R @ X2d.T + t).T
        except Exception:
            pass
        # Keep exact coordinates for fixed_clusters
        for c in fixed_list:
            X2d[cluster_to_idx[c]] = np.asarray(layout[c], dtype=float)
        return X2d
    elif len(fixed_list) == 1:
        c0 = fixed_list[0]
        idx0 = cluster_to_idx[c0]
        shift = np.asarray(layout[c0], dtype=float) - X2d[idx0]
        X2d = X2d + shift[None, :]

    # 2. If start_layout is provided, align onto start_layout via 2D Procrustes
    if start_layout:
        ref_indices = []
        ref_targets = []
        for i, c in enumerate(clusters):
            pts = [start_layout[s] for s in c.split(":") if s in start_layout]
            if pts:
                ref_indices.append(i)
                ref_targets.append(np.mean(pts, axis=0))
        if len(ref_indices) >= 2:
            A = X2d[ref_indices].T
            B = np.array(ref_targets, dtype=float).T
            try:
                R, t = superimpose2d.optimal_rotation_translation(A, B, allow_mirror=True)
                X2d = (R @ X2d.T + t).T
            except Exception:
                pass
        elif len(ref_indices) == 1:
            idx0 = ref_indices[0]
            shift = np.asarray(ref_targets[0], dtype=float) - X2d[idx0]
            X2d = X2d + shift[None, :]

    return X2d


def compute_dimred_layout(
    layout,
    a_layout,
    c_layout,
    ca_layout,
    sses_distance,
    sses_sizes,
    movable_clusters,
    fixed_clusters,
    start_layout,
    sses_positions=None,
    flags=None,
    method="mds",
):
    """
    Drop-in replacement / alternative to `modify_layout_via_evolution` in `utils_matrix2layout.py`.

    Computes 2D cluster coordinates `layout` using the specified dimensionality reduction
    algorithm (`pca`, `pcoa`, `mds`, `nmmds`, `sammon`, `isomap`, `spectral`, `lle`,
    `kamada_kawai`, `tsne`, `umap`), aligns to `fixed_clusters` / `start_layout`,
    and optimizes non-trivial beta-sheet cluster orientations (`a_layout`, `c_layout`, `ca_layout`).
    """
    import utils_matrix2layout

    movable_list = list(movable_clusters)
    fixed_list = [c for c in fixed_clusters if c not in movable_clusters]
    all_clusters = fixed_list + movable_list

    if not movable_list:
        return layout, a_layout, c_layout, ca_layout

    movable_sses = [s for cl in movable_list for s in cl.split(":")]
    relevant_sses = [s for cl in all_clusters for s in cl.split(":")]

    D, weights, X3d = _build_cluster_distance_matrix(
        all_clusters, sses_distance, sses_sizes, sses_positions=sses_positions, c_layout=c_layout
    )

    X2d = embed_distance_matrix(D, weights=weights, X3d=X3d, method=method)
    X2d = align_embedding_to_references(all_clusters, X2d, layout, fixed_list, start_layout)

    # Optional local Sammon/gradient refinement when fixed_clusters or start_layout are active
    # or when "+refine" is requested
    base_method = method.split("+")[0].strip().lower()
    wants_refine = ("+refine" in method) or (base_method in ("sammon", "wmds") and (fixed_list or start_layout))
    if wants_refine and len(all_clusters) >= 2:
        fixed_indices = [i for i, c in enumerate(all_clusters) if c in fixed_clusters]
        fixed_coords = [layout[all_clusters[i]] for i in fixed_indices]
        start_indices = []
        start_coords = []
        if start_layout:
            for i, c in enumerate(all_clusters):
                pts = [start_layout[s] for s in c.split(":") if s in start_layout]
                if pts:
                    start_indices.append(i)
                    start_coords.append(np.mean(pts, axis=0))
        X2d = _sammon_2dprots_optimize(
            D,
            weights,
            X2d,
            fixed_indices=fixed_indices,
            fixed_coords=fixed_coords,
            start_indices=start_indices,
            start_coords=start_coords,
            max_iter=200,
        )

    for i, c in enumerate(all_clusters):
        if c in movable_clusters:
            layout[c] = (float(X2d[i, 0]), float(X2d[i, 1]))

    # Optimize rotation and reflection of non-trivial beta-sheet clusters
    a_layout, c_layout, ca_layout = utils_matrix2layout.find_best_rotation_reflection(
        movable_clusters,
        layout,
        a_layout,
        c_layout,
        ca_layout,
        sses_sizes,
        relevant_sses,
        movable_sses,
        sses_distance,
        start_layout,
        flags or {},
    )

    return layout, a_layout, c_layout, ca_layout


def compute_projected_angles(
    a_layout,
    ca_layout,
    sses_angles,
    sses_sizes,
    movable_clusters,
    fixed_clusters,
    start_angles,
    layout,
    c_layout,
    sses_positions,
    method="mds",
):
    """
    Compute 2D orientation angles `a_layout` for trivial clusters (individual helices/strands)
    by projecting their 3D direction vectors (start -> end) into the 2D embedding via
    locally-weighted 3D->2D Procrustes/SVD Jacobian, followed by fast pairwise angle refinement.
    """
    import utils_matrix2layout
    import utils_membrane

    sses_layout, _ = utils_matrix2layout.count_sse_coords_based_on_cluster_coords(
        layout, a_layout, c_layout, ca_layout
    )

    valid_sses = [
        s for s in sses_layout
        if sses_positions and s in sses_positions and len(sses_positions[s]) >= 2
    ]

    if len(valid_sses) >= 2:
        mids_3d = np.array(
            [0.5 * (np.asarray(sses_positions[s][0], dtype=float) + np.asarray(sses_positions[s][1], dtype=float)) for s in valid_sses],
            dtype=float,
        )
        pts_2d = np.array([sses_layout[s] for s in valid_sses], dtype=float)
        sse_to_idx = {s: i for i, s in enumerate(valid_sses)}

        # Pairwise 3D distances between midpoints for local neighborhood weighting
        diff_3d = mids_3d[:, None, :] - mids_3d[None, :, :]
        dist_3d = np.sqrt(np.sum(diff_3d ** 2, axis=-1))

        base_method = (method or "mds").split("+")[0].strip().lower()
        use_global = base_method in ("pca", "pcoa", "cmds")

        for cl in movable_clusters:
            parts = cl.split(":")
            if len(parts) != 1:
                continue
            s = parts[0]
            if start_angles and s in start_angles:
                a_layout[cl] = float(start_angles[s])
                continue
            if s not in sse_to_idx:
                continue

            idx_s = sse_to_idx[s]
            if use_global or len(valid_sses) <= 4:
                w = np.ones(len(valid_sses), dtype=float)
            else:
                # Gaussian neighborhood weights (sigma ~ 18 A) so local domain rotation is captured
                w = np.exp(-(dist_3d[idx_s] ** 2) / (2.0 * (18.0 ** 2)))
                w[idx_s] = 1.0

            w_sum = np.sum(w)
            if w_sum < 1e-6:
                w = np.ones(len(valid_sses), dtype=float)
                w_sum = float(len(valid_sses))

            c3 = np.sum(mids_3d * w[:, None], axis=0) / w_sum
            c2 = np.sum(pts_2d * w[:, None], axis=0) / w_sum
            M3 = (mids_3d - c3) * np.sqrt(w)[:, None]
            P2 = (pts_2d - c2) * np.sqrt(w)[:, None]

            # 3x2 cross-covariance matrix H = M3^T @ P2
            H = M3.T @ P2
            U, _, Vt = np.linalg.svd(H, full_matrices=False)
            proj_3to2 = U @ Vt  # 3x2 orthogonal projection matrix

            vec_3d = np.asarray(sses_positions[s][1], dtype=float) - np.asarray(sses_positions[s][0], dtype=float)
            vec_2d = vec_3d @ proj_3to2
            if np.linalg.norm(vec_2d) > 1e-6:
                ang = math.degrees(math.atan2(float(vec_2d[0]), float(vec_2d[1]))) % 360.0
                a_layout[cl] = utils_membrane.snap_angle(ang)

    # Snap all angles cleanly
    for k in a_layout:
        a_layout[k] = utils_membrane.snap_angle(a_layout[k])

    return a_layout, ca_layout



def compute_protein_dimred_layout(
    layout,
    a_layout,
    c_layout,
    ca_layout,
    sses_distance,
    sses_sizes,
    method="mds",
):
    """
    Multi-chain / multi-domain protein assembly using dimensionality reduction
    (alternative to `modify_protein_layout_via_evolution` in `utils_matrix2layout.py`).

    1. Builds the inter-domain/inter-chain 3D distance matrix from `sses_distance[d1][c1][d2][c2]`.
    2. Embeds the domain/chain centers into 2D using `method`.
    3. Optimizes the 2D rotation and reflection of each domain/chain so inter-domain SSE
       contacts match their 3D distances.
    """
    import utils_matrix2layout

    pairs = []
    for d in layout:
        for c in layout[d]:
            if c_layout.get(d, {}).get(c):
                pairs.append((d, c))

    n = len(pairs)
    if n <= 1:
        return layout, a_layout, c_layout, ca_layout

    D = np.zeros((n, n), dtype=float)
    weights = np.ones(n, dtype=float)
    radii_2d = np.zeros(n, dtype=float)

    for i, (d1, c1) in enumerate(pairs):
        w_i = sum(max(1.0, float(sses_sizes[d1][c1].get(s, 5.0))) for s in c_layout[d1][c1])
        weights[i] = w_i
        pts_2d = np.array(list(c_layout[d1][c1].values()), dtype=float)
        if len(pts_2d) > 1:
            center_2d = np.mean(pts_2d, axis=0)
            radii_2d[i] = float(np.percentile(np.linalg.norm(pts_2d - center_2d, axis=1), 75))
        else:
            radii_2d[i] = 3.0

    for i in range(n):
        d1, c1 = pairs[i]
        for j in range(i + 1, n):
            d2, c2 = pairs[j]
            dist_dict = sses_distance.get(d1, {}).get(c1, {}).get(d2, {}).get(c2, {})
            all_dists = []
            all_w = []
            for s1, row in dist_dict.items():
                w1 = max(1.0, float(sses_sizes[d1][c1].get(s1, 5.0)))
                for s2, val in row.items():
                    w2 = max(1.0, float(sses_sizes[d2][c2].get(s2, 5.0)))
                    all_dists.append(float(val))
                    all_w.append(w1 + w2)
            if all_dists:
                all_dists_arr = np.array(all_dists, dtype=float)
                all_w_arr = np.array(all_w, dtype=float)
                mean_d = float(np.sum(all_dists_arr * all_w_arr) / np.sum(all_w_arr))
            else:
                mean_d = 35.0
            # Ensure domain centers are at least separated by their combined internal radii
            min_sep = max(8.0, (radii_2d[i] + radii_2d[j]) * 1.15)
            eff_d = max(min_sep, mean_d)
            D[i, j] = eff_d
            D[j, i] = eff_d

    X2d = embed_distance_matrix(D, weights=weights, X3d=None, method=method)

    for i, (d, c) in enumerate(pairs):
        layout[d][c] = (float(X2d[i, 0]), float(X2d[i, 1]))

    # Optimize rotation (in 30-deg steps) and reflection of each domain/chain
    for _pass in range(2):
        for d, c in pairs:
            if len(c_layout[d][c]) <= 1:
                continue
            best_err = None
            best_ang = a_layout[d][c]
            best_refl = False

            for refl in (False, True):
                if refl:
                    for s in c_layout[d][c]:
                        c_layout[d][c][s] = (c_layout[d][c][s][0], -c_layout[d][c][s][1])
                        ca_layout[d][c][s] = (360.0 - ca_layout[d][c][s]) % 360.0

                for step_i in range(12):
                    ang = step_i * 30.0
                    a_layout[d][c] = ang
                    err = utils_matrix2layout.count_protein_error(
                        layout, a_layout, c_layout, ca_layout, d, c, sses_distance, sses_sizes
                    )
                    if best_err is None or err < best_err:
                        best_err = err
                        best_ang = ang
                        best_refl = refl

                if refl:
                    # Revert reflection for now
                    for s in c_layout[d][c]:
                        c_layout[d][c][s] = (c_layout[d][c][s][0], -c_layout[d][c][s][1])
                        ca_layout[d][c][s] = (360.0 - ca_layout[d][c][s]) % 360.0

            if best_refl:
                for s in c_layout[d][c]:
                    c_layout[d][c][s] = (c_layout[d][c][s][0], -c_layout[d][c][s][1])
                    ca_layout[d][c][s] = (360.0 - ca_layout[d][c][s]) % 360.0
            a_layout[d][c] = best_ang

    return layout, a_layout, c_layout, ca_layout
