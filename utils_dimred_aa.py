#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""
Amino-Acid (C-alpha Residue) Level Dimensionality Reduction Core for 2DProts.

Projects protein 3D structures into 2D directly by the positions of individual
amino acids (C-alpha atoms) rather than pre-clustered secondary structure elements,
supporting all 12 benchmarked dimensionality reduction algorithms:
  - "evolution"    : Multi-scale stochastic + local-contact evolutionary optimization on residues
  - "pca"          : Orthogonal 3D->2D Principal Component Analysis on C-alpha coordinates
  - "pcoa"         : Classical MDS (Principal Coordinates Analysis) on C-alpha distance matrix
  - "mds"          : Metric Multidimensional Scaling (SMACOF, with Nyström/Landmark acceleration)
  - "nmmds"        : Non-Metric Multidimensional Scaling (Kruskal rank-order MDS)
  - "sammon"       : Locally-Weighted MDS / Sammon Mapping (<40 Å cutoff + peptide backbone constraints)
  - "isomap"       : Isometric Feature Mapping (geodesic k-NN + backbone graph + Classical MDS)
  - "spectral"     : Laplacian Eigenmaps / Spectral Embedding on C-alpha contact + backbone graph
  - "lle"          : Locally Linear Embedding on 3D C-alpha neighborhood
  - "kamada_kawai" : Force-directed spring layout on the C-alpha distance & backbone graph
  - "tsne"         : t-Distributed Stochastic Neighbor Embedding on C-alpha distances
  - "umap"         : Uniform Manifold Approximation and Projection on C-alpha distances

Also derives induced 2D SSE centers and angles directly from the 2D positions of each
SSE's constituent amino acids for direct comparison with traditional 2DProts layouts.
"""

import math
import logging
import numpy as np
from scipy.optimize import minimize
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
from sklearn.decomposition import PCA
from sklearn.manifold import MDS, TSNE, SpectralEmbedding, LocallyLinearEmbedding
import networkx as nx

import superimpose2d

try:
    import umap
    HAS_UMAP = True
except ImportError:
    HAS_UMAP = False


AVAILABLE_AA_METHODS = [
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


def _classical_mds(D):
    """
    Deterministic Classical MDS (PCoA / Torgerson) from an NxN distance matrix D.
    """
    n = D.shape[0]
    if n <= 1:
        return np.zeros((n, 2), dtype=float)
    if n == 2:
        d = max(2.0, float(D[0, 1]))
        return np.array([[-0.5 * d, 0.0], [0.5 * d, 0.0]], dtype=float)

    D2 = D ** 2
    row_mean = np.mean(D2, axis=1, keepdims=True)
    col_mean = np.mean(D2, axis=0, keepdims=True)
    tot_mean = np.mean(D2)
    B = -0.5 * (D2 - row_mean - col_mean + tot_mean)
    B = 0.5 * (B + B.T)

    eigvals, eigvecs = np.linalg.eigh(B)
    idx = np.argsort(eigvals)[::-1]
    eigvals = np.maximum(eigvals[idx[:2]], 1e-6)
    eigvecs = eigvecs[:, idx[:2]]
    return eigvecs * np.sqrt(eigvals)[None, :]


def _rescale_and_orient_aa(X2d, D, X_ref):
    """
    Rescale 2D amino-acid coordinates X2d so their local (<40 Å) RMS distance matches
    the 3D distance matrix D, and align orientation to reference X_ref via 2D Procrustes.
    """
    n = X2d.shape[0]
    if n <= 1:
        return np.zeros((n, 2), dtype=float)

    X2d = np.asarray(X2d, dtype=float)
    X2d = X2d - np.mean(X2d, axis=0, keepdims=True)

    # Subsample up to 400 residues if N is very large for fast RMS scale estimation
    if n > 400:
        step = max(1, n // 400)
        sub = np.arange(0, n, step)
        D_sub = D[np.ix_(sub, sub)]
        X_sub = X2d[sub]
    else:
        D_sub = D
        X_sub = X2d

    m = D_sub.shape[0]
    diff = X_sub[:, None, :] - X_sub[None, :, :]
    d2d = np.sqrt(np.sum(diff ** 2, axis=-1))

    mask = np.triu(np.ones((m, m), dtype=bool), k=1)
    local_mask = mask & (D_sub < 40.0) & (D_sub > 1.0)
    if not np.any(local_mask):
        local_mask = mask

    rms_3d = float(np.sqrt(np.mean(D_sub[local_mask] ** 2)))
    rms_2d = float(np.sqrt(np.mean(d2d[local_mask] ** 2)))
    if rms_2d > 1e-6 and rms_3d > 1e-6:
        X2d = X2d * (rms_3d / rms_2d)

    if X_ref is not None and X_ref.shape == X2d.shape and np.linalg.norm(X_ref) > 1e-6:
        try:
            R, t = superimpose2d.optimal_rotation_translation(X2d.T, X_ref.T, allow_mirror=True)
            X2d = (R @ X2d.T + t).T
        except Exception:
            pass

    return X2d


def _select_landmarks(n, max_landmarks=180):
    """Select evenly-spaced landmark indices covering N residues (always including endpoints)."""
    if n <= max_landmarks:
        return np.arange(n, dtype=int)
    idx = np.unique(np.round(np.linspace(0, n - 1, max_landmarks)).astype(int))
    return idx


def _interpolate_from_landmarks(X_lm, lm_indices, D, X3d=None, pca_ref=None):
    """
    Smoothly extend a 2D landmark embedding `X_lm` (of shape M x 2) to all N amino acids.
    First rescales & aligns `X_lm` to Ångström scale on `D[lm_indices, lm_indices]` so local
    3D->2D affine interpolation has matching units across all methods.
    """
    n = D.shape[0]
    m = len(lm_indices)
    D_lm = D[np.ix_(lm_indices, lm_indices)]
    ref_lm = pca_ref[lm_indices] if (pca_ref is not None and len(pca_ref) == n) else None
    X_lm_scaled = _rescale_and_orient_aa(X_lm, D_lm, ref_lm)

    if m == n:
        return X_lm_scaled.copy()

    X_full = np.zeros((n, 2), dtype=float)
    X_full[lm_indices] = X_lm_scaled

    D_to_lm = D[:, lm_indices]  # N x M
    k = min(8, m)

    for i in range(n):
        d_row = D_to_lm[i]
        nn_idx = np.argsort(d_row)[:k]
        d_nn = d_row[nn_idx]
        if d_nn[0] < 1e-5:
            X_full[i] = X_lm_scaled[nn_idx[0]]
            continue
        w = 1.0 / ((d_nn + 1.5) ** 2) * np.exp(-((d_nn / 14.0) ** 2))
        w_sum = np.sum(w)
        if w_sum < 1e-9:
            w = 1.0 / (d_nn + 1.0)
            w_sum = np.sum(w)
        w = w / w_sum

        if X3d is not None:
            P3 = X3d[lm_indices[nn_idx]]
            P2 = X_lm_scaled[nn_idx]
            c3 = np.sum(P3 * w[:, None], axis=0)
            c2 = np.sum(P2 * w[:, None], axis=0)
            M3 = (P3 - c3) * np.sqrt(w)[:, None]
            M2 = (P2 - c2) * np.sqrt(w)[:, None]
            H = M3.T @ M2
            U, _, Vt = np.linalg.svd(H, full_matrices=False)
            proj = U @ Vt
            # Local scale factor of the 2D landmark patch vs 3D landmark patch
            norm3 = float(np.sqrt(np.sum(M3 ** 2)))
            norm2 = float(np.sqrt(np.sum(M2 ** 2)))
            loc_scale = np.clip(norm2 / max(norm3, 1e-4), 0.4, 1.8)
            X_full[i] = c2 + loc_scale * ((X3d[i] - c3) @ proj)
        else:
            X_full[i] = np.sum(X_lm_scaled[nn_idx] * w[:, None], axis=0)

    return X_full


def _sammon_aa_optimize(D, X_init, backbone_bonds=None, max_iter=180, pca_ref=None):
    """
    Quasi-Newton (L-BFGS-B) optimization of the 2DProts local-contact objective (<40 Å)
    on amino acids, with covalent peptide-bond springs and steric repulsion.
    """
    n = D.shape[0]
    if n <= 2:
        return X_init.copy()

    if n > 190:
        lm_idx = _select_landmarks(n, max_landmarks=170)
        D_lm = D[np.ix_(lm_idx, lm_idx)]
        X_lm_init = X_init[lm_idx]
        lm_bonds = []
        for a in range(len(lm_idx) - 1):
            i1, i2 = lm_idx[a], lm_idx[a + 1]
            if D[i1, i2] < 25.0:
                lm_bonds.append((a, a + 1, float(D[i1, i2])))
        X_lm_opt = _sammon_aa_optimize_core(D_lm, X_lm_init, lm_bonds, max_iter=max_iter)
        X_full = _interpolate_from_landmarks(X_lm_opt, lm_idx, D, X3d=None, pca_ref=pca_ref)
        return _relax_backbone_2d(X_full, D, backbone_bonds, steps=20)
    else:
        return _sammon_aa_optimize_core(D, X_init, backbone_bonds, max_iter=max_iter)


def _sammon_aa_optimize_core(D, X_init, backbone_bonds=None, max_iter=180):
    """Vectorized L-BFGS-B solver for N <= 360 amino acids / landmarks."""
    n = D.shape[0]
    W = np.zeros((n, n), dtype=float)
    mask_local = (D > 0.5) & (D < 42.0)
    d_safe = np.maximum(D, 3.2)
    W[mask_local] = (1.0 - (D[mask_local] / 45.0) ** 2) / d_safe[mask_local]
    mask_close = (D > 0.5) & (D < 9.0)
    W[mask_close] *= 2.2
    W[~mask_local & (D > 0.5)] = 0.015 / d_safe[~mask_local & (D > 0.5)]
    np.fill_diagonal(W, 0.0)

    if backbone_bonds:
        for i, j, d0 in backbone_bonds:
            if 0 <= i < n and 0 <= j < n:
                W[i, j] = max(W[i, j], 4.5 / max(d0, 3.0))
                W[j, i] = W[i, j]

    W_triu = np.triu(W, k=1)

    def obj_and_grad(x_flat):
        X = x_flat.reshape(n, 2)
        diff = X[:, None, :] - X[None, :, :]
        dist2 = np.sum(diff ** 2, axis=-1)
        dist = np.sqrt(dist2 + 1e-4)

        res = dist - D
        loss = 0.5 * np.sum(W_triu * (res ** 2))

        coeff = W * res / dist
        np.fill_diagonal(coeff, 0.0)

        rep_mask = (dist < 2.0) & (np.eye(n) == 0)
        if np.any(rep_mask):
            pen = 2.0 - dist
            loss += 2.0 * np.sum(np.triu(rep_mask * (pen ** 2), k=1))
            coeff = coeff - 4.0 * rep_mask * pen / dist

        grad = np.sum(coeff[:, :, None] * diff, axis=1)
        return float(loss), grad.ravel()

    res = minimize(
        obj_and_grad,
        X_init.ravel(),
        method="L-BFGS-B",
        jac=True,
        options={"maxiter": max_iter, "ftol": 1e-5},
    )
    return res.x.reshape(n, 2)


def _relax_backbone_2d(X2d, D, backbone_bonds, steps=20):
    """Fast vectorized spring relaxation along covalent backbone and local contacts."""
    if not backbone_bonds:
        return X2d
    X = X2d.copy()
    idx_i = np.array([b[0] for b in backbone_bonds], dtype=int)
    idx_j = np.array([b[1] for b in backbone_bonds], dtype=int)
    d_target = np.array([b[2] for b in backbone_bonds], dtype=float)

    for _ in range(steps):
        diff = X[idx_j] - X[idx_i]
        dist = np.sqrt(np.sum(diff ** 2, axis=1) + 1e-6)
        err = (dist - d_target) / dist
        step_vec = 0.18 * err[:, None] * diff
        np.add.at(X, idx_i, step_vec)
        np.add.at(X, idx_j, -step_vec)
        X = 0.92 * X + 0.08 * X2d
    return X


def _evolution_aa_optimize(D, X3d, X_init, backbone_bonds=None, random_state=42, pca_ref=None):
    """
    Amino-acid evolutionary / stochastic hill-climbing optimizer inspired by 2DProts.
    """
    rng = np.random.default_rng(random_state)
    n = D.shape[0]
    if n <= 2:
        return X_init.copy()

    if n > 190:
        lm_idx = _select_landmarks(n, max_landmarks=170)
        D_lm = D[np.ix_(lm_idx, lm_idx)]
        X_lm = _evolution_aa_core(D_lm, X_init[lm_idx], rng, n_steps=1500)
        X_full = _interpolate_from_landmarks(X_lm, lm_idx, D, X3d=X3d, pca_ref=pca_ref)
        return _relax_backbone_2d(X_full, D, backbone_bonds, steps=20)
    else:
        X_opt = _evolution_aa_core(D, X_init, rng, n_steps=1800)
        return _relax_backbone_2d(X_opt, D, backbone_bonds, steps=15)


def _evolution_aa_core(D, X_init, rng, n_steps=2000):
    """Fast NumPy stochastic hill-climbing on residue/landmark positions."""
    n = D.shape[0]
    X = X_init.copy()
    mask_local = D < 40.0
    np.fill_diagonal(mask_local, False)

    diff = X[:, None, :] - X[None, :, :]
    D2 = np.sqrt(np.sum(diff ** 2, axis=-1))
    abs_err = np.abs(D2 - D) * mask_local

    for step in range(n_steps):
        temp = max(0.15, 3.5 * (1.0 - step / float(n_steps)) ** 1.6)
        i = int(rng.integers(0, n))
        delta = rng.normal(0.0, temp, size=2)
        new_pos = X[i] + delta

        row_mask = mask_local[i]
        if not np.any(row_mask):
            continue
        old_row_err = float(np.sum(abs_err[i, row_mask]))
        d_new = np.sqrt(np.sum((X[row_mask] - new_pos[None, :]) ** 2, axis=1))
        new_row_err = float(np.sum(np.abs(d_new - D[i, row_mask])))

        seq_pen_old = 0.0
        seq_pen_new = 0.0
        for nb in (i - 1, i + 1):
            if 0 <= nb < n and D[i, nb] < 12.0:
                seq_pen_old += 2.5 * abs(float(D2[i, nb]) - float(D[i, nb]))
                d_nb_new = float(np.linalg.norm(X[nb] - new_pos))
                seq_pen_new += 2.5 * abs(d_nb_new - float(D[i, nb]))

        if (new_row_err + seq_pen_new) < (old_row_err + seq_pen_old):
            X[i] = new_pos
            all_d_new = np.sqrt(np.sum((X - new_pos[None, :]) ** 2, axis=1))
            D2[i, :] = all_d_new
            D2[:, i] = all_d_new
            new_abs = np.abs(all_d_new - D[i, :]) * row_mask
            abs_err[i, :] = new_abs
            abs_err[:, i] = new_abs

    return X


def _stabilize_nonlinear_embedding(X2d, D, pca_ref):
    """
    Ensure non-linear eigen-embeddings (Spectral, LLE) do not collapse bulk coordinates
    into a 1D line or outlier spike when a protein has flexible terminal tails.
    """
    X2d = _rescale_and_orient_aa(X2d, D, pca_ref)
    radii = np.linalg.norm(X2d - np.mean(X2d, axis=0, keepdims=True), axis=1)
    med_r = float(np.median(radii))
    p98_r = float(np.percentile(radii, 98))
    if med_r < 1e-3 or (p98_r / max(med_r, 0.5)) > 3.6:
        # Softly clip extreme tail outliers and regularize with 20% PCoA/PCA geometry
        cap = max(med_r * 2.8, 15.0)
        scale_fac = np.minimum(1.0, cap / np.maximum(radii, 1e-6))
        X2d = X2d * scale_fac[:, None]
        X2d = _rescale_and_orient_aa(X2d, D, pca_ref)
        X2d = 0.80 * X2d + 0.20 * pca_ref
        X2d = _rescale_and_orient_aa(X2d, D, pca_ref)
    return X2d


def embed_amino_acids(D, X3d, X3d_raw=None, backbone_bonds=None, method="mds", random_state=42):
    """
    Project N amino acids from 3D into 2D using one of the 12 benchmark algorithms.
    Returns an Nx2 numpy array of 2D amino-acid coordinates in Ångström scale.
    """
    n = D.shape[0]
    if n == 0:
        return np.zeros((0, 2), dtype=float)
    if n == 1:
        return np.zeros((1, 2), dtype=float)

    base_method = method.split("+")[0].strip().lower()

    # Canonical reference orientation via PCA on X3d
    X_centered = X3d - np.mean(X3d, axis=0, keepdims=True)
    pca_ref = PCA(n_components=2, random_state=random_state).fit_transform(X_centered)
    pca_ref = _rescale_and_orient_aa(pca_ref, D, None)

    try:
        if base_method == "pca":
            pts = X3d_raw if (X3d_raw is not None and X3d_raw.shape[0] == n) else X3d
            pts_c = pts - np.mean(pts, axis=0, keepdims=True)
            X2d = PCA(n_components=2, random_state=random_state).fit_transform(pts_c)
            return _rescale_and_orient_aa(X2d, D, pca_ref)

        elif base_method in ("pcoa", "cmds", "classical_mds"):
            D_mod = D * (1.0 + 0.12 * np.exp(-((D / 14.0) ** 2)))
            np.fill_diagonal(D_mod, 0.0)
            X2d = _classical_mds(D_mod)
            return _rescale_and_orient_aa(X2d, D, pca_ref)

        elif base_method == "evolution":
            X_init = _classical_mds(D)
            X_init = _rescale_and_orient_aa(X_init, D, pca_ref)
            X2d = _evolution_aa_optimize(
                D, X3d, X_init, backbone_bonds=backbone_bonds, random_state=random_state, pca_ref=pca_ref
            )
            return _rescale_and_orient_aa(X2d, D, pca_ref)

        elif base_method == "mds":
            if n <= 190:
                X_init = _classical_mds(D)
                mds = MDS(
                    n_components=2,
                    metric=True,
                    dissimilarity="precomputed",
                    n_init=1,
                    max_iter=180,
                    random_state=random_state,
                )
                X2d = mds.fit_transform(D, init=X_init)
            else:
                lm_idx = _select_landmarks(n, max_landmarks=170)
                D_lm = D[np.ix_(lm_idx, lm_idx)]
                X_lm_init = _classical_mds(D_lm)
                mds = MDS(
                    n_components=2,
                    metric=True,
                    dissimilarity="precomputed",
                    n_init=1,
                    max_iter=180,
                    random_state=random_state,
                )
                X_lm = mds.fit_transform(D_lm, init=X_lm_init)
                X2d = _interpolate_from_landmarks(X_lm, lm_idx, D, X3d=X3d, pca_ref=pca_ref)
                X2d = _relax_backbone_2d(X2d, D, backbone_bonds, steps=15)
            return _rescale_and_orient_aa(X2d, D, pca_ref)

        elif base_method == "nmmds":
            if n <= 190:
                X_init = _classical_mds(D)
                mds = MDS(
                    n_components=2,
                    metric=False,
                    dissimilarity="precomputed",
                    n_init=1,
                    max_iter=180,
                    random_state=random_state,
                )
                X2d = mds.fit_transform(D, init=X_init)
            else:
                lm_idx = _select_landmarks(n, max_landmarks=170)
                D_lm = D[np.ix_(lm_idx, lm_idx)]
                X_lm_init = _classical_mds(D_lm)
                mds = MDS(
                    n_components=2,
                    metric=False,
                    dissimilarity="precomputed",
                    n_init=1,
                    max_iter=180,
                    random_state=random_state,
                )
                X_lm = mds.fit_transform(D_lm, init=X_lm_init)
                X2d = _interpolate_from_landmarks(X_lm, lm_idx, D, X3d=X3d, pca_ref=pca_ref)
                X2d = _relax_backbone_2d(X2d, D, backbone_bonds, steps=15)
            return _rescale_and_orient_aa(X2d, D, pca_ref)

        elif base_method in ("sammon", "wmds"):
            X_init = _classical_mds(D) if n <= 240 else pca_ref
            X_init = _rescale_and_orient_aa(X_init, D, pca_ref)
            X2d = _sammon_aa_optimize(D, X_init, backbone_bonds=backbone_bonds, max_iter=140, pca_ref=pca_ref)
            return _rescale_and_orient_aa(X2d, D, pca_ref)

        elif base_method == "isomap":
            k = max(4, min(12, n - 1))
            rows, cols, vals = [], [], []
            for i in range(n):
                nbrs = np.argsort(D[i])[1 : k + 1]
                for j in nbrs:
                    d_ij = float(D[i, j])
                    rows.extend([i, j])
                    cols.extend([j, i])
                    vals.extend([d_ij, d_ij])
            if backbone_bonds:
                for i, j, d0 in backbone_bonds:
                    rows.extend([i, j])
                    cols.extend([j, i])
                    vals.extend([float(d0), float(d0)])
            G_csr = csr_matrix((vals, (rows, cols)), shape=(n, n))
            if n <= 220:
                geo_D = shortest_path(G_csr, directed=False)
                inf_mask = ~np.isfinite(geo_D)
                if np.any(inf_mask):
                    geo_D[inf_mask] = D[inf_mask] * 1.4
                X2d = _classical_mds(geo_D)
            else:
                lm_idx = _select_landmarks(n, max_landmarks=180)
                geo_lm_to_all = shortest_path(G_csr, directed=False, indices=lm_idx)
                geo_lm = geo_lm_to_all[:, lm_idx]
                inf_mask = ~np.isfinite(geo_lm)
                if np.any(inf_mask):
                    D_lm = D[np.ix_(lm_idx, lm_idx)]
                    geo_lm[inf_mask] = D_lm[inf_mask] * 1.4
                geo_lm = 0.5 * (geo_lm + geo_lm.T)
                X_lm = _classical_mds(geo_lm)
                X2d = _interpolate_from_landmarks(X_lm, lm_idx, D, X3d=X3d, pca_ref=pca_ref)
            return _rescale_and_orient_aa(X2d, D, pca_ref)

        elif base_method == "spectral":
            if n <= 4:
                return pca_ref
            if n <= 220:
                med_d = float(np.median(D[D > 0])) if np.any(D > 0) else 15.0
                sigma = max(6.0, med_d * 0.75)
                affinity = 0.85 * np.exp(-(D ** 2) / (2.0 * (sigma ** 2))) + 0.15 * np.exp(
                    -(D ** 2) / (2.0 * ((2.5 * sigma) ** 2))
                )
                if backbone_bonds:
                    for i, j, _ in backbone_bonds:
                        affinity[i, j] = max(affinity[i, j], 0.92)
                        affinity[j, i] = affinity[i, j]
                np.fill_diagonal(affinity, 1.0)
                spec = SpectralEmbedding(n_components=2, affinity="precomputed", random_state=random_state)
                X2d = spec.fit_transform(affinity)
            else:
                lm_idx = _select_landmarks(n, max_landmarks=180)
                D_lm = D[np.ix_(lm_idx, lm_idx)]
                med_d = float(np.median(D_lm[D_lm > 0])) if np.any(D_lm > 0) else 15.0
                sigma = max(6.0, med_d * 0.75)
                affinity = 0.85 * np.exp(-(D_lm ** 2) / (2.0 * (sigma ** 2))) + 0.15 * np.exp(
                    -(D_lm ** 2) / (2.0 * ((2.5 * sigma) ** 2))
                )
                np.fill_diagonal(affinity, 1.0)
                spec = SpectralEmbedding(n_components=2, affinity="precomputed", random_state=random_state)
                X_lm = spec.fit_transform(affinity)
                X_lm = _stabilize_nonlinear_embedding(X_lm, D_lm, pca_ref[lm_idx])
                X2d = _interpolate_from_landmarks(X_lm, lm_idx, D, X3d=X3d, pca_ref=pca_ref)
            return _stabilize_nonlinear_embedding(X2d, D, pca_ref)

        elif base_method == "lle":
            if n <= 5:
                return pca_ref
            if n <= 220:
                k = max(8, min(22, n // 10))
                lle = LocallyLinearEmbedding(
                    n_neighbors=k,
                    n_components=2,
                    method="standard",
                    reg=5e-2,
                    random_state=random_state,
                )
                X2d = lle.fit_transform(X3d)
            else:
                lm_idx = _select_landmarks(n, max_landmarks=180)
                k = max(8, min(22, len(lm_idx) // 10))
                lle = LocallyLinearEmbedding(
                    n_neighbors=k,
                    n_components=2,
                    method="standard",
                    reg=5e-2,
                    random_state=random_state,
                )
                X_lm = lle.fit_transform(X3d[lm_idx])
                X_lm = _stabilize_nonlinear_embedding(X_lm, D[np.ix_(lm_idx, lm_idx)], pca_ref[lm_idx])
                X2d = _interpolate_from_landmarks(X_lm, lm_idx, D, X3d=X3d, pca_ref=pca_ref)
            return _stabilize_nonlinear_embedding(X2d, D, pca_ref)

        elif base_method in ("kamada_kawai", "kk"):
            if n <= 160:
                G = nx.Graph()
                for i in range(n):
                    G.add_node(i)
                dist_dict = {i: {j: float(D[i, j]) for j in range(n)} for i in range(n)}
                init_pos = {i: pca_ref[i] for i in range(n)}
                pos = nx.kamada_kawai_layout(G, dist=dist_dict, pos=init_pos)
                X2d = np.array([pos[i] for i in range(n)], dtype=float)
            else:
                lm_idx = _select_landmarks(n, max_landmarks=140)
                m = len(lm_idx)
                D_lm = D[np.ix_(lm_idx, lm_idx)]
                G = nx.Graph()
                for i in range(m):
                    G.add_node(i)
                dist_dict = {i: {j: float(D_lm[i, j]) for j in range(m)} for i in range(m)}
                init_pos = {i: pca_ref[lm_idx[i]] for i in range(m)}
                pos = nx.kamada_kawai_layout(G, dist=dist_dict, pos=init_pos)
                X_lm = np.array([pos[i] for i in range(m)], dtype=float)
                X2d = _interpolate_from_landmarks(X_lm, lm_idx, D, X3d=X3d, pca_ref=pca_ref)
                X2d = _relax_backbone_2d(X2d, D, backbone_bonds, steps=15)
            return _rescale_and_orient_aa(X2d, D, pca_ref)

        elif base_method == "tsne":
            if n <= 5:
                return pca_ref
            if n <= 180:
                perp = max(5.0, min(22.0, float(n - 1) / 3.0))
                tsne = TSNE(
                    n_components=2,
                    perplexity=perp,
                    metric="precomputed",
                    init=pca_ref,
                    learning_rate="auto",
                    max_iter=300,
                    random_state=random_state,
                )
                X2d = tsne.fit_transform(D)
            else:
                lm_idx = _select_landmarks(n, max_landmarks=160)
                D_lm = D[np.ix_(lm_idx, lm_idx)]
                perp = max(5.0, min(22.0, float(len(lm_idx) - 1) / 3.0))
                tsne = TSNE(
                    n_components=2,
                    perplexity=perp,
                    metric="precomputed",
                    init=pca_ref[lm_idx],
                    learning_rate="auto",
                    max_iter=300,
                    random_state=random_state,
                )
                X_lm = tsne.fit_transform(D_lm)
                X2d = _interpolate_from_landmarks(X_lm, lm_idx, D, X3d=X3d, pca_ref=pca_ref)
            return _rescale_and_orient_aa(X2d, D, pca_ref)

        elif base_method == "umap":
            if not HAS_UMAP or n <= 6:
                return _sammon_aa_optimize(D, pca_ref, backbone_bonds=backbone_bonds, max_iter=120, pca_ref=pca_ref)
            if n <= 180:
                n_neighbors = max(5, min(18, n - 1))
                reducer = umap.UMAP(
                    n_components=2,
                    n_neighbors=n_neighbors,
                    min_dist=0.25,
                    metric="precomputed",
                    init=pca_ref,
                    random_state=random_state,
                    n_epochs=180,
                )
                X2d = reducer.fit_transform(D)
            else:
                lm_idx = _select_landmarks(n, max_landmarks=160)
                D_lm = D[np.ix_(lm_idx, lm_idx)]
                n_neighbors = max(5, min(18, len(lm_idx) - 1))
                reducer = umap.UMAP(
                    n_components=2,
                    n_neighbors=n_neighbors,
                    min_dist=0.25,
                    metric="precomputed",
                    init=pca_ref[lm_idx],
                    random_state=random_state,
                    n_epochs=180,
                )
                X_lm = reducer.fit_transform(D_lm)
                X2d = _interpolate_from_landmarks(X_lm, lm_idx, D, X3d=X3d, pca_ref=pca_ref)
            return _rescale_and_orient_aa(X2d, D, pca_ref)

        else:
            logging.warning(f"Unknown AA dimred method '{method}', falling back to 'mds'")
            return embed_amino_acids(D, X3d, X3d_raw=X3d_raw, backbone_bonds=backbone_bonds, method="mds")

    except Exception as exc:
        logging.warning(f"AA dimensionality reduction '{method}' failed ({exc}), using PCoA fallback.")
        return pca_ref


def derive_sse_layout_from_aa(aa_records, X2d_aa, all_sses_3d):
    """
    Given the 2D amino-acid coordinates `X2d_aa` (Nx2) and `aa_records`,
    derive the 2D center `(x, y)`, 2D orientation angle `deg` (in [0, 360)),
    and mapping of SSE label -> constituent residue indices.
    """
    sse_to_indices = {}
    for idx, rec in enumerate(aa_records):
        s_lbl = rec.get("sse_label")
        if s_lbl:
            sse_to_indices.setdefault(s_lbl, []).append(idx)

    sse_layout = {}
    sse_angles = {}

    for s_lbl, indices in sse_to_indices.items():
        pts = X2d_aa[indices]
        center = np.mean(pts, axis=0)
        sse_layout[s_lbl] = (float(center[0]), float(center[1]))

        if len(pts) >= 2:
            centered = pts - center[None, :]
            _, _, vt = np.linalg.svd(centered, full_matrices=False)
            direction = vt[0]
            end_to_start = pts[-1] - pts[0]
            if float(np.dot(end_to_start, direction)) < 0:
                direction = -direction
            ang_deg = math.degrees(math.atan2(float(direction[0]), float(direction[1]))) % 360.0
            sse_angles[s_lbl] = ang_deg
        else:
            sse_angles[s_lbl] = 0.0

    for s_lbl in all_sses_3d:
        if s_lbl not in sse_layout:
            sse_layout[s_lbl] = (0.0, 0.0)
            sse_angles[s_lbl] = 0.0

    return sse_layout, sse_angles, sse_to_indices
