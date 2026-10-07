#!/usr/bin/env python3
"""
aSAM Ensemble Clustering and Analysis Script.
Clusters the SAM2 conformational ensemble using SVD PCA + Scipy KMeans clustering on C-alpha coordinates,
extracts one representative structure (frame) per cluster into a representative DCD trajectory,
and writes metadata for downstream template extraction.
"""

import sys
import os
import numpy as np
import scipy.cluster.vq as vq
import MDAnalysis as mda
from MDAnalysis.analysis import align

def main():
    if len(sys.argv) < 5:
        print("Usage: python3 asam_cluster_analysis.py <top_pdb> <traj_dcd> <ref_structure> <outdir> [n_clusters] [pca_dim]")
        sys.exit(1)

    top_path = sys.argv[1]
    traj_path = sys.argv[2]
    ref_path = sys.argv[3]
    outdir = sys.argv[4]
    
    n_clusters_target = int(sys.argv[5]) if len(sys.argv) > 5 else 10
    pca_dim_target = int(sys.argv[6]) if len(sys.argv) > 6 else 10

    os.makedirs(os.path.join(outdir, "clustering"), exist_ok=True)

    print(f"Loading ensemble: topology={top_path}, trajectory={traj_path}")
    u = mda.Universe(top_path, traj_path)
    n_frames = len(u.trajectory)
    print(f"Total ensemble frames: {n_frames}")

    if n_frames < 2:
        print("Too few frames for clustering.")
        meta_path = os.path.join(outdir, "cluster_meta.sh")
        with open(meta_path, "w") as f:
            f.write(f"ASAM_N_FRAMES={n_frames}\n")
            f.write("ASAM_N_CLUSTERS_USED=0\n")
            f.write("ASAM_USED_CLUSTERING=false\n")
        sys.exit(0)

    # Select CA atoms for alignment and clustering
    ca_atoms = u.select_atoms("name CA")
    if len(ca_atoms) == 0:
        ca_atoms = u.atoms

    # Align all frames to the first frame
    aligner = align.AlignTraj(u, u, select="name CA" if len(u.select_atoms("name CA")) > 0 else "all", in_memory=True)
    aligner.run()

    # Extract CA coordinate vectors for each frame: shape (n_frames, n_ca * 3)
    coords = []
    for ts in u.trajectory:
        coords.append(ca_atoms.positions.flatten())
    coords = np.array(coords)

    # PCA reduction using Numpy SVD
    coords_centered = coords - np.mean(coords, axis=0)
    n_components = min(pca_dim_target, coords_centered.shape[0], coords_centered.shape[1])
    _, _, vt = np.linalg.svd(coords_centered, full_matrices=False)
    coords_pca = np.dot(coords_centered, vt.T[:, :n_components])

    # KMeans clustering using Scipy
    actual_k = min(n_clusters_target, n_frames)
    centroids, cluster_labels = vq.kmeans2(coords_pca, actual_k, minit='points', seed=42)

    # Find representative frame for each cluster (closest frame to centroid)
    rep_frames = []
    for i in range(actual_k):
        cluster_indices = np.where(cluster_labels == i)[0]
        if len(cluster_indices) == 0:
            continue
        centroid = centroids[i]
        distances = np.linalg.norm(coords_pca[cluster_indices] - centroid, axis=1)
        rep_idx = cluster_indices[np.argmin(distances)]
        rep_frames.append(rep_idx)

    rep_frames = sorted(list(set(rep_frames)))
    print(f"Extracted {len(rep_frames)} representative frames from {actual_k} clusters: {rep_frames}")

    # Write representative DCD
    rep_dcd_path = os.path.join(outdir, "clustering", "cluster_representatives.dcd")
    with mda.Writer(rep_dcd_path, u.atoms.n_atoms) as W:
        for f_idx in rep_frames:
            u.trajectory[f_idx]
            W.write(u.atoms)

    # Write metadata script
    meta_path = os.path.join(outdir, "cluster_meta.sh")
    with open(meta_path, "w") as f:
        f.write(f"ASAM_N_FRAMES={n_frames}\n")
        f.write(f"ASAM_N_CLUSTERS_USED={len(rep_frames)}\n")
        f.write("ASAM_USED_CLUSTERING=true\n")

    print(f"Clustering complete. Meta script written to {meta_path}")

if __name__ == "__main__":
    main()
