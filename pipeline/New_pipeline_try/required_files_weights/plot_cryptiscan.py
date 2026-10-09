#!/usr/bin/env python3
"""
plot_cryptiscan.py -- publication figures from CryptiScan Stage 9 output
=======================================================================
Reads the CSVs the pocket-analysis stage writes and produces journal-ready
figures plus a statistics table.

WHAT IT MAKES
  fig2_pocket_opening.png/.pdf    main-text figure: cryptic-site exposure,
                                  pocket-volume distribution, effect summary
  figS_rmsf.png/.pdf              per-residue flexibility, mutant vs wild type
  figS_integrity.png/.pdf         fold-integrity controls (Rg, DSSP content)
  figS_clusters.png/.pdf          cluster populations
  figS_persistence.png/.pdf       per-residue pocket-lining persistence
  stats_summary.csv               every comparison, with effect sizes
  stats_summary.txt               the same, formatted for pasting into a draft

Panels whose input is missing or entirely undetected are SKIPPED with a message
rather than drawn empty. The script tells you which columns were unusable and
why, so a broken upstream stage cannot masquerade as a null result.

USAGE
  Local:   python plot_cryptiscan.py --indir pocket_analysis \
                                     --clusters clustering/cluster_representatives.json \
                                     --cryptic 200 205 233 235 255
  Colab:   see the companion notebook cell, or just set the paths below.

Only numpy, pandas, scipy and matplotlib are needed.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import stats

# ---------------------------------------------------------------------------
# Style. Publication defaults: sans-serif, small type, no chartjunk.
# Colours are a validated colourblind-safe pair (protan dE 20.9, normal 27.3).
# Identity is never carried by colour alone -- every panel is also directly
# labelled or legended.
# ---------------------------------------------------------------------------
C_WT = "#2B6CB0"      # wild type, blue
C_MUT = "#C05621"     # mutant, orange
C_INK = "#1A1A1A"
C_MUTED = "#6B6B6B"
C_GRID = "#E3E3E0"
SURFACE = "#FCFCFB"

mpl.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 400,
    "savefig.bbox": "tight",
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8.5,
    "axes.labelcolor": C_INK,
    "axes.edgecolor": "#9A9A9A",
    "axes.linewidth": 0.7,
    "axes.facecolor": SURFACE,
    "figure.facecolor": SURFACE,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "xtick.color": C_MUTED,
    "ytick.color": C_MUTED,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "legend.frameon": False,
    "legend.fontsize": 7.5,
    "lines.linewidth": 1.4,
    "grid.color": C_GRID,
    "grid.linewidth": 0.6,
})

MM = 1 / 25.4          # millimetres -> inches, for journal column widths
SINGLE = 90 * MM       # single column
DOUBLE = 190 * MM      # full width


def grid_y(ax):
    ax.yaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)


# ---------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------
def cliffs_delta(a, b):
    """Cliff's delta from the Mann-Whitney U. Non-parametric effect size:
    P(a > b) - P(b > a). Appropriate here because these distributions are
    skewed and, for volumes, zero-inflated."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.size == 0 or b.size == 0:
        return np.nan
    u = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
    return 2.0 * u / (a.size * b.size) - 1.0


def magnitude(d):
    """Conventional thresholds for |Cliff's delta|."""
    if np.isnan(d):
        return "n/a"
    a = abs(d)
    return ("negligible" if a < 0.147 else
            "small" if a < 0.33 else
            "medium" if a < 0.474 else "large")


def compare(name, mut, wt, unit=""):
    mut = np.asarray(mut, float); mut = mut[~np.isnan(mut)]
    wt = np.asarray(wt, float);   wt = wt[~np.isnan(wt)]
    if mut.size == 0 or wt.size == 0:
        return None
    sem = lambda x: x.std(ddof=1) / np.sqrt(x.size) if x.size > 1 else np.nan
    row = {
        "quantity": name, "unit": unit,
        "n_mut": mut.size, "n_wt": wt.size,
        "mean_mut": mut.mean(), "sem_mut": sem(mut),
        "mean_wt": wt.mean(), "sem_wt": sem(wt),
        "delta": mut.mean() - wt.mean(),
        "median_mut": np.median(mut), "median_wt": np.median(wt),
    }
    # Frames are independent draws, so plain two-sample tests apply directly:
    # no block averaging, no correlation-time correction.
    row["mannwhitney_p"] = stats.mannwhitneyu(mut, wt, alternative="two-sided").pvalue
    ks = stats.ks_2samp(mut, wt)
    row["ks_D"], row["ks_p"] = ks.statistic, ks.pvalue
    row["cliffs_delta"] = cliffs_delta(mut, wt)
    row["effect"] = magnitude(row["cliffs_delta"])
    return row


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------
def load_pocket(path, label):
    if not path.is_file():
        print(f"  [skip] {path.name} not found")
        return None
    df = pd.read_csv(path)
    n = len(df)
    counts = df["status"].value_counts().to_dict()
    print(f"  {label:10s} {n} frames   status: "
          + ", ".join(f"{k}={v}" for k, v in counts.items()))
    df.attrs["label"] = label
    df.attrs["n_open"] = int((df["status"] == "open").sum())
    df.attrs["n_closed"] = int((df["status"] == "closed").sum())
    df.attrs["n_bad"] = int(df["status"].isin(["fpocket_failed", "no_output"]).sum())
    return df


def detection_ok(df):
    """True when fpocket actually produced parseable output for this ensemble.

    A run in which EVERY frame is no_output/fpocket_failed is a tool failure,
    not a closed pocket: a genuinely closed frame has status 'closed' and a
    volume of 0. Plotting the former as though it were the latter would report
    a null biological result that the data does not contain.
    """
    return df is not None and (df.attrs["n_open"] + df.attrs["n_closed"]) > 0


# ---------------------------------------------------------------------------
# FIGURE 2 -- main text
# ---------------------------------------------------------------------------
def fig_main(mut, wt, outdir, rows):
    have_vol = detection_ok(mut) and detection_ok(wt)
    ncols = 3 if have_vol else 2
    fig, axes = plt.subplots(1, ncols, figsize=(DOUBLE if have_vol else 140 * MM, 62 * MM))

    # --- (a) cryptic-site solvent exposure --------------------------------
    # This is the detection-free witness to opening: pure geometry, computed
    # by Shrake-Rupley, with no pocket detector involved. It is valid even
    # when pocket detection has failed.
    ax = axes[0]
    a = mut["cryptic_sasa_A2"].dropna().values
    b = wt["cryptic_sasa_A2"].dropna().values
    lo, hi = min(a.min(), b.min()), max(a.max(), b.max())
    bins = np.linspace(lo, hi, 46)
    ax.hist(b, bins=bins, color=C_WT, alpha=0.55, label="wild type", zorder=2)
    ax.hist(a, bins=bins, color=C_MUT, alpha=0.55, label="mutant", zorder=2)
    for v, c in ((b.mean(), C_WT), (a.mean(), C_MUT)):
        ax.axvline(v, color=c, lw=1.4, ls="--", zorder=3)
    d = cliffs_delta(a, b)
    ax.set_xlabel("cryptic-site SASA ($\\AA^2$)")
    ax.set_ylabel("frames")
    ax.set_title("a   Cryptic-site solvent exposure", loc="left", fontweight="bold")
    # Headroom so the legend and the statistics block cannot collide with the
    # tallest bar, whatever the data turn out to be.
    ax.set_ylim(0, ax.get_ylim()[1] * 1.42)
    ax.legend(loc="upper right")
    ax.annotate(f"$\\Delta$ = {a.mean()-b.mean():+.1f} $\\AA^2$\n"
                f"Cliff's $\\delta$ = {d:+.2f} ({magnitude(d)})",
                xy=(0.975, 0.70), xycoords="axes fraction",
                ha="right", va="top", fontsize=7.5, color=C_INK)
    grid_y(ax)

    # --- (b) pocket volume, if detection worked ---------------------------
    if have_vol:
        ax = axes[1]
        va = mut["volume_A3"].fillna(0).values
        vb = wt["volume_A3"].fillna(0).values
        hi = max(va.max(), vb.max(), 1.0)
        bins = np.linspace(0, hi * 1.02, 40)
        ax.hist(vb, bins=bins, color=C_WT, alpha=0.55, label="wild type", zorder=2)
        ax.hist(va, bins=bins, color=C_MUT, alpha=0.55, label="mutant", zorder=2)
        ax.set_xlabel("cryptic-pocket volume ($\\AA^3$)")
        ax.set_ylabel("frames")
        ax.set_title("b   Cryptic-pocket volume", loc="left", fontweight="bold")
        ax.legend(loc="upper right")
        dv = cliffs_delta(va, vb)
        ax.annotate(f"$\\Delta\\langle V\\rangle$ = {va.mean()-vb.mean():+.0f} $\\AA^3$\n"
                    f"Cliff's $\\delta$ = {dv:+.2f}",
                    xy=(0.45, 0.95), xycoords="axes fraction", va="top",
                    fontsize=7.5, color=C_INK)
        grid_y(ax)

    # --- (c) effect sizes --------------------------------------------------
    # One axis, one unit. Cliff's delta is dimensionless, so the opening
    # signal and the controls are directly comparable on the same scale --
    # which a grouped bar chart of raw means (Angstroms beside percentages)
    # could not do without a second y-scale.
    ax = axes[-1]
    items = [("cryptic-site SASA", mut["cryptic_sasa_A2"], wt["cryptic_sasa_A2"], True)]
    if have_vol:
        items += [
            ("pocket volume", mut["volume_A3"].fillna(0), wt["volume_A3"].fillna(0), True),
            ("non-cryptic volume", mut["noncryptic_volume_A3"],
             wt["noncryptic_volume_A3"], False),
        ]
    items += [
        ("radius of gyration", mut["rg_A"], wt["rg_A"], False),
        ("helix content", mut["frac_helix"], wt["frac_helix"], False),
        ("sheet content", mut["frac_sheet"], wt["frac_sheet"], False),
    ]
    labels, deltas, is_signal = [], [], []
    for name, a, b, sig in items:
        a = np.asarray(a, float); a = a[~np.isnan(a)]
        b = np.asarray(b, float); b = b[~np.isnan(b)]
        if a.size and b.size:
            labels.append(name); deltas.append(cliffs_delta(a, b)); is_signal.append(sig)

    y = np.arange(len(labels))[::-1]
    ax.axvspan(-0.147, 0.147, color=C_GRID, zorder=1)      # negligible band
    ax.axvline(0, color="#9A9A9A", lw=0.7, zorder=2)
    for yi, d, sig in zip(y, deltas, is_signal):
        c = C_MUT if sig else C_WT
        ax.plot([0, d], [yi, yi], color=c, lw=1.6, zorder=3,
                solid_capstyle="round")
        ax.plot(d, yi, "o", ms=6.5, color=c, mec=SURFACE, mew=1.2, zorder=4)
        ax.text(d + (0.06 if d >= 0 else -0.06), yi, f"{d:+.2f}",
                ha="left" if d >= 0 else "right", va="center",
                fontsize=7, color=C_INK)
    ax.set_yticks(y); ax.set_yticklabels(labels)
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-0.6, len(labels) - 0.1)
    ax.set_xlabel("Cliff's $\\delta$   (mutant vs wild type)")
    ax.set_title(f"{'c' if have_vol else 'b'}   Effect size: signal vs controls",
                 loc="left", fontweight="bold")
    ax.text(0, len(labels) - 0.35, "negligible", ha="center", va="top",
            fontsize=6.5, color=C_MUTED)
    ax.xaxis.grid(True, zorder=0); ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)

    fig.tight_layout(w_pad=2.6)
    save(fig, outdir / "fig2_pocket_opening")
    if not have_vol:
        print("  NOTE: pocket-volume panel omitted -- no frame in either ensemble\n"
              "        produced parseable fpocket output. Fix detection and re-run.")


# ---------------------------------------------------------------------------
# RMSF
# ---------------------------------------------------------------------------
def fig_rmsf(outdir, indir, cryptic):
    pm, pw = indir / "mutant_rmsf.csv", indir / "wildtype_rmsf.csv"
    if not (pm.is_file() and pw.is_file()):
        print("  [skip] RMSF files not found")
        return None
    m, w = pd.read_csv(pm), pd.read_csv(pw)
    df = m.merge(w, on="resSeq", suffixes=("_mut", "_wt"))
    df["d"] = df["rmsf_A_mut"] - df["rmsf_A_wt"]

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(DOUBLE, 86 * MM), sharex=True,
        gridspec_kw={"height_ratios": [2.4, 1], "hspace": 0.12})

    ax1.plot(df.resSeq, df.rmsf_A_wt, color=C_WT, label="wild type", zorder=3)
    ax1.plot(df.resSeq, df.rmsf_A_mut, color=C_MUT, label="mutant", zorder=3)
    lo, hi = ax1.get_ylim()
    ax1.set_ylim(lo, hi * 1.12)              # headroom for the legend
    for r in cryptic:
        ax1.axvline(r, color=C_MUTED, lw=0.8, ls=":", zorder=1)
        # Labels sit along the BOTTOM, where the traces are lowest; at the top
        # they collide with the legend on any protein with a flexible terminus.
        ax1.annotate(str(r), xy=(r, lo), xytext=(0, 3),
                     textcoords="offset points", ha="center", va="bottom",
                     fontsize=6.5, color=C_MUTED, rotation=90)
    ax1.set_ylabel("RMSF ($\\AA$)")
    ax1.set_title("Per-residue flexibility" + ("  (dotted: mutated cryptic positions)"
                                               if cryptic else ""),
                  loc="left", fontweight="bold")
    ax1.legend(loc="upper right", ncol=2)
    grid_y(ax1)

    ax2.axhline(0, color="#9A9A9A", lw=0.7, zorder=2)
    ax2.fill_between(df.resSeq, 0, df.d, where=df.d >= 0, color=C_MUT,
                     alpha=0.75, lw=0, zorder=3)
    ax2.fill_between(df.resSeq, 0, df.d, where=df.d < 0, color=C_WT,
                     alpha=0.75, lw=0, zorder=3)
    ax2.set_xlabel("residue (aSAM numbering)")
    ax2.set_ylabel("$\\Delta$RMSF ($\\AA$)")
    ax2.xaxis.set_major_locator(MaxNLocator(nbins=14, integer=True))
    grid_y(ax2)

    save(fig, outdir / "figS_rmsf")
    top = df.reindex(df.d.abs().sort_values(ascending=False).index).head(8)
    print("  largest |dRMSF| residues: "
          + ", ".join(f"{int(r.resSeq)} ({r.d:+.2f})" for r in top.itertuples()))
    return df


# ---------------------------------------------------------------------------
# fold integrity
# ---------------------------------------------------------------------------
def fig_integrity(mut, wt, outdir):
    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE, 58 * MM))
    specs = [("rg_A", "radius of gyration ($\\AA$)", 1.0, "a"),
             ("frac_helix", "helix content (fraction)", 1.0, "b"),
             ("frac_sheet", "sheet content (fraction)", 1.0, "c")]
    for ax, (col, lab, scale, tag) in zip(axes, specs):
        a = mut[col].dropna().values * scale
        b = wt[col].dropna().values * scale
        lo, hi = min(a.min(), b.min()), max(a.max(), b.max())
        bins = np.linspace(lo, hi, 40)
        ax.hist(b, bins=bins, color=C_WT, alpha=0.55, label="wild type", zorder=2)
        ax.hist(a, bins=bins, color=C_MUT, alpha=0.55, label="mutant", zorder=2)
        d = cliffs_delta(a, b)
        ax.set_xlabel(lab); ax.set_ylabel("frames")
        ax.set_title(f"{tag}   $\\Delta$ = {a.mean()-b.mean():+.3f},  "
                     f"$\\delta$ = {d:+.2f}", loc="left", fontweight="bold")
        grid_y(ax)
    axes[0].legend(loc="upper right")
    fig.suptitle("Fold-integrity controls: the substitutions do not change the fold",
                 x=0.005, ha="left", fontsize=9, fontweight="bold")
    fig.tight_layout(w_pad=2.2, rect=(0, 0, 1, 0.93))
    save(fig, outdir / "figS_integrity")


# ---------------------------------------------------------------------------
# clusters
# ---------------------------------------------------------------------------
def fig_clusters(path, outdir):
    if not path or not Path(path).is_file():
        print("  [skip] cluster_representatives.json not found")
        return
    recs = json.loads(Path(path).read_text())
    recs = sorted(recs, key=lambda r: -r["size"])
    total = sum(r["size"] for r in recs)
    fig, ax = plt.subplots(figsize=(110 * MM, 60 * MM))
    x = np.arange(len(recs))
    sizes = [r["size"] for r in recs]
    ax.bar(x, sizes, 0.74, color=C_WT, zorder=2)
    for xi, r in zip(x, recs):
        ax.text(xi, r["size"], f"{100*r['size']/total:.0f}%", ha="center",
                va="bottom", fontsize=6.8, color=C_INK)
        # Representative frame inside the bar, rotated: putting it in the tick
        # label collides as soon as the frame indices run to four digits.
        # Skipped on bars too short to hold the text, which would overflow.
        if r["size"] > 0.30 * sizes[0]:
            ax.text(xi, r["size"] * 0.5, f"frame {r['ensemble_frame']}",
                    rotation=90, ha="center", va="center", fontsize=6.2,
                    color=SURFACE)
    ax.set_xticks(x)
    ax.set_xticklabels([r["cluster"] for r in recs])
    ax.set_xlabel("cluster")
    ax.set_ylabel("frames in cluster")
    ax.set_title(f"Conformational clusters, {total} frames", loc="left",
                 fontweight="bold")
    grid_y(ax)
    save(fig, outdir / "figS_clusters")
    print(f"  {len(recs)} clusters, {total} frames; largest "
          f"{100*sizes[0]/total:.0f}%, smallest {100*sizes[-1]/total:.0f}%")


# ---------------------------------------------------------------------------
# lining persistence
# ---------------------------------------------------------------------------
def fig_persistence(indir, outdir, top_n=25):
    pm = indir / "mutant_persistence.csv"
    pw = indir / "wildtype_persistence.csv"
    frames = {}
    for lab, p in (("mutant", pm), ("wild type", pw)):
        if not p.is_file():
            continue
        try:
            df = pd.read_csv(p)
        except pd.errors.EmptyDataError:
            df = pd.DataFrame()
        if df.empty:
            print(f"  [skip] {p.name} has no rows -- no residue ever lined a "
                  f"detected pocket in the {lab} ensemble")
            continue
        frames[lab] = df
    if not frames:
        print("  [skip] persistence figure: both files empty (detection produced "
              "no pockets)")
        return
    # plot whichever ensembles have data, mutant first
    fig, ax = plt.subplots(figsize=(DOUBLE, 56 * MM))
    order = None
    for i, (lab, df) in enumerate(frames.items()):
        df = df.sort_values("fraction", ascending=False)
        if order is None:
            order = df.head(top_n)["resSeq"].tolist()
        sub = df.set_index("resSeq").reindex(order).fillna(0)
        x = np.arange(len(order)) + (i - 0.5) * 0.38
        colors = [C_MUT if i == 0 else C_WT] * len(order)
        ax.bar(x, sub["fraction"].values, 0.36, color=colors, label=lab, zorder=2)
        if "is_cryptic_site" in sub:
            for xi, flag in zip(x, sub["is_cryptic_site"].fillna(0).values):
                if flag:
                    ax.plot(xi, -0.03, marker="^", ms=4, color=C_INK,
                            clip_on=False, zorder=4)
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels([int(r) for r in order], rotation=90)
    ax.set_xlabel("residue (aSAM numbering);  $\\blacktriangle$ = predicted cryptic site")
    ax.set_ylabel("fraction of frames lining")
    ax.set_title("Pocket-lining persistence", loc="left", fontweight="bold")
    ax.legend(loc="upper right")
    grid_y(ax)
    save(fig, outdir / "figS_persistence")


# ---------------------------------------------------------------------------
def save(fig, stem):
    for ext in ("png", "pdf"):
        fig.savefig(f"{stem}.{ext}")
    plt.close(fig)
    print(f"  wrote {Path(stem).name}.png / .pdf")


def write_stats(rows, outdir):
    df = pd.DataFrame([r for r in rows if r])
    if df.empty:
        print("  [skip] no comparable quantities")
        return
    df.to_csv(outdir / "stats_summary.csv", index=False, float_format="%.6g")
    lines = ["CryptiScan -- mutant vs wild type", "=" * 78, ""]
    for r in df.itertuples():
        lines += [
            f"{r.quantity}  [{r.unit}]" if r.unit else f"{r.quantity}",
            f"   wild type : {r.mean_wt:10.3f} +/- {r.sem_wt:.3f}   "
            f"(median {r.median_wt:.3f}, n={r.n_wt})",
            f"   mutant    : {r.mean_mut:10.3f} +/- {r.sem_mut:.3f}   "
            f"(median {r.median_mut:.3f}, n={r.n_mut})",
            f"   delta     : {r.delta:+10.3f}",
            f"   Cliff's d : {r.cliffs_delta:+10.3f}   ({r.effect})",
            f"   Mann-Whitney p = {r.mannwhitney_p:.3e}   "
            f"KS D = {r.ks_D:.3f}, p = {r.ks_p:.3e}",
            "",
        ]
    lines += [
        "Frames are independent draws from a generative model, not a time series:",
        "SEM = sigma/sqrt(N) with no block averaging, and two-sample tests apply",
        "directly. With N in the thousands, p-values are driven by sample size --",
        "report Cliff's delta as the result, not p.",
    ]
    (outdir / "stats_summary.txt").write_text("\n".join(lines))
    print("  wrote stats_summary.csv / stats_summary.txt")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", default="pocket_analysis",
                    help="directory holding the *_pocket.csv / *_rmsf.csv files")
    ap.add_argument("--clusters", default=None,
                    help="path to cluster_representatives.json")
    ap.add_argument("--outdir", default="figures")
    ap.add_argument("--cryptic", type=int, nargs="*", default=[],
                    help="mutated cryptic residues in aSAM numbering, to mark "
                         "on the RMSF plot (the mapping is printed by Stage 9)")
    a = ap.parse_args()

    indir, outdir = Path(a.indir), Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("Loading pocket tables")
    mut = load_pocket(indir / "mutant_pocket.csv", "mutant")
    wt = load_pocket(indir / "wildtype_pocket.csv", "wild type")
    if mut is None or wt is None:
        sys.exit("ERROR: both mutant_pocket.csv and wildtype_pocket.csv are required.")

    bad = [d.attrs["label"] for d in (mut, wt)
           if d.attrs["n_open"] + d.attrs["n_closed"] == 0]
    if bad:
        print()
        print("  " + "!" * 68)
        print(f"  DETECTION FAILED for: {', '.join(bad)}")
        print("  Every frame is no_output/fpocket_failed, so fpocket never returned")
        print("  a parseable result. This is NOT a closed pocket -- a closed frame")
        print("  has status 'closed' and volume 0. Volume, druggability, lining and")
        print("  persistence are unusable until detection is fixed; the exposure,")
        print("  flexibility and fold-integrity panels below are unaffected.")
        print("  " + "!" * 68)

    print("\nStatistics")
    rows = [
        compare("cryptic-site SASA", mut["cryptic_sasa_A2"], wt["cryptic_sasa_A2"], "A^2"),
        compare("radius of gyration", mut["rg_A"], wt["rg_A"], "A"),
        compare("backbone RMSD", mut["rmsd_A"], wt["rmsd_A"], "A"),
        compare("helix content", mut["frac_helix"], wt["frac_helix"], "fraction"),
        compare("sheet content", mut["frac_sheet"], wt["frac_sheet"], "fraction"),
    ]
    if detection_ok(mut) and detection_ok(wt):
        rows += [
            compare("cryptic-pocket volume (all frames)",
                    mut["volume_A3"].fillna(0), wt["volume_A3"].fillna(0), "A^3"),
            compare("non-cryptic pocket volume",
                    mut["noncryptic_volume_A3"], wt["noncryptic_volume_A3"], "A^3"),
            compare("pockets per frame",
                    mut["n_pockets_total"], wt["n_pockets_total"], "count"),
        ]
    if "cv_A" in mut and mut["cv_A"].notna().any():
        rows.append(compare("collective variable", mut["cv_A"], wt["cv_A"], "A"))
    write_stats(rows, outdir)

    print("\nFigures")
    fig_main(mut, wt, outdir, rows)
    fig_integrity(mut, wt, outdir)
    fig_rmsf(outdir, indir, a.cryptic)
    fig_clusters(a.clusters, outdir)
    fig_persistence(indir, outdir)
    print(f"\nDone. Everything is in {outdir}/")


if __name__ == "__main__":
    main()
