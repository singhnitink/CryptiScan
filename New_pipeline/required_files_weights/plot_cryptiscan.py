#!/usr/bin/env python3
"""
plot_cryptiscan.py -- publication figures from CryptiScan Stage 9 output
=======================================================================
Reads the per-frame CSVs the quantification stage writes and produces
journal-ready figures plus a statistics table.

FIGURES
  fig2_opening.png/.pdf       main text, three panels:
                                (a) cryptic-site exposure distributions
                                (b) ECDFs -- makes the KS separation visible
                                (c) effect sizes, signal against every control
  figS1_integrity.png/.pdf    fold-integrity controls
  figS2_rmsf.png/.pdf         per-residue flexibility and its difference
  figS3_persistence.png/.pdf  pocket-lining persistence (needs fpocket output)
  figS4_clusters.png/.pdf     cluster populations (not cited in the main text)
  stats_summary.csv / .txt    every comparison, with effect sizes

Panels whose input is missing or entirely undetected are SKIPPED with a printed
reason rather than drawn empty, so a broken upstream stage cannot masquerade as
a null result.

USAGE
    python plot_cryptiscan.py --indir pocket_analysis \
        --clusters clustering/cluster_representatives.json \
        --cryptic 200 205 233 235 255 \
        --label "TEM-1 beta-lactamase"

Needs only numpy, pandas, scipy and matplotlib, so it runs in Colab unchanged.
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
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
from scipy import stats

# ---------------------------------------------------------------------------
# Visual system
#
# Two series, so two hues, assigned in fixed order and never cycled. The pair is
# validated colourblind-safe: protan dE 20.9, tritan 28.1, normal-vision 27.3,
# all above the dE >= 8 CVD floor and the 15 normal-vision floor. Identity is
# never carried by colour alone -- every panel is legended or directly labelled.
# ---------------------------------------------------------------------------
C_WT   = "#2B6CB0"      # wild type
C_MUT  = "#C05621"      # mutant
C_INK  = "#16181D"      # primary text
C_SUB  = "#5C6370"      # secondary text
C_LINE = "#B6BBC4"      # axis lines
C_GRID = "#E8EAEE"
SURF   = "#FFFFFF"
BAND   = "#EDEFF3"      # "negligible" shading

mpl.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42,          # editable text in Illustrator / Inkscape
    "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 7.5,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "axes.labelcolor": C_INK,
    "axes.edgecolor": C_LINE,
    "axes.linewidth": 0.6,
    "axes.facecolor": SURF,
    "figure.facecolor": SURF,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.labelpad": 3.0,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "xtick.color": C_SUB,
    "ytick.color": C_SUB,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "xtick.minor.size": 1.6,
    "ytick.minor.size": 1.6,
    "xtick.minor.width": 0.5,
    "ytick.minor.width": 0.5,
    "legend.frameon": False,
    "legend.fontsize": 7,
    "legend.handlelength": 1.1,
    "legend.handletextpad": 0.5,
    "legend.labelspacing": 0.3,
    "lines.linewidth": 1.3,
    "lines.solid_capstyle": "round",
    "grid.color": C_GRID,
    "grid.linewidth": 0.5,
    "mathtext.fontset": "dejavusans",
})

MM = 1 / 25.4
SINGLE = 88 * MM          # single journal column
ONEHALF = 130 * MM
DOUBLE = 180 * MM         # full width


def panel_label(ax, letter, dx=-0.16, dy=1.06):
    """Bold panel letter set outside the axes, journal style."""
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=9.5,
            fontweight="bold", va="top", ha="left", color=C_INK)


def fig_label(fig, x, y, letter):
    """Panel letter in FIGURE coordinates. Axes-relative placement fails on a
    panel whose tick labels are long -- the letter lands in its neighbour."""
    fig.text(x, y, letter, fontsize=9.5, fontweight="bold", va="top",
             ha="left", color=C_INK)


def tidy(ax, grid="y", minor=True):
    if grid in ("y", "both"):
        ax.yaxis.grid(True, zorder=0)
    if grid in ("x", "both"):
        ax.xaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)
    if minor:
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))


def note(ax, text, xy=(0.97, 0.97), ha="right", va="top", size=6.8):
    """Statistics block with a soft backing so it never fights the data."""
    ax.text(*xy, text, transform=ax.transAxes, ha=ha, va=va, fontsize=size,
            color=C_INK, linespacing=1.45,
            bbox=dict(boxstyle="round,pad=0.42", fc=SURF, ec=C_GRID, lw=0.5,
                      alpha=0.93))


# ---------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------
def cliffs_delta(a, b):
    """P(a > b) - P(b > a), from the Mann-Whitney U. Non-parametric, because
    these distributions are skewed and, for volumes, zero-inflated."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.size == 0 or b.size == 0:
        return np.nan
    u = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
    return 2.0 * u / (a.size * b.size) - 1.0


def magnitude(d):
    if np.isnan(d):
        return "n/a"
    a = abs(d)
    return ("negligible" if a < 0.147 else "small" if a < 0.33
            else "medium" if a < 0.474 else "large")


def compare(name, mut, wt, unit=""):
    mut = np.asarray(mut, float); mut = mut[~np.isnan(mut)]
    wt = np.asarray(wt, float);   wt = wt[~np.isnan(wt)]
    if mut.size == 0 or wt.size == 0:
        return None
    sem = lambda x: x.std(ddof=1) / np.sqrt(x.size) if x.size > 1 else np.nan
    ks = stats.ks_2samp(mut, wt)
    return {
        "quantity": name, "unit": unit, "n_mut": mut.size, "n_wt": wt.size,
        "mean_mut": mut.mean(), "sem_mut": sem(mut),
        "mean_wt": wt.mean(), "sem_wt": sem(wt),
        "delta": mut.mean() - wt.mean(),
        "median_mut": float(np.median(mut)), "median_wt": float(np.median(wt)),
        # Frames are independent draws: plain two-sample tests apply directly,
        # with no block averaging or correlation-time correction.
        "mannwhitney_p": stats.mannwhitneyu(mut, wt, alternative="two-sided").pvalue,
        "ks_D": ks.statistic, "ks_p": ks.pvalue,
        "cliffs_delta": cliffs_delta(mut, wt),
        "effect": magnitude(cliffs_delta(mut, wt)),
    }


def fmt_p(p):
    if p == 0 or p < 1e-250:
        return r"$p < 10^{-250}$"
    return rf"$p = {p:.0e}$".replace("e-0", r"\times10^{-").replace("e-", r"\times10^{-") + "}$" \
        if False else rf"$p = {p:.1e}$"


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------
def load_pocket(path, label):
    if not path.is_file():
        print(f"  [skip] {path.name} not found")
        return None
    df = pd.read_csv(path)
    counts = df["status"].value_counts().to_dict()
    print(f"  {label:10s} {len(df):5d} frames   "
          + ", ".join(f"{k}={v}" for k, v in counts.items()))
    df.attrs.update(label=label,
                    n_open=int((df.status == "open").sum()),
                    n_closed=int((df.status == "closed").sum()))
    return df


def detection_ok(df):
    """True when fpocket returned a parseable result for at least one frame.

    A run in which EVERY frame is no_output/fpocket_failed is a TOOL FAILURE,
    not a closed pocket: a genuinely closed frame has status 'closed' and a
    volume of 0. The two must never be conflated in a figure.
    """
    return df is not None and (df.attrs["n_open"] + df.attrs["n_closed"]) > 0


# ---------------------------------------------------------------------------
# FIGURE 2 -- main text
# ---------------------------------------------------------------------------
def fig_main(mut, wt, outdir, label):
    a = mut["cryptic_sasa_A2"].dropna().values
    b = wt["cryptic_sasa_A2"].dropna().values
    have_vol = detection_ok(mut) and detection_ok(wt)

    # With cavity data there are two independent opening observables --
    # exposure (detector-free) and volume (detector-based) -- and both
    # deserve a distribution panel, so the figure becomes 2x2. Without it
    # the old 1x3 layout is kept, with the ECDF standing in for volume.
    if have_vol:
        fig = plt.figure(figsize=(DOUBLE, 104 * MM))
        gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.0],
                              height_ratios=[1.0, 1.0],
                              wspace=0.30, hspace=0.52,
                              left=0.075, right=0.985,
                              top=0.92, bottom=0.085)
        slots = [gs[0, 0], gs[0, 1], gs[1, 0], gs[1, 1]]
    else:
        fig = plt.figure(figsize=(DOUBLE, 58 * MM))
        gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 1.15],
                              wspace=0.42, left=0.055, right=0.995,
                              top=0.88, bottom=0.17)
        slots = [gs[0], gs[1], gs[2]]

    # ---- (a) exposure distributions -------------------------------------
    ax = fig.add_subplot(slots[0])
    lo, hi = min(a.min(), b.min()), max(a.max(), b.max())
    pad = 0.04 * (hi - lo)
    grid = np.linspace(lo - pad, hi + pad, 400)
    bins = np.linspace(lo, hi, 44)

    for v, c, lab in ((b, C_WT, "wild type"), (a, C_MUT, "mutant")):
        ax.hist(v, bins=bins, density=True, color=c, alpha=0.22,
                zorder=2, lw=0)
        kde = stats.gaussian_kde(v)
        ax.plot(grid, kde(grid), color=c, lw=1.5, zorder=4, label=lab)
        ax.fill_between(grid, 0, kde(grid), color=c, alpha=0.14, lw=0, zorder=3)

    ymax = ax.get_ylim()[1]
    for v, c in ((b, C_WT), (a, C_MUT)):
        ax.plot([v.mean(), v.mean()], [0, ymax * 0.985], color=c, lw=1.0,
                ls=(0, (3, 2)), zorder=5)
    # The arrow and its label are lifted clear of the taller KDE peak and of
    # the legend: with well-separated distributions the span is wide, so a
    # label centred on it would otherwise sit on top of a curve.
    y_arrow = ymax * 1.07
    ax.annotate("", xy=(a.mean(), y_arrow), xytext=(b.mean(), y_arrow),
                arrowprops=dict(arrowstyle="-|>", color=C_INK, lw=0.9,
                                shrinkA=0, shrinkB=0, mutation_scale=8),
                zorder=6)
    ax.text((a.mean() + b.mean()) / 2, ymax * 1.10,
            rf"$\Delta = {a.mean()-b.mean():+.1f}\,\mathrm{{\AA}}^2$",
            ha="center", va="bottom", fontsize=7, color=C_INK, zorder=6)

    ax.set_xlabel(r"cryptic-site SASA  ($\mathrm{\AA}^2$)")
    ax.set_ylabel("probability density")
    ax.set_ylim(0, ymax * 1.26)
    ax.legend(loc="upper right", bbox_to_anchor=(1.02, 0.88))
    tidy(ax)

    # ---- (b) pocket volume, when a detector ran --------------------------
    # Volume is zero in a closed frame, so the distribution is a spike at
    # zero plus the open population. Plotting only the open frames would
    # hide exactly the thing that changed, so closed frames are kept at
    # zero and the open fraction is stated on the panel.
    if have_vol:
        ax = fig.add_subplot(slots[1])
        va = mut["volume_A3"].fillna(0.0).values
        vb = wt["volume_A3"].fillna(0.0).values
        vhi = max(np.percentile(va, 99.5), np.percentile(vb, 99.5))
        vbins = np.linspace(0, vhi, 46)
        for v, c, lab in ((vb, C_WT, "wild type"), (va, C_MUT, "mutant")):
            ax.hist(np.clip(v, 0, vhi), bins=vbins, density=True, color=c,
                    alpha=0.42, lw=0, zorder=2, label=lab)
        vy = ax.get_ylim()[1]
        for v, c in ((vb, C_WT), (va, C_MUT)):
            ax.plot([v.mean()] * 2, [0, vy * 0.97], color=c, lw=1.0,
                    ls=(0, (3, 2)), zorder=5)
        fo_m = float((mut["status"] == "open").mean())
        fo_w = float((wt["status"] == "open").mean())
        ax.set_ylim(0, vy * 1.30)
        ax.text(0.97, 0.96,
                f"open frames\nmutant {fo_m:.0%}\nwild type {fo_w:.0%}",
                transform=ax.transAxes, ha="right", va="top", fontsize=6.6,
                color=C_INK, linespacing=1.45)
        ax.set_xlabel(r"cryptic-pocket volume  ($\mathrm{\AA}^3$)")
        ax.set_ylabel("probability density")
        ax.legend(loc="upper left", fontsize=6.6)
        tidy(ax)

    # ---- ECDF -- makes the KS separation visible ------------------------
    ax = fig.add_subplot(slots[2] if have_vol else gs[1])
    for v, c, lab in ((b, C_WT, "wild type"), (a, C_MUT, "mutant")):
        x = np.sort(v)
        y = np.arange(1, x.size + 1) / x.size
        ax.step(x, y, where="post", color=c, lw=1.5, zorder=4, label=lab)

    # the KS statistic, drawn where it is attained
    xs = np.sort(np.concatenate([a, b]))
    fa = np.searchsorted(np.sort(a), xs, side="right") / a.size
    fb = np.searchsorted(np.sort(b), xs, side="right") / b.size
    k = int(np.argmax(np.abs(fa - fb)))
    ax.plot([xs[k], xs[k]], [fa[k], fb[k]], color=C_INK, lw=1.0, zorder=5)
    ax.plot([xs[k]] * 2, [fa[k], fb[k]], "o", ms=3.0, color=C_INK, zorder=6)
    ax.annotate(rf"$D = {abs(fa[k]-fb[k]):.2f}$",
                xy=(xs[k], (fa[k] + fb[k]) / 2), xytext=(7, 0),
                textcoords="offset points", fontsize=7, va="center",
                color=C_INK, zorder=6)

    ax.set_xlabel(r"cryptic-site SASA  ($\mathrm{\AA}^2$)")
    ax.set_ylabel("cumulative fraction")
    ax.set_ylim(-0.02, 1.02)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    tidy(ax)

    # ---- (c) effect sizes -------------------------------------------------
    # One axis, one unit. Cliff's delta is dimensionless, so the opening signal
    # and every control are directly comparable on a single scale -- which a
    # grouped bar chart of raw means (Angstroms beside percentages) could not
    # do without a second y-scale.
    ax = fig.add_subplot(slots[3] if have_vol else gs[2])
    items = [("cryptic-site SASA", mut["cryptic_sasa_A2"], wt["cryptic_sasa_A2"], True)]
    if have_vol:
        items += [
            ("pocket volume", mut["volume_A3"].fillna(0), wt["volume_A3"].fillna(0), True),
            # Druggability is shown alongside the opening signal, not hidden
            # among the controls, because for a charge-driven opening it can
            # move the OTHER way: fpocket rewards apolar lining, and the
            # arginines that force the site open make it polar. A negative
            # bar here is a real property of the mutant pocket, not a
            # failure of the opening.
            # OPEN FRAMES ONLY. Scoring a closed frame as druggability 0
            # would mix "how often is it open" into "how good is it when it
            # is open", and those are separate questions -- the open
            # fractions are already on panel (b).
            ("druggability (open frames)",
             mut.loc[mut["status"] == "open", "drug_score"],
             wt.loc[wt["status"] == "open", "drug_score"], True),
            ("non-cryptic volume", mut["noncryptic_volume_A3"],
             wt["noncryptic_volume_A3"], False)]
    items += [
        ("backbone RMSD", mut["rmsd_A"], wt["rmsd_A"], False),
        ("radius of gyration", mut["rg_A"], wt["rg_A"], False),
        ("helix content", mut["frac_helix"], wt["frac_helix"], False),
        ("sheet content", mut["frac_sheet"], wt["frac_sheet"], False)]

    labels, deltas, sig = [], [], []
    for nm, x, y, s in items:
        x = np.asarray(x, float); x = x[~np.isnan(x)]
        y = np.asarray(y, float); y = y[~np.isnan(y)]
        if x.size and y.size:
            labels.append(nm); deltas.append(cliffs_delta(x, y)); sig.append(s)

    ypos = np.arange(len(labels))[::-1]
    ax.axvspan(-0.147, 0.147, color=BAND, zorder=1, lw=0)
    ax.axvline(0, color=C_LINE, lw=0.7, zorder=2)
    for yi, d, s in zip(ypos, deltas, sig):
        c = C_MUT if s else C_WT
        ax.plot([0, d], [yi, yi], color=c, lw=1.5, zorder=3)
        ax.plot(d, yi, "o", ms=5.5, color=c, mec=SURF, mew=1.1, zorder=4)
        # A positive value's label sits outside the point, to the right, where
        # there is empty margin. A negative one cannot: the left margin holds
        # the y tick labels, so an outward label lands on top of them. Negative
        # labels are therefore placed on the inner side, toward zero.
        if d >= 0:
            tx, ha = d + 0.07, "left"
        else:
            tx, ha = d + 0.07, "left"
        ax.text(tx, yi, f"{d:+.2f}", ha=ha, va="center",
                fontsize=6.8, color=C_INK, zorder=5,
                bbox=dict(fc=SURF, ec="none", pad=0.8, alpha=0.85))
    ax.set_yticks(ypos); ax.set_yticklabels(labels)
    ax.set_xlim(-1.03, 1.03); ax.set_ylim(-0.75, len(labels) - 0.25)
    ax.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax.set_xlabel(r"Cliff's $\delta$   (mutant $-$ wild type)")
    ax.text(0, -0.70, "negligible", ha="center", va="bottom", fontsize=6.3,
            color=C_SUB, style="italic")
    ax.xaxis.grid(True, zorder=0); ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))

    # Panel letters follow the layout actually used: four on a 2x2 grid when
    # there is cavity data, three across when there is not.
    if have_vol:
        places = ((0.008, 0.975, "a"), (0.508, 0.975, "b"),
                  (0.008, 0.505, "c"), (0.508, 0.505, "d"))
    else:
        places = ((0.012, 0.955, "a"), (0.368, 0.955, "b"),
                  (0.648, 0.955, "c"))
    for xx, yy, ll in places:
        fig_label(fig, xx, yy, ll)

    if label:
        fig.suptitle(label, x=0.012, y=1.035, ha="left", fontsize=8,
                     fontweight="bold", color=C_INK)
    save(fig, outdir / "fig2_opening")
    if not have_vol:
        print("  note: no cavity data -- volume rows omitted from panel (c)")


# ---------------------------------------------------------------------------
def fig_rmsf(outdir, indir, cryptic):
    pm, pw = indir / "mutant_rmsf.csv", indir / "wildtype_rmsf.csv"
    if not (pm.is_file() and pw.is_file()):
        print("  [skip] RMSF files not found")
        return
    d = pd.read_csv(pm).merge(pd.read_csv(pw), on="resSeq",
                              suffixes=("_mut", "_wt"))
    d["diff"] = d.rmsf_A_mut - d.rmsf_A_wt

    fig, (a1, a2) = plt.subplots(
        2, 1, figsize=(DOUBLE, 72 * MM), sharex=True,
        gridspec_kw=dict(height_ratios=[2.5, 1], hspace=0.10,
                         left=0.065, right=0.995, top=0.93, bottom=0.12))

    for r in cryptic:
        for ax in (a1, a2):
            ax.axvline(r, color=C_SUB, lw=0.5, ls=(0, (1, 2.5)), zorder=1)
    a1.plot(d.resSeq, d.rmsf_A_wt, color=C_WT, lw=1.2, zorder=3, label="wild type")
    a1.plot(d.resSeq, d.rmsf_A_mut, color=C_MUT, lw=1.2, zorder=4, label="mutant")
    lo, hi = a1.get_ylim(); a1.set_ylim(lo, hi * 1.30)
    # Labels are staggered: mutated positions a few residues apart would
    # otherwise print on top of each other. Two levels are not enough -- a
    # run of three or more close positions (e.g. 140/141/144) needs the
    # level to keep climbing, so it advances per consecutive close neighbour
    # and resets once a gap is wide enough to clear the previous label.
    srt = sorted(cryptic)
    gap = 0.06 * (d.resSeq.max() or 1)
    levels, lvl = [], 0
    for j, r in enumerate(srt):
        lvl = lvl + 1 if (j > 0 and (r - srt[j - 1]) < gap) else 0
        levels.append(lvl % 3)
    for j, r in enumerate(srt):
        ty = hi * (1.045 + 0.070 * levels[j])
        a1.plot(r, hi * 1.02, marker="v", ms=3.4, color=C_INK, zorder=5,
                clip_on=False)
        a1.plot([r, r], [hi * 1.02, ty], color=C_SUB, lw=0.4, zorder=4,
                clip_on=False)
        a1.text(r, ty, str(r), ha="center", va="bottom", fontsize=6,
                color=C_SUB, clip_on=False)
    a1.set_ylabel(r"RMSF  ($\mathrm{\AA}$)")
    a1.legend(loc="upper left", ncol=2, bbox_to_anchor=(0.0, 1.0))
    tidy(a1, minor=False)
    panel_label(a1, "a", dx=-0.055, dy=1.10)

    a2.axhline(0, color=C_LINE, lw=0.6, zorder=2)
    a2.fill_between(d.resSeq, 0, d["diff"], where=d["diff"] >= 0, color=C_MUT,
                    alpha=0.85, lw=0, zorder=3, interpolate=True)
    a2.fill_between(d.resSeq, 0, d["diff"], where=d["diff"] < 0, color=C_WT,
                    alpha=0.85, lw=0, zorder=3, interpolate=True)
    a2.set_xlabel("residue (internal numbering)")
    a2.set_ylabel(r"$\Delta$RMSF  ($\mathrm{\AA}$)")
    a2.xaxis.set_major_locator(MaxNLocator(nbins=14, integer=True))
    tidy(a2, minor=False)
    panel_label(a2, "b", dx=-0.055, dy=1.22)

    save(fig, outdir / "figS2_rmsf")
    top = d.reindex(d["diff"].abs().sort_values(ascending=False).index).head(6)
    print("  largest |dRMSF|: "
          + ", ".join(f"{int(r.resSeq)} ({getattr(r,'diff'):+.2f})"
                      for r in top.itertuples()))


# ---------------------------------------------------------------------------
def fig_integrity(mut, wt, outdir):
    specs = [("rg_A", r"radius of gyration  ($\mathrm{\AA}$)", "a"),
             ("rmsd_A", r"backbone RMSD  ($\mathrm{\AA}$)", "b"),
             ("frac_helix", "helix content  (fraction)", "c"),
             ("frac_sheet", "sheet content  (fraction)", "d")]
    fig, axes = plt.subplots(1, 4, figsize=(DOUBLE, 46 * MM))
    fig.subplots_adjust(wspace=0.36, left=0.05, right=0.995, top=0.80, bottom=0.22)
    for ax, (col, lab, tag) in zip(axes, specs):
        x = mut[col].dropna().values
        y = wt[col].dropna().values
        lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
        bins = np.linspace(lo, hi, 42)
        for v, c in ((y, C_WT), (x, C_MUT)):
            ax.hist(v, bins=bins, density=True, color=c, alpha=0.5, lw=0, zorder=2)
        d = cliffs_delta(x, y)
        ax.set_xlabel(lab)
        ax.set_title(rf"$\Delta = {x.mean()-y.mean():+.3f}$   "
                     rf"$\delta = {d:+.2f}$", fontsize=7.2, color=C_INK, pad=4)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        tidy(ax, grid="x", minor=False)
        panel_label(ax, tag, dx=-0.10, dy=1.30)
    axes[0].set_ylabel("density")
    h = [plt.Rectangle((0, 0), 1, 1, fc=C_WT, alpha=0.5, lw=0),
         plt.Rectangle((0, 0), 1, 1, fc=C_MUT, alpha=0.5, lw=0)]
    axes[-1].legend(h, ["wild type", "mutant"], loc="upper right")
    save(fig, outdir / "figS1_integrity")


# ---------------------------------------------------------------------------
def fig_clusters(path, outdir):
    if not path or not Path(path).is_file():
        print("  [skip] cluster_representatives.json not found")
        return
    recs = sorted(json.loads(Path(path).read_text()), key=lambda r: -r["size"])
    total = sum(r["size"] for r in recs)
    fig, ax = plt.subplots(figsize=(ONEHALF, 50 * MM))
    fig.subplots_adjust(left=0.10, right=0.99, top=0.90, bottom=0.20)
    x = np.arange(len(recs))
    sizes = [r["size"] for r in recs]
    ax.bar(x, sizes, 0.70, color=C_WT, zorder=3, lw=0)
    for xi, r in zip(x, recs):
        ax.text(xi, r["size"] + total * 0.004, f"{100*r['size']/total:.1f}%",
                ha="center", va="bottom", fontsize=6.6, color=C_INK)
        if r["size"] > 0.33 * sizes[0]:
            ax.text(xi, r["size"] * 0.5, f"frame {r['ensemble_frame']}",
                    rotation=90, ha="center", va="center", fontsize=6,
                    color=SURF)
    ax.set_xticks(x); ax.set_xticklabels([r["cluster"] for r in recs])
    ax.set_xlabel("cluster, ordered by population")
    ax.set_ylabel("frames")
    ax.set_ylim(0, max(sizes) * 1.16)
    ax.set_title(f"{len(recs)} clusters from {total} frames", fontsize=8,
                 loc="left", color=C_INK, pad=6)
    tidy(ax, minor=False)
    save(fig, outdir / "figS4_clusters")
    print(f"  {len(recs)} clusters; largest {100*sizes[0]/total:.1f}%, "
          f"smallest {100*sizes[-1]/total:.1f}%")


# ---------------------------------------------------------------------------
def fig_persistence(indir, outdir, top_n=24):
    frames = {}
    for lab, fn in (("mutant", "mutant_persistence.csv"),
                    ("wild type", "wildtype_persistence.csv")):
        p = indir / fn
        if not p.is_file():
            continue
        try:
            df = pd.read_csv(p)
        except pd.errors.EmptyDataError:
            df = pd.DataFrame()
        if df.empty:
            print(f"  [skip] {fn}: no residue ever lined a detected pocket "
                  f"({lab})")
            continue
        frames[lab] = df
    if not frames:
        print("  [skip] persistence figure: no cavity data")
        return
    fig, ax = plt.subplots(figsize=(DOUBLE, 48 * MM))
    fig.subplots_adjust(left=0.07, right=0.99, top=0.90, bottom=0.26)
    order = None
    for i, (lab, df) in enumerate(frames.items()):
        df = df.sort_values("fraction", ascending=False)
        if order is None:
            order = df.head(top_n)["resSeq"].tolist()
        sub = df.set_index("resSeq").reindex(order).fillna(0)
        ax.bar(np.arange(len(order)) + (i - 0.5) * 0.40,
               sub["fraction"].values, 0.38,
               color=C_MUT if i == 0 else C_WT, label=lab, zorder=3, lw=0)
        if "is_cryptic_site" in sub:
            for xi, f in zip(np.arange(len(order)) + (i - 0.5) * 0.40,
                             sub["is_cryptic_site"].fillna(0).values):
                if f and i == 0:
                    ax.plot(xi, -0.035, marker="^", ms=3.6, color=C_INK,
                            clip_on=False, zorder=4)
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels([int(r) for r in order], rotation=90)
    ax.set_xlabel(r"residue;  $\blacktriangle$ = predicted cryptic site")
    ax.set_ylabel("fraction of frames lining")
    ax.legend(loc="upper right")
    tidy(ax, minor=False)
    save(fig, outdir / "figS3_persistence")


# ---------------------------------------------------------------------------
def save(fig, stem):
    for ext in ("png", "pdf"):
        fig.savefig(f"{stem}.{ext}")
    plt.close(fig)
    print(f"  wrote {Path(stem).name}.png / .pdf")


def write_stats(rows, outdir):
    df = pd.DataFrame([r for r in rows if r])
    if df.empty:
        print("  [skip] nothing comparable")
        return
    df.to_csv(outdir / "stats_summary.csv", index=False, float_format="%.6g")
    L = ["CryptiScan -- mutant versus wild type", "=" * 76, ""]
    for r in df.itertuples():
        L += [f"{r.quantity}" + (f"  [{r.unit}]" if r.unit else ""),
              f"   wild type : {r.mean_wt:10.3f} +/- {r.sem_wt:.3f}"
              f"   (median {r.median_wt:.3f}, n={r.n_wt})",
              f"   mutant    : {r.mean_mut:10.3f} +/- {r.sem_mut:.3f}"
              f"   (median {r.median_mut:.3f}, n={r.n_mut})",
              f"   delta     : {r.delta:+10.3f}",
              f"   Cliff's d : {r.cliffs_delta:+10.3f}   ({r.effect})",
              f"   Mann-Whitney p = {r.mannwhitney_p:.3e}"
              f"   KS D = {r.ks_D:.3f}, p = {r.ks_p:.3e}", ""]
    L += ["Frames are independent draws from a generative model, not a time",
          "series: SEM = sigma/sqrt(N) with no block averaging, and two-sample",
          "tests apply directly. With N in the thousands a p-value is driven by",
          "sample size -- report Cliff's delta as the result."]
    (outdir / "stats_summary.txt").write_text("\n".join(L))
    print("  wrote stats_summary.csv / .txt")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", default="pocket_analysis")
    ap.add_argument("--clusters", default=None)
    ap.add_argument("--outdir", default="figures")
    ap.add_argument("--cryptic", type=int, nargs="*", default=[],
                    help="mutated residues in internal (aSAM) numbering. Omit "
                         "it: they are read from pocket_summary.json, or from "
                         "mutant_persistence.csv, whichever Stage 10 wrote.")
    ap.add_argument("--label", default="",
                    help="target name, set as the figure-2 running head")
    args = ap.parse_args()

    indir, outdir = Path(args.indir), Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # The cryptic residues are recorded by Stage 10 in pocket_summary.json, so
    # the launcher never has to pass them and the user never has to know them.
    # An explicit --cryptic still wins, for re-plotting an older run by hand.
    if not args.cryptic:
        summ = indir / "pocket_summary.json"
        if summ.is_file():
            try:
                args.cryptic = [int(r) for r in
                                json.loads(summ.read_text()).get(
                                    "cryptic_residues", [])]
                if args.cryptic:
                    print(f"  cryptic residues from {summ.name}: "
                          + ", ".join(str(r) for r in args.cryptic))
            except Exception as exc:
                print(f"  [warn] could not read {summ.name}: {exc}")
        else:
            # Fall back to the persistence table, which flags them per residue.
            pers = indir / "mutant_persistence.csv"
            if pers.is_file():
                try:
                    d = pd.read_csv(pers)
                    if "is_cryptic_site" in d:
                        args.cryptic = sorted(
                            int(r) for r in
                            d.loc[d.is_cryptic_site == 1, "resSeq"])
                        if args.cryptic:
                            print("  cryptic residues from "
                                  f"{pers.name}: "
                                  + ", ".join(str(r) for r in args.cryptic))
                except Exception:
                    pass

    print("Loading")
    mut = load_pocket(indir / "mutant_pocket.csv", "mutant")
    wt = load_pocket(indir / "wildtype_pocket.csv", "wild type")
    if mut is None or wt is None:
        sys.exit("ERROR: both mutant_pocket.csv and wildtype_pocket.csv are required.")

    dead = [d.attrs["label"] for d in (mut, wt) if not detection_ok(d)]
    if dead:
        print("\n  " + "!" * 70)
        print(f"  CAVITY DETECTION PRODUCED NOTHING for: {', '.join(dead)}")
        print("  Every frame is no_output/fpocket_failed, so fpocket never returned")
        print("  a parseable result. This is a TOOL FAILURE, not a closed pocket:")
        print("  a genuinely closed frame has status 'closed' and volume 0.")
        print("  Volume, druggability and persistence are therefore omitted from")
        print("  the figures. Exposure, flexibility and fold integrity are")
        print("  unaffected and are plotted in full.")
        print("  " + "!" * 70)

    print("\nStatistics")
    rows = [compare("cryptic-site SASA", mut.cryptic_sasa_A2, wt.cryptic_sasa_A2, "A^2"),
            compare("radius of gyration", mut.rg_A, wt.rg_A, "A"),
            compare("backbone RMSD", mut.rmsd_A, wt.rmsd_A, "A"),
            compare("helix content", mut.frac_helix, wt.frac_helix, "fraction"),
            compare("sheet content", mut.frac_sheet, wt.frac_sheet, "fraction")]
    if detection_ok(mut) and detection_ok(wt):
        rows += [compare("cryptic-pocket volume (all frames)",
                         mut.volume_A3.fillna(0), wt.volume_A3.fillna(0), "A^3"),
                 compare("non-cryptic pocket volume",
                         mut.noncryptic_volume_A3, wt.noncryptic_volume_A3, "A^3")]
    if "cv_A" in mut and mut.cv_A.notna().any():
        rows.append(compare("collective variable", mut.cv_A, wt.cv_A, "A"))
    write_stats(rows, outdir)

    print("\nFigures")
    fig_main(mut, wt, outdir, args.label)
    fig_rmsf(outdir, indir, args.cryptic)
    fig_integrity(mut, wt, outdir)
    fig_clusters(args.clusters, outdir)
    fig_persistence(indir, outdir)
    print(f"\nDone -> {outdir}/")


if __name__ == "__main__":
    main()
