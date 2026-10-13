#!/usr/bin/env python3
"""
make_flowchart.py -- Figure 1 for the CryptiScan manuscript.

Twelve stages in three bands, in the visual language of the previous figure
(numbered cards, a coloured band header per phase, METHOD / OUTPUT footers).
Writes PDF for LaTeX and PNG for slides.

    python3 make_flowchart.py [--outdir Figs]
"""

import argparse
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

# --- palette -----------------------------------------------------------
# One accent per band, all from the same hue family as the original figure
# so the three phases read as one system rather than three charts.
ACCENT = ["#8C1D3F", "#1F5673", "#2F6B3C"]
TINT = ["#FBEFF3", "#EDF4F8", "#EEF5EF"]
INK = "#16181D"
SUB = "#5C6370"
EDGE = "#D9DCE1"
SURF = "#FFFFFF"

mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "savefig.bbox": "tight", "savefig.pad_inches": 0.04,
})

BANDS = [
    ("A", "TARGET ACQUISITION, STRUCTURAL CURATION & CRYPTIC SITE DISCOVERY", [
        ("1", "Target Ingestion",
         "Fetch the entry from RCSB, or read a user-supplied coordinate file.",
         "RCSB API / local file", "<pdb>_raw.pdb"),
        ("2", "Structure Curation",
         "Drop waters, separate ions and ligands, collapse alternate conformers.",
         "MDAnalysis 2.7", "<pdb>_protein.pdb"),
        ("3", "Chain Selection & Checks",
         "Isolate the target chain; verify standard residues and backbone continuity.",
         "MDAnalysis topology", "<pdb>.pdb (curated)"),
        ("4", "Cryptic Pocket Scoring",
         "Score per-residue crypticity from sequence; rank all positions.",
         "CryptoBank ProtT5", "ranked crypticity profile"),
    ]),
    ("B", "DESTABILISING VARIANT SELECTION, MUTAGENESIS & GENERATIVE SAMPLING", [
        ("5", "ESM-Scan Selection",
         "Score every substitution; choose the best within the charged/polar set.",
         "ESM-1v (650M)", r"$\Delta s$ fitness matrix"),
        ("6", "In Silico Mutagenesis",
         "Apply mutations sequentially under dihedral restraints, then relax.",
         "MODELLER 10.8", "<pdb>_mutant.pdb"),
        ("7", "Mutant Ensemble",
         "Latent diffusion draws independent conformations at 320 K.",
         "aSAM / SAM2 (mdCATH)", "2000 states (.dcd)"),
        ("8", "Wild-Type Ensemble",
         "The matched reference, generated under identical settings.",
         "aSAM / SAM2 (mdCATH)", "2000 states (.dcd)"),
    ]),
    ("C", "POCKET QUANTIFICATION, CONTROLS & REPORTING", [
        ("9", "Ensemble Clustering",
         "PCA and k-means on backbone coordinates; one frame per cluster.",
         "scipy + MDAnalysis", "10 representatives"),
        ("10", "Pocket Quantification",
         "Guided per-frame detection at the predicted site, in both ensembles.",
         "fpocket (-P, -u 2)", "per-frame volume, SASA"),
        ("11", "Statistics & Figures",
         "Effect sizes reported against specificity and fold controls.",
         "Cliff's delta, KS, MWU", "figures/, stats_summary"),
        ("12", "Report & Archive",
         "Run manifest, engineered FASTA, every table, compressed.",
         "make_report.py", "<job>.zip"),
    ]),
]


def wrap(text, width):
    words, lines, cur = text.split(), [], ""
    for w in words:
        trial = (cur + " " + w).strip()
        if len(trial) <= width:
            cur = trial
        else:
            lines.append(cur)
            cur = w
    if cur:
        lines.append(cur)
    return lines


def draw(outdir):
    # Geometry in figure units; one unit = one card width.
    ncol, cw, ch = 4, 1.0, 0.60
    gapx, gapy = 0.26, 0.26
    band_h = 0.17
    W = ncol * cw + (ncol - 1) * gapx
    H = len(BANDS) * (band_h + ch + gapy + 0.10) + 0.34

    fig = plt.figure(figsize=(W * 3.05, H * 2.55))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(-0.06, W + 0.06)
    ax.set_ylim(-0.04, H)
    ax.axis("off")

    ax.text(W / 2, H - 0.13, "CryptiScan Computational Pipeline",
            ha="center", va="top", fontsize=15, fontweight="bold", color=INK)

    y = H - 0.46
    for bi, (letter, title, cards) in enumerate(BANDS):
        acc, tint = ACCENT[bi], TINT[bi]

        # band header
        ax.add_patch(Rectangle((0, y - band_h), 0.30, band_h,
                               facecolor=acc, edgecolor="none", zorder=3))
        ax.text(0.15, y - band_h / 2, letter, ha="center", va="center",
                fontsize=11, fontweight="bold", color=SURF, zorder=4)
        ax.text(0.40, y - band_h / 2, title, ha="left", va="center",
                fontsize=9.2, fontweight="bold", color=acc, zorder=4)
        ax.plot([0, W], [y - band_h - 0.035, y - band_h - 0.035],
                color=acc, lw=1.1, zorder=2)

        top = y - band_h - 0.10
        for ci, (num, name, body, method, output) in enumerate(cards):
            x = ci * (cw + gapx)
            bot = top - ch

            ax.add_patch(FancyBboxPatch(
                (x, bot), cw, ch, boxstyle="round,pad=0,rounding_size=0.022",
                facecolor=SURF, edgecolor=EDGE, lw=0.9, zorder=3))
            # title strip
            ax.add_patch(FancyBboxPatch(
                (x, top - 0.165), cw, 0.165,
                boxstyle="round,pad=0,rounding_size=0.022",
                facecolor=tint, edgecolor=EDGE, lw=0.9, zorder=4))
            ax.add_patch(Rectangle((x, top - 0.175), cw, 0.012,
                                   facecolor=tint, edgecolor="none", zorder=4))

            ax.add_patch(plt.Circle((x + 0.075, top - 0.083), 0.036,
                                    facecolor=acc, edgecolor="none", zorder=6))
            ax.text(x + 0.075, top - 0.083, num, ha="center", va="center",
                    fontsize=6.6, fontweight="bold", color=SURF, zorder=7)
            ax.text(x + 0.135, top - 0.083, name, ha="left", va="center",
                    fontsize=8.1, fontweight="bold", color=INK, zorder=7)

            ty = top - 0.215
            for line in wrap(body, 40):
                ax.text(x + 0.055, ty, line, ha="left", va="top",
                        fontsize=6.6, color=SUB, zorder=5)
                ty -= 0.058

            ax.plot([x + 0.055, x + cw - 0.055], [bot + 0.175, bot + 0.175],
                    color=EDGE, lw=0.8, zorder=5)
            ax.text(x + 0.055, bot + 0.118, "METHOD:", ha="left", va="center",
                    fontsize=5.9, fontweight="bold", color=acc, zorder=5)
            ax.text(x + 0.245, bot + 0.118, method, ha="left", va="center",
                    fontsize=6.3, color=INK, zorder=5)
            ax.text(x + 0.055, bot + 0.050, "OUTPUT:", ha="left", va="center",
                    fontsize=5.9, fontweight="bold", color=acc, zorder=5)
            ax.text(x + 0.245, bot + 0.050, output, ha="left", va="center",
                    fontsize=6.3, color=INK, zorder=5)

            if ci < ncol - 1:
                ax.add_patch(FancyArrowPatch(
                    (x + cw + 0.055, bot + ch / 2),
                    (x + cw + gapx - 0.055, bot + ch / 2),
                    arrowstyle="-|>", mutation_scale=9, lw=1.0,
                    color=acc, zorder=5))

        # connector into the next band
        if bi < len(BANDS) - 1:
            ymid = top - ch - gapy / 2
            ax.plot([0.10, W - 0.10], [ymid, ymid], color=EDGE, lw=0.8,
                    zorder=1)
            ax.add_patch(FancyArrowPatch(
                (W / 2, ymid + 0.055), (W / 2, ymid - 0.055),
                arrowstyle="-|>", mutation_scale=10, lw=1.2,
                color=ACCENT[bi + 1], zorder=4))

        y = top - ch - gapy

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for ext, dpi in (("pdf", None), ("png", 400)):
        fig.savefig(outdir / f"flowchart.{ext}", **({"dpi": dpi} if dpi else {}))
        print(f"  wrote {outdir / ('flowchart.' + ext)}")
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="Figs")
    draw(ap.parse_args().outdir)
