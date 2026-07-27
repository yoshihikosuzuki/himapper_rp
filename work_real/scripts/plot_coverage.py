#!/usr/bin/env python3
"""Coverage-depth tracks over MAT/PAT from normalized mapping tables.

Depth is accumulated exactly (per read: +1 over [start,end) via a difference array +
cumsum), so it is a true depth track, not a midpoint histogram.  Same sample (HG002
reads on the HG002 assembly) => UNIFORM coverage is GOOD.

Usage:
  plot_coverage.py [--window-bp BP] [--xlo-mb MB] [--xhi-mb MB]
      <out.png> name=mapping.tsv [...]
    - --window-bp defaults to 1000.
    - The x range defaults to the full configured contig length.
    - One or more explicitly named mapping tables are required.
"""
import argparse
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

contigs = {"chr6_MATERNAL": 174742601, "chr6_PATERNAL": 174411882}
FIXED_COLORS = {"himapper": "#d62728", "minimap2": "#1f77b4", "winnowmap2": "#2ca02c"}
PALETTE = ["#9467bd", "#ff7f0e", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22"]

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--window-bp", type=int, default=1000,
                    help="coverage window size in bp [1000]")
parser.add_argument("--xlo-mb", type=float, default=0,
                    help="left edge of the plotted region in Mb [0]")
parser.add_argument("--xhi-mb", type=float,
                    help="right edge of the plotted region in Mb [full contig]")
parser.add_argument("output")
parser.add_argument("tracks", nargs="+", metavar="name=mapping.tsv")
args = parser.parse_intermixed_args()
if args.window_bp <= 0:
    parser.error("--window-bp must be positive")
if args.xlo_mb < 0:
    parser.error("--xlo-mb must be non-negative")

W = args.window_bp
OUT = args.output
XLO = args.xlo_mb * 1e6
XHI = args.xhi_mb * 1e6 if args.xhi_mb is not None else max(contigs.values())
if XHI <= XLO:
    parser.error("--xhi-mb must be greater than --xlo-mb")
if XHI > max(contigs.values()):
    parser.error("--xhi-mb exceeds the configured contig length")
specs = []
for spec in args.tracks:
    if "=" not in spec:
        parser.error(f"invalid track specification: {spec}")
    name, path = spec.split("=", 1)
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        parser.error(f"missing or empty track: {path}")
    specs.append((name, path))
tracks = specs


def build_depth(path):
    """True depth track: +1 at start bin, -1 after end bin, then cumsum."""
    cov = {c: np.zeros(L // W + 3, dtype=np.int64) for c, L in contigs.items()}
    with open(path) as fh:
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 5 or p[2] not in cov:
                continue
            s, e = int(p[3]), int(p[4])
            a = cov[p[2]]
            sb, eb = s // W, e // W
            if sb < len(a):
                a[sb] += 1
            if eb + 1 < len(a):
                a[eb + 1] -= 1
    for c in cov:
        cov[c] = np.cumsum(cov[c])
    return cov


data = []
for i, (name, path) in enumerate(tracks):
    data.append((name, build_depth(path), FIXED_COLORS.get(name, PALETTE[i % len(PALETTE)])))
    print(f"{name}: loaded {path}", file=sys.stderr)

fig, axes = plt.subplots(2, 1, figsize=(13, 6.5), sharex=True)
for ax, ctg in zip(axes, ["chr6_MATERNAL", "chr6_PATERNAL"]):
    for name, cov, color in data:
        a = cov[ctg]
        x = (np.arange(len(a)) * W + W / 2) / 1e6
        m = (x >= XLO / 1e6) & (x <= XHI / 1e6)
        ax.plot(x[m], a[m], color=color, lw=0.7, label=name, alpha=0.85)
    ax.axvspan(56, 66, color="grey", alpha=0.08, zorder=0)  # annotated chr6 centromere
    ax.set_ylabel(f"depth ({W}bp win)")
    ax.set_title(ctg, fontsize=10, loc="left")
    ax.set_xlim(XLO / 1e6, XHI / 1e6)
    ax.grid(alpha=0.15)
axes[0].legend(ncol=3, fontsize=9, loc="upper right")
axes[1].set_xlabel("chr6 position (Mb)   [grey band = centromere 56-66Mb; same sample => flat = good]")
fig.suptitle(f"chr6 coverage depth — {XLO/1e6:g}-{XHI/1e6:g}Mb ({W}bp window)", fontsize=11)
fig.tight_layout()
fig.savefig(OUT, dpi=140)
fig.savefig(OUT.replace(".png", ".pdf"))
print(f"wrote {OUT}", file=sys.stderr)
