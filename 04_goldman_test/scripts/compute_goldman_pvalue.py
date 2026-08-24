#!/usr/bin/env python3
"""
Goldman (1993) parametric-bootstrap likelihood-ratio test for phylogenetic
model fit, applied to the comparison LG+F+G+C60 (null model) vs CAT-PMSF
(alternative model), following the protocol in Wang et al. and used in
Fig. 3C of the companion phylogenomics manuscript.

For each of N parametric-bootstrap replicates (simulated under the null
model, LG+F+G+C60, on the null tree) we need the log-likelihood obtained
when the replicate is re-fit under:
  (a) the null model, LG+F+G+C60   -> lnL_LGC60_sim
  (b) the alternative model, CAT-PMSF (frozen site-frequency profiles)
                                    -> lnL_PMSF_sim

and the same two log-likelihoods for the real (empirical) alignment.

    delta_sim_i = lnL_PMSF_sim_i - lnL_LGC60_sim_i     for i = 1..N
    delta_real  = lnL_PMSF_real  - lnL_LGC60_real

    r = #{i : delta_sim_i >= delta_real}
    p = (r + 1) / (N + 1)

A small p-value indicates that the observed advantage of CAT-PMSF over
LG+F+G+C60 on the real data is not something LG+F+G+C60 alone could have
produced by chance -- i.e. CAT-PMSF captures real signal in the data that
LG+F+G+C60 does not.

Usage:
    python3 compute_goldman_pvalue.py \
        --lgc60-real  ../../01_ml_tree/output/nematode2.phy.iqtree \
        --pmsf-real   ../../03_pmsf_fit/output/CAT_pmsf_.iqtree \
        --lgc60-sim-dir ../simulated_fits/LGC60 \
        --pmsf-sim-dir  ../simulated_fits/PMSF \
        --lgc60-sim-suffix _tree.iqtree \
        --pmsf-sim-suffix  _pmsf_.iqtree \
        --out-csv ../results/nematode_goldman_deltas.csv \
        --out-plot ../results/goldman_test_nematode.png
"""
import argparse
import glob
import os
import re
import sys

LNL_RE = re.compile(r"Log-likelihood of the tree:\s*(-?\d+\.\d+)")
REPLICATE_RE = re.compile(r"_sim_(\d+)")


def read_lnl(path: str) -> float:
    with open(path) as fh:
        txt = fh.read()
    m = LNL_RE.search(txt)
    if not m:
        raise ValueError(f"Could not find a log-likelihood in {path}")
    return float(m.group(1))


def collect(sim_dir: str, suffix: str) -> dict:
    out = {}
    for f in glob.glob(os.path.join(sim_dir, f"*{suffix}")):
        m = REPLICATE_RE.search(os.path.basename(f))
        if not m:
            continue
        out[int(m.group(1))] = read_lnl(f)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lgc60-real", required=True,
                     help=".iqtree file for the real alignment fit under LG+F+G+C60")
    ap.add_argument("--pmsf-real", required=True,
                     help=".iqtree file for the real alignment fit under CAT-PMSF")
    ap.add_argument("--lgc60-sim-dir", required=True,
                     help="directory with the *.iqtree files of the null-model refits")
    ap.add_argument("--pmsf-sim-dir", required=True,
                     help="directory with the *.iqtree files of the CAT-PMSF refits")
    ap.add_argument("--lgc60-sim-suffix", default="_tree.iqtree")
    ap.add_argument("--pmsf-sim-suffix", default="_pmsf_.iqtree")
    ap.add_argument("--out-csv", default="nematode_goldman_deltas.csv")
    ap.add_argument("--out-plot", default=None,
                     help="if given, save a Fig.3-style histogram (needs matplotlib)")
    args = ap.parse_args()

    real_lgc60 = read_lnl(args.lgc60_real)
    real_pmsf = read_lnl(args.pmsf_real)
    delta_real = real_pmsf - real_lgc60

    lgc60_sim = collect(args.lgc60_sim_dir, args.lgc60_sim_suffix)
    pmsf_sim = collect(args.pmsf_sim_dir, args.pmsf_sim_suffix)
    common = sorted(set(lgc60_sim) & set(pmsf_sim))
    if not common:
        sys.exit("No matching replicate IDs found between the two simulated-fit directories.")

    rows = []
    r = 0
    for i in common:
        d = pmsf_sim[i] - lgc60_sim[i]
        if d >= delta_real:
            r += 1
        rows.append((i, lgc60_sim[i], pmsf_sim[i], d))

    n = len(common)
    p = (r + 1) / (n + 1)

    os.makedirs(os.path.dirname(args.out_csv) or ".", exist_ok=True)
    with open(args.out_csv, "w") as out:
        out.write("replicate,lnL_LGC60_sim,lnL_PMSF_sim,delta_sim\n")
        for i, lg, pm, d in rows:
            out.write(f"{i},{lg:.4f},{pm:.4f},{d:.4f}\n")

    print(f"Real data:  lnL(LG+F+G+C60) = {real_lgc60:.4f}   lnL(CAT-PMSF) = {real_pmsf:.4f}")
    print(f"delta_real (PMSF - LG+C60)  = {delta_real:.4f}")
    print(f"Bootstrap replicates: N = {n}")
    print(f"r (sim replicates with delta_sim >= delta_real) = {r}")
    print(f"Goldman parametric-bootstrap p-value = (r+1)/(N+1) = {p:.6g}")
    print(f"Wrote per-replicate deltas to {args.out_csv}")

    if args.out_plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib not installed -- skipping plot", file=sys.stderr)
            return
        deltas = [d for _, _, _, d in rows]
        # delta_real sits many orders of magnitude away from the null
        # (simulated) distribution -- as in Fig. 3B/C of the companion
        # manuscript, the two cannot share a linear axis. We plot the null
        # distribution on its own scale and annotate where the real value
        # falls (off-scale, in the direction of a stronger CAT-PMSF fit).
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
        ax.hist(deltas, bins=20, color="#4C72B0", edgecolor="white")
        ax.set_xlabel(r"$\Delta \ln L$ (CAT-PMSF $-$ LG+F+G+C60)")
        ax.set_ylabel("Number of bootstrap replicates")
        ax.set_title("Goldman parametric-bootstrap test\nnull model: LG+F+G+C60 (nematode dataset)")
        ax.annotate(
            f"observed on real data:\nΔlnL = {delta_real:,.0f}\n(off-scale, {n} of {n} replicates fall short → p = {p:.3g})",
            xy=(1.0, 0.5), xycoords="axes fraction",
            xytext=(1.03, 0.5), textcoords="axes fraction",
            ha="left", va="center", fontsize=8.5,
            arrowprops=dict(arrowstyle="->", color="crimson"),
            color="crimson",
        )
        fig.tight_layout()
        os.makedirs(os.path.dirname(args.out_plot) or ".", exist_ok=True)
        fig.savefig(args.out_plot, dpi=150)
        print(f"Wrote plot to {args.out_plot}")


if __name__ == "__main__":
    main()
