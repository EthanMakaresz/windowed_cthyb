#!/usr/bin/env python
"""Stock vs windowed mode on a single-orbital Hubbard impurity given by Delta(tau) and h_loc0
(Delta interface, as solve_generic uses it). The HDF5 file holds the keys Delta_tau, h_loc0, U and beta,
e.g. the Delta(tau) of a first DMFT iteration.

    python validation/real_problem.py DELTA.h5 [--seeds 4] [--moves 2e7] [--jobs 11] [--w 1.5]

Variants: stock; windowed with w (and the default w4 = w, d_min = beta/4); windowed with
w = beta/2, d_min = 0, which proposes like the stock moves. Reports, over the seeds, <k>, n and
double occupancy with z against stock, the move acceptance, the time per proposed move and
CTHYB's autocorrelation time. The last two give the cost per independent sample.
"""
import argparse, json, os, re, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
import numpy as np

STAT_RE = re.compile(r"Move set (.+?): Proposed = (\d+), Accepted = (\d+), Rate")
DUR_RE  = re.compile(r"Move set (.+?): Duration = ([0-9.eE+-]+)")


def worker(delta, variant, seed, moves, w):
    import triqs.operators  # noqa: F401
    from h5 import HDFArchive
    from triqs.operators import n as n_op
    from triqs_cthyb import Solver
    with HDFArchive(delta, "r") as ar:
        Delta_tau, h_loc0, U, beta = ar["Delta_tau"], ar["h_loc0"], float(ar["U"]), float(ar["beta"])
    gf_struct = [(name, g.target_shape[0]) for name, g in Delta_tau]
    S = Solver(beta=beta, gf_struct=gf_struct, n_iw=5000, n_tau=10001, delta_interface=True)
    S.Delta_tau << Delta_tau
    lc = 100
    p = dict(h_int=U * n_op("up", 0) * n_op("down", 0), h_loc0=h_loc0, n_cycles=int(moves / lc), length_cycle=lc,
             n_warmup_cycles=int(0.1 * moves / lc), random_seed=34788 + 928374 * seed, measure_density_matrix=True,
             use_norm_as_weight=True, verbosity=3, measure_G_tau=False, perform_post_proc=False)
    if variant.endswith("_nd"):  # variant without the four-operator moves
        variant = variant[:-3]
        p.update(move_double=False)
    if variant == "W": p.update(windowed_moves=True, window_half_width=w)
    if variant == "W_half": p.update(windowed_moves=True, window_half_width=beta / 2, window_double_min_distance=0.0)
    t0 = time.time()
    S.solve(**p)
    elapsed = time.time() - t0
    fs = S.h_loc_diagonalization.fock_states
    probs = [float(np.real(S.density_matrix[i][0, 0])) for i in range(len(fs))]
    pc = [bin(int(f[0])).count("1") for f in fs]
    out = dict(k=float(S.average_order), occ=sum(c * q for c, q in zip(pc, probs)), docc=sum((c == 2) * q for c, q in zip(pc, probs)),
               tau=float(S.auto_corr_time), elapsed=elapsed, moves=moves * 1.1)
    print("RESULT " + json.dumps(out), flush=True)


def run_one(delta, variant, seed, moves, w):
    # one BLAS/OpenMP thread: otherwise idle OpenBLAS threads spin and compete with the solvers for the cores
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    p = subprocess.run([sys.executable, "-W", "ignore", __file__, delta, "--worker", variant, str(seed), "--moves", str(moves), "--w", str(w)],
                       capture_output=True, text=True, env=env)
    res, stats, durs = None, {}, {}
    for line in p.stdout.splitlines():
        line = re.sub(r"^\[Rank \d+\]\s*", "", line)
        if line.startswith("RESULT "): res = json.loads(line[7:])
        if (m := STAT_RE.search(line)): stats[m.group(1)] = (int(m.group(2)), int(m.group(3)))
        if (m := DUR_RE.search(line)): durs[m.group(1)] = float(m.group(2))
    if res is None: raise RuntimeError(p.stdout[-2000:] + p.stderr[-2000:])
    res["stats"], res["durations"] = stats, durs
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("delta")
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--moves", type=float, default=2e7)
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--w", type=float, default=1.5)
    ap.add_argument("--seed0", type=int, default=1, help="seeds seed0 .. seed0 + seeds - 1 (random_seed = 34788 + 928374 * seed)")
    ap.add_argument("--worker", nargs=2)
    a = ap.parse_args()
    if a.worker: return worker(a.delta, a.worker[0], int(a.worker[1]), a.moves, a.w)

    variants = ["stock", "W", "W_half"]
    seeds = range(a.seed0, a.seed0 + a.seeds)
    jobs = [(v, s) for s in seeds for v in variants] # interleaved, so that all variants see the same machine load
    with ThreadPoolExecutor(a.jobs) as ex:
        res = {j: f.result() for j, f in {j: ex.submit(run_one, a.delta, *j, a.moves, a.w) for j in jobs}.items()}
    ref, t_set = {}, {}
    print(f"{a.delta}: {a.seeds} seeds x {a.moves:.1e} moves (+10% warmup); W: w = {a.w}; W_half: w = beta/2, d_min = 0")
    for v in variants:
        runs = [res[(v, s)] for s in seeds]
        line = f"  {v:7s}"
        for key in ("k", "occ", "docc"):
            x = np.array([r[key] for r in runs])
            mean, err = x.mean(), x.std(ddof=1) / np.sqrt(len(x))
            if v == "stock": ref[key] = (mean, err)
            z = (mean - ref[key][0]) / np.hypot(err, ref[key][1]) if v != "stock" else 0.0
            line += f"  {key} {mean:.5f}+-{err:.5f} (z {z:+.1f})"
        us = 1e6 * np.mean([r["elapsed"] / r["moves"] for r in runs])
        tau = np.mean([r["tau"] for r in runs])
        acc = {n: sum(r["stats"][n][1] for r in runs) / sum(r["stats"][n][0] for r in runs) for n in runs[0]["stats"]}
        line += f"\n           {us:.3f} us/move, autocorr {tau:.1f} cycles, cost/indep. sample {us * tau:.1f} (arb.), acceptance " + \
                ", ".join(f"{n.replace(' operators', '')}: {x:.4f}" for n, x in acc.items())
        print(line)
        t_set[v] = {n: 1e6 * sum(r["durations"][n] for r in runs) / sum(r["stats"][n][0] for r in runs) for n in runs[0]["stats"]}
    print("  time per proposal of each move set (us; mc_generic's move timers, measurements excluded), ratio to stock:")
    for n in t_set["stock"]:
        print(f"    {n:24s} " + "  ".join(f"{v} {t_set[v][n]:.3f} ({t_set[v][n] / t_set['stock'][n]:.3f})" for v in variants))


if __name__ == "__main__":
    main()
