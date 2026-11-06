#!/usr/bin/env python
"""Validation of the windowed moves: longer runs than the ctest, several move-set variants.

For each variant and each model of test/python/windowed_ed.py (Kanamori + bath, "diag" and "offdiag"
block structures), N independent runs (seeds spaced by 2, see windowed_ed.compare) are done in
parallel and compared with exact diagonalisation:
  - z of every static observable, chi^2/dof over all observables;
  - flux balance: accepted insertions minus accepted removals of each sub-move (block, or block pair),
    summed over the runs, relative to the Poisson scale sqrt(sum). Detailed balance makes every sub-move's
    expected net flux zero; global balance alone would not (a 2-operator insertion could be undone by
    a 4-operator removal);
  - time per proposed move of each move set, and the acceptance rates.
The variant "W_half" (w = beta/2, d_min = 0) proposes like the stock moves, so its acceptance rates must
equal stock's, and its time per move measures the overhead of the windowed bookkeeping.

    python validation/validate_windowed.py [--seeds 8] [--cycles 100000] [--jobs 8] [--variants ...]

Run from the repository with the fork's triqs_cthyb on PYTHONPATH (e.g. build/python).
"""
import argparse, json, os, re, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "test", "python"))

BETA = 10.0
VARIANTS = {
    "stock":      dict(),
    "stock_nodouble": dict(move_double=False),                                             # control for W_nodouble
    "W":          dict(windowed_moves=True, window_half_width=BETA / 10),                 # defaults: w4 = w, d_min = beta/4
    "W_small":    dict(windowed_moves=True, window_half_width=BETA / 40),
    "W_4paths":   dict(windowed_moves=True, window_half_width=BETA / 10, window_double_half_width=1.5,
                       window_double_min_distance=0.5),                                    # d_min < 2 w4: all 4 labelled paths
    "W_nodouble": dict(windowed_moves=True, window_half_width=BETA / 10, move_double=False),
    "W_noshift":  dict(windowed_moves=True, window_half_width=BETA / 10, move_shift=False),
    "W_half":     dict(windowed_moves=True, window_half_width=BETA / 2, window_double_min_distance=0.0),
}

STAT_RE = re.compile(r"Move (set )?(.+?): Proposed = (\d+), Accepted = (\d+), Rate")
DUR_RE  = re.compile(r"Move (set )?(.+?): Duration = ([0-9.eE+-]+)")


def worker(kind, variant, seed, n_cycles):
    """One run; prints a JSON line with the observables and the move statistics"""
    import windowed_ed as m
    t0 = time.time()
    S = m.solve(kind, seed=seed, n_cycles=n_cycles, verbosity=3, **VARIANTS[variant])
    out = dict(obs=m.qmc_averages(S, kind), wall=time.time() - t0, order=S.average_order, sign=S.average_sign)
    print("RESULT " + json.dumps(out), flush=True)


def run_one(kind, variant, seed, n_cycles):
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")  # no spinning BLAS threads
    p = subprocess.run([sys.executable, "-W", "ignore", __file__, "--worker", kind, variant, str(seed), str(n_cycles)],
                       capture_output=True, text=True, env=env)
    res = None
    stats, durs = {}, {}
    for line in p.stdout.splitlines():
        line = re.sub(r"^\[Rank \d+\]\s*", "", line)
        if line.startswith("RESULT "): res = json.loads(line[7:])
        if (mt := STAT_RE.search(line)): stats[mt.group(2)] = (int(mt.group(3)), int(mt.group(4)))
        if (mt := DUR_RE.search(line)): durs[mt.group(2)] = float(mt.group(3))
    if res is None: raise RuntimeError(f"run {kind} {variant} {seed} failed:\n{p.stdout[-2000:]}\n{p.stderr[-2000:]}")
    res["stats"], res["durations"] = stats, durs
    return res


def analyse(kind, variant, runs, exact):
    lines = []
    keys = list(exact)
    zs = []
    for k in keys:
        v = np.array([r["obs"][k] for r in runs])
        mean, err = v.mean(), v.std(ddof=1) / np.sqrt(len(v))
        z = (mean - exact[k]) / err
        zs.append(z)
        lines.append(f"    {k:10s} {mean: .5f} +- {err:.5f}  ED {exact[k]: .5f}  z {z: .2f}")
    chi2 = float(np.mean(np.square(zs)))
    # flux balance per sub-move: insertion "Insert window Delta_X" vs removal "Remove window Delta_X"
    flux = []
    for name in runs[0]["stats"]:
        if not name.startswith("Insert") or name.endswith("operators"): continue
        rname = "Remove" + name[len("Insert"):]
        a_ins = sum(r["stats"][name][1] for r in runs)
        a_rem = sum(r["stats"][rname][1] for r in runs)
        flux.append((name.split("Delta_")[-1], a_ins, a_rem, (a_ins - a_rem) / np.sqrt(max(a_ins + a_rem, 1))))
    # acceptance and time per proposed move of each move set
    sets = {}
    for name in runs[0]["stats"]:
        if name.endswith("operators") or name == "Shift one operator":
            prop = sum(r["stats"][name][0] for r in runs)
            acc = sum(r["stats"][name][1] for r in runs)
            dur = sum(r["durations"].get(name, 0) for r in runs)
            sets[name] = (acc / prop, 1e6 * dur / prop)
    wall = np.mean([r["wall"] for r in runs])
    return chi2, max(abs(z) for z in zs), lines, flux, sets, wall, np.mean([r["order"] for r in runs]), np.mean([r["sign"] for r in runs])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--cycles", type=int, default=100000)
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--variants", nargs="*", default=list(VARIANTS))
    ap.add_argument("--kinds", nargs="*", default=["diag", "offdiag"])
    ap.add_argument("--seed0", type=int, default=1001, help="seeds are seed0 + 2k (RandMT: seeds 2k, 2k+1 coincide)")
    ap.add_argument("--worker", nargs=4)
    a = ap.parse_args()
    if a.worker:
        kind, variant, seed, n_cycles = a.worker
        return worker(kind, variant, int(seed), int(n_cycles))

    import windowed_ed as m
    assert m.beta == BETA
    jobs = [(k, v, a.seed0 + 2 * s) for k in a.kinds for s in range(a.seeds) for v in a.variants]  # variants interleaved
    with ThreadPoolExecutor(a.jobs) as ex:
        futs = {j: ex.submit(run_one, *j, a.cycles) for j in jobs}
        results = {j: f.result() for j, f in futs.items()}

    summary = {}
    for kind in a.kinds:
        exact = m.ed_averages(kind)
        for v in a.variants:
            runs = [results[(kind, v, a.seed0 + 2 * s)] for s in range(a.seeds)]
            chi2, zmax, lines, flux, sets, wall, order, sign = analyse(kind, v, runs, exact)
            summary[(kind, v)] = sets
            print(f"\n== {kind} / {v}: {VARIANTS[v]}")
            print(f"  {a.seeds} runs x {a.cycles} cycles, wall {wall:.1f} s/run, <order> {order:.2f}, <sign> {sign:.3f}")
            print(f"  chi2/dof vs ED = {chi2:.2f} (dof {len(exact)}), max |z| = {zmax:.2f}")
            print("\n".join(lines))
            print("  flux balance (accepted insert / remove, z):  " +
                  ",  ".join(f"{n} {i}/{r} {z:+.1f}" for n, i, r, z in flux))
            print("  move sets (acceptance, us per proposal):  " + ",  ".join(f"{n}: {acc:.4f}, {t:.2f}" for n, (acc, t) in sets.items()))
    if "stock" in a.variants and "W_half" in a.variants:
        print("\n== W_half (w = beta/2, d_min = 0) vs stock: acceptance ratio, time-per-proposal ratio")
        for kind in a.kinds:
            st, wh = summary[(kind, "stock")], summary[(kind, "W_half")]
            for n in st:
                if n in wh: print(f"  {kind:8s} {n:24s} acc {wh[n][0] / st[n][0]:.3f}   time {wh[n][1] / st[n][1]:.3f}")


if __name__ == "__main__":
    main()
