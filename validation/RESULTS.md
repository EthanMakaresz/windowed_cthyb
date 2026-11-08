# Validation results

All on the code of commit `2f5702b` (on `main`). Raw outputs: `validate_windowed.out.txt`
and `real_problem_mu7.out.txt`.

## Exact diagonalisation (`validate_windowed.py --seeds 32 --seed0 9001`)

Model of `test/python/windowed_ed.py`: a two-orbital Kanamori impurity (U = 2, J = 0.3, with
spin-flip and pair hopping) plus one bath site per spin-orbital, at β = 10, ⟨k⟩ ≈ 10.

- *diag*: four 1×1 blocks.
- *offdiag*: two 2×2 blocks with an orbital-off-diagonal hybridisation. Its average sign is ≈ 0.63.

Each move set gets 32 independent runs of 10⁵ cycles. The nine static observables are compared
with ED: occupations, double occupancies, the spin-flip and pair-hopping correlators, and ⟨H_int⟩.
The 9 observables are strongly correlated, so χ²/dof fluctuates more than for independent ones.

"Flux" is the per-sub-move z of (accepted insertions − accepted removals) / √(sum). It is only
informative when 4-operator moves are present; with 2-operator moves alone, it is fixed by the
change in perturbation order.

| move set | diag χ²/dof | diag max \|z\| | offdiag χ²/dof | offdiag max \|z\| | flux max \|z\| (diag / offdiag) |
|---|---|---|---|---|---|
| stock | 1.22 | 1.78 | 0.68 | 1.84 | 2.3 / 1.5 |
| W: w = β/10 (defaults w4 = w, d_min = β/4) | 0.66 | 1.45 | 0.54 | 1.42 | 2.1 / 1.5 |
| W_small: w = β/40 | 0.91 | 1.73 | 1.37 | 2.37 | 1.6 / 1.6 |
| W_4paths: w = β/10, w4 = 1.5, d_min = 0.5 (all 4 labelled paths) | 0.41 | 1.42 | 0.91 | 1.68 | 2.4 / 1.1 |
| W_noshift: W without shift moves | 1.62 | 2.78 | 1.54 | 2.46 | 2.5 / 2.4 |
| W_half: w = β/2, d_min = 0 (stock-like proposals) | 1.62 | 2.54 | 0.42 | 1.63 | 2.0 / 1.1 |
| stock, `move_double=False` | 0.86 | 1.80 | **1030** | **84** | – |
| W, `move_double=False` | 0.37 | 1.06 | **634** | **61** | – |

All windowed variants with 4-operator moves agree with ED. That includes small windows, all four
labelled paths, and no shift moves.

**Without 4-operator moves the off-diagonal model is not ergodic, for stock cthyb and windowed mode
alike.** Both converge to the same wrong distribution:

| | ED | stock, no 4-op | windowed, no 4-op |
|---|---|---|---|
| spin-flip ⟨c†↑₀c↓₀c†↓₁c↑₁⟩ | 0.05042 | 0.04260(18) | 0.04257(19) |
| pair-hop ⟨c†↑₀c†↓₀c↓₁c↑₁⟩ | −0.02058 | −0.01455(7) | −0.01460(10) |
| average sign | 0.629 (with 4-op) | 0.853 | 0.851 |

This is cthyb's known need for 4-operator moves. Here an orbital-off-diagonal pair has zero weight
unless a pair in the other spin block compensates it. The windowed cross-block 4-operator moves
reach those configurations, as the stock ones do.

### W_half against stock: same proposals, same acceptance

| | diag | offdiag |
|---|---|---|
| 2-operator acceptance, W_half / stock | 1.001 (insert), 1.001 (remove) | 0.999, 1.000 |
| 4-operator acceptance, W_half / stock | 1.093, 1.092 | 1.073, 1.073 |

The 2-operator moves match. The 4-operator acceptance is higher by design: for same-block removals,
stock may pick the same c† or c twice and then rejects. Its proposal ratio is therefore
(sβ/(n+2))⁴, while the windowed one is s⁴β⁴/((n+2)²(n+1)²). Both are correct, and the windowed one
accepts more often.

## Real problem (`real_problem.py`, single-orbital Hubbard impurity, U = 15, β = 40, μ = 7)

The Δ(τ) and h_loc0 are those of a first DMFT iteration (`diag_mu7_delta.h5`, not included).

24 seeds per move set, 2·10⁷ moves each. Observables are compared with stock;
`auto_corr_time` is cthyb's (sign and perturbation order, in cycles of 100 moves).

| move set | ⟨k⟩ | n | double occ. | acceptance (2-op / 4-op) | autocorr. (cycles) |
|---|---|---|---|---|---|
| stock | 11.340(17) | 0.99971(47) | 0.00962(33) | 0.0094 / 0.0002 | 82.8 |
| W: w = 1.5 | 11.326(6), z −0.8 | 0.99972(15), z +0.0 | 0.00959(11), z −0.1 | 0.102 / 0.017 | 6.9 |
| W_half: w = β/2, d_min = 0 | 11.321(23), z −0.7 | 1.00032(56), z +0.8 | 0.00960(45), z −0.0 | 0.0094 / 0.0002 | 83.1 |

Cost, from `perf stat` (user-space instructions and cycles per move, one solve of 10⁷ moves each, one BLAS thread):

| move set | instructions / move | cycles / move | cycles × autocorr., relative to stock |
|---|---|---|---|
| stock | 6582 | 3225 | 1 |
| W_half | 7027 (+6.8%) | 3563 (+10.5%) | 1.11 |
| W, w = 1.5 | 8247 (+25%) | 4327 (+34%) | **0.11** (8.9× cheaper per independent sample) |
| stock, no 4-op | 6103 | 2769 | – |
| W_half, no 4-op | 6173 (+1.1%) | 2828 (+2.2%) | – |

The 8.9× is per independent sample of the sign and the perturbation order, the quantities that
`auto_corr_time` covers. G(τ) was not measured here. On the Sr2RuO4 problem
(`speedup_sr2ruo4/`), the gain on G(β/2) was two to four times smaller than on the sign and the
order at β ≥ 70.

With identical proposals, the windowed 2-operator moves cost 1–2% more per move. The rest of
W_half's overhead is in the 4-operator moves, explained in the README's Overhead section. Wall-clock
times per move set from mc_generic's timers vary by ±20% between identical runs on this machine,
so the instruction and cycle counts are the numbers to use.
