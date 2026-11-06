# triqs_cthyb with windowed insertion moves

This is a fork of [TRIQS/cthyb](https://github.com/TRIQS/cthyb), branch `4.0.x` (base commit
`a6af361`, version 4.0.1). It adds an optional **windowed mode**, in which the insert and remove moves
propose the c and c† of each pair close together in imaginary time. With the default parameters the
solver is stock cthyb 4.0.1: same move set, same random-number stream, bit-identical results.

Two commits on branch `windowed-moves`, on top of upstream `4.0.x`:

1. **Σ₁ tail fix**, from [TRIQS/cthyb#197](https://github.com/TRIQS/cthyb/pull/197) (not merged
   upstream when this was written). `sigma_high_frequency_moments` subtracted the element-wise
   square of Σ₀ from Σ₁; the correct term is the matrix product Σ₀Σ₀. This matters for blocks where
   Σ₀ has off-diagonal elements, used as known moments in `solve_generic` tail fits. Its ED test is
   `test/python/sigma_moments_offdiag.py`. Revert this commit once a 4.0.x release contains the fix.
2. **Windowed moves**, described below.

**How much faster than stock?** On a realistic three-orbital Kanamori problem (the t2g bands of
Sr2RuO4, with spin-flip and pair hopping), windowed mode cuts the CPU time per independent sample by
1.3× at β = 10 and 17× at β = 200, the gain growing as ≈ β^0.9, with the cost per move included.
See [`validation/speedup_sr2ruo4/`](validation/speedup_sr2ruo4/README.md) for the plot, the model,
the method and the caveats; measured on G(β/2) directly, the gain is smaller (≈ 5–12× at β = 200).

---

## Windowed mode

### Why

At low temperature, the stock moves insert a c†–c pair at two independent uniform times on
[0, β). Most such proposals are far apart and rejected, since the hybridisation and the local
trace favour short pairs. Proposing the c within a window around its c† gives higher acceptance per
proposal. The proposal ratio corrects exactly for the bias, so results are unchanged.

### Usage

```python
S.solve(h_int=h_int, n_cycles=..., ...,
        windowed_moves=True,          # replace the stock insert/remove moves
        window_half_width=2.0)        # w: |tau_c - tau_cdag| < w within each pair
```

Windowed mode is all or nothing. With `windowed_moves=True`, only the windowed insert/remove moves
are in the move set; the stock ones are not added:

| stock move set | windowed mode |
|---|---|
| "Insert two operators": `Insert Delta_<b>` | `Insert window Delta_<b>` |
| "Remove two operators": `Remove Delta_<b>` | `Remove window Delta_<b>` |
| "Insert four operators": `Insert Delta_<b1>_<b2>` (if `move_double`) | `Insert window Delta_<b1>_<b2>` (if `move_double`) |
| "Remove four operators": `Remove Delta_<b1>_<b2>` (if `move_double`) | `Remove window Delta_<b1>_<b2>` (if `move_double`) |
| shift (`move_shift`), global moves (`move_global`) | unchanged, controlled by their own parameters |

Each move set keeps its weight. `proposal_prob` keeps weighting the blocks, as in stock cthyb.

### Parameters

| parameter | default | meaning |
|---|---|---|
| `windowed_moves` | `False` | Turn windowed mode on. |
| `window_half_width` | `None` → β/20 | w, the half-width of the window of the two-operator moves: the c is proposed uniformly in (τ_c† − w, τ_c† + w). Must be in (0, β/2]. If unset, β/20 is used and a **warning is printed**. |
| `window_double_half_width` | `None` → `window_half_width` | w4, the same half-width, within each of the two pairs of the four-operator moves. Must be in (0, β/2]. |
| `window_double_min_distance` | `None` → β/4 | d_min, which only applies to four-operator moves within **one block**: the two c† are proposed more than d_min apart. Must be in [0, β/2). |

All distances are cyclic on the imaginary-time circle. If any `window_*` parameter is set while
`windowed_moves=False`, a warning says it is ignored. With `verbosity >= 2`, the values used are
printed. The parameters are stored in `last_solve_parameters` and round-trip through HDF5.

Why d_min applies within one block only: with two pairs in the same block, the det is updated by
one rank-2 update (`try_insert2`/`try_remove2`). If the two c† are close together, the two new rows
are nearly linearly dependent and the fast update is poorly conditioned. Pairs in different blocks
are two independent rank-1 updates of different matrices, which have no such problem. Their second
c† is drawn uniformly on [0, β), so a single move can still insert correlated cross-block
operators, such as the spin-flip and pair-hopping vertex pairs of a Kanamori interaction.

### Choosing w

w should cover the c†–c distances that are actually accepted. A window that is too small proposes
only short pairs. It stays correct (see Ergodicity below), but long pairs must then be built up from
several moves, and the autocorrelation time grows. A window that is too large wastes proposals, as
the stock moves do. With w = β/2 the two-operator moves propose exactly as the stock ones do.

To estimate w, run stock cthyb once with `performance_analysis=True` and look at the histograms
`insert_length_accepted_<block>` (τ_c − τ_c† mod β, so short pairs sit near both 0 and β).
Then compare the autocorrelation time per CPU second (`S.auto_corr_time`, wall time) for a few
values of w. Windowed mode has no analysis tool of its own, and its moves fill no histograms.

---

## The moves

Notation: s = size of the block, n = number of c† in the block, all counts taken in the larger of the
two configurations (C' = C plus the inserted operators). W(C) is the configuration weight
(det × trace). "near X" means within the window of X; "far from X" means farther than d_min from X.

**Two operators** (`moves/insert_window.cpp`, `moves/remove_window.cpp`)
- Insert: c† at X, uniform on [0, β); c at x, uniform in the window of X; inner indices uniform.
- Remove: X uniformly among the n c† of the block; then x uniformly among the m c of the block near
  X (m ≥ 1 is required).

**Four operators, different blocks b1 ≠ b2** (`moves/double_insert_window.cpp`, `double_remove_window.cpp`)
- Insert: pair (X, x) in b1 as above; pair (Y, y) in b2 with Y uniform on [0, β) and y near Y.
- Remove: X among the c† of b1, x among the c of b1 near X, Y among the c† of b2, y among the c of
  b2 near Y.

**Four operators, same block**
- Insert: X uniform; x near X; Y uniform among the points farther than d_min from X; y near Y.
- Remove: X among the n c†; x among the c near X; Y among the c† far from X (X itself is never
  far); y among the c near Y, excluding x.

## Detailed balance

Each insertion sub-move (block b, or block pair (b1, b2)) is paired with the removal of the same b or
(b1, b2). Both are proposed with the same probability: the insert and remove move sets have equal
weights, and the same `proposal_prob` weights. Balance therefore holds separately for each pair:

    P(C→C') / P(C'→C) = [Σ_P f(P)] / [Σ_P r(P)],     R_insert = W(C')/W(C) · Σ_P r(P) / Σ_P f(P),

with R_remove = 1/R_insert evaluated for the same pair (C, C'). The sum runs over the labelled
paths P, the ways the move's choices can name the same set of new operators. f(P) is the proposal
density of the insertion along P, and r(P) the proposal probability of the removal along P.
Neither the det nor the trace depends on how the operators are paired, so only the set of
operators matters.

- **Two operators.** There is one path. f = 1/(β · L_w · s²) and r = 1/(n · m), so

      R_insert = W(C')/W(C) · s² β L_w / (n m),

  where m is the number of c near X in C'. L_w ≈ 2w is the measure of the window (below). At
  w = β/2 this is the stock move, R = W'/W · (sβ/n)².
- **Four operators, different blocks.** There is one path (the blocks tell the pairs apart), and the move
  factorises into two windowed pairs:

      R_insert = W(C')/W(C) · s1² s2² β² L_w4² / (n1 m_x · n2 m_y).

  (b1, b2) and (b2, b1) are distinct sub-moves, each balanced with its own reverse.
- **Four operators, same block.** Four labelled paths name the same four new operators: which new
  c† is X, times which new c is paired with it. With I(P) = [x near X]·[Y far X]·[y near Y]:

      f(P) = I(P) / (β · L_w4 · L_far · L_w4 · s⁴),
      r(P) = I(P) / (n · m_x(X) · m_Y(X) · (m_x(Y) − [x near Y])),

  where m_x(X) is the number of c near X, m_Y(X) the number of c† far from X, and the bracket
  removes x from y's candidates. [Y far X] is symmetric in X and Y, so it is common to all four
  paths. If d_min ≥ 2·w4, only the two X↔Y relabellings can be valid. The implementation is
  `window::same_block_proposal_ratio` in `moves/window_util.hpp`.

Proposals that collide with an existing operator time are rejected (as in stock). The reverse
move never proposes such a configuration, so this does not bias the chain.

**Exactness on the time grid.** cthyb stores times on a uint64 grid of Nmax points, where time
arithmetic is cyclic. Windows are defined on this grid: a window of half-width w has
Wn = ⌊Nmax·w/β⌋ grid steps, and "near" means a cyclic grid distance below Wn, which is exactly
2Wn − 1 points. "Far" means Dn < e < Nmax − Dn, which is exactly Nmax − 2Dn − 1 points. The
insertions draw uniformly from exactly these point sets. The removals count with exactly the same
integer predicates, and the ratios use L = (number of points) · β/Nmax. There is no rounding at the
window edges, and no proposal ever lands outside its window.

## Ergodicity

The windowed moves restrict proposals, not configurations: every configuration remains reachable.

- **Transport.** In a configuration containing c†(X), insert a pair (c†(X'), c(s)) with s near X
  and X' near s, then remove the pair (c†(X), c(s)). The net effect is that c†(X) has moved by
  less than 2w, and its inner index can have changed. The same works for a c. Repeating this, any
  operator can be carried anywhere, so any configuration can be reduced to the empty one and built
  from it.
- **Nonzero intermediate weights.** The extra short pair can always be chosen with nonzero trace. To
  first order, a same-flavour pair at nearly equal times inserts n_a or 1 − n_a (depending on which
  operator comes first). Since Tr[…n_a…] + Tr[…(1 − n_a)…] = Tr[…] ≠ 0, at least one of the two
  orderings has nonzero weight.
- **Four-operator moves.** Cross-block windowed moves still propose a second pair anywhere. So, as
  in stock, one move can place correlated operators of two blocks close together (spin-flip or
  pair-hopping vertices) when the 2-operator intermediate has zero weight. The vertices must then
  be within w4 of each other. Same-block moves never put the two c† within d_min; those
  configurations are reached with two-operator moves.
- **Keep `move_double=True`** (the default) when the hybridisation is off-diagonal within a block and
  H_loc has spin-flip or pair-hopping terms. In the off-diagonal ED test, stock cthyb *and* windowed
  mode without 4-operator moves are both badly non-ergodic, and both give the same wrong answer (the
  spin-flip and pair-hopping correlators are off by 15–30%). With 4-operator moves, both are
  correct.
- **Limits.** As with stock cthyb, ergodicity for an arbitrary local Hamiltonian is not provable,
  because intermediate configurations can have exactly zero trace. Ergodicity has therefore also
  been checked numerically: see the ED tests and validation below, on models with spin-flip and
  pair-hopping terms, diagonal and multi-orbital blocks, and shift moves on and off.

## Overhead

The windowed moves follow the stock moves line for line: tree insert, det update, quick-abandon
trace with `rng.preview()`, accept and reject. The window bookkeeping is designed to cost little:

- The dets keep their c† and c sorted by time. The near points of t0 form one cyclic arc of the
  grid, [t0 − (Wn − 1), t0 + (Wn − 1)], and the far points another, [t0 + Dn + 1, t0 − Dn − 1]. Each
  arc's two end points are computed once per query. The operators inside the arc are then a single
  cyclically contiguous index range of the det. It is found by two binary searches that only
  compare time points, O(log n), with no allocation and no time arithmetic per probe. Counting the
  range or drawing its k-th element is constant time after that. Det insertion positions also use
  binary search; stock uses a linear scan.
- Window sizes, lengths and arc offsets are computed once, when the moves are constructed.
- No histograms are filled. The same-block 4-path sum reuses the counts and arcs that the move
  already has.

**Measured** on a single-orbital Hubbard impurity, U = 15, β = 40, μ = 7, ⟨k⟩ ≈ 11 (the Δ(τ) of a
first DMFT iteration).
Instructions and cycles per move are from `perf stat`, with one BLAS thread:

| move set | instructions / move | cycles / move | vs stock |
|---|---|---|---|
| stock, no 4-op moves | 6103 | 2769 | — |
| windowed, w = β/2 (stock proposals), no 4-op moves | 6173 | 2828 | +1.1% instr., +2.2% cycles |
| stock | 6582 | 3225 | — |
| windowed, w = β/2, d_min = 0 (stock proposals) | 7027 | 3563 | +6.8% instr., +10.5% cycles |
| windowed, w = 1.5 | 8247 | 4327 | +25% instr., +34% cycles; 10× the acceptance |

The windowed 2-operator moves cost about 1–2% over stock for the same proposals. The 4-operator
moves add more, for two reasons. First, their counts (four arcs, same block). Second, stock's
same-block 4-op removal returns at once, with no work, when it picks the same c† or c twice
(about a third of the time at n ≈ 6), while the windowed one never wastes a proposal that way.
With a real window (w = 1.5), a move costs 34% more cycles, because 10× more proposals pass the
quick-abandon test and get full det and trace updates. The cost per independent sample is still
**8× lower than stock** (see `validation/RESULTS.md`).

---

## Tests and validation

ctest, run with `ctest --test-dir build`:
- All of upstream's tests pass with their stock reference files, which shows that the default path
  is unchanged. The same holds for the Σ₁ test of #197.
- `window_geometry` (C++): the near/far counts, k-th selection and insertion positions are compared
  with brute force on random sorted grids. These include wrap-around, t0 in and out of the array,
  operators exactly at the window edges, w = β/2 and d_min = 0. The draws are checked to land exactly
  on the window sets, including their first and last points. `same_block_proposal_ratio` is compared
  with Σr/Σf from brute-force enumeration of every choice sequence of the removal, and of the
  labelled insertion paths. The cases include d_min < 2w4, so that all four paths occur. Deliberate
  bugs each make it fail: no x exclusion, a wrong path normalisation, off-by-ones at the window and
  far-set edges, in the wrap-around count, in the rank map and in the draw offset.
- `Py_windowed_ed`: a two-orbital Kanamori impurity (U = 2, J = 0.3, spin-flip and pair-hopping)
  with one bath site per spin-orbital, β = 10, w = β/10. It has 8 modes, so it is solved exactly
  with `AtomDiag`. Two block structures are used: diagonal (4 blocks of size 1), and per spin with an
  orbital-off-diagonal hybridisation (2 blocks of size 2). The test compares, over 8 runs:
  occupations, double occupancies, the ⟨c†↑₀c↓₀c†↓₁c↑₁⟩ and ⟨c†↑₀c†↓₀c↓₁c↑₁⟩ correlators, and ⟨H_int⟩.
  These come from the impurity density matrix, which does not use M⁻¹.

`validation/validate_windowed.py` is a longer run with several move-set variants. It checks χ² against
ED, the per-sub-move flux balance and the acceptance, and measures the time per proposal.
`validation/real_problem.py` runs stock and windowed mode side by side on a Δ(τ) file, here a
single-orbital Hubbard impurity with U = 15, β = 40. The results of both are in
[`validation/RESULTS.md`](validation/RESULTS.md). In short: every windowed variant agrees with ED,
and on the U = 15 problem, w = 1.5 is 8.9× cheaper per independent sample than stock.
[`validation/speedup_sr2ruo4/`](validation/speedup_sr2ruo4/README.md) measures the speedup over
stock against temperature, β = 10 … 200, on a three-orbital Sr2RuO4 impurity problem: 1.3× at
β = 10, 17× at β = 200.

Notes for running tests:
- The G(τ) estimator built from M⁻¹ has a heavy-tailed distribution in these models, for stock and
  windowed alike. Compare G(τ) between runs only with many independent seeds and robust error bars.
- Run each solver process with `OPENBLAS_NUM_THREADS=1` (and `OMP_NUM_THREADS=1`) unless BLAS
  threads are wanted. Otherwise, idle OpenBLAS threads of numpy and scipy spin and take CPU from the
  solvers when several run on one node: they were about 23% of all CPU samples in a profile of one
  solve here. This applies to stock cthyb as well.
- The default RNG (`random_name=""`, i.e. RandMT) forces the seed odd, so seeds 2k and 2k + 1 give
  identical streams. Space seeds by 2, or use `random_name="mt19937"`. cthyb's usual
  `34788 + 928374*rank` is unaffected.

## Build and install

**Requirements:** the same as upstream cthyb 4.0. You need an installed TRIQS 4.0.x, CMake ≥ 3.22,
the C++20 compiler and Python that TRIQS was built with, and network access at configure time.
The configure step downloads two build dependencies, c2py and GoogleTest (see Offline builds below).

Throughout, `$TRIQS_PREFIX` is the TRIQS installation prefix and `$PREFIX` is where this fork goes.
To keep a stock cthyb next to it, use a separate `$PREFIX`. If you leave out
`CMAKE_INSTALL_PREFIX`, the fork is installed into the TRIQS prefix and replaces any cthyb there.

```bash
source $TRIQS_PREFIX/share/triqs/triqsvars.sh          # sets TRIQS_ROOT, CMAKE_PREFIX_PATH, PYTHONPATH

git clone -b windowed-moves <url of this repository> cthyb.src
cmake -S cthyb.src -B cthyb.build \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER=<the compiler TRIQS was built with> \
      -DCMAKE_INSTALL_PREFIX=$PREFIX
cmake --build cthyb.build -j <N>                        # about 2 GB of memory per compile job
ctest --test-dir cthyb.build -j <N> --output-on-failure
cmake --install cthyb.build
```

To find the compiler TRIQS was built with, run
`grep TRIQS_CXX_COMPILER $TRIQS_PREFIX/lib/cmake/triqs/TRIQSConfig.cmake`.

**Using it.** With a separate `$PREFIX`, source the environment script that the install writes. Do
this *after* `triqsvars.sh`, so that the fork comes first on `PYTHONPATH`:

```bash
source $TRIQS_PREFIX/share/triqs/triqsvars.sh
source $PREFIX/share/triqs_cthyb/triqs_cthybvars.sh
python -c "import triqs_cthyb; print(triqs_cthyb.__file__)"   # check which cthyb is imported
```

In batch jobs, put the same two lines in the job script. MPI launchers normally pass the
environment on to all ranks; if a rank still imports another cthyb, export `PYTHONPATH` explicitly
(e.g. `mpirun -x PYTHONPATH` with Open MPI). Removing the second line switches back to the stock
cthyb.

**CMake options.** These are upstream's, unchanged:

| option | default | meaning |
|---|---|---|
| `-DHybridisation_is_complex=ON` | OFF | allow a complex Δ(τ) |
| `-DLocal_hamiltonian_is_complex=ON` | OFF | allow a complex H_loc |
| `-DMeasureG2=ON` | OFF | two-particle measurements (needs NFFT) |
| `-DBuild_Tests=OFF` | ON | skip the tests |

If you replace an existing cthyb install, use the same complex options it was built with: scripts
that pass complex Δ or H_loc need them. To read them back from that install, run
`grep -o 'BOOL:[^>]*>:[A-Z_]*IS_COMPLEX' <its prefix>/lib/cmake/triqs_cthyb/triqs_cthyb-targets.cmake`,
where `BOOL:>` means OFF and `BOOL:ON>` means ON. Build cthyb with the same compiler and flags as
TRIQS (for example the same `-march`), and on hardware that can run them.

**Python bindings.** `python/triqs_cthyb/solver_core.wrap.cxx` was edited by hand for the new
parameters. Do not configure with `-DUpdate_Python_Bindings=ON`: that regenerates the bindings from
the C++ headers, needs the clair tool, and is not required.

**Offline builds** (e.g. compute nodes without internet). Clone the two dependencies into `deps/`
beforehand; configure then uses them instead of downloading. The directory names must be exactly these:

```bash
git clone -b 0.9.x https://github.com/flatironinstitute/c2py cthyb.src/deps/c2py
git clone https://github.com/google/googletest cthyb.src/deps/GTest
```

---

## Upstream README

[![build](https://github.com/TRIQS/cthyb/workflows/build/badge.svg)](https://github.com/TRIQS/cthyb/actions?query=workflow%3Abuild)
[![conda-forge](https://img.shields.io/conda/vn/conda-forge/triqs_cthyb.svg)](https://anaconda.org/conda-forge/triqs_cthyb)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21690767-blue.svg)](https://doi.org/10.5281/zenodo.21690767)

### TRIQS continuous-time hybridization-expansion solver

Copyright (C) 2014, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet \
Copyright (C) 2017, H. U.R. Strand, N. Wentzell \
Copyright (C) 2018, The Simons Foundation

#### Documentation

For documentation and installation instructions see [triqs.github.io/cthyb](https://triqs.github.io/cthyb).

#### License

This application is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version (see http://www.gnu.org/licenses/).

It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

#### Support

<picture>
  <source media="(prefers-color-scheme: dark)" width="20%" srcset="doc/_static/CCQ-dark.png">
  <img alt="Flatiron Center for Computational Quantum Physics logo." width="20%" src="doc/_static/CCQ.png">
</picture>

TRIQS/cthyb is supported by the Flatiron Institute, a division of the Simons Foundation.
