# Test of measure_G_auto_corr_time: the autocorrelation time of G(tau) averaged over a window around beta/2.
#
# A two-orbital Kanamori impurity (spin-flip and pair hopping) with one bath site per spin-orbital, one 2x2 block
# per spin with an orbital-off-diagonal hybridisation, solved with stock and with windowed moves:
#   - turning the measure on leaves the Markov chain unchanged: G_tau, the average order and the sign are
#     bit-identical with the same seed;
#   - the results have one entry per block and diagonal element, tau finite and non-negative, the saturation
#     flags 0 or 1;
#   - the window average of G(tau) from the measure's own samples agrees with G_tau averaged over the same
#     window (they differ only by the bins at the window edges);
#   - the results and the two parameters round-trip through HDF5 with the solver;
#   - a window outside (0, 1] is refused.
import numpy as np
import triqs.utility.mpi as mpi
from h5 import HDFArchive
from triqs.gfs import GfImFreq, BlockGf
from triqs.operators.util.hamiltonians import h_int_kanamori
from triqs_cthyb import Solver

beta, n_orb = 10.0, 2
U, J = 2.0, 0.3
eps_d = np.array([-2.0, -1.8])
E = np.array([-0.4, 0.5])
V = np.array([[0.7, 0.3], [0.2, 0.6]])
spins = ("up", "down")
window = 0.1
n_tau = 1001
n_iw = 400

gf_struct = [(s, n_orb) for s in spins]
h_int = h_int_kanamori(spins, n_orb, np.array([[0, U - 3 * J], [U - 3 * J, 0]]), np.array([[U, U - 2 * J], [U - 2 * J, U]]), J,
                       off_diag=True)


def g0():
    g = GfImFreq(beta=beta, n_points=n_iw, target_shape=(n_orb, n_orb))
    iw = np.array([complex(x) for x in g.mesh])
    delta = sum(np.einsum("i,j,w->wij", V[:, k], V[:, k], 1 / (iw - E[k])) for k in range(n_orb))
    g.data[:] = np.linalg.inv(iw[:, None, None] * np.eye(n_orb) - np.diag(eps_d) - delta)
    return g


def solve(windowed, measure, seed=34788 + 928374 * mpi.rank, **extra):
    S = Solver(beta=beta, gf_struct=gf_struct, n_tau=n_tau, n_iw=n_iw)
    for s in spins:
        S.G0_iw[s] << g0()
    params = dict(h_int=h_int, n_cycles=20000, length_cycle=50, n_warmup_cycles=2000, random_seed=seed, random_name="mt19937",
                  verbosity=0, measure_G_auto_corr_time=measure, G_auto_corr_window=window, **extra)
    if windowed:
        params.update(windowed_moves=True, window_half_width=beta / 10)
    S.solve(**params)
    return S


def window_average(S):
    """G_tau averaged over the mesh points within the window, per block and diagonal element"""
    out = {}
    for name, g in S.G_tau:
        tau = np.array([float(t) for t in g.mesh])
        mid = np.abs(tau - 0.5 * beta) <= 0.5 * window * beta
        out[name] = np.diagonal(g.data[mid], axis1=1, axis2=2).real.mean(axis=0)
    return out


for windowed in (False, True):
    label = "windowed" if windowed else "stock"
    off = solve(windowed, False)
    on = solve(windowed, True)

    # the measure only reads the configuration: same chain, same results
    for name, g in off.G_tau:
        assert np.array_equal(g.data, on.G_tau[name].data), f"{label}: G_tau changed by the measure"
    assert off.average_order == on.average_order and off.average_sign == on.average_sign, f"{label}: chain changed"
    assert off.G_auto_corr_time == {}, f"{label}: G_auto_corr_time filled with the measure off"

    tau, conv, mean = on.G_auto_corr_time, on.G_auto_corr_time_converged, on.G_auto_corr_mean
    assert sorted(tau) == sorted(spins) and all(len(tau[s]) == n_orb for s in spins), tau
    for s in spins:
        assert all(np.isfinite(t) and t >= 0 for t in tau[s]), tau
        assert all(c in (0, 1) for c in conv[s]), conv

    # the measure's window average against G_tau over the same window
    ref = window_average(on)
    for s in spins:
        diff = np.abs(np.array(mean[s]) - ref[s])
        assert np.all(diff <= 0.03 * np.abs(ref[s]) + 2e-3), f"{label} {s}: {mean[s]} against G_tau {ref[s]}"

    # the results round-trip through HDF5 with the solver
    if mpi.is_master_node():
        with HDFArchive(f"G_auto_corr_time_{label}.h5", "w") as A:
            A["solver"] = on
        with HDFArchive(f"G_auto_corr_time_{label}.h5", "r") as A:
            back = A["solver"]
        for key in ("G_auto_corr_time", "G_auto_corr_time_converged", "G_auto_corr_mean"):
            a, b = getattr(on, key), getattr(back, key)
            assert sorted(a) == sorted(b) and all(np.array_equal(a[k], b[k]) for k in a), f"{label}: {key} lost in HDF5"
        assert back.last_solve_parameters.measure_G_auto_corr_time
        assert back.last_solve_parameters.G_auto_corr_window == window

    if mpi.is_master_node():
        print(f"{label}: tau_G {tau} (saturated {conv}), auto_corr_time {on.auto_corr_time:.3f}; "
              f"window mean {mean} against G_tau {ref}")

# a window outside (0, 1] is refused
try:
    S = Solver(beta=beta, gf_struct=gf_struct, n_tau=n_tau, n_iw=n_iw)
    for s in spins:
        S.G0_iw[s] << g0()
    S.solve(h_int=h_int, n_cycles=10, length_cycle=10, n_warmup_cycles=0, verbosity=0, measure_G_auto_corr_time=True,
            G_auto_corr_window=0.0)
    raise AssertionError("G_auto_corr_window = 0 was accepted")
except RuntimeError:
    pass
