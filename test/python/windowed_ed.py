# Test of the windowed moves (windowed_moves=True) against exact diagonalisation.
#
# A two-orbital Kanamori impurity (with spin-flip and pair-hopping) coupled to one bath site per
# spin-orbital is small enough (8 modes, 256 states) to be solved exactly with AtomDiag. Thermal
# averages of local operators (occupations, double occupancies, spin-flip and pair-hopping
# correlators, interaction energy) are compared with CTHYB in windowed mode, from the impurity
# density matrix (which, unlike G(tau), does not use M^-1 and has no heavy-tailed noise). The window
# is small (w = beta/10) so that the windowed proposals differ substantially from the stock ones:
#   - "diag": one block per spin-orbital, diagonal hybridisation. Exercises the 2-operator moves
#     and the different-block 4-operator moves;
#   - "offdiag": one 2x2 block per spin with an orbital-off-diagonal hybridisation. Exercises the
#     multi-orbital blocks and the same-block 4-operator moves (rank-2 det updates, d_min).
# Agreement within statistical errors requires both detailed balance and ergodicity of the
# windowed move set.
import sys
import numpy as np
from itertools import product
import triqs.utility.mpi as mpi
from triqs.gfs import GfImFreq, BlockGf
from triqs.operators import c, c_dag, n, Operator
from triqs.operators.util.hamiltonians import h_int_kanamori
from triqs.atom_diag import AtomDiag

beta, n_orb = 10.0, 2
U, J = 2.0, 0.3
eps_d = np.array([-2.0, -1.8])     # impurity levels, chemical potential included (near half filling)
spins = ('up', 'down')


def model(kind):
    """gf_struct, h_int, the impurity + bath Hamiltonian, and G0_iw of the impurity"""
    if kind == "diag":
        E = np.diag([-0.4, 0.5])                      # bath levels
        V = np.diag([0.7, 0.6])                       # V[orb, bath]
        blocks = [(f"{s}_{o}", s, [o]) for s in spins for o in range(n_orb)]
        off_diag = False
    else:
        E = np.diag([-0.4, 0.5])
        V = np.array([[0.7, 0.3], [0.2, 0.6]])
        blocks = [(s, s, list(range(n_orb))) for s in spins]
        off_diag = True
    gf_struct = [(name, len(orbs)) for name, _, orbs in blocks]
    # Kanamori operators carry the block names of set_operator_structure(spins, n_orb, off_diag)
    h_int = h_int_kanamori(spins, n_orb, np.array([[0, U - 3 * J], [U - 3 * J, 0]]), np.array([[U, U - 2 * J], [U - 2 * J, U]]), J,
                           off_diag=off_diag)

    def d(s, o):    # impurity operator (block name, index) as in gf_struct
        return (f"{s}_{o}", 0) if kind == "diag" else (s, o)

    H = Operator()
    H += h_int  # (H = h_int would alias: += acts in place)
    fops = []
    for s in spins:
        for o in range(n_orb):
            fops.append(d(s, o))
            H += eps_d[o] * n(*d(s, o))
        for k in range(n_orb):
            fops.append((f"bath_{s}", k))
            H += E[k, k] * n(f"bath_{s}", k)
            for o in range(n_orb):
                if V[o, k] != 0:
                    H += V[o, k] * (c_dag(*d(s, o)) * c(f"bath_{s}", k) + c_dag(f"bath_{s}", k) * c(*d(s, o)))

    def g0(name, orbs):
        g = GfImFreq(beta=beta, n_points=1025, target_shape=(len(orbs), len(orbs)))
        iw = np.array([complex(x) for x in g.mesh])
        Vb = V[orbs, :]
        delta = sum(np.einsum('i,j,w->wij', Vb[:, k], Vb[:, k], 1 / (iw - E[k, k])) for k in range(n_orb))
        g.data[:] = np.linalg.inv(iw[:, None, None] * np.eye(len(orbs)) - np.diag(eps_d[orbs]) - delta)
        return g

    G0 = BlockGf(name_list=[b[0] for b in blocks], block_list=[g0(name, orbs) for name, _, orbs in blocks], make_copies=False)
    return gf_struct, h_int, H, fops, G0


def observables(kind):
    """Local operators whose thermal averages are compared: occupations, double occupancies, the
    spin-flip and pair-hopping correlators between the orbitals, and the interaction energy"""
    gf_struct, h_int, H, fops, G0 = model(kind)
    d = (lambda s, o: (f"{s}_{o}", 0)) if kind == "diag" else (lambda s, o: (s, o))
    ops = {}
    for s in spins:
        for o in range(n_orb):
            ops[f"n_{s}{o}"] = n(*d(s, o))
    for o in range(n_orb):
        ops[f"docc_{o}"] = n(*d("up", o)) * n(*d("down", o))
    ops["spin_flip"] = c_dag(*d("up", 0)) * c(*d("down", 0)) * c_dag(*d("down", 1)) * c(*d("up", 1))
    ops["pair_hop"]  = c_dag(*d("up", 0)) * c_dag(*d("down", 0)) * c(*d("down", 1)) * c(*d("up", 1))
    ops["h_int"]     = h_int
    return ops


def ed_averages(kind):
    from triqs.atom_diag import atomic_density_matrix, trace_rho_op
    gf_struct, h_int, H, fops, G0 = model(kind)
    ad  = AtomDiag(H, fops)
    rho = atomic_density_matrix(ad, beta)
    return {k: trace_rho_op(rho, O, ad).real for k, O in observables(kind).items()}


def solve(kind, seed, **params):
    from triqs_cthyb import Solver
    gf_struct, h_int, H, fops, G0 = model(kind)
    S = Solver(beta=beta, gf_struct=gf_struct, n_tau=2051, n_iw=1025)
    S.G0_iw << G0
    p = dict(h_int=h_int, n_cycles=20000, length_cycle=50, n_warmup_cycles=2000, random_seed=seed, verbosity=0, measure_G_tau=False,
             measure_density_matrix=True, use_norm_as_weight=True, perform_post_proc=False)
    p.update(params)
    S.solve(**p)
    return S


def qmc_averages(S, kind):
    from triqs.atom_diag import trace_rho_op
    return {k: trace_rho_op(S.density_matrix, O, S.h_loc_diagonalization).real for k, O in observables(kind).items()}


def compare(kind, n_seeds, **params):
    """Mean and standard error over n_seeds independent runs, and z = (mean - exact) / error per observable.
    The seeds are spaced by 2: the default RNG (RandMT) gives the same stream for seeds 2k and 2k+1."""
    exact = ed_averages(kind)
    runs  = [qmc_averages(solve(kind, seed=1001 + 2 * k, **params), kind) for k in range(n_seeds)]
    res   = {}
    for key, ex in exact.items():
        v         = np.array([r[key] for r in runs])
        mean, err = v.mean(), v.std(ddof=1) / np.sqrt(n_seeds)
        res[key]  = (mean, err, ex, (mean - ex) / err)
    return res


if __name__ == "__main__":
    # 8 independent runs per case; |z| < 4.5 for every observable (t distribution with 7 degrees
    # of freedom: a false alarm has probability ~3e-3 per observable for an unbiased sampler)
    z_max = 4.5
    ok    = True
    for kind in ("diag", "offdiag"):
        res = compare(kind, 8, windowed_moves=True, window_half_width=beta / 10)
        for key, (mean, err, ex, z) in res.items():
            print(f"{kind:8s} {key:10s} windowed {mean: .5f} +- {err:.5f}   ED {ex: .5f}   z = {z: .2f}")
            ok = ok and abs(z) < z_max
    if not ok:
        print("FAILED: windowed CTHYB differs from exact diagonalisation")
        sys.exit(1)
