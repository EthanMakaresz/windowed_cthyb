/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 ******************************************************************************/
#pragma once
#include <triqs/utility/time_pt.hpp>
#include <triqs/utility/exceptions.hpp>

#include <algorithm>
#include <cstdint>

// Geometry of the windowed moves (moves/insert_window, remove_window, double_insert_window,
// double_remove_window), on the integer time grid of time_pt.
//
// The imaginary-time circle is the grid Z_Nmax (time_pt arithmetic is cyclic modulo Nmax). For two
// points t, t0 let e = (t0 - t) mod Nmax. A window of half-width w has Wn = make_time_pt(w) grid steps and
//     near(t, t0)  <=>  e < Wn  or  e > Nmax - Wn                 (2 Wn - 1 points; t0 itself included),
// and a minimum distance d_min has Dn grid steps and
//     far(t, t0)   <=>  Dn < e < Nmax - Dn                        (Nmax - 2 Dn - 1 points; t0 never included).
// Both relations are symmetric in (t, t0). Equivalently, the near points of t0 form the cyclic arc
// [t0 - (Wn - 1), t0 + (Wn - 1)] and the far points the arc [t0 + Dn + 1, t0 - Dn - 1] (arc_t). The proposals draw exactly uniformly from these sets, and the
// reverse moves count with exactly the same predicates, so there is no rounding at the window edge.
// The measure of a set in imaginary time is (number of grid points) * beta / Nmax, as for the stock
// moves, which draw a time uniformly among the Nmax points (density 1/beta).
//
// The dets store the c^dagger (x) and the c (y) in decreasing time order: get_x(i).first, get_y(i).first.
// The helpers below take an accessor get(i) -> time_pt of such a decreasing sequence, of length n.
// The near set of t0 is the cyclic arc of grid points [t0 - (Wn - 1), t0 + (Wn - 1)], the far set the arc
// [t0 + Dn + 1, t0 - Dn - 1]. The end points are computed once per query; the elements of the sequence in an
// arc then form one cyclically contiguous index range, found with two binary searches that only compare
// time points (no time arithmetic per probe): O(log n) per count or draw, no allocation.

namespace triqs_cthyb::window {

  using triqs::utility::time_pt;
  constexpr uint64_t Nmax = time_pt::Nmax;

  // Grid index of a time point (time_pt::n is private; floor_div by one grid step returns it)
  inline uint64_t grid_index(time_pt const &t, double beta) { return floor_div(t, time_pt{1, beta}); }

  // Draw uniformly in [0, n), n >= 1. Clamped: rng(n) = floor(n u) can round up to n for n ~ 2^64.
  template <typename RNG> inline uint64_t draw_index(RNG &rng, uint64_t n) { return std::min<uint64_t>(rng(n), n - 1); }

  // The cyclic arc of grid points [lo, hi], from lo forward in time to hi (lo > hi: it wraps through 0).
  // Membership is two comparisons, without time arithmetic.
  struct arc_t {
    time_pt lo, hi;
    bool contains(time_pt const &t) const { return lo <= hi ? (lo <= t && t <= hi) : (t >= lo || t <= hi); }
  };

  // Window of half-width w: the near relation, its draw and its measure
  struct window_t {
    double beta = 1;
    uint64_t Wn = 1;
    time_pt lo, hi; // e < lo or e > hi
    time_pt offset; // Wn - 1: the near arc of t0 is [t0 - offset, t0 + offset]
    uint64_t n_pts = 1; // 2 Wn - 1
    double length  = 0; // n_pts beta / Nmax, ~ 2 w

    window_t() = default;
    window_t(double w, double beta_) : beta(beta_) {
      if (!(w > 0 && w <= 0.5 * beta)) TRIQS_RUNTIME_ERROR << "window half-width " << w << " must be in (0, beta/2], beta = " << beta;
      Wn = grid_index(triqs::utility::time_segment{beta}.make_time_pt(w), beta);
      Wn = std::min<uint64_t>(Wn, (Nmax - 1) / 2 + 1); // 2 Wn - 1 <= Nmax
      if (Wn == 0) TRIQS_RUNTIME_ERROR << "window half-width " << w << " is below the time grid resolution";
      lo     = time_pt{Wn, beta};
      hi     = time_pt{Nmax - Wn, beta};
      offset = time_pt{Wn - 1, beta};
      n_pts  = 2 * Wn - 1;
      length = double(n_pts) * (beta / double(Nmax));
    }
    bool near_e(time_pt const &e) const { return e < lo || e > hi; }
    bool near(time_pt const &t, time_pt const &t0) const { return near_e(t0 - t); }
    // The near points of t0 as an arc: near(t, t0) <=> arc(t0).contains(t)
    arc_t arc(time_pt const &t0) const { return {t0 - offset, t0 + offset}; }
    // A point uniformly in the window around t0: t0 + (k - (Wn - 1)), k uniform in [0, 2 Wn - 1)
    template <typename RNG> time_pt draw(time_pt const &t0, RNG &rng) const {
      uint64_t k = draw_index(rng, n_pts), o = Wn - 1;
      return t0 + time_pt{k >= o ? k - o : Nmax - (o - k), beta};
    }
  };

  // Minimum distance d_min: the far relation, its draw and its measure
  struct distance_t {
    double beta = 1;
    uint64_t Dn = 0;
    time_pt lo, hi; // lo < e < hi
    time_pt step;   // Dn + 1: the far arc of t0 is [t0 + step, t0 - step]
    uint64_t n_pts = 1; // Nmax - 2 Dn - 1
    double length  = 0; // n_pts beta / Nmax, ~ beta - 2 d_min

    distance_t() = default;
    distance_t(double d_min, double beta_) : beta(beta_) {
      if (!(d_min >= 0 && d_min < 0.5 * beta)) TRIQS_RUNTIME_ERROR << "minimum distance " << d_min << " must be in [0, beta/2), beta = " << beta;
      Dn = grid_index(triqs::utility::time_segment{beta}.make_time_pt(d_min), beta);
      if (Dn >= (Nmax - 1) / 2) TRIQS_RUNTIME_ERROR << "minimum distance " << d_min << " too close to beta/2";
      lo     = time_pt{Dn, beta};
      hi     = time_pt{Nmax - Dn, beta};
      step   = time_pt{Dn + 1, beta};
      n_pts  = Nmax - 2 * Dn - 1;
      length = double(n_pts) * (beta / double(Nmax));
    }
    bool far_e(time_pt const &e) const { return e > lo && e < hi; }
    bool far(time_pt const &t, time_pt const &t0) const { return far_e(t0 - t); }
    // The far points of t0 as an arc: far(t, t0) <=> arc(t0).contains(t)
    arc_t arc(time_pt const &t0) const { return {t0 + step, t0 - step}; }
    // A point uniformly among those far from t0
    template <typename RNG> time_pt draw(time_pt const &t0, RNG &rng) const { return t0 + time_pt{Dn + 1 + draw_index(rng, n_pts), beta}; }
  };

  // First index i in [0, n) with get(i) < t (n if none). get is decreasing.
  // Also the insertion position of a new time t in the det.
  template <typename Get> long first_below(Get const &get, long n, time_pt const &t) {
    long lo = 0, hi = n;
    while (lo < hi) {
      long mid = lo + (hi - lo) / 2;
      if (get(mid) < t)
        hi = mid;
      else
        lo = mid + 1;
    }
    return lo;
  }

  // First index i in [0, n) with get(i) <= t (n if none)
  template <typename Get> long first_le(Get const &get, long n, time_pt const &t) {
    long lo = 0, hi = n;
    while (lo < hi) {
      long mid = lo + (hi - lo) / 2;
      if (get(mid) <= t)
        hi = mid;
      else
        lo = mid + 1;
    }
    return lo;
  }

  // The elements of a decreasing sequence that lie in the cyclic arc [lo, hi] of grid points (from lo forward in
  // time to hi; lo > hi means the arc wraps through 0). They are the cyclically contiguous indices
  // start, start + 1, ... (mod n), count of them.
  struct arc_set {
    long n = 0, start = 0, count = 0;
    template <typename Get> arc_set(Get const &get, long n_, arc_t const &arc) : n(n_) {
      auto const &[lo, hi] = arc;
      long a = first_le(get, n, hi);    // first index with t <= hi
      long b = first_below(get, n, lo); // first index with t < lo
      start  = a; // may be n (wrap case with no element <= hi): index and rank wrap it
      count  = (lo <= hi ? b - a : (n - a) + b); // [a, b) ; or [a, n) then [0, b)
    }
    // Index of the k-th element, k in [0, count)
    long index(long k) const {
      long i = start + k;
      return i < n ? i : i - n;
    }
    // Rank of the element at index i, or -1 if it is not in the arc
    long rank(long i) const {
      long k = (i >= start ? i - start : i + n - start);
      return k < count ? k : -1;
    }
  };

  // The elements near t0, and those far from t0
  template <typename Get> arc_set near_set(Get const &get, long n, time_pt const &t0, window_t const &win) {
    return arc_set(get, n, win.arc(t0));
  }
  template <typename Get> arc_set far_set(Get const &get, long n, time_pt const &t0, distance_t const &dist) {
    return arc_set(get, n, dist.arc(t0));
  }

  template <typename Get> long count_near(Get const &get, long n, time_pt const &t0, window_t const &win) {
    return near_set(get, n, t0, win).count;
  }
  template <typename Get> long count_far(Get const &get, long n, time_pt const &t0, distance_t const &dist) {
    return far_set(get, n, t0, dist).count;
  }

  // ---------------------------------------------------------------------------------------------
  // Proposal ratio of the windowed 4-operator moves within one block (double_insert_window /
  // double_remove_window with block1 = block2), s = block size.
  //
  // Insert, C -> C': X = c^dag uniform on the circle, x = c uniform in the window w4 of X,
  // Y = c^dag uniform among the points far (> d_min) from X, y = c uniform in the window w4 of Y,
  // inner indices uniform. A labelled path P = (X, x, Y, y) has the density
  //     f(P) = I(P) / (beta L_w4 L_far L_w4 s^4),   I(P) = [x near X] [Y far X] [y near Y].
  // Remove, C' -> C: X uniform among the n c^dag, x uniform among the m_x(X) c near X, Y uniform among
  // the m_Y(X) c^dag far from X (X itself is never far), y uniform among the c near Y other than x:
  //     r(P) = I(P) / (n m_x(X) m_Y(X) (m_x(Y) - [x near Y])),   all counts in C'.
  // The same C' is reached by the 4 labelled paths that name the same 4 operators (which new c^dag is X,
  // times which new c goes with it). The kernels are sums over these paths, hence
  //     R_insert = W(C')/W(C) sum_P r(P) / sum_P f(P),   R_remove = 1 / (that proposal ratio),
  // and detailed balance holds between the insert and the remove of the same block.
  // [Y far X] is symmetric, so common to the 4 paths; with d_min >= 2 w4 only (X1,x1,X2,x2) and
  // (X2,x2,X1,x1) can satisfy I.
  //
  // Arguments: n_cdag = n in C'; x[k] the c paired with Xk in C' (k = 0, 1); near_X[k] = the near arc of Xk;
  // c_near[k] = m_x(Xk), the number of c near Xk in C'; cdag_far[k] = m_Y(Xk), the number of c^dag far
  // from Xk in C'; X_far = [X2 far X1]. Returns sum_P r / sum_P f, or 0 if no path is valid.
  inline double same_block_proposal_ratio(long n_cdag, time_pt const (&x)[2], arc_t const (&near_X)[2], long const (&c_near)[2],
                                          long const (&cdag_far)[2], bool X_far, window_t const &w4, distance_t const &dist, int block_size) {
    if (!X_far) return 0;
    double sum_r = 0;
    int n_valid  = 0;
    for (int a = 0; a < 2; ++a)   // X = X[a], Y = X[1-a]
      for (int p = 0; p < 2; ++p) { // x = x[p], y = x[1-p]
        int b = 1 - a, q = 1 - p;
        if (!(near_X[a].contains(x[p]) && near_X[b].contains(x[q]))) continue;
        ++n_valid;
        long m_y = c_near[b] - (near_X[b].contains(x[p]) ? 1 : 0);
        sum_r += 1.0 / (double(n_cdag) * double(c_near[a]) * double(cdag_far[a]) * double(m_y));
      }
    if (n_valid == 0) return 0;
    double s2  = double(block_size) * block_size;
    double f_1 = 1.0 / (w4.beta * w4.length * dist.length * w4.length * s2 * s2); // f(P) of one valid path
    return sum_r / (n_valid * f_1);
  }

} // namespace triqs_cthyb::window
