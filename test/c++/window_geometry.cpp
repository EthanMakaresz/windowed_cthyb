// Tests of the geometry of the windowed moves (c++/triqs_cthyb/moves/window_util.hpp) against brute force:
// - near/far sets (arc_set: counts, k-th element, rank) and insertion positions on decreasing time sequences, with
//   wrap-around, t0 inside and outside the sequence, points exactly at the window edges;
// - the draws land exactly on the window sets, including their first and last points;
// - same_block_proposal_ratio equals sum_P r(P) / sum_P f(P), with r obtained by enumerating every choice
//   sequence of the remove move that removes a given set of 4 operators, and f by enumerating the labelled
//   insertion paths. The counts are also computed the way the insert move does (det + new operators).
#include <triqs/test_tools/arrays.hpp>
#include <triqs_cthyb/moves/window_util.hpp>

#include <algorithm>
#include <cmath>
#include <random>
#include <set>
#include <vector>

using namespace triqs_cthyb::window;
using tp = triqs::utility::time_pt;

constexpr double beta = 2.5;

struct fixed_rng { // returns a fixed value, to probe the edges of the draws
  uint64_t v;
  uint64_t operator()(uint64_t) const { return v; }
};

struct test_rng { // uniform integers, like mc_tools::random_generator
  std::mt19937_64 gen;
  template <typename T> T operator()(T n) { return std::uniform_int_distribution<T>(0, n - 1)(gen); }
};

std::vector<tp> random_times(std::mt19937_64 &gen, int n, std::vector<tp> extra = {}) {
  std::set<uint64_t> s;
  auto e = extra;
  for (auto &t : e) s.insert(grid_index(t, beta));
  while (int(s.size()) < n) s.insert(std::uniform_int_distribution<uint64_t>(0, Nmax - 1)(gen));
  std::vector<tp> v;
  for (auto i : s) v.emplace_back(i, beta);
  std::sort(v.begin(), v.end(), std::greater<>()); // decreasing, as in the dets
  return v;
}

tp shift(tp t, long long d) { return d >= 0 ? t + tp{uint64_t(d), beta} : t - tp{uint64_t(-d), beta}; }

TEST(window, counts_selection_positions) {
  std::mt19937_64 gen(1);
  for (double w : {0.01, 0.3, 1.0, 1.25}) {
    window_t win(w, beta);
    distance_t dist(w == 1.25 ? 0.0 : 0.4 * w, beta); // includes w = beta/2 and d_min = 0
    for (int trial = 0; trial < 400; ++trial) {
      tp t0{std::uniform_int_distribution<uint64_t>(0, Nmax - 1)(gen), beta};
      if (trial % 4 == 1) t0 = tp{std::uniform_int_distribution<uint64_t>(0, 1000)(gen), beta};               // near 0: wrap
      if (trial % 4 == 2) t0 = tp{Nmax - 1 - std::uniform_int_distribution<uint64_t>(0, 1000)(gen), beta}; // near beta: wrap
      // operators at the window edges and at t0 itself (trial odd)
      long long W = (long long)win.Wn, D = (long long)dist.Dn;
      std::vector<tp> extra = {shift(t0, W - 1), shift(t0, -(W - 1)), shift(t0, W), shift(t0, -W),
                               shift(t0, D), shift(t0, -D), shift(t0, D + 1), shift(t0, -(D + 1))};
      if (trial % 2) extra.push_back(t0);
      int n  = int(extra.size()) + int(gen() % 12);
      auto v = random_times(gen, n, extra);
      n      = v.size();
      auto get = [&v](long i) -> tp const & { return v[i]; };

      // brute force
      std::set<long> near_bf, far_bf;
      for (long i = 0; i < n; ++i) {
        if (win.near(v[i], t0)) near_bf.insert(i);
        if (dist.far(v[i], t0)) far_bf.insert(i);
      }
      long below_bf = std::find_if(v.begin(), v.end(), [&](tp const &t) { return t < t0; }) - v.begin();

      EXPECT_EQ(count_near(get, n, t0, win), long(near_bf.size()));
      EXPECT_EQ(count_far(get, n, t0, dist), long(far_bf.size()));
      EXPECT_EQ(first_below(get, n, t0), below_bf);

      auto ns = near_set(get, n, t0, win);
      auto fs = far_set(get, n, t0, dist);
      std::set<long> near_sel, far_sel;
      for (long k = 0; k < ns.count; ++k) {
        near_sel.insert(ns.index(k));
        EXPECT_EQ(ns.rank(ns.index(k)), k);
      }
      for (long k = 0; k < fs.count; ++k) {
        far_sel.insert(fs.index(k));
        EXPECT_EQ(fs.rank(fs.index(k)), k);
      }
      EXPECT_EQ(near_sel, near_bf);
      EXPECT_EQ(far_sel, far_bf);
      for (long i = 0; i < n; ++i) {
        if (!near_bf.count(i)) { EXPECT_EQ(ns.rank(i), -1); }
        if (!far_bf.count(i)) { EXPECT_EQ(fs.rank(i), -1); }
        // the arcs are the predicates
        EXPECT_EQ(win.arc(t0).contains(v[i]), win.near(v[i], t0));
        EXPECT_EQ(dist.arc(t0).contains(v[i]), dist.far(v[i], t0));
      }
    }
  }
}

TEST(window, draws_and_edges) {
  std::mt19937_64 gen(2);
  for (double w : {1e-6, 0.3, 1.25}) {
    window_t win(w, beta);
    EXPECT_NEAR(win.length, 2 * w, 1e-12 * beta);
    for (uint64_t t0n : {uint64_t(0), uint64_t(5), Nmax - 3, Nmax / 3}) {
      tp t0{t0n, beta};
      // first and last points of the draw are near; one step further is not (unless w = beta/2)
      fixed_rng r_first{0}, r_last{win.n_pts - 1};
      tp first = win.draw(t0, r_first), last = win.draw(t0, r_last);
      EXPECT_TRUE(win.near(first, t0));
      EXPECT_TRUE(win.near(last, t0));
      if (w < 0.5 * beta) {
        EXPECT_FALSE(win.near(shift(first, -1), t0));
        EXPECT_FALSE(win.near(shift(last, 1), t0));
      }
      test_rng rng{std::mt19937_64(gen())};
      for (int k = 0; k < 200; ++k) EXPECT_TRUE(win.near(win.draw(t0, rng), t0));
    }
  }
  window_t half(0.5 * beta, beta);
  EXPECT_EQ(half.n_pts, Nmax); // the whole circle
  for (double d : {0.0, 0.2, 1.2}) {
    distance_t dist(d, beta);
    EXPECT_NEAR(dist.length, beta - 2 * d, 1e-12 * beta);
    for (uint64_t t0n : {uint64_t(0), Nmax - 3, Nmax / 3}) {
      tp t0{t0n, beta};
      fixed_rng r_first{0}, r_last{dist.n_pts - 1};
      tp first = dist.draw(t0, r_first), last = dist.draw(t0, r_last);
      EXPECT_TRUE(dist.far(first, t0));
      EXPECT_TRUE(dist.far(last, t0));
      EXPECT_FALSE(dist.far(shift(first, -1), t0));
      EXPECT_FALSE(dist.far(shift(last, 1), t0));
      EXPECT_FALSE(dist.far(t0, t0));
      test_rng rng{std::mt19937_64(gen())};
      for (int k = 0; k < 200; ++k) EXPECT_TRUE(dist.far(dist.draw(t0, rng), t0));
    }
  }
}

// sum over the remove move's choice sequences that remove {cdag a, b} and {c p, q} from (xs, ys)
double brute_sum_r(std::vector<tp> const &xs, std::vector<tp> const &ys, std::set<long> const &cd, std::set<long> const &cc, window_t const &win,
                   distance_t const &dist) {
  long n     = xs.size();
  double sum = 0;
  for (long X = 0; X < n; ++X) {
    std::vector<long> xc, Yc;
    for (long i = 0; i < n; ++i)
      if (win.near(ys[i], xs[X])) xc.push_back(i);
    for (long i = 0; i < n; ++i)
      if (dist.far(xs[i], xs[X])) Yc.push_back(i);
    for (long x : xc)
      for (long Y : Yc) {
        std::vector<long> yc;
        for (long i = 0; i < n; ++i)
          if (i != x && win.near(ys[i], xs[Y])) yc.push_back(i);
        for (long y : yc)
          if (std::set<long>{X, Y} == cd && std::set<long>{x, y} == cc)
            sum += 1.0 / (double(n) * double(xc.size()) * double(Yc.size()) * double(yc.size()));
      }
  }
  return sum;
}

TEST(window, same_block_proposal_ratio) {
  std::mt19937_64 gen(3);
  test_rng rng{std::mt19937_64(4)};
  int n_multi_path = 0;
  // d_min < 2 w4 (all 4 paths possible), and the defaults' regime d_min >= 2 w4
  for (auto [w4, dmin] : std::vector<std::pair<double, double>>{{0.6, 0.3}, {0.9, 0.1}, {0.2, 0.6}, {1.25, 0.0}}) {
    window_t win(w4, beta);
    distance_t dist(dmin, beta);
    for (int trial = 0; trial < 300; ++trial) {
      int s = 1 + gen() % 3;
      // C: n random c^dag and c; the inserted pairs drawn as the insert move does
      int n   = gen() % 7;
      auto xs = random_times(gen, n), ys = random_times(gen, n);
      tp X1 = tp{std::uniform_int_distribution<uint64_t>(0, Nmax - 1)(gen), beta};
      tp x1 = win.draw(X1, rng), X2 = dist.draw(X1, rng), x2 = win.draw(X2, rng);
      if (X1 == X2 || x1 == x2) continue;
      auto getx = [&xs](long i) -> tp const & { return xs[i]; };
      auto gety = [&ys](long i) -> tp const & { return ys[i]; };

      // counts as in move_insert_window_double: det of C + the new operators
      long c_near_ins[2]   = {count_near(gety, n, X1, win) + win.near(x1, X1) + win.near(x2, X1),
                              count_near(gety, n, X2, win) + win.near(x1, X2) + win.near(x2, X2)};
      long cdag_far_ins[2] = {count_far(getx, n, X1, dist) + dist.far(X2, X1), count_far(getx, n, X2, dist) + dist.far(X1, X2)};
      double ratio_ins     = same_block_proposal_ratio(n + 2, {x1, x2}, {win.arc(X1), win.arc(X2)}, c_near_ins, cdag_far_ins, dist.far(X2, X1), win,
                                                       dist, s);

      // C' and the counts as in move_remove_window_double: on C' itself
      auto xs2 = xs, ys2 = ys;
      xs2.push_back(X1), xs2.push_back(X2), ys2.push_back(x1), ys2.push_back(x2);
      std::sort(xs2.begin(), xs2.end(), std::greater<>());
      std::sort(ys2.begin(), ys2.end(), std::greater<>());
      long n2    = n + 2;
      auto getx2 = [&xs2](long i) -> tp const & { return xs2[i]; };
      auto gety2 = [&ys2](long i) -> tp const & { return ys2[i]; };
      long c_near_rem[2]   = {count_near(gety2, n2, X1, win), count_near(gety2, n2, X2, win)};
      long cdag_far_rem[2] = {count_far(getx2, n2, X1, dist), count_far(getx2, n2, X2, dist)};
      EXPECT_EQ(c_near_rem[0], c_near_ins[0]);
      EXPECT_EQ(c_near_rem[1], c_near_ins[1]);
      EXPECT_EQ(cdag_far_rem[0], cdag_far_ins[0]);
      EXPECT_EQ(cdag_far_rem[1], cdag_far_ins[1]);
      // the remove move may name the pairs either way round
      double ratio_rem = same_block_proposal_ratio(n2, {x2, x1}, {win.arc(X2), win.arc(X1)}, {c_near_rem[1], c_near_rem[0]},
                                                   {cdag_far_rem[1], cdag_far_rem[0]}, dist.far(X1, X2), win, dist, s);

      // brute force
      auto idx = [](std::vector<tp> const &v, tp t) { return long(std::find(v.begin(), v.end(), t) - v.begin()); };
      double sum_r = brute_sum_r(xs2, ys2, {idx(xs2, X1), idx(xs2, X2)}, {idx(ys2, x1), idx(ys2, x2)}, win, dist);
      tp const Xs[2] = {X1, X2}, xs_[2] = {x1, x2};
      int n_paths = 0;
      for (int a = 0; a < 2; ++a)
        for (int p = 0; p < 2; ++p)
          if (win.near(xs_[p], Xs[a]) && dist.far(Xs[1 - a], Xs[a]) && win.near(xs_[1 - p], Xs[1 - a])) ++n_paths;
      double f_1   = 1.0 / (beta * win.length * dist.length * win.length * std::pow(double(s), 4));
      double sum_f = n_paths * f_1;
      if (n_paths > 2) ++n_multi_path;

      ASSERT_GT(sum_r, 0);
      EXPECT_NEAR(ratio_ins / (sum_r / sum_f), 1.0, 1e-12);
      EXPECT_NEAR(ratio_rem / (sum_r / sum_f), 1.0, 1e-12);
    }
  }
  EXPECT_GT(n_multi_path, 0); // the cross-paired paths were exercised
}


