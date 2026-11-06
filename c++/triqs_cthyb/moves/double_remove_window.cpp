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
#include "./double_remove_window.hpp"

namespace triqs_cthyb {

  move_remove_window_double::move_remove_window_double(int block_index1, int block_index2, int block_size1, int block_size2, double w4,
                                                       double d_min, qmc_data &data, mc_tools::random_generator &rng)
     : data(data),
       config(data.config),
       rng(rng),
       block_index1(block_index1),
       block_index2(block_index2),
       block_size1(block_size1),
       block_size2(block_size2),
       win(w4, data.config.beta()),
       dist(d_min, data.config.beta()) {}

  mc_weight_t move_remove_window_double::attempt() {

    bool same_block = (block_index1 == block_index2);
    double beta     = config.beta();
    auto &det1      = data.dets[block_index1];
    auto &det2      = data.dets[block_index2];
    long det1_size  = det1.size();
    long det2_size  = det2.size();
    if (same_block) {
      if (det1_size < 2) return 0; // det does not contain two operators
    } else {
      if ((det1_size == 0) || (det2_size == 0)) return 0; // one of two dets is empty
    }
    auto get_x1 = [&det1](long i) -> time_pt const & { return det1.get_x(i).first; };
    auto get_y1 = [&det1](long i) -> time_pt const & { return det1.get_y(i).first; };
    auto get_x2 = [&det2](long i) -> time_pt const & { return det2.get_x(i).first; };
    auto get_y2 = [&det2](long i) -> time_pt const & { return det2.get_y(i).first; };
    det_scalar_t det_ratio;

    // X: any C^dagger of block 1; x: a C of block 1 near X
    long num_c_dag1 = rng(det1_size);
    auto X          = get_x1(num_c_dag1);
    auto const near_arc_X = win.arc(X);
    window::arc_set near_x(get_y1, det1_size, near_arc_X);
    if (near_x.count == 0) return 0;
    long num_c1 = near_x.index(rng(near_x.count));
    auto x      = get_y1(num_c1);

    // Y: any C^dagger of block 2 (same block: one far from X)
    long num_c_dag2, m_Y = 0;
    if (same_block) {
      auto far_Y = window::far_set(get_x1, det1_size, X, dist);
      m_Y        = far_Y.count;
      if (m_Y == 0) return 0;
      num_c_dag2 = far_Y.index(rng(m_Y));
    } else {
      num_c_dag2 = rng(det2_size);
    }
    auto Y = get_x2(num_c_dag2);

    // y: a C of block 2 near Y (same block: other than x)
    auto const near_arc_Y = win.arc(Y);
    window::arc_set near_y(get_y2, det2_size, near_arc_Y);
    long rank_x = (same_block ? near_y.rank(num_c1) : -1);
    long m_y    = near_y.count - (rank_x >= 0 ? 1 : 0);
    if (m_y == 0) return 0;
    long k = rng(m_y);
    if (rank_x >= 0 && k >= rank_x) ++k;
    long num_c2 = near_y.index(k);

    // Proposal ratio of the insertion that leads to the current configuration (counts before the removal)
    mc_weight_t t_ratio;
    if (same_block) {
      // same det: near_y counts the c near Y; X is far from Y by construction
      window::arc_t const near_X[2] = {near_arc_X, near_arc_Y};
      long const c_near[2]          = {near_x.count, near_y.count};
      long const cdag_far[2]        = {m_Y, window::count_far(get_x1, det1_size, Y, dist)};
      time_pt const xs[2]           = {x, get_y1(num_c2)};
      t_ratio = window::same_block_proposal_ratio(det1_size, xs, near_X, c_near, cdag_far, true, win, dist, block_size1);
      if (t_ratio == 0) TRIQS_RUNTIME_ERROR << "(remove window double) the chosen path is outside the windows";
    } else {
      double s = double(block_size1) * block_size2;
      t_ratio  = s * s * beta * beta * win.length * win.length / (double(det1_size) * double(near_x.count) * double(det2_size) * double(near_y.count));
    }

    // now mark 4 nodes for deletion
    tau1 = data.imp_trace.try_delete(num_c1, block_index1, false);
    tau2 = data.imp_trace.try_delete(num_c_dag1, block_index1, true);
    tau3 = data.imp_trace.try_delete(num_c2, block_index2, false);
    tau4 = data.imp_trace.try_delete(num_c_dag2, block_index2, true);

    if (same_block) {
      det_ratio = det1.try_remove2(num_c_dag1, num_c_dag2, num_c1, num_c2);
    } else {
      auto det_ratio1 = det1.try_remove(num_c_dag1, num_c1);
      auto det_ratio2 = det2.try_remove(num_c_dag2, num_c2);
      det_ratio       = det_ratio1 * det_ratio2;
    }

    // For quick abandon
    double random_number = rng.preview();
    if (random_number == 0.0) return 0;
    double p_yee = std::abs(det_ratio / t_ratio / data.atomic_weight);

    // recompute the trace
    std::tie(new_atomic_weight, new_atomic_reweighting) = data.imp_trace.compute(p_yee, random_number);
    if (new_atomic_weight == 0.0) return 0;
    auto atomic_weight_ratio = new_atomic_weight / data.atomic_weight;
    if (!isfinite(atomic_weight_ratio))
      TRIQS_RUNTIME_ERROR << "(remove window double) atomic_weight_ratio not finite " << new_atomic_weight << " " << data.atomic_weight
                          << " in config " << config.get_id();

    mc_weight_t p = atomic_weight_ratio * det_ratio;
    if (!isfinite(p / t_ratio))
      TRIQS_RUNTIME_ERROR << "(remove window double) p / t_ratio not finite p : " << p << " t_ratio :  " << t_ratio << " in config "
                          << config.get_id();
    return p / t_ratio;
  }

  mc_weight_t move_remove_window_double::accept() {

    // remove from the tree
    data.imp_trace.confirm_delete();

    // remove from the configuration
    config.erase(tau1);
    config.erase(tau2);
    config.erase(tau3);
    config.erase(tau4);
    config.finalize();

    // remove from the determinants
    data.dets[block_index1].complete_operation();
    if (block_index1 != block_index2) data.dets[block_index2].complete_operation();
    data.update_sign();
    data.atomic_weight      = new_atomic_weight;
    data.atomic_reweighting = new_atomic_reweighting;
    return data.current_sign / data.old_sign;
  }

  void move_remove_window_double::reject() {
    config.finalize();
    data.imp_trace.cancel_delete();
    data.dets[block_index1].reject_last_try();
    if (block_index1 != block_index2) data.dets[block_index2].reject_last_try();
  }
} // namespace triqs_cthyb
