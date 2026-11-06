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
#include "./double_insert_window.hpp"

namespace triqs_cthyb {

  move_insert_window_double::move_insert_window_double(int block_index1, int block_index2, int block_size1, int block_size2, double w4,
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

  mc_weight_t move_insert_window_double::attempt() {

    bool same_block = (block_index1 == block_index2);
    double beta     = config.beta();

    // Pick up the value of alpha and choose the operators
    auto rs1 = rng(block_size1), rs2 = rng(block_size1), rs3 = rng(block_size2), rs4 = rng(block_size2);
    op1 = op_desc{block_index1, rs1, true, data.linindex[std::make_pair(block_index1, rs1)]};
    op2 = op_desc{block_index1, rs2, false, data.linindex[std::make_pair(block_index1, rs2)]};
    op3 = op_desc{block_index2, rs3, true, data.linindex[std::make_pair(block_index2, rs3)]};
    op4 = op_desc{block_index2, rs4, false, data.linindex[std::make_pair(block_index2, rs4)]};

    // X = tau1 anywhere, x = tau2 near X, Y = tau3 anywhere (same block: far from X), y = tau4 near Y
    tau1 = data.tau_seg.get_random_pt(rng);
    tau2 = win.draw(tau1, rng);
    tau3 = (same_block ? dist.draw(tau1, rng) : data.tau_seg.get_random_pt(rng));
    tau4 = win.draw(tau3, rng);
    if ((tau1 == tau3) or (tau2 == tau4)) return 0; // trying to insert two operators at exactly the same time

    try {
      data.imp_trace.try_insert(tau1, op1);
      data.imp_trace.try_insert(tau2, op2);
      data.imp_trace.try_insert(tau3, op3);
      data.imp_trace.try_insert(tau4, op4);
    } catch (rbt_insert_error const &) {
      data.imp_trace.cancel_insert();
      return 0;
    }

    // Computation of det ratio
    auto &det1     = data.dets[block_index1];
    auto &det2     = data.dets[block_index2];
    long det1_size = det1.size();
    long det2_size = det2.size();
    auto get_x1    = [&det1](long i) -> time_pt const & { return det1.get_x(i).first; };
    auto get_y1    = [&det1](long i) -> time_pt const & { return det1.get_y(i).first; };
    auto get_x2    = [&det2](long i) -> time_pt const & { return det2.get_x(i).first; };
    auto get_y2    = [&det2](long i) -> time_pt const & { return det2.get_y(i).first; };
    det_scalar_t det_ratio;

    // Proposal ratio, with the counts of the configuration after the insertion
    mc_weight_t t_ratio;
    if (same_block) {
      // arcs of the near and far points of X = tau1 and Y = tau3, computed once
      window::arc_t const near_X[2] = {win.arc(tau1), win.arc(tau3)};
      window::arc_t const far_1 = dist.arc(tau1), far_3 = dist.arc(tau3);
      long const c_near[2]   = {window::arc_set(get_y1, det1_size, near_X[0]).count + near_X[0].contains(tau2) + near_X[0].contains(tau4),
                                window::arc_set(get_y1, det1_size, near_X[1]).count + near_X[1].contains(tau2) + near_X[1].contains(tau4)};
      long const cdag_far[2] = {window::arc_set(get_x1, det1_size, far_1).count + far_1.contains(tau3),
                                window::arc_set(get_x1, det1_size, far_3).count + far_3.contains(tau1)};
      time_pt const x[2]     = {tau2, tau4};
      t_ratio = window::same_block_proposal_ratio(det1_size + 2, x, near_X, c_near, cdag_far, far_1.contains(tau3), win, dist, block_size1);
      if (t_ratio == 0) TRIQS_RUNTIME_ERROR << "(insert window double) the proposed path is outside the windows";
    } else {
      long m_x = window::count_near(get_y1, det1_size, tau1, win) + 1;
      long m_y = window::count_near(get_y2, det2_size, tau3, win) + 1;
      double s = double(block_size1) * block_size2;
      t_ratio  = s * s * beta * beta * win.length * win.length / (double(det1_size + 1) * double(m_x) * double(det2_size + 1) * double(m_y));
    }

    // Position for insertion in the determinant (it stores the operators in decreasing time order)
    long num_c_dag1 = window::first_below(get_x1, det1_size, tau1);
    long num_c1     = window::first_below(get_y1, det1_size, tau2);
    long num_c_dag2 = window::first_below(get_x2, det2_size, tau3);
    long num_c2     = window::first_below(get_y2, det2_size, tau4);

    // Insert in the det. Returns the ratio of dets (Cf det_manip doc).
    if (same_block) {
      // The determinant positions that need to be passed to det_manip are those in the *final* det of size N+2.
      // Shift the operator at the smaller time one step further in the determinant to account for the larger operator.
      if (tau1 < tau3)
        num_c_dag1++;
      else
        num_c_dag2++;
      if (tau2 < tau4)
        num_c1++;
      else
        num_c2++;
      det_ratio = det1.try_insert2(num_c_dag1, num_c_dag2, num_c1, num_c2, {tau1, op1.inner_index}, {tau3, op3.inner_index}, {tau2, op2.inner_index},
                                   {tau4, op4.inner_index});
    } else {
      auto det_ratio1 = det1.try_insert(num_c_dag1, num_c1, {tau1, op1.inner_index}, {tau2, op2.inner_index});
      auto det_ratio2 = det2.try_insert(num_c_dag2, num_c2, {tau3, op3.inner_index}, {tau4, op4.inner_index});
      det_ratio       = det_ratio1 * det_ratio2;
    }

    // For quick abandon
    double random_number = rng.preview();
    if (random_number == 0.0) return 0;
    double p_yee = std::abs(t_ratio * det_ratio / data.atomic_weight);

    // computation of the new atomic_weight after insertion
    std::tie(new_atomic_weight, new_atomic_reweighting) = data.imp_trace.compute(p_yee, random_number);
    if (new_atomic_weight == 0.0) return 0;
    auto atomic_weight_ratio = new_atomic_weight / data.atomic_weight;
    if (!isfinite(atomic_weight_ratio))
      TRIQS_RUNTIME_ERROR << "(insert window double) atomic_weight_ratio not finite " << new_atomic_weight << " " << data.atomic_weight
                          << " in config " << config.get_id();

    mc_weight_t p = atomic_weight_ratio * det_ratio;
    if (!isfinite(p * t_ratio))
      TRIQS_RUNTIME_ERROR << "(insert window double) p * t_ratio not finite p : " << p << " t_ratio :  " << t_ratio << " in config "
                          << config.get_id();
    return p * t_ratio;
  }

  mc_weight_t move_insert_window_double::accept() {

    // insert in the tree
    data.imp_trace.confirm_insert();

    // insert in the configuration
    config.insert(tau1, op1);
    config.insert(tau2, op2);
    config.insert(tau3, op3);
    config.insert(tau4, op4);
    config.finalize();

    // insert in the determinant
    data.dets[block_index1].complete_operation();
    if (block_index1 != block_index2) data.dets[block_index2].complete_operation();
    data.update_sign();
    data.atomic_weight      = new_atomic_weight;
    data.atomic_reweighting = new_atomic_reweighting;
    return data.current_sign / data.old_sign;
  }

  void move_insert_window_double::reject() {
    config.finalize();
    data.imp_trace.cancel_insert();
    data.dets[block_index1].reject_last_try();
    if (block_index1 != block_index2) data.dets[block_index2].reject_last_try();
  }
} // namespace triqs_cthyb
