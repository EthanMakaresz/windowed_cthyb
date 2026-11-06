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
#include "./insert_window.hpp"

namespace triqs_cthyb {

  move_insert_window::move_insert_window(int block_index, int block_size, double w, qmc_data &data, mc_tools::random_generator &rng)
     : data(data), config(data.config), rng(rng), block_index(block_index), block_size(block_size), win(w, data.config.beta()) {}

  mc_weight_t move_insert_window::attempt() {

    // Pick up the value of alpha and choose the operators
    auto rs1 = rng(block_size), rs2 = rng(block_size);
    op1 = op_desc{block_index, rs1, true, data.linindex[std::make_pair(block_index, rs1)]};
    op2 = op_desc{block_index, rs2, false, data.linindex[std::make_pair(block_index, rs2)]};

    // The C^dagger anywhere, the C in the window around it
    tau1 = data.tau_seg.get_random_pt(rng);
    tau2 = win.draw(tau1, rng);

    try {
      data.imp_trace.try_insert(tau1, op1);
      data.imp_trace.try_insert(tau2, op2);
    } catch (rbt_insert_error const &) {
      data.imp_trace.cancel_insert();
      return 0;
    }

    // Computation of det ratio
    auto &det    = data.dets[block_index];
    long det_size = det.size();
    auto get_x   = [&det](long i) -> time_pt const & { return det.get_x(i).first; };
    auto get_y   = [&det](long i) -> time_pt const & { return det.get_y(i).first; };

    // Position for insertion in the determinant (it stores the operators in decreasing time order)
    long num_c_dag = window::first_below(get_x, det_size, tau1);
    long num_c     = window::first_below(get_y, det_size, tau2);

    // m': the C near the new C^dagger after the insertion, the new C included
    long m_new = window::count_near(get_y, det_size, tau1, win) + 1;

    // Insert in the det. Returns the ratio of dets (Cf det_manip doc).
    auto det_ratio = det.try_insert(num_c_dag, num_c, {tau1, op1.inner_index}, {tau2, op2.inner_index});

    // proposition probability
    mc_weight_t t_ratio = double(block_size) * block_size * config.beta() * win.length / (double(det_size + 1) * double(m_new));

    // For quick abandon
    double random_number = rng.preview();
    if (random_number == 0.0) return 0;
    double p_yee = std::abs(t_ratio * det_ratio / data.atomic_weight);

    // computation of the new trace after insertion
    std::tie(new_atomic_weight, new_atomic_reweighting) = data.imp_trace.compute(p_yee, random_number);
    if (new_atomic_weight == 0.0) return 0;
    auto atomic_weight_ratio = new_atomic_weight / data.atomic_weight;
    if (!isfinite(atomic_weight_ratio))
      TRIQS_RUNTIME_ERROR << "(insert window) trace_ratio not finite " << new_atomic_weight << " " << data.atomic_weight << " in config "
                          << config.get_id();

    mc_weight_t p = atomic_weight_ratio * det_ratio;
    if (!isfinite(p * t_ratio))
      TRIQS_RUNTIME_ERROR << "(insert window) p * t_ratio not finite p : " << p << " t_ratio : " << t_ratio << " in config " << config.get_id();
    return p * t_ratio;
  }

  mc_weight_t move_insert_window::accept() {

    // insert in the tree
    data.imp_trace.confirm_insert();

    // insert in the configuration
    config.insert(tau1, op1);
    config.insert(tau2, op2);
    config.finalize();

    // insert in the determinant
    data.dets[block_index].complete_operation();
    data.update_sign();
    data.atomic_weight      = new_atomic_weight;
    data.atomic_reweighting = new_atomic_reweighting;
    return data.current_sign / data.old_sign;
  }

  void move_insert_window::reject() {
    config.finalize();
    data.imp_trace.cancel_insert();
    data.dets[block_index].reject_last_try();
  }
} // namespace triqs_cthyb
