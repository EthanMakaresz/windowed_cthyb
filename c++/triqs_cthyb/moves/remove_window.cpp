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
#include "./remove_window.hpp"

namespace triqs_cthyb {

  move_remove_window::move_remove_window(int block_index, int block_size, double w, qmc_data &data, mc_tools::random_generator &rng)
     : data(data), config(data.config), rng(rng), block_index(block_index), block_size(block_size), win(w, data.config.beta()) {}

  mc_weight_t move_remove_window::attempt() {

    auto &det     = data.dets[block_index];
    long det_size = det.size();
    if (det_size == 0) return 0; // nothing to remove
    auto get_y = [&det](long i) -> time_pt const & { return det.get_y(i).first; };

    // A C^dagger at random, then a C at random among those near it
    long num_c_dag = rng(det_size);
    auto near = window::near_set(get_y, det_size, det.get_x(num_c_dag).first, win);
    if (near.count == 0) return 0;
    long num_c = near.index(rng(near.count));

    // now mark 2 nodes for deletion
    tau1 = data.imp_trace.try_delete(num_c, block_index, false);
    tau2 = data.imp_trace.try_delete(num_c_dag, block_index, true);

    auto det_ratio = det.try_remove(num_c_dag, num_c);

    // proposition probability of the insertion from the smaller configuration (n, m: before the removal)
    mc_weight_t t_ratio = double(block_size) * block_size * config.beta() * win.length / (double(det_size) * double(near.count));

    // For quick abandon
    double random_number = rng.preview();
    if (random_number == 0.0) return 0;
    double p_yee = std::abs(det_ratio / t_ratio / data.atomic_weight);

    // recompute the atomic_weight
    std::tie(new_atomic_weight, new_atomic_reweighting) = data.imp_trace.compute(p_yee, random_number);
    if (new_atomic_weight == 0.0) return 0;
    auto atomic_weight_ratio = new_atomic_weight / data.atomic_weight;
    if (!isfinite(atomic_weight_ratio))
      TRIQS_RUNTIME_ERROR << "(remove window) atomic_weight_ratio not finite " << new_atomic_weight << " " << data.atomic_weight << " in config "
                          << config.get_id();

    mc_weight_t p = atomic_weight_ratio * det_ratio;
    if (!isfinite(p / t_ratio))
      TRIQS_RUNTIME_ERROR << "(remove window) p / t_ratio not finite p : " << p << " t_ratio : " << t_ratio << " in config " << config.get_id();
    return p / t_ratio;
  }

  mc_weight_t move_remove_window::accept() {

    // remove from the tree
    data.imp_trace.confirm_delete();

    // remove from the configuration
    config.erase(tau1);
    config.erase(tau2);
    config.finalize();

    // remove from the determinants
    data.dets[block_index].complete_operation();
    data.update_sign();
    data.atomic_weight      = new_atomic_weight;
    data.atomic_reweighting = new_atomic_reweighting;
    return data.current_sign / data.old_sign;
  }

  void move_remove_window::reject() {
    config.finalize();
    data.imp_trace.cancel_delete();
    data.dets[block_index].reject_last_try();
  }
} // namespace triqs_cthyb
