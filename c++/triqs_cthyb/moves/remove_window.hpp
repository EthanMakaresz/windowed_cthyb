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
#include <triqs/mc_tools.hpp>
#include "../qmc_data.hpp"
#include "./window_util.hpp"

namespace triqs_cthyb {

  // Removal of a C^dagger, C pair, the reverse of move_insert_window: X uniformly among the n
  // C^dagger of the block, x uniformly among the m C near X (n, m in the current configuration),
  //     R = W(C)/W(C')  n m / (s^2 beta L_w),
  // C being the configuration after the removal.
  class move_remove_window {

    qmc_data &data;
    configuration &config;
    mc_tools::random_generator &rng;
    int block_index, block_size;
    window::window_t win;
    h_scalar_t new_atomic_weight, new_atomic_reweighting;
    time_pt tau1, tau2;

    public:
    move_remove_window(int block_index, int block_size, double w, qmc_data &data, mc_tools::random_generator &rng);

    mc_weight_t attempt();
    mc_weight_t accept();
    void reject();
  };
} // namespace triqs_cthyb
