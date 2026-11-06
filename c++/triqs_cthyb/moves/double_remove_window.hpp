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

  // Removal of two C^dagger, C pairs, the reverse of move_insert_window_double (windowed mode):
  // X uniformly among the C^dagger of block 1, x uniformly among the C of block 1 near X (window w4),
  // Y uniformly among the C^dagger of block 2 (same block: among those farther than d_min from X),
  // y uniformly among the C of block 2 near Y (same block: x excluded).
  class move_remove_window_double {

    qmc_data &data;
    configuration &config;
    mc_tools::random_generator &rng;
    int block_index1, block_index2, block_size1, block_size2;
    window::window_t win;
    window::distance_t dist;
    h_scalar_t new_atomic_weight, new_atomic_reweighting;
    time_pt tau1, tau2, tau3, tau4;

    public:
    move_remove_window_double(int block_index1, int block_index2, int block_size1, int block_size2, double w4, double d_min, qmc_data &data,
                              mc_tools::random_generator &rng);

    mc_weight_t attempt();
    mc_weight_t accept();
    void reject();
  };
} // namespace triqs_cthyb
