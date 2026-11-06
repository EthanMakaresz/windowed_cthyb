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

  // Insertion of a C^dagger, C pair with the C within the window w of the C^dagger (windowed mode).
  //
  // C^dagger X uniform on the circle (density 1/beta), C x uniform in the window of X (1/L_w,
  // L_w ~ 2w, see window_util.hpp), inner indices uniform (1/s^2). The reverse, move_remove_window,
  // picks X uniformly among the n' C^dagger of the block and x uniformly among the m' C near X:
  // 1/(n' m'), with n', m' counted after the insertion. One labelled path connects C and C' (the
  // pairing of X with x is a property of the proposal only), so
  //     R = W(C')/W(C)  s^2 beta L_w / (n' m').
  // With w = beta/2 this is the stock move_insert_c_cdag.
  class move_insert_window {

    qmc_data &data;
    configuration &config;
    mc_tools::random_generator &rng;
    int block_index, block_size;
    window::window_t win;
    h_scalar_t new_atomic_weight, new_atomic_reweighting;
    time_pt tau1, tau2;
    op_desc op1, op2;

    public:
    move_insert_window(int block_index, int block_size, double w, qmc_data &data, mc_tools::random_generator &rng);

    mc_weight_t attempt();
    mc_weight_t accept();
    void reject();
  };
} // namespace triqs_cthyb
