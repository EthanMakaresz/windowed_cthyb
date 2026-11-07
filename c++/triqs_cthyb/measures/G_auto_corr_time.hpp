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

#include "../qmc_data.hpp"

#include <triqs/stat/log_binning.hpp>

#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace triqs_cthyb {

  /// Integrated autocorrelation time (in samples) of a log-binning accumulator, by the rule of measure_auto_corr_time:
  /// tau at the deepest bin level with >= 64 effective samples, flagged as a lower bound (converged = false) if it
  /// exceeds its finer neighbour (half the bin size) by more than 2 errors. Returns (tau, converged).
  inline std::pair<double, bool> log_binning_tau(triqs::stat::log_binning<dcomplex> const &acc, mpi::communicator const &comm) {
    constexpr int min_samples = 64;
    constexpr double n_sigma  = 2.0;
    auto d_tau                = [](double tau, long M) { return (tau + 0.5) * std::sqrt(2.0 / static_cast<double>(M - 1)); };

    auto [mean, errs, taus, effs] = acc.mean_errors_and_taus(comm, min_samples);
    if (taus.empty()) return {0.0, false};           // too short to estimate tau at all
    if (!std::isfinite(taus.back())) return {0.0, true}; // zero variance: tau ~ 0
    double const tau_a = taus.back();
    auto const n       = taus.size();
    if (n < 2 || !std::isfinite(taus[n - 2])) return {tau_a, false};
    return {tau_a, !(tau_a - taus[n - 2] > n_sigma * d_tau(tau_a, effs[n - 1]))};
  }

  /// Autocorrelation time of G(tau) averaged over a window around beta/2, for each block and diagonal element.
  ///
  /// One sample per cycle (tau in units of MC cycles): the instantaneous estimator of measure_G_tau restricted to the
  /// pairs whose time difference lies within the window, -sum_pairs s M / (beta * window width), the same quantity
  /// whose run average divided by the average sign is G(tau) averaged over the window. Each sample goes into its own
  /// log-binning accumulator, and tau is obtained as for auto_corr_time (log_binning_tau). The measure only reads the
  /// configuration and draws no random numbers, so it leaves the Markov chain unchanged.
  struct measure_G_auto_corr_time {

    /// window: the width of the window as a fraction of beta, centred on beta/2.
    measure_G_auto_corr_time(qmc_data const &_data, gf_struct_t const &gf_struct, double window,
                             std::map<std::string, std::vector<double>> &_tau, std::map<std::string, std::vector<int>> &_converged,
                             std::map<std::string, std::vector<double>> &_mean)
       : data(_data), tau(_tau), converged(_converged), mean(_mean) {
      double const beta = data.config.beta();
      half_beta         = 0.5 * beta;
      half_width        = 0.5 * window * beta;
      norm              = -1.0 / (beta * window * beta);
      long n            = 0;
      for (auto const &[name, size] : gf_struct) {
        names.push_back(name);
        sizes.push_back(size);
        offsets.push_back(n);
        n += size;
      }
      sample.assign(n, dcomplex{0.0});
      accs.assign(n, triqs::stat::log_binning<dcomplex>{dcomplex{0.0}, -1});
    }

    void accumulate(mc_weight_t s) {
      s *= data.atomic_reweighting;
      sum_sign += dcomplex(s);
      ++n_samples;
      std::fill(sample.begin(), sample.end(), dcomplex{0.0});
      for (long b = 0; b < long(names.size()); ++b) {
        foreach (data.dets[b], [this, s, b](op_t const &x, op_t const &y, det_scalar_t M) {
          if (x.second != y.second) return;
          double const dtau = double(y.first - x.first);
          if (std::abs(dtau - half_beta) > half_width) return;
          sample[offsets[b] + x.second] += dcomplex((y.first >= x.first ? s : -s) * M);
        })
          ;
      }
      for (size_t k = 0; k < accs.size(); ++k) accs[k] << sample[k] * norm;
    }

    void collect_results(mpi::communicator const &comm) {
      dcomplex const total_sign = mpi::all_reduce(sum_sign, comm);
      long const total_n        = mpi::all_reduce(n_samples, comm);
      dcomplex const mean_sign  = total_sign / double(std::max(total_n, 1L));
      tau.clear();
      converged.clear();
      mean.clear();
      for (size_t b = 0; b < names.size(); ++b) {
        for (long a = 0; a < sizes[b]; ++a) {
          auto const &acc       = accs[offsets[b] + a];
          auto [t, conv]        = log_binning_tau(acc, comm);
          auto [m, e, ts, effs] = acc.mean_errors_and_taus(comm, 2);
          tau[names[b]].push_back(t);
          converged[names[b]].push_back(conv ? 1 : 0);
          mean[names[b]].push_back(std::real(m / mean_sign));
        }
      }
    }

    private:
    qmc_data const &data;
    std::map<std::string, std::vector<double>> &tau;
    std::map<std::string, std::vector<int>> &converged;
    std::map<std::string, std::vector<double>> &mean;
    double half_beta = 0, half_width = 0, norm = 0;
    std::vector<std::string> names;
    std::vector<long> sizes, offsets;
    std::vector<dcomplex> sample;
    std::vector<triqs::stat::log_binning<dcomplex>> accs;
    dcomplex sum_sign = 0;
    long n_samples    = 0;
  };

} // namespace triqs_cthyb
