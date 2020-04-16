/* Copyright 2017-2018 PaGMO development team

This file is part of the PaGMO library.

The PaGMO library is free software; you can redistribute it and/or modify
it under the terms of either:

  * the GNU Lesser General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your
    option) any later version.

or

  * the GNU General Public License as published by the Free Software
    Foundation; either version 3 of the License, or (at your option) any
    later version.

or both in parallel, as here.

The PaGMO library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received copies of the GNU General Public License and the
GNU Lesser General Public License along with the PaGMO library.  If not,
see https://www.gnu.org/licenses/. */

#ifndef PAGMO_ALGORITHMS_TS_HPP
#define PAGMO_ALGORITHMS_TS_HPP

#include <vector>
#include <math.h>
#include <random>

#include <pagmo/rng.hpp>

namespace pagmo
{

std::vector<std::vector<double>> generate_sigma_points(int dim,
                                                       double factor = 1) {
  std::vector<std::vector<double>> result;
  result.push_back(std::vector<double>(dim, 0));
  for (int i = 0; i < dim; i++) {
    std::vector<double> pos(dim), neg(dim);
    pos[i] = std::sqrt(dim * factor);
    neg[i] = -std::sqrt(dim * factor);
    result.push_back(pos);
    result.push_back(neg);
  }
  return result;
}

std::vector<double> transform_gaussian(std::vector<double> gaussian,
                                       double interval = 1) {
  std::vector<double> result = gaussian;
  int dim = gaussian.size();

  if (interval == 1) {
    return result;
  }

  double r = 0;
  for (int d = 0; d < dim; d++) {
    r += gaussian[d] * gaussian[d];
  }
  r = std::sqrt(r);

  const double r_transformed = std::exp(-(r * r) * 0.5);
  const double new_r = std::sqrt(-2 * std::log(interval * r_transformed));

  for (int d = 0; d < dim; d++) {
    result[d] *= new_r / r;
  }

  return result;
}

std::vector<std::vector<double>> generate_tail_scented(int dim, int num_total,
                                                       detail::random_engine_type && g,
                                                       double ratio_tail = 0.3,
                                                       bool verbose = false) {
  std::vector<std::vector<double>> result;
  if (ratio_tail > 1 || ratio_tail < 0) {
    throw std::invalid_argument("ratio_tail must be in [0,1], but is" +
                                std::to_string(ratio_tail));
  }

  if (num_total == 0) {
    return result;
  }
  result.reserve(num_total);

  int num_tail = int(num_total * ratio_tail);
  int num_sigma = num_total - num_tail;
  double sigma_scaling = 1;

  std::vector<std::vector<double>> tail_points(num_tail);

  if (ratio_tail > 0) {
    double min_r = std::sqrt(-2 * std::log(ratio_tail));
    double var_tail = -std::log(ratio_tail) + 1;
    double var_total = 1;
    double var_sigma =
        (var_total * (num_total - 1) - var_tail * (std::max(num_tail - 1, 0))) /
        (num_sigma - 1);
    sigma_scaling = var_sigma;

    // generate gaussian numbers
    std::normal_distribution<double> normally_distributed_number(0., 1.);
    for (int i = 0; i < num_tail; i++) {
      tail_points[i].resize(dim);
      for (int d = 0; d < dim; d++) {
        tail_points[i][d] = normally_distributed_number(g);
      }
      tail_points[i] = transform_gaussian(tail_points[i], ratio_tail);
    }
  }

  // generate sigma points
  std::vector<std::vector<double>> sigma_points =
      generate_sigma_points(dim, sigma_scaling);
  std::shuffle(sigma_points.begin(), sigma_points.end(), g);

  if (num_sigma < sigma_points.size()) {
    // truncate
    sigma_points.erase(sigma_points.begin() + num_sigma, sigma_points.end());
  }

  if (num_sigma > sigma_points.size()) {
    // extend
    throw std::invalid_argument("Asked for " + std::to_string(num_sigma) +
                           " sigma points, but only " +
                           std::to_string(sigma_points.size()) + " are available.");
  }

  if (num_sigma != sigma_points.size()) {
    throw std::logic_error("Expected" + std::to_string(num_sigma) +
                           " sigma points, but found " +
                           std::to_string(sigma_points.size()));
  }

  std::merge(sigma_points.begin(), sigma_points.end(), tail_points.begin(),
             tail_points.end(), std::back_inserter(result));

  if (num_total != result.size()) {
    throw std::logic_error("Expected" + std::to_string(num_total) +
                           " points, but found " +
                           std::to_string(result.size()));
  }

  return result;
}
}

#endif