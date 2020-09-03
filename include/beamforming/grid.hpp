/*
 * grid.hpp
 *
 * Copyright 2020 Florian Thomas <>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 *
 *
 */

#pragma once

#include <armadillo>

template <typename value_t>
class Grid {

    public:

        arma::Col<value_t> coords;

        value_t R;

        //Grid(int grid_size, value_t radius);
        Grid(int grid_size);
        void define_grid(value_t radius);

        std::vector<arma::Mat<value_t>> get_dists_to_points(arma::Mat<value_t> points);
        std::vector<arma::Mat<value_t>> get_grid_time_delay(arma::Mat<value_t> points);
        std::vector<arma::Mat<value_t>> get_phis_for_points(arma::Mat<value_t> points);

    private:
        int grid_size;


};
