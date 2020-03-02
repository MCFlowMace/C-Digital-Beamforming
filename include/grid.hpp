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

        //Grid(int grid_size, value_t radius);
        Grid(int grid_size);
        void define_grid(value_t radius);

        std::vector<arma::Mat<value_t>> get_dists_to_points(arma::Mat<value_t> points);
        std::vector<arma::Mat<value_t>> get_grid_time_delay(arma::Mat<value_t> points);
        std::vector<arma::Mat<value_t>> get_phis_for_points(arma::Mat<value_t> points);

    private:
        int grid_size;


};

template <typename value_t>
Grid<value_t>::Grid(int grid_size):
grid_size(grid_size),
coords(grid_size)
{

}

template <typename value_t>
void Grid<value_t>::define_grid(value_t radius)
{

    for(int i=0; i<grid_size; ++i) {
        value_t val = (-1+(value_t)(2*i+1)/grid_size)*radius;
        coords(i) = val;
    }

}

template <typename value_t>
std::vector<arma::Mat<value_t>> Grid<value_t>::get_dists_to_points(arma::Mat<value_t> points)
{

    int N=points.n_cols;

    std::vector<arma::Mat<value_t>> grid_dists(N);

    for(int i=0; i<N; ++i) {
        arma::Mat<value_t> dists(grid_size, grid_size);
        for(int j=0; j<grid_size; ++j) {
            for(int k=0; k<grid_size; ++k) {
                //std::cout << i << " " << j << " " << k << std::endl;
                value_t dx = points(0,i)-coords(k);
                value_t dy = points(1,i)-coords(j);
                dists(k,j) = sqrt(dx*dx + dy*dy);
            }
        }
        grid_dists[i] = dists;
    }

    return grid_dists;
}

template <typename value_t>
std::vector<arma::Mat<value_t>> Grid<value_t>::get_grid_time_delay(arma::Mat<value_t> points)
{
    std::vector<arma::Mat<value_t>> grid_dists = this->get_dists_to_points(points);

    int N=points.n_cols;
    std::vector<arma::Mat<value_t>> time_delays(N);


    for(int i=0; i<N; ++i) {
        time_delays[i] = grid_dists[i]/(arma::datum::c_0*100.0f);
    }

    return time_delays;
}

template <typename value_t>
std::vector<arma::Mat<value_t>> Grid<value_t>::get_phis_for_points(arma::Mat<value_t> points)
{
    std::vector<arma::Mat<value_t>> phis(points.n_cols);

   // points.print("points");

    //#grid_center=points_ - self.coords
    //#return np.arctan2(-grid_center[:,:,:,1],-grid_center[:,:,:,0])+np.pi

    int N=points.n_cols;
    for(int i=0; i<N; ++i) {
        arma::Mat<value_t> grid_phi(grid_size, grid_size);
        for(int j=0; j<grid_size; ++j) {
            for(int k=0; k<grid_size; ++k) {
                value_t dx = points(0,i)-coords(k);
                value_t dy = points(1,i)-coords(j);
                grid_phi(k,j) = atan2(-dy, -dx) + M_PI;
            }
        }
        phis[i] = grid_phi;
    }

    return phis;
}
