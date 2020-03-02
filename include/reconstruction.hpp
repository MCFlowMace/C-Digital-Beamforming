/*
 * reconstruction.hpp
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
#include "antenna.hpp"
#include "grid.hpp"
#include "hpc_helpers.hpp"

template <typename value_t>
class Reconstruction {


    public:

        arma::Mat<value_t> img;
        arma::Col<value_t> frequency;

        //arma::field<arma::Cube<std::complex<value_t>>> grid_phase;

        Reconstruction(int grid_size, Sample<value_t> sample);

        void set_antenna_array(const Antenna_Array<value_t>& array);
        void run(std::vector<Sample<value_t>> samples);

        std::complex<value_t>* grid_phase;

    private:

        int grid_size;
        Grid<value_t> grid;
        int N;
        value_t R;
};

template <typename value_t>
Reconstruction<value_t>::Reconstruction(int grid_size, Sample<value_t> sample):
grid_size(grid_size),
img(grid_size,grid_size),
grid(grid_size)
{
    frequency = sample.frequency;

    grid_phase=nullptr;
}

template <typename value_t>
void Reconstruction<value_t>::set_antenna_array(const Antenna_Array<value_t>& array)
{

    //this->array=array;
    N = array.N;
    R = array.R;

    arma::Mat<value_t> coords(2,N);


    for(int i=0; i<N; ++i) {
        coords(0,i) = array.antennas[i].x;
        coords(1,i) = array.antennas[i].y;
    }

    grid.define_grid(R);


    std::vector<arma::Mat<value_t>> grid_time_delays=grid.get_grid_time_delay(coords);

    std::vector<arma::Mat<value_t>> grid_phis = grid.get_phis_for_points(coords);

    grid_phase = (std::complex<value_t>*)calloc(N*grid_size*grid_size*frequency.n_elem,sizeof(std::complex<value_t>));
    //grid_phase = (std::complex<value_t>*)calloc(30*100*100*499,sizeof(std::complex<value_t>));

    value_t wmix = array.wmix;

    for(int j=0; j<grid_size; ++j) {
        for(int k=0; k<grid_size; ++k) {
            for(int l=0; l<frequency.n_elem; ++l) {
                for(int i=0; i<N; ++i) {
                    value_t phi = grid_time_delays[i](k,j)*(2*M_PI*frequency(l)+wmix);
                    phi += grid_phis[i](k,j);
                    grid_phase[((grid_size*j+k)*frequency.n_elem+l)*N+i] = std::complex<value_t>(cos(phi), sin(phi));

                }
            }
        }
    }

}

template <typename value_t>
void Reconstruction<value_t>::run(std::vector<Sample<value_t>> samples)
{
    arma::Cube<value_t> reconstructed(frequency.n_elem, grid_size, grid_size);


    TIMERSTART(REC)

#ifdef PARALLEL
    #pragma omp parallel num_threads(4)
    {
        std::cout << "threads: " << omp_get_num_threads() << std::endl;
        #pragma omp for collapse(3)
#endif
        for(int j=0; j<grid_size; ++j) {
            for(int k=0; k<grid_size; ++k) {
                for(int l=0; l<frequency.n_elem; ++l) {
                    std::complex<value_t> accum(0);
                    for(int i=0; i<N; ++i) {
                        accum += samples[i].frequency_data(l)*grid_phase[((grid_size*j+k)*frequency.n_elem+l)*N+i];
                    }
                    reconstructed(l, k, j) = std::abs(accum);
                }
            }
        }
#ifdef PARALLEL
    }
#endif

    TIMERSTOP(REC)

    int channel_max=arma::index_max(samples[0].frequency_data);


    if(grid_phase)
        free(grid_phase);

    for(int i=0; i<grid_size; ++i) {
        value_t y = grid.coords(i);

        for(int j=0; j<grid_size; ++j) {
            value_t x = grid.coords(j);
            value_t d = x*x+y*y;

            if(d>R*R) {
                img(j,i)=-1;
            } else {
                img(j,i) = reconstructed(channel_max, j, i);
            }
        }
    }
}
