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
grid(grid_size,0)
{
    frequency = sample.frequency;
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

    grid = Grid<value_t>(grid_size,R);
    std::vector<arma::Mat<value_t>> grid_time_delays=grid.get_grid_time_delay(coords);
    std::vector<arma::Mat<value_t>> grid_phis = grid.get_phis_for_points(coords);

    //std::vector<arma::cube<std::complex<value_t>>> grid_phase(N);
    grid_phase = (std::complex<value_t>*)malloc(N*grid_size*grid_size*frequency.n_elem*sizeof(std::complex<value_t>));
    value_t * grid_phi2 = (value_t*)malloc(N*grid_size*grid_size*frequency.n_elem*sizeof(value_t));

    value_t wmix = array.wmix;

    for(int j=0; j<grid_size; ++j) {
        for(int k=0; k<grid_size; ++k) {
            for(int l=0; l<frequency.n_elem; ++l) {
                for(int i=0; i<N; ++i) {
                    //std::cout << j << " " << k << " " << l << " " << i << std::endl;
                    value_t phi = grid_time_delays[i](k,j)*(2*M_PI*frequency(l)+wmix);
                    phi += grid_phis[i](k,j);
                    grid_phase[((grid_size*j+k)*frequency.n_elem+l)*N+i] = std::complex<value_t>(cos(phi), sin(phi));
                    grid_phi2[((grid_size*j+k)*frequency.n_elem+l)*N+i] = phi;
                }
            }
        }
    }

    //~ for(int j=0; j<grid_size; ++j) {
        //~ int l=437;
        //~ int i=0;
        //~ for(int k=0; k<grid_size; ++k) {
            //~ std::cout << std::arg(grid_phase[((grid_size*j+k)*frequency.n_elem+l)*N+i]) << " " ;
            //~ //std::cout << grid_time_delays[0](k,j) << " " ;
        //~ }
        //~ std::cout << std::endl;
    //~ }

    //~ for(int j=0; j<grid_size; ++j) {
        //~ int l=437;
        //~ int i=0;
        //~ for(int k=0; k<grid_size; ++k) {
            //~ std::cout << grid_phi2[((grid_size*j+k)*frequency.n_elem+l)*N+i] << " " ;
            //~ //std::cout << grid_time_delays[25](i,j) << " " ;
        //~ }
        //~ std::cout << std::endl;
    //~ }

    //~ for(int i=0; i<grid_size; ++i) {
        //~ int l=437;
        //~ int k=0;
        //~ for(int j=0; j<grid_size; ++j) {
            //~ std::cout << grid_phis[0](i,j) << " " ;
        //~ }
        //~ std::cout << std::endl;
    //~ }

}

template <typename value_t>
void Reconstruction<value_t>::run(std::vector<Sample<value_t>> samples)
{
    arma::Cube<value_t> reconstructed(frequency.n_elem, grid_size, grid_size);

    std::cout << "pixels: " << grid_size*grid_size << " f_bins: " << frequency.n_elem << std::endl;


    TIMERSTART(REC)
    for(int j=0; j<grid_size; ++j) {
        for(int k=0; k<grid_size; ++k) {
            for(int l=0; l<frequency.n_elem; ++l) {
                std::complex<value_t> accum(0);
                for(int i=0; i<N; ++i) {
                   // std::cout << j << " " << k << " " << l << " " << i << std::endl;
                    //if(l==437)
                    //    std::cout << grid_phase[grid_size*j+frequency.n_elem*k+N*l+i] << std::endl;
                    accum += samples[i].frequency_data(l)*grid_phase[((grid_size*j+k)*frequency.n_elem+l)*N+i];
                    //if(l==437 && ((j==0 && k==0)||(j==70 && k==49)))
                     //   std::cout << "a " << samples[i].frequency_data(l) << " b " << grid_phase[grid_size*j+frequency.n_elem*k+N*l+i] << " a*b "  << samples[i].frequency_data(l)*grid_phase[grid_size*j+frequency.n_elem*k+N*l+i] << std::endl;
                }
                reconstructed(l, k, j) = std::abs(accum);
                //if(l==437 && j==0 && k==0)
                //    std::cout << accum << std::endl;
            }
        }
    }
    TIMERSTOP(REC)

    int channel_max=arma::index_max(samples[0].frequency_data);

   // std::cout << "channel " << channel_max << std::endl;

    for(int i=0; i<grid_size; ++i) {
        value_t y = grid.coords(i);
        //std::cout << y << std::endl;
        for(int j=0; j<grid_size; ++j) {
            value_t x = grid.coords(j);
            value_t d = x*x+y*y;
            //std::cout << d << std::endl;
            if(d>R*R) {
                img(j,i)=-1;
            } else {
                img(j,i) = reconstructed(channel_max, j, i);
            }
        }
    }
}
