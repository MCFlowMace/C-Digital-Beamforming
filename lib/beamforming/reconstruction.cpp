/*
 * reconstruction.cpp
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
 
#include <cfloat>

#include "beamforming/reconstruction.hpp"
#include "utility/utility_macros.hpp"

#include <stdexcept>
#include <string>
#include <limits>

template <typename value_t>
Reconstruction<value_t>::Reconstruction(int grid_size, int n_packets,
					arma::Col<value_t> frequency,
				    const Antenna_Array<value_t>& array, 
				    bool weighted, value_t r_grid):
grid_size(grid_size),
grid(grid_size),
frequency(frequency),
n_packets(n_packets),
bins(frequency.n_elem),
R(r_grid),
weighted(weighted)
{
    //frequency = sample.frequency;
    //reconstructed(frequency.n_elem, grid_size, grid_size);
    set_antenna_array(array);
}

template <typename value_t>
Grid<value_t> Reconstruction<value_t>::get_grid() const {
	return this->grid;
}

template <typename value_t>
void Reconstruction<value_t>::set_antenna_array(const Antenna_Array<value_t>& array)
{

    //this->array=array;
    N = array.N;
    //R = array.R;
    wmix = array.wmix;

    arma::Mat<value_t> coords(2,N);


    for(int i=0; i<N; ++i) {
        coords(0,i) = array.antennas[i].x;
        coords(1,i) = array.antennas[i].y;
    }

    grid.define_grid(R);


    //std::vector<arma::Mat<value_t>> grid_time_delays=grid.get_grid_time_delay(coords);
    
    grid_time_delays=grid.get_grid_time_delay(coords);

    //std::vector<arma::Mat<value_t>> grid_phis = grid.get_phis_for_points(coords);
    
    grid_phis = grid.get_phis_for_points(coords);

}

template <typename value_t>
value_t Reconstruction<value_t>::get_mean_val(unsigned int bin, 
												unsigned int packet)
{
    arma::Mat<value_t> img = get_img(packet, bin);

    value_t mean {0};
    unsigned int count=0;
    for(int i=0; i<grid_size; ++i) {
        for(int j=0; j<grid_size; ++j) {
            if(img(j, i) ==-1)
                continue;

            count++;
            mean+=img(j,i);
        }
    }

    return mean/=count;
}

template <typename value_t>
value_t Reconstruction<value_t>::get_max_val(unsigned int bin,
												unsigned int packet)
{
    arma::Mat<value_t> img = get_img(packet, bin);

    value_t max_val = std::numeric_limits<value_t>::min();
    unsigned int x;
    unsigned int y;

    for(int i=0; i<grid_size; ++i) {
        for(int j=0; j<grid_size; ++j) {
			
			value_t val = img(j,i);
            if(val ==-1)
                continue;

            if(val>max_val && !std::isinf(val)) {
                max_val=val;
                x=i;
                y=j;
            }
        }
    }

    //std::cerr << "max at: " << this->grid.coords(x) << " " << this->grid.coords(y) << std::endl;
    return max_val;
}

template <typename value_t>
arma::Mat<value_t> Reconstruction<value_t>::get_max_vals()
{
	
	arma::Mat<value_t> max_vals(frequency.n_elem, n_packets);
	
	for(int j=0; j<n_packets; ++j) {
		
#if defined PARALLEL || defined USE_GPU
		#pragma omp parallel for schedule(static)
#endif
		for(int i=0; i<frequency.n_elem; ++i)
			max_vals(i,j) = this->get_max_val(i, j);
		
	}
	
	return max_vals;
}

DEFINE_TEMPLATES(Reconstruction)
