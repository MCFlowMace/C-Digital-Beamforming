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

/*

#include "reconstruction.hpp"
#include "utility_macros.hpp"
#include "reconstruction_cpu.hpp"

#include <stdexcept>
#include <string>
#include <limits>

template <typename value_t>
Reconstruction_CPU<value_t>::Reconstruction_CPU(int grid_size, int buffer_size,
						arma::Col<value_t> frequency,
                        const Antenna_Array<value_t>& array):
grid_size(grid_size),
grid(grid_size),
frequency(frequency),
reconstructed(frequency.n_elem*buffer_size, grid_size, grid_size)
{

    grid_phase=nullptr;
    set_antenna_array(array);

}

template <typename value_t>
Reconstruction_CPU<value_t>::~Reconstruction_CPU()
{
    free_grid_phase();
}

template <typename value_t>
void Reconstruction_CPU<value_t>::set_antenna_array(const Antenna_Array<value_t>& array)
{

    //this->array=array;
    N = array.N;
    R = array.R;
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

    //TODO needs a "free" somewhere + class needs better design with respect to this calloc
    std::complex<value_t>* grid_phase_local = (std::complex<value_t>*)calloc(N*grid_size*grid_size*frequency.n_elem,sizeof(std::complex<value_t>));

	TIMERSTART(CALC_PHASE)
    calc_phase(grid_time_delays, grid_phis, grid_phase_local);
	TIMERSTOP(CALC_PHASE)
	
    set_grid_phase(&grid_phase_local);

    //data is on GPU now
    if(grid_phase_local) {
        free(grid_phase_local);
        grid_phase_local=nullptr;
    }

}


template <typename value_t>
void Reconstruction_CPU<value_t>::calc_phase(
                        const std::vector<arma::Mat<value_t>>& grid_time_delays,
                        const std::vector<arma::Mat<value_t>>& grid_phis,
                        std::complex<value_t>* const grid_phase)
{

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
void Reconstruction_CPU<value_t>::set_grid_phase(std::complex<value_t>** grid_phase)
{
    this->grid_phase= *grid_phase;
    *grid_phase=nullptr;
}

template <typename value_t>
void Reconstruction_CPU<value_t>::run(const std::vector<Data_Packet<value_t>>& samples)
{
    //arma::Cube<value_t> reconstructed(frequency.n_elem, grid_size, grid_size);

    reconstructed.ones();
    reconstructed*=-1;

    TIMERSTART(REC)

#ifdef PARALLEL
    #pragma omp parallel
    {
        std::cerr << "threads: " << omp_get_num_threads() << std::endl;
        #pragma omp for
#endif
        for(int j=0; j<grid_size; ++j) {
            value_t y = grid.coords(j);
            for(int k=0; k<grid_size; ++k) {
                value_t x = grid.coords(k);
                value_t d = x*x+y*y;

                if(d<=R*R) {
                    for(int l=0; l<frequency.n_elem; ++l) {
                        std::complex<value_t> accum(0);
                        for(int i=0; i<N; ++i) {
                            accum += samples[i].frequency_data(l)*grid_phase[((grid_size*j+k)*frequency.n_elem+l)*N+i];
                        }
                        reconstructed(l, k, j) = std::abs(accum);
                    }
                }
            }
        }
#ifdef PARALLEL
    }
#endif

    TIMERSTOP(REC)

}

template <typename value_t>
void Reconstruction_CPU<value_t>::free_grid_phase()
{
    if(grid_phase) {
        free(grid_phase);
        grid_phase=nullptr;
    }
}

template <typename value_t>
unsigned int Reconstruction_CPU<value_t>::get_max_bin()
{

    value_t max_val = std::numeric_limits<value_t>::min();
    unsigned int index;
    for(unsigned int i=0; i<frequency.n_elem; ++i) {
        for(unsigned int j=0; j<grid_size; ++j) {
            for(unsigned int l=0; l<grid_size; ++l) {

                if(reconstructed(i,j,l)>max_val) {
                    max_val=reconstructed(i,j,l);
                    index = i;
                }
            }
        }
    }

    std::cerr << "Max frequency: " << frequency[index] << std::endl;

    return index;
}

template <typename value_t>
arma::Mat<value_t> Reconstruction_CPU<value_t>::get_img(unsigned int bin)
{
    if(bin >= frequency.n_elem) {
        std::string err = "Bin " + std::to_string(bin) + " is not a valid frequency bin!";
        throw std::out_of_range(err);
    }

    arma::Mat<value_t> img(grid_size, grid_size);

    for(int i=0; i<grid_size; ++i) {
        //value_t y = grid.coords(i);

        for(int j=0; j<grid_size; ++j) {
            //value_t x = grid.coords(j);
            //~ value_t d = x*x+y*y;

            //~ if(d>R*R) {
                //~ img(j,i)=-1;
            //~ } else {
                //~ img(j,i) = reconstructed(bin, j, i);
            //~ }
            img(j,i) = reconstructed(bin, j, i);
        }
    }

    return img;
}

template <typename value_t>
value_t Reconstruction_CPU<value_t>::get_mean_val(unsigned int bin)
{
    arma::Mat<value_t> img = get_img(bin);

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
value_t Reconstruction_CPU<value_t>::get_max_val(unsigned int bin)
{
    arma::Mat<value_t> img = get_img(bin);

    value_t max_val = std::numeric_limits<value_t>::min();
    unsigned int x;
    unsigned int y;

    for(int i=0; i<grid_size; ++i) {
        for(int j=0; j<grid_size; ++j) {
            if(img(j, i) ==-1)
                continue;

            if(img(j, i) > max_val) {
                max_val=img(j,i);
                x=i;
                y=j;
            }
        }
    }

    std::cerr << "max at: " << x << " " << y << std::endl;
    return max_val;
}

DEFINE_TEMPLATES(Reconstruction_CPU)

*/
