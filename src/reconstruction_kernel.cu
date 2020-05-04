/*
 * reconstruction_kernel.cu
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

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/complex.h>
#include <thrust/device_ptr.h>

#include "hpc_helpers.hpp"
#include "data_packet.hpp"
#include "reconstruction.hpp"


#ifdef USE_GPU

//not very nice solution
#include "../src/reconstruction.cpp"

template <typename value_t>
__global__ void reconstruction(thrust::complex<value_t>* samples,
                                thrust::complex<value_t>* grid_phase,
                                value_t* coords, value_t* rec, value_t R,
                                int bins, int grid_size, int N)
{
    int tid = blockDim.x*blockIdx.x+threadIdx.x;

    if(tid<grid_size*grid_size*bins) {

        int i = tid%grid_size;
        int j = (tid/grid_size)%grid_size;
        int k = tid/(grid_size*grid_size);

        value_t x = coords[j];
        value_t y = coords[i];
        value_t d = x*x+y*y;

        if(d<=R*R) {

            thrust::complex<value_t> accum(0);

            for(int l=0; l<N; ++l) {
                //accum += samples[l*bins+k]*grid_phase[((grid_size*i+j)*bins+k)*N+l]; //not coalesced
                accum += samples[l*bins+k]*grid_phase[((l*bins+k)*grid_size+j)*grid_size+i];
            }
            rec[(grid_size*i+j)*bins+k] = thrust::abs(accum);

        }

    }
}

template <typename value_t>
void Reconstruction<value_t>::calc_phase(
                        const std::vector<arma::Mat<value_t>>& grid_time_delays,
                        const std::vector<arma::Mat<value_t>>& grid_phis,
                        std::complex<value_t>* const grid_phase)
{
	
	int bins = frequency.n_elem;

    for(int i=0; i<N; ++i) {
        for(int l=0; l<bins; ++l) {
            for(int k=0; k<grid_size; ++k) {
                for(int j=0; j<grid_size; ++j) {
                    value_t phi = grid_time_delays[i](k,j)*(2*M_PI*frequency(l)+wmix);
                    phi += grid_phis[i](k,j);
                    //grid_phase[((grid_size*j+k)*frequency.n_elem+l)*N+i] = std::complex<value_t>(cos(phi), sin(phi));
					grid_phase[((i*bins+l)*grid_size+k)*grid_size+j] = std::complex<value_t>(cos(phi), sin(phi));
                }
            }
        }
    }

}

template <typename value_t>
void Reconstruction<value_t>::set_grid_phase(std::complex<value_t>** grid_phase)
{
    cudaMalloc(&(this->grid_phase),
        N*grid_size*grid_size*frequency.n_elem*sizeof(std::complex<value_t>)); CUERR

    cudaMemcpy(this->grid_phase, *grid_phase,
        N*grid_size*grid_size*frequency.n_elem*sizeof(std::complex<value_t>),
        H2D);                           CUERR
}

template <typename value_t>
void Reconstruction<value_t>::free_grid_phase()
{
    if(grid_phase) {
        cudaFree(grid_phase);   CUERR
        grid_phase=nullptr;
    }
}

template <typename value_t>
void Reconstruction<value_t>::run(const std::vector<Data_Packet<value_t>>& samples)
{


    TIMERSTART(REC)

    int bins = samples[0].frequency.n_elem;

    //reorder data
    thrust::host_vector<thrust::complex<value_t> > samples_H(N*bins, thrust::complex<value_t>(0,1));   CUERR

    TIMERSTART(COPY_CPU)
    for(int i=0; i<samples.size(); ++i) {
        thrust::copy(samples[i].frequency_data.begin(),
                        samples[i].frequency_data.end(),
                        samples_H.begin()+i*bins);                      CUERR
    }
    TIMERSTOP(COPY_CPU)

    //copy data to GPU
    TIMERSTART(COPY_GPU)
    thrust::device_vector<thrust::complex<value_t> > samples_D = samples_H;  CUERR
                                                            //(
                                                            //samples_H.begin(),
                                                            //samples_H.end()); CUERR

    thrust::device_vector<value_t> coords_D(grid.coords.begin(),
                                            grid.coords.end());         CUERR
    TIMERSTOP(COPY_GPU)

    thrust::device_vector<value_t> reconstructed_D(grid_size*grid_size*bins);
    thrust::fill(reconstructed_D.begin(), reconstructed_D.end(), value_t(-1));

    thrust::complex<value_t>* samples_dev = thrust::raw_pointer_cast(samples_D.data());
    value_t* rec_dev = thrust::raw_pointer_cast(reconstructed_D.data());
    value_t* coords_dev = thrust::raw_pointer_cast(coords_D.data());

    int threads = 512;
    int tasks=grid_size*grid_size*bins;
    int blocks = SDIV(tasks, threads);

    std::cerr << bins << " " << grid_size << " " << N << std::endl;

    TIMERSTART(KERNEL)
    reconstruction<<<blocks, threads>>>(samples_dev,
                                        (thrust::complex<value_t>*)grid_phase,
                                        coords_dev, rec_dev, grid.R,
                                        bins, grid_size, N);   CUERR
    TIMERSTOP(KERNEL)

    //~ thrust::host_vector<value_t> reconstructed_H(reconstructed_D);

    //~ int count=0;
    //~ for(auto i = reconstructed_H.begin(); i!=reconstructed_H.end(); i++) {

        //~ if(*i==0)
            //~ count++;
    //~ }
    //~ std::cout << "zeros: " << count << std::endl;

    //copy back result
    TIMERSTART(COPY_BACK)
    thrust::copy(reconstructed_D.begin(), reconstructed_D.end(),
                    reconstructed.begin());                             CUERR
	TIMERSTOP(COPY_BACK)

    TIMERSTOP(REC)
}

DEFINE_TEMPLATES(Reconstruction)

#endif
