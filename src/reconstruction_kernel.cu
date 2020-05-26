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

//speed of light in cm/s
#define SPEED_OF_LIGHT 29979245800 


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
                accum += samples[l*bins+k]*grid_phase[((l*bins+k)*grid_size+j)*grid_size+i];
            }
            rec[(grid_size*i+j)*bins+k] = thrust::abs(accum);

        }

    }
}

#define NANTENNAS 30
 /*
template <typename value_t>
__global__ void reconstruction_red(thrust::complex<value_t>* samples,
                                value_t* frequencies, value_t* time_delays,
                                value_t* phis, value_t* rec, value_t R,
                                value_t wmix, int bins, int grid_size, int N)
{
    int tid = blockDim.x*blockIdx.x+threadIdx.x;

    if(tid<grid_size*grid_size*bins) {

        int i = tid%grid_size;
        int j = (tid/grid_size)%grid_size;
        int k = tid/(grid_size*grid_size);

        value_t x = (-1+(value_t)(2*j+1)/grid_size)*R; //coords[j];
        value_t y = (-1+(value_t)(2*i+1)/grid_size)*R; //coords[i];
        value_t d = x*x+y*y;

        if(d<=R*R) {
			
			value_t fc = 2*M_PI*frequencies[k] + wmix; //maybe to shared memory
			//value_t f = wmix/2;

            thrust::complex<value_t> accum(0);
            int index_0 = grid_size*j+i;
            int grid_size2 = grid_size*grid_size;

			#pragma unroll
            for(int l=0; l<NANTENNAS; ++l) {
				//~ value_t phi = time_delays[(l*grid_size+j)*grid_size+i]*fc;
				//~ phi += phis[(l*grid_size+j)*grid_size+i];
				
				int index_l = grid_size2*l+index_0;
				value_t phi = time_delays[index_l]*fc + phis[index_l];
				
				//value_t phi = ((l*grid_size+j)*grid_size+i)*(2*M_PI*f+wmix);
				//phi += (l*grid_size+j)*grid_size+i;
				//value_t phi2 = phi*phi;
				//value_t phi3 = phi2*phi;
				
				//~ value_t real;
				//~ value_t imag;
				//~ __sincosf(phi, &imag, &real);
				//~ thrust::complex<value_t> phase(real, imag);
				thrust::complex<value_t> phase(__cosf(phi), __sinf(phi));
				
				//thrust::complex<value_t> phase(time_delays[index_l], phis[index_l]);

                accum += samples[l*bins+k]*phase;
            }
            rec[(grid_size*i+j)*bins+k] = thrust::norm(accum);
            //rec[(grid_size*i+j)*bins+k] = accum.real()*accum.real() + accum.imag()*accum.imag();

        }

    }
} */

template <typename value_t>
__global__ void reconstruction_red(thrust::complex<value_t>* samples,
                                value_t* frequencies, value_t* time_delays,
                                value_t* phis, value_t* rec, value_t R,
                                value_t wmix, int bins, int grid_size, int N,
                                int packet)
{
	
	/*
	 * Fastest so far
	 * 
	 * */
	
    int tid = blockDim.x*blockIdx.x+threadIdx.x;

    if(tid<grid_size*grid_size*bins) {
		
		int k = tid%bins;
        int i = (tid/bins)%grid_size;
        int j = tid/(bins*grid_size);

        value_t x = (-1+(value_t)(2*j+1)/grid_size)*R; //coords[j];
        value_t y = (-1+(value_t)(2*i+1)/grid_size)*R; //coords[i];
        value_t d = x*x+y*y;

        if(d<=R*R) {
			
			value_t fc = 2*M_PI*frequencies[k] + wmix; //maybe to shared memory
			//value_t f = wmix/2;

            thrust::complex<value_t> accum(0);

            for(int l=0; l<N; ++l) {
				value_t phi = time_delays[(j*grid_size+i)*N+l]*fc + phis[(j*grid_size+i)*N+l];
				
				thrust::complex<value_t> phase(__cosf(phi), __sinf(phi));
				
				//thrust::complex<value_t> phase(time_delays[(j*grid_size+i)*N+l]*fc, phis[(j*grid_size+i)*N+l]);

                accum += samples[(packet*N+l)*bins+k]*phase;
            }
            rec[((packet*grid_size+i)*grid_size+j)*bins+k] = thrust::norm(accum);

        }

    }
}

template <typename value_t>
__global__ void reconstruction_red2(thrust::complex<value_t>* samples,
                                value_t* frequencies,
                                value_t* coords, value_t* rec, value_t R,
                                value_t wmix, int bins, int grid_size, int N)
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
			
			value_t f = frequencies[k]; //maybe to shared memory

            thrust::complex<value_t> accum(0);

            for(int l=0; l<N; ++l) {
				value_t angle = frequencies[l]; //(value_t)l/N*2*M_PI;
				value_t xi = R*cos(angle);
				value_t yi = R*sin(angle);
				
				value_t dx = xi-x;
				value_t dy = yi-y;
				value_t dist = sqrt(dx*dx + dy*dy);
				
				value_t time_delay = dist/SPEED_OF_LIGHT;
				value_t phi_c = atan2(-dy, -dx) + M_PI;
				
				//value_t phi = time_delays[(l*grid_size+j)*grid_size+i]*(2*M_PI*f+wmix);
				//phi += phis[(l*grid_size+j)*grid_size+i];
				
				value_t phi = time_delay*(2*M_PI*f+wmix)+phi_c;
				thrust::complex<value_t> phase(cos(phi), sin(phi));
                accum += samples[l*bins+k]*phase;
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
		arma::Mat<value_t> time_delay = grid_time_delays[i].t();
		arma::Mat<value_t> phi_c = grid_phis[i].t();
        for(int l=0; l<bins; ++l) {
			value_t f = frequency(l);
            for(int k=0; k<grid_size; ++k) {
                for(int j=0; j<grid_size; ++j) {
                    //~ value_t phi = grid_time_delays[i](j,k)*(2*M_PI*frequency(l)+wmix);
                    //~ phi += grid_phis[i](j,k);
                    value_t phi = time_delay(j,k)*(2*M_PI*f+wmix);
                    phi += phi_c(j,k);
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
    //~ cudaMalloc(&(this->grid_phase),
        //~ N*grid_size*grid_size*frequency.n_elem*sizeof(std::complex<value_t>)); CUERR

    //~ cudaMemcpy(this->grid_phase, *grid_phase,
        //~ N*grid_size*grid_size*frequency.n_elem*sizeof(std::complex<value_t>),
        //~ H2D);                           CUERR
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
void Reconstruction<value_t>::init_gpu()
{

    thrust::host_vector<value_t> time_delays_H(N*grid_size*grid_size);	CUERR
    thrust::host_vector<value_t> phis_H(N*grid_size*grid_size);			CUERR

	for(int i=0; i<N; ++i) {
		arma::Mat<value_t> mat_delay= grid_time_delays[i].t();
		arma::Mat<value_t> mat_phi = grid_phis[i].t();
		thrust::copy(mat_delay.begin(),mat_delay.end(),
					time_delays_H.begin()+i*grid_size*grid_size);		CUERR
		thrust::copy(mat_phi.begin(),mat_phi.end(),
					phis_H.begin()+i*grid_size*grid_size);		CUERR
	}
	
								
	cudaMalloc(&(this->frequencies_dev), frequency.n_elem*sizeof(value_t)); CUERR

	cudaMalloc(&(this->time_delays_dev), N*grid_size*grid_size*sizeof(value_t)); CUERR

	cudaMalloc(&(this->phis_dev), N*grid_size*grid_size*sizeof(value_t)); CUERR

	
	//cudaMemcpy(dev_test, test.data(), test.size(), H2D);		CUERR
	
	std::cerr << "copy2" << std::endl;
	//~ thrust::copy(time_delays_H.begin(), time_delays_H.end(), 
											//~ this->time_delays_dev);		CUERR
	cudaMemcpy(this->time_delays_dev, time_delays_H.data(), 
										time_delays_H.size(), H2D);		CUERR
	
	std::cerr << "copy3" << std::endl;				
	//~ thrust::copy(phis_H.begin(), phis_H.end(), 
											//~ this->phis_dev);		CUERR
											
	cudaMemcpy(this->phis_dev, phis_H.data(), 
										phis_H.size(), H2D);		CUERR
	
	std::cerr << "copy1" << std::endl;
	//~ thrust::copy(this->frequency.begin(),this->frequency.end(), 
											//~ this->frequencies_dev);		CUERR
											
	cudaMemcpy(this->frequencies_dev, this->frequency.memptr(), 
										this->frequency.size(), H2D);		CUERR
}

#define FREE(pointer) \
    do { \
		if(pointer) { \
			cudaFree(pointer); \
			CUERR \
			pointer=nullptr; \
		} \
    } while (false)  

template <typename value_t>
void Reconstruction<value_t>::free_gpu()
{
	FREE(grid_phase);
	FREE(frequencies_dev);
	FREE(time_delays_dev);
	FREE(phis_dev);
}

template <typename value_t>
void Reconstruction<value_t>::run(const std::vector<std::vector<Data_Packet<value_t>>>& samples)
{
	std::cerr << "run" << std::endl;

    TIMERSTART(REC)

    int bins = samples[0][0].frequency.n_elem;
    int n_packets = samples[0].size();
    std::cerr << "packets: " << n_packets << " antennas: " << samples.size() << std::endl;
    std::cerr << samples.size() << " " << samples[0].size() << std::endl;
    
    
    //~ for(int j=0; j<n_packets; ++j) {

        //~ std::vector<Data_Packet<float>> data_in;

        //~ for(int i=0; i<data_out.size(); ++i)
            //~ data_in.push_back(data_out[i][j]);

    //~ }

    //reorder data
    thrust::host_vector<thrust::complex<value_t> > samples_H(n_packets*N*bins, 
										thrust::complex<value_t>(0,1));   CUERR
										
	std::cerr << "Copy CPU" << std::endl;
	size_t memSize=sizeof(thrust::complex<value_t>)*n_packets*N*bins;
    TIMERSTART(COPY_CPU)
    for(int j=0; j<n_packets; ++j) {
		for(int i=0; i<N; ++i) {
			//std::cerr << i << " " << j << " " << samples[i][j].frequency_data[0] << std::endl;
			thrust::copy(samples[i][j].frequency_data.begin(),
							samples[i][j].frequency_data.end(),
							samples_H.begin()+(j*N+i)*bins);                      CUERR
		}
	}
    //time_delays[(l*grid_size+j)*grid_size+i]
    //TIMERSTOP(COPY_CPU)
    TIMERBW(memSize, COPY_CPU)

    //copy data to GPU
    
    std::cerr << "Copy GPU" << std::endl;
    TIMERSTART(COPY_GPU)
    thrust::device_vector<thrust::complex<value_t> > samples_D = samples_H;  CUERR
                                                            //(
                                                            //samples_H.begin(),
                                                            //samples_H.end()); CUERR

    //thrust::device_vector<value_t> coords_D(grid.coords.begin(),
    //                                        grid.coords.end());         CUERR
                                            

    //TIMERSTOP(COPY_GPU)
    TIMERBW(memSize, COPY_GPU)

    thrust::device_vector<value_t> reconstructed_D(grid_size*grid_size*bins*n_packets);
    thrust::fill(reconstructed_D.begin(), reconstructed_D.end(), value_t(-1));

    thrust::complex<value_t>* samples_dev = thrust::raw_pointer_cast(samples_D.data());
    value_t* rec_dev = thrust::raw_pointer_cast(reconstructed_D.data());
    //value_t* coords_dev = thrust::raw_pointer_cast(coords_D.data());
    //~ value_t* time_delays_dev = thrust::raw_pointer_cast(time_delays_D.data());
	//~ value_t* phis_dev = thrust::raw_pointer_cast(phis_D.data());
	//~ value_t* frequencies_dev = thrust::raw_pointer_cast(frequencies_D.data());
	
    int threads = 512;
    int tasks=grid_size*grid_size*bins;
    int blocks = SDIV(tasks, threads);


    //~ TIMERSTART(KERNEL)
    //~ reconstruction<<<blocks, threads>>>(samples_dev,
                                        //~ (thrust::complex<value_t>*)grid_phase,
                                        //~ coords_dev, rec_dev, grid.R,
                                        //~ bins, grid_size, N);   CUERR
    //~ TIMERSTOP(KERNEL)
    
    TIMERSTART(KERNELS)
    for(int packet=0; packet<n_packets; ++packet) {
		std::cerr << packet << std::endl;
		TIMERSTART(KERNEL_RED)
		reconstruction_red<<<blocks, threads>>>(samples_dev,
											frequencies_dev, time_delays_dev,
											phis_dev, rec_dev, grid.R,
											wmix, bins, grid_size, N, packet);   CUERR
		TIMERSTOP(KERNEL_RED)
	}
	TIMERSTOP(KERNELS)

    //copy back result
    size_t memsize_res = grid_size*grid_size*bins*n_packets*sizeof(value_t);
    TIMERSTART(COPY_BACK)
    thrust::copy(reconstructed_D.begin(), reconstructed_D.end(),
                    reconstructed.begin());                             CUERR
	//TIMERSTOP(COPY_BACK)
	TIMERBW(memsize_res, COPY_BACK)

    TIMERSTOP(REC)
}

//DEFINE_TEMPLATES(Reconstruction)
template class Reconstruction<float>;

#endif
