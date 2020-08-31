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
#include <stdexcept>
#include <cfloat>

#include <thrust/complex.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include "hpc_helpers.hpp"
#include "data_packet.hpp"
#include "reconstruction_gpu.hpp"
#include "utility_macros.hpp"

//speed of light in cm/s
#define SPEED_OF_LIGHT 29979245800


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
__device__ inline value_t weighted_beamforming(value_t const fc, int const N, 
												int const grid_size, 
												int const packet, int const bins,
												int const i, int const j, int const k,
												value_t const * const time_delays,
												value_t const * const phis,
												thrust::complex<value_t> const * const samples)
{
	thrust::complex<value_t> accum(0);
	value_t A{0};
	for(int l=0; l<N; ++l) {
		value_t t = time_delays[(j*grid_size+i)*N+l];
		value_t phi = t*fc + phis[(j*grid_size+i)*N+l];
		thrust::complex<value_t> phase(__cosf(phi), __sinf(phi));
		
		A += 1/(t*t);

		accum += samples[(packet*N+l)*bins+k]*phase/t;
		
	}
	return thrust::norm(accum)*SPEED_OF_LIGHT/A;
}

template <typename value_t>
__device__ inline value_t beamforming(value_t const fc, int const N, 
										int const grid_size, 
										int const packet, int const bins,
										int const i, int const j, int const k,
										value_t const * const time_delays,
										value_t const * const phis,
										thrust::complex<value_t> const * const samples)
{
	
	thrust::complex<value_t> accum(0);
	for(int l=0; l<N; ++l) {
		value_t phi = time_delays[(j*grid_size+i)*N+l]*fc + phis[(j*grid_size+i)*N+l];
		
		thrust::complex<value_t> phase(__cosf(phi), __sinf(phi));
		
		//thrust::complex<value_t> phase(time_delays[(j*grid_size+i)*N+l]*fc, phis[(j*grid_size+i)*N+l]);

		accum += samples[(packet*N+l)*bins+k]*phase;
		
	}
	return thrust::norm(accum);			
}

template <typename value_t, bool weighted>
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
    
    //printf("id: %d\n", tid);

    if(tid<grid_size*grid_size*bins) {
		
		int k = tid%bins; //frequency
        int i = (tid/bins)%grid_size; //y
        int j = tid/(bins*grid_size); //x

        value_t x = (-1+(value_t)(2*j+1)/grid_size)*R; //coords[j];
        value_t y = (-1+(value_t)(2*i+1)/grid_size)*R; //coords[i];
        value_t d = x*x+y*y;

        if(d<=R*R) {
			
			value_t fc = 2*M_PI*frequencies[k] + wmix; //maybe to shared memory
			//value_t f = wmix/2;
			value_t result;
            
            
            if(weighted) {
				result = weighted_beamforming(fc, N, grid_size, packet, bins,
												i, j, k, time_delays, phis, samples);
			} else {
				result = beamforming(fc, N, grid_size, packet, bins,
												i, j, k, time_delays, phis, samples);
			}
			rec[((packet*grid_size+i)*grid_size+j)*bins+k] = result;
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
Reconstruction_GPU<value_t>::Reconstruction_GPU(int grid_size, int n_packets,
						arma::Col<value_t> frequency,
						const Antenna_Array<value_t>& array, bool weighted):
Reconstruction<value_t>(grid_size, n_packets, frequency, array, weighted)
{
	init_gpu();
}

template <typename value_t>
Reconstruction_GPU<value_t>::~Reconstruction_GPU()
{
	FREE_CUDA(frequencies_D);
	FREE_CUDA(time_delays_D);
	FREE_CUDA(phis_D);
	FREE_CUDA(samples_D);
	FREE_CUDA(reconstructed_D);
	
	FREE_HOST(samples_H);
	//FREE_HOST(reconstructed_H);
	free(reconstructed_H);
}

/*
template <typename value_t>
void Reconstruction_GPU<value_t>::calc_phase(
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
} */

template <typename value_t>
unsigned int Reconstruction_GPU<value_t>::get_max_bin(unsigned int packet)
{
	
	int bins = this->bins;
	int grid_size = this->grid_size;

	value_t max_val = std::numeric_limits<value_t>::min();

	unsigned int index;
	//int packet = 0;
	for(unsigned int i=0; i<bins; ++i) {
		for(unsigned int j=0; j<grid_size; ++j) {
			for(unsigned int l=0; l<grid_size; ++l) {
				
				value_t val = reconstructed_H[((packet*grid_size+j)*grid_size+l)*bins+i];
				
				//std::cerr << i << " " << val << " " << max_val << std::endl;
				if(val>max_val && !std::isinf(val)) {
					max_val=val;
					index = i;
				}
			}
		}
	}

	//std::cerr << "Max frequency: " << this->frequency[index] 
	//			<< " (bin: " << index << ") val: " << max_val << " packet: "
	//			<< packet << std::endl;
				
	/*if(index>bins) {
		for(unsigned int i=0; i<bins; ++i) {
			for(unsigned int j=0; j<grid_size; ++j) {
				for(unsigned int l=0; l<grid_size; ++l) {
					
					value_t val = reconstructed_H[((packet*grid_size+j)*grid_size+l)*bins+i];
					
					std::cerr << i << " " << j << " " << l << " "
								<< val << " " << max_val << std::endl;
				}
			}
		}
	}*/

	return index;
}

template <typename value_t>
void Reconstruction_GPU<value_t>::init_gpu()
{
	
	std::cerr << "init" << std::endl;
	
	int N = this->N;
	int bins = this->bins;
	int n_packets = this->n_packets;
	int grid_size = this->grid_size;

    thrust::host_vector<value_t> time_delays_H(N*grid_size*grid_size);	CUERR
    thrust::host_vector<value_t> phis_H(N*grid_size*grid_size);			CUERR

	for(int l=0; l<N; ++l) {
		/*arma::Mat<value_t> mat_delay= this->grid_time_delays[i].t();
		arma::Mat<value_t> mat_phi = this->grid_phis[i].t();
		thrust::copy(mat_delay.begin(),mat_delay.end(),
					time_delays_H.begin()+i*grid_size*grid_size);		CUERR
		thrust::copy(mat_phi.begin(),mat_phi.end(),
					phis_H.begin()+i*grid_size*grid_size);				CUERR*/
		for(int i=0; i<grid_size; ++i) {
			for(int j=0; j<grid_size; ++j) {
				value_t time_delay = this->grid_time_delays[l](j,i);
				value_t phi = this->grid_phis[l](j,i);
				time_delays_H[(j*grid_size+i)*N+l] = time_delay;
				phis_H[(j*grid_size+i)*N+l] = phi;
			}
		}
	}
	
	cudaMallocHost(&(this->samples_H), 
				n_packets*N*bins*sizeof(thrust::complex<value_t>));		CUERR
	std::cerr << "allocating " 
				<< (grid_size*grid_size*bins*n_packets*sizeof(value_t)/1e9) 
				<< "GB of unpinned host memory" << std::endl;
	//using pinned memory results in frequent errors here for whatever reason ...
	this->reconstructed_H=(value_t*) malloc(grid_size*grid_size*bins*n_packets*sizeof(value_t));
	//cudaMallocHost(&(this->reconstructed_H), 
	//			grid_size*grid_size*bins*n_packets*sizeof(value_t));	CUERR
				
	cudaMalloc(&(this->samples_D), 
				n_packets*N*bins*sizeof(thrust::complex<value_t>));		CUERR
	cudaMalloc(&(this->reconstructed_D), 
				grid_size*grid_size*bins*n_packets*sizeof(value_t));	CUERR
				
	cudaMalloc(&(this->frequencies_D), bins*sizeof(value_t)); 			CUERR
	cudaMalloc(&(this->time_delays_D), 
					N*grid_size*grid_size*sizeof(value_t)); 			CUERR
	cudaMalloc(&(this->phis_D), N*grid_size*grid_size*sizeof(value_t)); CUERR


	cudaMemcpy(this->time_delays_D, time_delays_H.data(), 
						time_delays_H.size()*sizeof(value_t), H2D);		CUERR									
	cudaMemcpy(this->phis_D, phis_H.data(), 
						phis_H.size()*sizeof(value_t), H2D);			CUERR						
	cudaMemcpy(this->frequencies_D, this->frequency.memptr(), 
						this->frequency.size()*sizeof(value_t), H2D);	CUERR
}

//~ template <typename value_t>
//~ void Reconstruction_GPU<value_t>::free_memory()
//~ {
	//~ //FREE_CUDA(grid_phase_D);
	
	//~ FREE_CUDA(frequencies_D);
	//~ FREE_CUDA(time_delays_D);
	//~ FREE_CUDA(phis_D);
	//~ FREE_CUDA(samples_D);
	//~ FREE_CUDA(reconstructed_D);
	
	//~ FREE_HOST(samples_H);
	//~ FREE_HOST(reconstructed_H);
//~ }

template <typename value_t>
arma::Mat<value_t> Reconstruction_GPU<value_t>::get_img(unsigned int packet, 
															unsigned int bin)
{
	//int packet=0;
	
	if(bin >= this->frequency.n_elem) {
		std::string err = "Bin " + std::to_string(bin) 
								+ " is not a valid frequency bin!";
		throw std::out_of_range(err);
	}
	
	//std::cerr << "Fetching image for packet " << packet << " and bin " << bin << std::endl;

	int grid_size = this->grid_size;
	int bins = this->bins;
	arma::Mat<value_t> img(grid_size, grid_size);

	for(int i=0; i<grid_size; ++i) {
		for(int j=0; j<grid_size; ++j) {
			img(j,i) = reconstructed_H[((packet*grid_size+i)*grid_size+j)*bins+bin];

		}
	}

	return img;
}

template <typename value_t>
void Reconstruction_GPU<value_t>::print(unsigned int packet)
{
	int grid_size = this->grid_size;
	int bins = this->bins;
	arma::Mat<value_t> img(grid_size, grid_size);
	
	for(int k=0; k<bins; ++k) {
		for(int i=0; i<grid_size; ++i) {
			for(int j=0; j<grid_size; ++j) {
				std::cout << reconstructed_H[((packet*grid_size+i)*grid_size+j)*bins+k] << ' ';
			}
			std::cout << '\n';
		}
	}

}

template <typename value_t>
void Reconstruction_GPU<value_t>::run(const std::vector<std::complex<value_t>>& samples)
{
	std::cerr << "run" << std::endl;
		
	int N = this->N;
	int bins = this->bins;
	int n_packets = this->n_packets;
	int grid_size = this->grid_size;
	
	//~ int bins_ = samples[0][0].frequency.n_elem;
    //~ int n_packets_ = samples[0].size();
	//~ int N_ = samples.size();
	
	//~ if(N!=N_ || bins != bins_ || n_packets != n_packets_)
        //~ throw std::invalid_argument(
				//~ "'samples' input array dimension is (" + std::to_string(N_) 
				//~ + "," + std::to_string(n_packets_) + "," 
				//~ + std::to_string(bins_) + ") but expected dimension was ("
				//~ + std::to_string(N) + "," + std::to_string(n_packets) + ","
				//~ + std::to_string(bins) + ")" );

	size_t rec_size = grid_size*grid_size*bins*n_packets;
	size_t samples_size = n_packets*N*bins;
	
    TIMERSTART(REC)
    								
	std::cerr << "Copy CPU" << std::endl;
	size_t memsize_samples=sizeof(thrust::complex<value_t>)*samples_size;
	

    TIMERSTART(COPY_CPU)
    //~ for(int j=0; j<n_packets; ++j) {
		//~ for(int i=0; i<N; ++i) {
			//~ //thrust::copy(samples[i][j].frequency_data.begin(),
			//~ //				samples[i][j].frequency_data.end(),
			//~ //				samples_H+(j*N+i)*bins);					CUERR
			//~ cudaMemcpy(this->samples_H, samples[i][j].frequency_data.memptr(), 
							//~ bins*sizeof(thrust::complex<value_t>), H2H);	CUERR
		//~ }
	//~ }
	cudaMemcpy(this->samples_H, samples.data(), 
											memsize_samples, H2H);			CUERR	
    //TIMERSTOP(COPY_CPU)
    TIMERBW(memsize_samples, COPY_CPU)

    //copy data to GPU
    
    std::cerr << "Copy GPU" << std::endl;
    TIMERSTART(COPY_GPU)
    
	cudaMemcpy(this->samples_D, this->samples_H, 
											memsize_samples, H2D);			CUERR	


    //TIMERSTOP(COPY_GPU)
    TIMERBW(memsize_samples, COPY_GPU)
	
	thrust::device_ptr<value_t> rec_thrust =
									thrust::device_pointer_cast(reconstructed_D);  
    thrust::fill(rec_thrust, rec_thrust+rec_size, value_t(-1));			CUERR
	
    int threads = 512;
    int tasks=grid_size*grid_size*bins;
    int blocks = SDIV(tasks, threads);

    //~ TIMERSTART(KERNEL)
    //~ reconstruction<<<blocks, threads>>>(samples_dev,
                                        //~ (thrust::complex<value_t>*)grid_phase,
                                        //~ coords_dev, rec_dev, grid.R,
                                        //~ bins, grid_size, N);   CUERR
    //~ TIMERSTOP(KERNEL)
    
    std::cerr << "start kernel" << std::endl;
    TIMERSTART(KERNELS)
    if(this->weighted) {
		for(int packet=0; packet<n_packets; ++packet) {
			reconstruction_red<value_t, true><<<blocks, threads>>>((thrust::complex<value_t>*)samples_D,
												frequencies_D, time_delays_D,
												phis_D, reconstructed_D, this->R,
												this->wmix, bins, grid_size, N, packet);   CUERR
		}
	} else {
		for(int packet=0; packet<n_packets; ++packet) {
			reconstruction_red<value_t, false><<<blocks, threads>>>((thrust::complex<value_t>*)samples_D,
												frequencies_D, time_delays_D,
												phis_D, reconstructed_D, this->R,
												this->wmix, bins, grid_size, N, packet);   CUERR
		}
	}
	TIMERSTOP(KERNELS)

    //copy back result
    size_t memsize_res = rec_size*sizeof(value_t);
    TIMERSTART(COPY_BACK)
    //thrust::copy(reconstructed_D.begin(), reconstructed_D.end(),
    //                reconstructed.begin());                             CUERR
	cudaMemcpy(this->reconstructed_H, this->reconstructed_D, 
											memsize_res, D2H);			CUERR	
	//TIMERSTOP(COPY_BACK)
	TIMERBW(memsize_res, COPY_BACK)

    TIMERSTOP(REC)
}

DEFINE_TEMPLATES(Reconstruction_GPU)

