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

#include "reconstruction.hpp"
#include "reconstruction_gpu.hpp"
 
#ifdef PRECISION
	typedef double value_t;
#else
	typedef float value_t;
#endif
 
#ifdef USE_GPU
	typedef Reconstruction_GPU<value_t> Rec_Type;
#else
	typedef Reconstruction_CPU<value_t> Rec_Type;
#endif

#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdlib>

#include "electron.hpp"
#include "antenna.hpp"
#include "antenna_array.hpp"
#include "simulation.hpp"
#include "hpc_helpers.hpp"

int main(int argc, char **argv)
{

    if(argc!=5) {
        std::cerr << "args: [grid_size] [N_samples] [snr] [seed]!" << std::endl;
        exit(0);
    }

    Simulation_Settings<value_t> settings;

    int grid_size = std::atoi(argv[1]);

    settings.n_events = 1;
    settings.w_min = 2*M_PI*24.6*1e9;
    settings.w_max = 2*M_PI*26.2*1e9;
    settings.R = 5.0f;
    settings.mean_event_lifetime = 1/(2*1e-4f);
    settings.trap_efficiency = 0.5f;
    settings.run_duration = 0.005f;
    settings.seed = std::atoi(argv[4]);
    
    if(settings.seed==-1)
		settings.manual = true;
	else
		settings.manual= false;
		
    settings.e_r = 3.5f;
    settings.e_phi = 0.0f;
    settings.w0 = 2*M_PI*26*1e9;

    //event observation and data generation
    settings.N = 30; //antennas
    settings.snr = std::atoi(argv[3]);
    settings.sample_rate = 3.2*1e9;
    settings.w_mix = 2*M_PI*24.6*1e9;
    settings.n_samples = std::atoi(argv[2]); //for fourier transform

    int n_packets = 1;
    settings.run_duration = n_packets*settings.n_samples/settings.sample_rate; //5 data packets
	value_t dt = 1/settings.sample_rate;
	
    Simulation<value_t> sim(settings);
    
    Antenna_Array<value_t> array(settings.N, settings.R, settings.snr, 
									settings.w_mix, settings.sample_rate);
    arma::Col<value_t> frequency = Data_Packet<value_t>::get_frequency(
														settings.n_samples, dt);

    Rec_Type rec(grid_size, n_packets, frequency, array);
    
	value_t t_start{0};

	TIMERSTART(SAMPLE)

	std::vector<std::complex<value_t>> data = sim.observation_flat(
													t_start, n_packets);
	TIMERSTOP(SAMPLE)

	rec.run(data);
	
    unsigned int index_max = rec.get_max_bin(0);

    float max_val = rec.get_max_val(index_max, 0);

    std::cerr << "val: " << max_val << std::endl;
   
	std::cerr << "frequency: " << settings.w_mix/(2*M_PI)+frequency[index_max]
				<< " deltaf: " << frequency[1]-frequency[0] << std::endl;

    auto img = rec.get_img(0, index_max);

    img.print();

	std::cerr << "done!" << std::endl;
}

