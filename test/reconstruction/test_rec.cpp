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

#include <iostream>

#include "utility/hpc_helpers.hpp"

#include "beamforming/beamformer.hpp"
#include "beamforming/simulation.hpp"

//~ template <typename value_t>
//~ void run_test(int grid_size, int n_samples,value_t snr, int seed, bool weighted, 
				//~ value_t e_r, value_t e_phi, value_t f0, int N)
//~ {
	//~ Simulation_Settings<value_t> settings;

    //~ std::cerr << weighted << std::endl;

    //~ settings.n_events = 1;
    //~ settings.w_min = 2*M_PI*24.6*1e9;
    //~ settings.w_max = 2*M_PI*26.2*1e9;
    //~ settings.R = 5.0f;
    //~ settings.mean_event_lifetime = 1/(2*1e-4f);
    //~ settings.trap_efficiency = 0.5f;
    //~ settings.run_duration = 0.005f;
    //~ settings.seed = seed;
		
	//~ settings.e_r = e_r;
    //~ settings.e_phi = e_phi;
    //~ settings.w0 = 2*M_PI*f0*1e9;
    
    //~ if(settings.e_r >=0)
		//~ settings.manual = true;
	//~ else
		//~ settings.manual = false;

    //~ //event observation and data generation
    //~ settings.N = N; //30; //antennas
    //~ settings.snr = snr;
    //~ settings.sample_rate = 3.2*1e9;
    //~ settings.w_mix = 2*M_PI*24.6*1e9;
    //~ settings.n_samples = n_samples; //for fourier transform

    //~ int n_packets = 2; //use two packets because in manual simulation mode half of the packets will be noise and we want 1 full packet of signal to test the reconstruction
    //~ settings.run_duration = n_packets*settings.n_samples/settings.sample_rate; //5 data packets
	//~ value_t dt = 1/settings.sample_rate;
	
    //~ Simulation<value_t> sim(settings);
    
    //~ Antenna_Array<value_t> array(settings.N, settings.R, settings.snr, 
									//~ settings.w_mix, settings.sample_rate);
    //~ arma::Col<value_t> frequency = Data_Packet<value_t>::get_frequency(
														//~ settings.n_samples, dt);

    //~ Reconstruction_GPU<value_t> rec(grid_size, n_packets, frequency, array, weighted);
    
	//~ value_t t_start{0};

	//~ TIMERSTART(SAMPLE)

	//~ std::vector<std::complex<value_t>> data = sim.observation_flat(
													//~ t_start, n_packets);
	//~ TIMERSTOP(SAMPLE)

	//~ rec.run(data);
	
    //~ unsigned int index_max = rec.get_max_bin(0);
	//~ float delta_f = frequency[1]-frequency[0];
	
   //~ /*
	//~ unsigned int index_max = 0;
	
	
	//~ for(int i=0; i<frequency.n_elem; ++i) {
		//~ value_t w =	settings.w_mix+(frequency[i]+delta_f/2)*2*M_PI;
		//~ if(w>settings.w0) {
			//~ index_max = i;
			//~ break;
		//~ }
	//~ } */
		
    //~ float max_val = rec.get_max_val(index_max, 0);
    
    //~ std::cerr << " max ind: " <<  index_max << std::endl;

    //~ std::cerr << "val: " << max_val << std::endl;
   
	//~ std::cerr << "frequency: " << (settings.w_mix/(2*M_PI)+frequency[index_max])/1e9
				//~ << " deltaf: " << delta_f << std::endl;

    //~ auto img = rec.get_img(0, index_max);

    //~ //img.print();
    //~ rec.print(0);
    //~ rec.print(1);

	//~ std::cerr << "done!" << std::endl;
//~ }


int main(int argc, char **argv)
{
	
	Simulation_Settings<float> settings{};
	
	int n_packets = 2;
	int grid_size = 100;
	bool weighted = true;
	settings.manual=true;
	settings.run_duration = n_packets*settings.n_samples/settings.sample_rate;
	
	Beamformer<float> bf{settings, grid_size, n_packets, weighted};
	
	auto img = bf.get_next_img();
	
	img.print();
}

