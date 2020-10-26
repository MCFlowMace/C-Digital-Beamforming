/*
 * beamformer.cpp
 * 
 * Copyright 2020 Florian Thomas <flthomas@students.uni-mainz.de>
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


#include "beamforming/beamformer.hpp"
#include "utility/utility_macros.hpp"


template <typename value_t>
Beamformer<value_t>::Beamformer(Simulation_Settings<value_t> settings, 
				int grid_size, int n_packets, 
				bool weighted, bool full_frequency,
				value_t r_grid):
sim{settings},
rec{grid_size, 
	n_packets, 
	Data_Packet<value_t>::get_frequency(settings.n_samples, 
						1/settings.sample_rate,
						full_frequency), 
	Antenna_Array<value_t>{settings.N, 
							settings.R, 
							settings.snr, 
							settings.w_mix, 
							settings.sample_rate}, 
	weighted,
	r_grid},
n_packets{n_packets},
packet_counter{0}
{
}


template <typename value_t>
void Beamformer<value_t>::run_next()
{
    std::cout << "running run_next" << std::endl;
	Simulation_Settings<value_t> settings{sim.get_settings()};
    int n_samples = settings.n_samples;
    value_t dt = 1/settings.sample_rate;
	value_t t_start = packet_counter*n_samples*dt;
	packet_counter += n_packets;

	TIMERSTART(SAMPLE)
	std::vector<std::complex<value_t>> data = sim.observation_flat(
													t_start, n_packets);
	TIMERSTOP(SAMPLE)

	this->run(data);
	
}

template <typename value_t>
void Beamformer<value_t>::run(std::vector<std::complex<value_t>> data)
{
	rec.run(data);	
}

template <typename value_t>
void Beamformer<value_t>::get_next(value_t *dest)
{
	this->run_next();
	this->rec.copy_res(dest);
}

template <typename value_t>
void Beamformer<value_t>::get_result(std::complex<value_t> *src, 
					value_t *dest)
{
        Simulation_Settings<value_t> settings{sim.get_settings()};
    
	size_t len = this->n_packets*settings.N*settings.n_samples;
	std::vector<std::complex<value_t>> data;
	data.assign(src, src+len);
	this->run(data);
	this->rec.copy_res(dest);
}

template <typename value_t>
arma::Mat<value_t> Beamformer<value_t>::get_next_img()
{
	this->run_next();

    unsigned int index_max = rec.get_max_bin(0);
    Simulation_Settings<value_t> settings{sim.get_settings()};
    
    arma::Col<value_t> frequency{Data_Packet<value_t>::get_frequency(settings.n_samples, 
										1/settings.sample_rate, false)};
										
	value_t delta_f = frequency[1]-frequency[0];
    
    std::cerr << " max ind: " <<  index_max << std::endl;
   
	std::cerr << "frequency: " << (settings.w_mix/(2*M_PI)+frequency[index_max])/1e9
				<< " deltaf: " << delta_f << std::endl;

    return rec.get_img(0, index_max);     
}

template <typename value_t>
arma::Mat<value_t> Beamformer<value_t>::get_next_max_vals()
{
	this->run_next();
	
	return rec.get_max_vals();
}
    
DEFINE_TEMPLATES(Beamformer)
