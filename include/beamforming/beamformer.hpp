/*
 * beamformer.hpp
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


#pragma once

#include "beamforming/simulation.hpp"
#include "beamforming/antenna_array.hpp"
#include "beamforming/reconstruction_gpu.hpp"

#include <armadillo>

template <typename value_t>
class Beamformer
{
	
	public:
	
		Beamformer(Simulation_Settings<value_t> settings, int grid_size, 
					int n_packets, bool weighted, bool full_frequency,
					value_t r_grid);

		void get_next(value_t* dest);
		void get_result(std::complex<value_t> *src, 
					value_t *dest);
		arma::Mat<value_t> get_next_img();
		arma::Mat<value_t> get_next_max_vals();
	
	private:
	
		Reconstruction_GPU<value_t> rec;
		Simulation<value_t> sim;
		
		void run_next();
		void run(std::vector<std::complex<value_t>> data);
		
		uint32_t packet_counter;
		int n_packets;
};


typedef Beamformer<float> Beamformerf;
typedef Beamformer<double> Beamformerd;
