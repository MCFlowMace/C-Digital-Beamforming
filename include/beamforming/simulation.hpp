/*
 * simulation.hpp
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

#include <vector>
#include <armadillo>

#include "beamforming/event.hpp"
#include "beamforming/data_packet.hpp"

template<typename value_t>
class Simulation_Settings
{

	public:
		//event generation
		int n_events;
		value_t w_min;
		value_t w_max;
		value_t R;
		value_t run_duration;
		value_t mean_event_lifetime;
		value_t trap_efficiency;
		long seed;
		
		//for non-random events
		bool manual;
		value_t e_r;
		value_t e_phi;
		value_t w0;

		//event observation and data generation
		int N; //antennas
		value_t snr;
		value_t sample_rate;
		value_t w_mix;
		int n_samples; //for fourier transform
		
		Simulation_Settings();
};

template <typename value_t>
class Simulation{

    public:

        //Simulation(value_t run_duration);
        Simulation(Simulation_Settings<value_t> settings);

        std::vector<std::vector<Data_Packet<value_t>>> observation(value_t t_start,
                                                        int n_packets);
        std::vector<std::complex<value_t>> observation_flat(value_t t_start,
                                                        int n_packets);

        arma::Mat<value_t> w_mat;
        
        Simulation_Settings<value_t> get_settings();

    private:

        Simulation_Settings<value_t> settings;
        std::vector<Event<value_t>> events;

        void generation(long seed=-1);
        void manual_event();
        void fill_w_mat();
};

typedef Simulation_Settings<float> Simulation_Settingsf;
