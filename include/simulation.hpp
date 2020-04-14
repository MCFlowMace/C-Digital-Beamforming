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

#include "event.hpp"

template<typename value_t>
struct Simulation_Settings
{

    //event generation
    int n_events;
    value_t w_min;
    value_t w_max;
    value_t R;
    value_t run_duration;
    value_t mean_event_lifetime;
    value_t trap_efficiency;

    //event observation and data generation
    int N; //antennas
    value_t snr;
    value_t sample_rate;
    value_t w_mix;
    int n_samples; //for fourier transform
};

template <typename value_t>
class Simulation{

    public:

        //Simulation(value_t run_duration);
        Simulation(Simulation_Settings<value_t> settings);

    private:

        Simulation_Settings<value_t> settings;
        std::vector<Event<value_t>> events;

        void event_generation();
        void event_observation();
};

