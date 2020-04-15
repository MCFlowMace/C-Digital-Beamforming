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
#include <cmath>
#include <armadillo>
#include <cstdlib>

#include "electron.hpp"
#include "antenna.hpp"
#include "antenna_array.hpp"
#include "reconstruction.hpp"
#include "simulation.hpp"

int main(int argc, char **argv)
{

    if(argc!=3) {
        std::cerr << "args: [grid_size] [N_samples] !" << std::endl;
        exit(0);
    }

    Simulation_Settings<float> settings;


    int grid_size = std::atoi(argv[1]);

    settings.n_events = 1;
    settings.w_min = 2*M_PI*24.6*1e9;
    settings.w_max = 2*M_PI*26.2*1e9;
    settings.R = 5.0f;
    settings.mean_event_lifetime = 1/(2*1e-4f);
    settings.trap_efficiency = 0.5f;
    settings.run_duration = 0.005f;

    //event observation and data generation
    settings.N = 30; //antennas
    settings.snr = 1000.0f;
    settings.sample_rate = 3.2*1e9;
    settings.w_mix = 2*M_PI*24.6*1e9;
    settings.n_samples = std::atoi(argv[2]); //for fourier transform

    Simulation<float> sim(settings);

    TIMERSTART(SAMPLE)
    std::vector<std::vector<Sample<float>>> data_out = sim.observation(0.0f, 2*1000*settings.n_samples/settings.sample_rate);
    TIMERSTOP(SAMPLE)

    std::vector<Sample<float>> data_in;

    for(int i=0; i<data_out.size(); ++i)
        data_in.push_back(data_out[i][0]);

    Reconstruction<float> rec(grid_size, data_in[0]);

    Antenna_Array<float> array(settings.N, settings.R, settings.snr, settings.w_mix, settings.sample_rate);

    rec.set_antenna_array(array);
    rec.run(data_in);
    auto img = rec.img;

    img.print();

}

