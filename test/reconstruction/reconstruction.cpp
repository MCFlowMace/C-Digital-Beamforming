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
#include "hpc_helpers.hpp"

int main(int argc, char **argv)
{

    if(argc!=3) {
        std::cerr << "args: [grid_size] [N_samples] !" << std::endl;
        exit(0);
    }

    float e_r = 3.0f;
    float e_phi = 0.0f;
    float w0 = 2*M_PI*26*1e9;
    float snr = 0.5f;
    float sample_rate = 3.2*1e9;
    float wmix = 2*M_PI*24.6*1e9;
    int N= 30;
    float R = 5.0f;
    int grid_size = std::atoi(argv[1]);
    int n_samples= std::atoi(argv[2]);

    auto emission = [](float t, float w, float phi) -> float {
                        float dw = 300*1e6;
                        return cosf((w+dw*t)*t+phi);
                    };

    Electron<float> e {w0, e_r*cosf(e_phi), e_r*sinf(e_phi), emission};


    Antenna_Array<float> array(N, R, snr, wmix, sample_rate);

    std::vector<Sample<float>> data;


    for(int i=0; i<N; ++i) {
        data.push_back(array.antennas[i].sample_data(n_samples, 0.0f, e));
    }

    Reconstruction<float> rec(grid_size, data[0].frequency, array);

    //rec.set_antenna_array(array);

    rec.run(data);

    unsigned int index_max = rec.get_max_bin();

    auto img = rec.get_img(index_max);

    img.print();

}

