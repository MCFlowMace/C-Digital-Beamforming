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

 /** currently not working, needs overhaul! **/


#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdlib>

#include "electron.hpp"
#include "antenna.hpp"
#include "antenna_array.hpp"
#include "reconstruction.hpp"
#include "hpc_helpers.hpp"

float calc_snr(float A, float N)
{
    return A*A/(N*N);
}

float calc_snr_err(float A, float N, float dA, float dN)
{
    float N2 = N*N;
    float val1 = dA*2*A/N2;
    float val2 = dN*2*A*A/(N2*N);
    return sqrt(val1*val1 + val2*val2);
}

int main(int argc, char **argv)
{

    if(argc!=6) {
        std::cerr << "args: [grid_size] [n_min] [n_max] [N_min] [N_max] !" << std::endl;
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

    int n_min=std::atoi(argv[2]);
    int n_max=std::atoi(argv[3]);
    int n_skip=100;

    int n_repeat = 50;

    int N_min=std::atoi(argv[4]);
    int N_skip=5;
    int N_max = std::atoi(argv[5]);

    for(int N=N_min; N<N_max; N+= N_skip) {
        for(int n_samples=n_min; n_samples<n_max; n_samples+=n_skip) {

            arma::Col<float> amp_max (n_repeat);
            arma::Col<float> amp_noise (n_repeat);

            for(int i=0; i<n_repeat; ++i) {

                std::cerr << N << " " << n_samples << std::endl;

                auto emission = [](float t, float w, float phi) -> float {
                                    float dw = 300*1e6;
                                    return cosf((w+dw*t)*t+phi);
                                };

                Electron<float> e {w0, e_r*cosf(e_phi), e_r*sinf(e_phi), emission};


                Antenna_Array<float> array(N, R, snr, wmix, sample_rate);

                std::vector<Data_Packet<float>> data;


                for(int i=0; i<N; ++i) {
                    data.push_back(array.antennas[i].sample_data(n_samples, 0.0f, e));
                }

                Reconstruction<float> rec(grid_size, data[0].frequency, array);

                //rec.set_antenna_array(array);

                rec.run(data);

                unsigned int index_max = rec.get_max_bin();

                amp_max[i] = rec.get_max_val(index_max);

                Electron<float> e2 {w0, e_r*cosf(e_phi), e_r*sinf(e_phi), [](float t, float w, float phi) -> float {return 0; }};

                data.clear();

                for(int i=0; i<N; ++i) {
                    data.push_back(array.antennas[i].sample_data(n_samples, 0.0f, e2));
                }

                rec.run(data);

                amp_noise[i] = rec.get_mean_val(index_max);

                //std::cout << n_samples << " " << snr << std::endl;
            }

            float signal_mean = arma::mean(amp_max);
            float noise_mean = arma::mean(amp_noise);

            float signal_err = arma::stddev(amp_max)/sqrt(n_repeat);
            float noise_err = arma::stddev(amp_noise)/sqrt(n_repeat);

            float snr =calc_snr(signal_mean, noise_mean);
            float dsnr = calc_snr_err(signal_mean, noise_mean, signal_err, noise_err);

            std::cout << N << " " << n_samples << " " << snr << " " << signal_err << " " << noise_err << " " << dsnr << std::endl;
        }
    }
}

