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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */

#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdlib>

#include "beamforming/antenna_array.hpp"
#include "beamforming/reconstruction_gpu.hpp"
#include "beamforming/simulation.hpp"

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
    int N = 30;
    float R = 5.0f;
    int grid_size = std::atoi(argv[1]);

    int n_min = std::atoi(argv[2]);
    int n_max = std::atoi(argv[3]);
    int n_skip = 100;

    int n_repeat = 50;

    int N_min = std::atoi(argv[4]);
    int N_skip = 5;
    int N_max = std::atoi(argv[5]);

    for(int N=N_min; N<N_max; N+=N_skip) {
        for(int n_samples=n_min; n_samples<n_max; n_samples+=n_skip) {

            arma::Col<float> amp_max(n_repeat);
            arma::Col<float> amp_noise(n_repeat);

            for(int i=0; i<n_repeat; ++i) {

                std::cerr << N << " " << n_samples << std::endl;

                Simulation_Settings<float> settings;
                settings.N = N;
                settings.R = R;
                settings.snr = snr;
                settings.sample_rate = sample_rate;
                settings.w_mix = wmix;
                settings.n_samples = n_samples;
                settings.manual = true;
                settings.e_r = e_r;
                settings.e_phi = e_phi;
                settings.w0 = w0;
                settings.run_duration = 2.0f*n_samples/sample_rate;

                Simulation<float> sim(settings);

                Antenna_Array<float> array(N, R, snr, wmix, sample_rate);

                auto data = sim.observation_flat(0.0f, 1);

                auto frequency = Data_Packet<float>::get_frequency(
                    n_samples, 1.0f/sample_rate, false);

                Reconstruction_GPU<float> rec(
                    grid_size, 1, frequency, array, false, R);

                rec.run(data);

                unsigned int index_max = rec.get_max_bin(0);

                amp_max[i] = rec.get_max_val(index_max);

                settings.w0 = 0.0f;
                Simulation<float> noise_sim(settings);
                data = noise_sim.observation_flat(0.0f, 1);

                rec.run(data);

                amp_noise[i] = rec.get_mean_val(index_max);
            }

            float signal_mean = arma::mean(amp_max);
            float noise_mean = arma::mean(amp_noise);

            float signal_err = arma::stddev(amp_max)/sqrt(n_repeat);
            float noise_err = arma::stddev(amp_noise)/sqrt(n_repeat);

            float snr = calc_snr(signal_mean, noise_mean);
            float dsnr = calc_snr_err(
                signal_mean, noise_mean, signal_err, noise_err);

            std::cout << N << " " << n_samples << " " << snr << " "
                      << signal_err << " " << noise_err << " " << dsnr
                      << std::endl;
        }
    }
}