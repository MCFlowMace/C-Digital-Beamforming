/*
 * simulation.cpp
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


#include "simulation.hpp"
#include "event_generator.hpp"
#include "antenna_array.hpp"
#include "hpc_helpers.hpp"

template <typename value_t>
Simulation<value_t>::Simulation(Simulation_Settings<value_t> settings):
settings(settings)
{
    this->generation();
}

template <typename value_t>
void Simulation<value_t>::generation()
{
    Event_Generator<value_t> gen(settings.mean_event_lifetime,
                                    settings.trap_efficiency);

    for(int i=0; i<settings.n_events; ++i)
        events.push_back(gen.generate(value_t{0}, settings.run_duration,
                                        settings.w_min, settings.w_max,
                                        settings.R));
}

template <typename value_t>
std::vector<std::vector<Sample<value_t>>> Simulation<value_t>::observation(
                                                value_t t_start, value_t t_end)
{

    Antenna_Array<value_t> array(settings.N, settings.R, settings.snr,
                                settings.w_mix, settings.sample_rate);

    std::vector<std::vector<Sample<value_t>>> data(settings.N);

    //std::cout << "sample_rate: " << settings.sample_rate << " t_end " << t_end << std::endl;

    value_t dt = 1/settings.sample_rate;
    value_t delta_t = t_end-t_start;
    int samples = (int) (delta_t*settings.sample_rate);
    int n_packets = samples/settings.n_samples; //only take full packets of data
    //SDIV(samples, settings.n_samples);

    //std::cerr << "delta_t " << delta_t << " run time " << settings.run_duration << std::endl;
    //std::cerr << "samples: " << samples << " n_packets: " << n_packets << std::endl;

#ifdef PARALLEL
    #pragma omp parallel
    {
        //std::cout << "threads: " << omp_get_num_threads() << std::endl;
    #pragma omp for
#endif
    for(int i=0; i<settings.N; ++i) {
        //std::cout << i << " " << omp_get_thread_num() << std::endl;
        value_t t {t_start};
        //std::vector<Sample<value_t>> data_i(n_packets);
        std::vector<Sample<value_t>> data_i;
        for(int j=0; j<n_packets; ++j) {
            //data_i[j] = array.antennas[i].sample_data(samples, t, this->events);
            data_i.push_back(array.antennas[i].sample_data(settings.n_samples, t, this->events));
            t+=dt*samples;

            //std::cout << "i, j " << i << ", " << j << std::endl;
        }
        //data.push_back(data_i);
        data[i] = data_i;
    }

#ifdef PARALLEL
    }
#endif

    return data;
}

template class Simulation<float>;
template class Simulation<double>;
