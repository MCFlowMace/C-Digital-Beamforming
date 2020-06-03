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
settings(settings),
w_mat((int) (settings.sample_rate*settings.run_duration),settings.n_events)
{
    this->generation();
    this->fill_w_mat();
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
void Simulation<value_t>::fill_w_mat()
{
    int n_samples = (int) (settings.sample_rate*settings.run_duration);
    value_t dt = 1/settings.sample_rate;
    for(int i=0; i<settings.n_events; ++i) {
        for(int j=0; j<n_samples; ++j) {
            w_mat(j,i) = events[i].get_w(dt*j);
        }
    }
}

template <typename value_t>
std::vector<std::vector<Data_Packet<value_t>>> Simulation<value_t>::observation(
                                                value_t t_start, int n_packets)
{

    Antenna_Array<value_t> array(settings.N, settings.R, settings.snr,
                                settings.w_mix, settings.sample_rate);

    std::vector<std::vector<Data_Packet<value_t>>> data(settings.N);

    value_t dt = 1/settings.sample_rate;
    //value_t delta_t = t_end-t_start;
    //int samples = (int) (delta_t*settings.sample_rate);
    //int n_packets = samples/settings.n_samples; //only take full packets of data
    //SDIV(samples, settings.n_samples);

#if defined PARALLEL || defined USE_GPU
    #pragma omp parallel
    {
        //std::cout << "threads: " << omp_get_num_threads() << std::endl;
    #pragma omp for
#endif
    for(int i=0; i<settings.N; ++i) {
        //std::cout << i << " " << omp_get_thread_num() << std::endl;
        value_t t {t_start};
        std::vector<Data_Packet<value_t>> data_i(n_packets);
        //std::vector<Sample<value_t>> data_i;
        for(int j=0; j<n_packets; ++j) {
            //data_i[j] = array.antennas[i].sample_data(samples, t, this->events);
            //data_i.push_back(std::move(array.antennas[i].sample_data(settings.n_samples, t, this->events)));
            data_i[j] = std::move(array.antennas[i].sample_packet(settings.n_samples, t, this->events));
            //std::cout << "i, j " << i << ", " << j << std::endl;
            t=t_start+(j+1)*dt*settings.n_samples;
        }
        //std::cout << "addr observation: " << &data_i[0].time_data[0] << std::endl;
        //data.push_back(data_i);
        data[i] = std::move(data_i);
        //std::cout << "addr observation 2: " << &data[i][0].time_data[0] << std::endl;
        
        //if(i==0)
		//	fprintf(stderr,"last t: %18.15f next t: %18.15f\n", t, t+dt);
    }

#if defined PARALLEL || defined USE_GPU
    }
#endif

    return data;
}

template <typename value_t>
std::vector<std::complex<value_t>> Simulation<value_t>::observation_flat(
                                                value_t t_start, int n_packets)
{

    Antenna_Array<value_t> array(settings.N, settings.R, settings.snr,
                                settings.w_mix, settings.sample_rate);
                                

    value_t dt = 1/settings.sample_rate;
    
    int bins = Data_Packet<value_t>::get_frequency(settings.n_samples, dt).n_elem;
    std::vector<std::complex<value_t>> data(settings.N*n_packets*bins);
    
    //std::cerr << "t start: " << t_start << std::endl;
    
    //value_t delta_t = t_end-t_start;
    //int samples = (int) (delta_t*settings.sample_rate);
    //int n_packets = samples/settings.n_samples; //only take full packets of data
    //SDIV(samples, settings.n_samples);

#if defined PARALLEL || defined USE_GPU
    #pragma omp parallel
    {
        //std::cout << "threads: " << omp_get_num_threads() << std::endl;
    #pragma omp for
#endif
    for(int i=0; i<n_packets; ++i) {

        value_t t = t_start+i*dt*settings.n_samples;

        for(int j=0; j<settings.N; ++j) {
			
			//std::cerr << "t: " << t << std::endl;

            arma::Col<std::complex<value_t>> data_j = std::move(array.antennas[j].sample_data(settings.n_samples, t, this->events));
            
            //improve this stupid copy later (as if I'll ever get around doing that ...)
            for(int k=0; k<bins; ++k)
				data[(i*settings.N+j)*bins+k] = data_j[k];
				

        }

    }

#if defined PARALLEL || defined USE_GPU
    }
#endif

    return data;
}

template class Simulation<float>;
template class Simulation<double>;
