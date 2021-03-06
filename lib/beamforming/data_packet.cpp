/*
 * data_packet.cpp
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

#include "beamforming/data_packet.hpp"
#include "utility/utility_macros.hpp"

template <typename value_t>
Data_Packet<value_t>::Data_Packet(arma::Col<value_t>&& time,
                                    arma::Col<value_t>&& time_data):
time(std::move(time)),
time_data(std::move(time_data))
{

   // this->time=time;
   // this->time_data=time_data;

    n_samples=this->time.n_elem;
    timestep=this->time[1]-this->time[0];

	frequency=get_frequency(n_samples, timestep, false);
	
	int upper = frequency.n_elem;

    frequency_data = arma::fft(this->time_data);
    
}

template <typename value_t>
arma::Col<value_t> Data_Packet<value_t>::get_frequency(int n_samples, 
							value_t dt,
							bool full_frequency) {
	
    if(full_frequency)
	return get_frequency_new(n_samples, dt);
    else
	return get_frequency_old(n_samples, dt);
}

template <typename value_t>
arma::Col<value_t> Data_Packet<value_t>::get_frequency_old(int n_samples, 
							value_t dt) {
	
    int n = n_samples/2 + 1; //integer division intended

	//see scipy.fft.rfftfreq    
    

    arma::Col<value_t> _frequency(n);

    for(int i=0; i<n; ++i)
        _frequency[i] = i/(dt*n_samples);

	return _frequency;
}

template <typename value_t>
arma::Col<value_t> Data_Packet<value_t>::get_frequency_new(int n_samples, 
							value_t dt) {

	//see scipy.fft.rfftfreq   
    int start, end;
    
    if(n_samples%2 == 0) {
	start = -n_samples/2;
	end = n_samples/2;
    } else {
	start = -(n_samples-1)/2;
	end = (n_samples-1)/2+1;
    }
    
   // std::cout << start << " " << end << " " << n_samples << std::endl;

    arma::Col<value_t> _frequency(n_samples);

    for(int i=start; i<end; ++i)
        _frequency[i-start] = i/(dt*n_samples);

    //_frequency.print();
    
    return _frequency;
}

DEFINE_TEMPLATES(Data_Packet)
