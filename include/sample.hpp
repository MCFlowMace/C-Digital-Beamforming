/*
 * sample.hpp
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

#include <armadillo>

template <typename value_t>
class Sample {

    public:
        Sample(arma::Col<value_t> time, arma::Col<value_t> data_time);

        arma::Col<std::complex<value_t>> frequency_data;
        arma::Col<value_t> frequency;

        arma::Col<value_t> time;

        arma::Col<value_t> time_data;

        int n_samples;
        value_t timestep;
};


template <typename value_t>
Sample<value_t>::Sample(arma::Col<value_t> time, arma::Col<value_t> time_data)
{

    this->time=time;
    this->time_data=time_data;
    n_samples=time.n_elem;
    timestep=time[1]-time[0];

    int upper, lower;

    if(!(n_samples%2)) {
        lower=-n_samples/2;
        upper=n_samples/2-1;
    } else {
        lower=-(n_samples-1)/2;
        upper=(n_samples-1)/2;
    }

    frequency = arma::Col<value_t>(upper);

    for(int i=0; i<upper; ++i) {

        frequency[i] = (i+1)/(timestep*n_samples);

    }

    arma::Col<std::complex<value_t>> fft = arma::fft(time_data);

    frequency_data = fft.subvec(1,upper+1);
    //frequency_data = fft;

}
