/*
 * data_packet.hpp
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
class Data_Packet {

    public:

        Data_Packet(arma::Col<value_t>&& time, arma::Col<value_t>&& data_time);
        Data_Packet() = default;

        arma::Col<std::complex<value_t>> frequency_data;
        arma::Col<value_t> frequency;
        
        static arma::Col<value_t> get_frequency(int n_samples, value_t dt,
                                                bool full_frequency);

        arma::Col<value_t> time;

        arma::Col<value_t> time_data;

        int n_samples;
        value_t timestep;
        
    private:
        static arma::Col<value_t> get_frequency_old(int n_samples, value_t dt);
        static arma::Col<value_t> get_frequency_new(int n_samples, value_t dt);
};

