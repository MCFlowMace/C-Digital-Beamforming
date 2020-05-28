/*
 * antenna.hpp
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
#include <cmath>
#include <vector>

#include "electron.hpp"
#include "data_packet.hpp"
#include "event.hpp"

template <typename value_t>
class Antenna {

    public:

        Antenna(value_t snr, value_t sample_rate, value_t wmix,
                value_t x, value_t y);

        Data_Packet<value_t> sample_packet(int n_samples, value_t t0,
                                    const std::vector<Event<value_t>>& events);
		arma::Col<std::complex<value_t>> sample_data(int n_samples, value_t t0,
                                    const std::vector<Event<value_t>>& events);

        value_t x;
        value_t y;


    private:

        value_t snr;
        value_t sample_rate;
        value_t w_mix;

        value_t time_delay(const Event<value_t>& event, value_t t);
        value_t spiral_phase(const Event<value_t>& event, value_t t);
        value_t distance_to_electron(const Event<value_t>& event, value_t t);

        arma::Col<value_t> sample_noise(int n_samples);
        arma::Col<value_t> get_mixed_sample(
                                    const arma::Col<value_t>& t,
                                    const std::vector<Event<value_t>>& events);

};

