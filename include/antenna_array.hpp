/*
 * antenna_array.hpp
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

#include "antenna.hpp"
#include <cmath>

template <typename value_t>
class Antenna_Array {

    public:
        int N;
        value_t R;
        std::vector<Antenna<value_t>> antennas;
        value_t wmix;

        Antenna_Array(int N, value_t R, value_t snr,
                                        value_t wmix, value_t sample_rate);
};


