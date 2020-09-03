/*
 * antenna_array.cpp
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

#include "antenna_array.hpp"
#include "utility_macros.hpp"

template <typename value_t>
Antenna_Array<value_t>::Antenna_Array(int N, value_t R, value_t snr,
                                        value_t wmix, value_t sample_rate,
                                        long seed):
N(N),
R(R),
wmix(wmix)
{
	if (seed<0)
		arma::arma_rng::set_seed_random();
	else
		arma::arma_rng::set_seed(seed);

    for(int i=0; i<N; ++i) {
        value_t phi=(value_t)i/N*2*M_PI;
        value_t x=R*cos(phi);
        value_t y=R*sin(phi);
        antennas.push_back(Antenna<value_t>(snr, sample_rate, wmix, x, y));
    }
}

DEFINE_TEMPLATES(Antenna_Array)
