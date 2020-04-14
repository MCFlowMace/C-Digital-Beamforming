/*
 * event_generator.hpp
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

#include <random>
#include "event.hpp"

template <typename value_t>
class Event_Generator
{

    public:

        Event_Generator(value_t lambda, value_t trap_efficiency, long seed=-1);
        Event<value_t> generate(value_t t_min, value_t t_max, value_t w_min,
                                        value_t w_max, value_t R);

    private:

        std::mt19937 generator;

        value_t lambda;
        value_t trap_efficiency; //pitch angles

        value_t generate_r0(value_t R);
        value_t generate_phi0();
        value_t get_E_loss();
        value_t generate_t0(value_t t_min, value_t t_max);
        value_t generate_w(value_t w_min, value_t w_max);

        value_t get_BW_val();
        value_t inv_BW(value_t x, value_t t, value_t s);

        value_t next_timestamp(value_t t_old);
        value_t new_frequency(value_t t, value_t t_old, value_t w_old, value_t w_max);
        bool has_left_trap();

};
