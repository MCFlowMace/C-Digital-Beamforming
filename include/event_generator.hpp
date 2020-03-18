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

#include <random>
#include "event.hpp"

template <typename value_t>
class Event_Generator
{

    public:

    Event_Generator();

    Event<value_t> generate(value_t t_max, value_t w_max);

    private:

    std::mt19937 generator;

    value_t lambda;
    value_t trap_efficiency; //pitch angles

    value_t generate_t0(value_t t_max);
    value_t generate_w0(value_t w_max);
    value_t next_timestamp(value_t t_old);
    value_t new_frequency(value_t w_old);
    bool has_left_trap();

};
