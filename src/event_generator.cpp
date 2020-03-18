/*
 * event_generator.cpp
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

#include "event_generator.hpp"


template <typename value_t>
Event_Generator<value_t>::Event_Generator():
generator(42)
{
}

//TODO generate position
template <typename value_t>
Event<value_t> Event_Generator<value_t>::generate(value_t t_max, value_t w_max)
{
    std::vector<value_t> timestamps;
    std::vector<value_t> w_vals;

    timestamps.push_back(generate_t0(t_max));
    w_vals.push_back(value_t(0));
    w_vals.push_back(generate_w0(w_max));

    do {
        int N=timestamps.size();
        timestamps.push_back(next_timestamp(timestamps[N-1]));

        if(!has_left_trap())
            break;

        w_vals.push_back(new_frequency(w_vals[N]));
    } while(true);

    w_vals.push_back(value_t(0));

    return Event<value_t>(value_t(0), value_t(0), timestamps, w_vals);
}

template <typename value_t>
value_t Event_Generator<value_t>::generate_t0(value_t t_max)
{
    std::uniform_real_distribution<value_t> dis(0.0, 1.0);
    value_t rand_val = dis(generator);

    return rand_val*t_max;
}

template <typename value_t>
value_t Event_Generator<value_t>::generate_w0(value_t w_max)
{
    std::uniform_real_distribution<value_t> dis(0.0, 1.0);
    value_t rand_val = dis(generator);

    return rand_val*w_max;
}

template <typename value_t>
bool Event_Generator<value_t>::has_left_trap()
{
    std::uniform_real_distribution<value_t> dis(0.0, 1.0);
    value_t rand_val = dis(generator);

    //assuming isotropic scattering
    return rand_val<trap_efficiency;
}

template <typename value_t>
value_t Event_Generator<value_t>::next_timestamp(value_t t_old)
{
    std::exponential_distribution<value_t> dis(lambda);

    return t_old + dis(generator);
}

template <typename value_t>
value_t Event_Generator<value_t>::new_frequency(value_t w_old)
{
    std::exponential_distribution<value_t> dis(lambda);

    return w_old + dis(generator);
}

template class Event_Generator<float>;
template class Event_Generator<double>;
