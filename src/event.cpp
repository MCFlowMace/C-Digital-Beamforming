/*
 * event.cpp
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

#define DW 300*1e6*M_PI //300 MHz/s
#define WK 1e-7*M_PI   //1 circumvolution every 100ns (magnetron motion)

#include "event.hpp"

template <typename value_t>
Event<value_t>::Event(value_t x0, value_t y0,
                        const std::vector<value_t>& timestamps,
                        const std::vector<value_t>& w_vals):
x0(x0),
y0(y0),
timestamps(timestamps),
w_vals(w_vals)
{

}

template <typename value_t>
value_t Event<value_t>::get_x(value_t t)
{
    return cos(WK*t)*x0 - sin(WK*t)*y0;
}

template <typename value_t>
value_t Event<value_t>::get_y(value_t t)
{
    return sin(WK*t)*x0 + cos(WK*t)*y0;
}

template <typename value_t>
value_t Event<value_t>::get_w(value_t t)
{
    value_t w0 {0};
    value_t t0 {0};

    for(int i=0; i<timestamps.size(); ++i) {
        if(t<timestamps[i]) {
            w0 = w_vals[i];
            break;
        }
        t0 = timestamps[i];
    }
    //if t>timestamps then the electron has left the trap -> w0=0
    return (t-t0)*DW + w0;
}

template class Event<float>;
template class Event<double>;
