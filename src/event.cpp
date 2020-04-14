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

#define DW 0.3*2*M_PI //300 MHz/s
#define WK 2*M_PI/(4*1e-3)   //1 circumvolution every 4ms (magnetron motion)

#include "event.hpp"

#include <utility>
#include <iostream>

template <typename value_t>
Event<value_t>::Event(value_t r0, value_t phi0,
                        std::vector<value_t>&& timestamps,
                        std::vector<value_t>&& w_vals):
x0(cos(phi0)*r0),
y0(sin(phi0)*r0),
n_scatter(timestamps.size()-1),
timestamps(std::move(timestamps)),
w_vals(std::move(w_vals))
{
}

template <typename value_t>
value_t Event<value_t>::get_x(value_t t) const
{
    return cos(WK*t)*x0 - sin(WK*t)*y0;
}

template <typename value_t>
value_t Event<value_t>::get_y(value_t t) const
{
    return sin(WK*t)*x0 + cos(WK*t)*y0;
}

template <typename value_t>
value_t Event<value_t>::get_w(value_t t) const
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

    value_t w {0};

    if(w0!=0) {
        w=calc_w(t, t0, w0);
        //std::cout << "calcw!" << std::endl;
    }

    return w;
}

template <typename value_t>
value_t Event<value_t>::calc_w(value_t t, value_t t0, value_t w0)
{
    return (t-t0)*DW + w0;
}

template <typename value_t>
value_t Event<value_t>::get_n_scatter() const
{
    return n_scatter;
}

template class Event<float>;
template class Event<double>;
