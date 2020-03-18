/*
 * event.hpp
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

#include <vector>
#include <cmath>


template <typename value_t>
class Event
{

    public:

        Event(value_t x0, value_t y0, const std::vector<value_t>& timestamps,
                const std::vector<value_t>& w_vals);

        value_t get_x(value_t t);
        value_t get_y(value_t t);
        value_t get_w(value_t t);

    private:

        value_t x0;
        value_t y0;
        //stores the times of the subevents [emerged, scatter_0, scatter_1 ..., dissappeared]
        std::vector<value_t> timestamps;
        std::vector<value_t> w_vals; //corresponding frequencies

};
