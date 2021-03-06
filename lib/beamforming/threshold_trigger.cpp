/*
 * threshold_trigger.cpp
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


#include "beamforming/threshold_trigger.hpp"
#include <iostream>

template <typename value_t>
Threshold_Trigger<value_t>::Threshold_Trigger(value_t threshold):
threshold {threshold}
{

}

template <typename value_t>
bool Threshold_Trigger<value_t>::classify(value_t x) const
{
	//std::cerr << "val: " << x << "threshold: " << threshold << std::endl;
    return x>threshold;
}

template class Threshold_Trigger<float>;
template class Threshold_Trigger<double>;

