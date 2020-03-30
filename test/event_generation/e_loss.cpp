/*
 * generate_events.cpp
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


#include <iostream>
#include "event.hpp"
#include "event_generator.hpp"
#include <cmath>
#include <vector>
#include <stdio.h>

#define float double

int main(int argc, char **argv)
{

    //float lambda {1/(5*1e-4f)}; //mean lifetime of 500us
    float lambda {1/(2*1e-4f)}; //mean lifetime of 0.2ms
    float trap_efficiency {0.5f}; //trapping probability

    //Event_Generator<float> gen(lambda, trap_efficiency, 12351);
    Event_Generator<float> gen(lambda, trap_efficiency);

    int N=100000;


    for(int i=0; i<N; ++i) {
        //fprintf(pos_data, "%20.10f %20.10f\n",event0.get_x(t), event0.get_y(t));
        float E_loss = gen.get_E_loss();
        std::cout << E_loss << "\n";
    }
}

