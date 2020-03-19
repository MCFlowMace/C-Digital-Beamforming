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

int main(int argc, char **argv)
{
    //float lambda {1/(5*1e-4f)}; //mean lifetime of 500us
    float lambda {1/(5*1e-3f)}; //mean lifetime of 5ms
    float trap_efficiency {0.75f}; //trapping probability

    //Event_Generator<float> gen(lambda, trap_efficiency, 12351);
    Event_Generator<float> gen(lambda, trap_efficiency);

    float t_max = 1.0f;
    Event<float> event0 = gen.generate(t_max, 2*M_PI*1.6);

    std::cout << "#scattered: " << event0.get_n_scatter() << " times!" << std::endl;

    int N=100000;
    float dt = t_max/N;
    for(int i=0; i<N; ++i) {
        float t = i*dt;
        //std::cout << t << " " << event0.get_w(t) << "\n";
        printf("%20.10f %20.10f\n",t, event0.get_w(t));
    }
}

