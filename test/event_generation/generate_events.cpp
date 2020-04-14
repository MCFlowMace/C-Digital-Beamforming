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
#include "hpc_helpers.hpp"

//#define float double

int main(int argc, char **argv)
{
    if(argc<6) {
        std::cerr << "Usage: [tmin] [tmax] [fmin] [fmax] [R]" << std::endl;
        return -1;
    }

    //float lambda {1/(5*1e-4f)}; //mean lifetime of 500us
    float lambda {1/(2*1e-4f)}; //mean lifetime of 0.2ms
    float trap_efficiency {0.5f}; //trapping probability

    //Event_Generator<float> gen(lambda, trap_efficiency, 12351);
    Event_Generator<float> gen(lambda, trap_efficiency);

    float t_min = std::atof(argv[1]); //0.05f;
    float t_max = std::atof(argv[2]); //0.1f;
    float w_min = 2*M_PI*std::atof(argv[3])*1e9; //24.6;
    float w_max = 2*M_PI*std::atof(argv[4])*1e9; //26.2;
    float R = std::atof(argv[5]); //5.0f;

    int N=100000;

    //~ TIMERSTART(GEN)
    //~ for(int i=0; i<N; ++i)
        //~ gen.generate(t_min, t_max, w_min, w_max, R);
    //~ TIMERSTOP(GEN)

    Event<float> event0 = gen.generate(t_min, t_max, w_min, w_max, R);

    std::cerr << "#scattered: " << event0.get_n_scatter() << " times!" << std::endl;
    FILE* freq_data;

    freq_data = fopen("event0_freq.dat", "w+");

    float dt = t_max/N; //dt=10us
    for(int i=0; i<N; ++i) {
        float t = i*dt;
        //std::cout << t << " " << event0.get_w(t) << "\n";
        fprintf(freq_data, "%20.10f %20.10f\n",t, event0.get_w(t));
    }
    fclose(freq_data);

    FILE* pos_data;

    pos_data = fopen("event0_pos.dat", "w+");

    for(int i=0; i<N; ++i) {
        float t = i*dt;
        /*float r = gen.generate_r0(R);
        float phi = gen.generate_phi0();
        float x = cos(phi)*r;
        float y = sin(phi)*r;*/
        if(event0.get_w(t)!=0.0)
            fprintf(pos_data, "%20.10f %20.10f\n",event0.get_x(t), event0.get_y(t));
        //fprintf(pos_data, "%20.10f %20.10f\n",x, y);
    }
    fclose(pos_data);

}

