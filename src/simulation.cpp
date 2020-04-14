/*
 * simulation.cpp
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


#include "simulation.hpp"
#include "event_generator.hpp"

template <typename value_t>
Simulation<value_t>::Simulation(Simulation_Settings<value_t> settings):
settings(settings)
{

}

template <typename value_t>
void Simulation<value_t>::event_generation()
{
    Event_Generator<value_t> gen(settings.mean_event_lifetime,
                                    settings.trap_efficiency);

    for(int i=0; i<settings.n_events; ++i)
        events.push_back(gen.generate(value_t{0}, settings.run_duration,
                                        settings.w_min, settings.w_max,
                                        settings.R));
}

template <typename value_t>
void Simulation<value_t>::event_observation()
{
    //~ float dt = t_max/N; //dt=10us
    //~ for(int i=0; i<N; ++i) {
        //~ float t = i*dt;
        //~ event0.get_w(t);
    //~ }

    //~ for(int i=0; i<N; ++i) {
        //~ float t = i*dt;
        //~ if(event0.get_w(t)!=0.0)
            //~ event0.get_y(t);
    //~ }
}

template class Simulation<float>;
template class Simulation<double>;
