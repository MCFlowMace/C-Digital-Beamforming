/*
 * reconstruction.hpp
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

#pragma once

#include <armadillo>
#include "beamforming/antenna_array.hpp"
#include "beamforming/grid.hpp"
#include "utility/hpc_helpers.hpp"
#include "beamforming/reconstruction.hpp"

template <typename value_t>
class Reconstruction_GPU : public Reconstruction<value_t> {


    public:

        Reconstruction_GPU(int grid_size, int n_packets,
						arma::Col<value_t> frequency,
                        const Antenna_Array<value_t>& array, bool weighted);

        ~Reconstruction_GPU();

        Reconstruction_GPU(const Reconstruction_GPU& temp_obj) = delete;
        Reconstruction_GPU& operator=(const Reconstruction_GPU& temp_obj) = delete;

        //virtual void run(const std::vector<std::vector<Data_Packet<value_t>>>& samples)=0;
        virtual void run(const std::vector<std::complex<value_t>>& samples);
		virtual arma::Mat<value_t> get_img(unsigned int packet, 
											unsigned int bin);
        virtual unsigned int get_max_bin(unsigned int packet);
        
        void print(unsigned int packet);
        void copy_res(value_t* dest);

    private:

        value_t* time_delays_D;
        value_t* phis_D;
        value_t* frequencies_D;
        value_t* reconstructed_D;
        std::complex<value_t>* samples_D;
        
        std::complex<value_t>* samples_H;
        value_t* reconstructed_H;
        
        void init_gpu();
};

