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
#include "antenna_array.hpp"
#include "grid.hpp"
#include "hpc_helpers.hpp"

template <typename value_t>
class Reconstruction {


    public:

        //arma::Mat<value_t> img;

        //arma::field<arma::Cube<std::complex<value_t>>> grid_phase;

        Reconstruction(int grid_size, arma::Col<value_t> frequency,
                        const Antenna_Array<value_t>& array);

        ~Reconstruction();

        Reconstruction(const Reconstruction& temp_obj) = delete;
        Reconstruction& operator=(const Reconstruction& temp_obj) = delete;

        void run(const std::vector<Data_Packet<value_t>>& samples);

        arma::Mat<value_t> get_img(unsigned int bin);

        unsigned int get_max_bin();
        value_t get_max_val(unsigned int bin);
        value_t get_mean_val(unsigned int bin);

    private:

        void set_antenna_array(const Antenna_Array<value_t>& array);

        int grid_size;
        Grid<value_t> grid;
        int N;
        value_t R;

        std::complex<value_t>* grid_phase;
        arma::Col<value_t> frequency;
        arma::Cube<value_t> reconstructed;
};

