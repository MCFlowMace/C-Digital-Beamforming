/*
 * ROC_evaluator.hpp
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

#include <vector>
#include "binary_classifier.hpp"

class ROC_Evaluator
{
    public:

        double get_FPR();
        double get_TPR();
        const std::vector<bool>& get_inference();

        template <typename input_t>
        void evaluate(const Binary_Classifier<input_t>& classifier,
                        const std::vector<input_t>& input,
                        const std::vector<bool>& label);

    private:

        double FPR;
        double TPR;
        std::vector<bool> inference;
};
