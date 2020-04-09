/*
 * ROC_evaluator.cpp
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

#include <stdexcept>

#include "ROC_evaluator.hpp"


double evaluate_TPR(const std::vector<bool>& test,
                    const std::vector<bool>& truth)
{
    if(test.size()!=truth.size())
        throw std::invalid_argument("test and truth must have same length!");

    int count_TP {0};
    int count_P {0};

    for(int i=0; i<test.size(); ++i) {
        if(truth[i]) {
            count_P++;

            if(test[i])
                count_TP++;
        }
    }

    return (double) count_TP/count_P;
}

double evaluate_FPR(const std::vector<bool>& test,
                    const std::vector<bool>& truth)
{
    if(test.size()!=truth.size())
        throw std::invalid_argument("test and truth must have same length!");

    std::vector<bool> inv_truth(test.size());

    for(int i=0; i<test.size(); ++i)
        inv_truth[i] = !truth[i];

    return evaluate_TPR(test, inv_truth);
}

double ROC_Evaluator::get_FPR()
{
    return FPR;
}

double ROC_Evaluator::get_TPR()
{
    return TPR;
}

const std::vector<bool>& ROC_Evaluator::get_inference()
{
    return this->inference;
}

template <typename input_t>
void ROC_Evaluator::evaluate(const Binary_Classifier<input_t>& classifier,
                                const std::vector<input_t>& input,
                                const std::vector<bool>& label)
{
    this->inference = classifier.classify(input);
    this->FPR = evaluate_FPR(this->inference, label);
    this->TPR = evaluate_TPR(this->inference, label);
}

template void ROC_Evaluator::evaluate<double>(const Binary_Classifier<double>&,
                                                const std::vector<double>&,
                                                const std::vector<bool>&);

template void ROC_Evaluator::evaluate<float>(const Binary_Classifier<float>&,
                                                const std::vector<float>&,
                                                const std::vector<bool>&);
