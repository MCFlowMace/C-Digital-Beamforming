/*
 * binary_classifier.cpp
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

#define DEFINE_TEMPLATE_CLASS(name, type) template class name<type>
#define DEFINE_TEMPLATES(name)\
    DEFINE_TEMPLATE_CLASS(name, float);\
    DEFINE_TEMPLATE_CLASS(name, double);

#include "binary_classifier.hpp"

template <typename input_t>
std::vector<bool> Binary_Classifier<input_t>::classify(
                                                const std::vector<input_t>& x) const
{
    std::vector<bool> label(x.size());

    for(int i=0; i<x.size(); ++i) {
        label[i] = this->classify(x[i]);
    }

    return label;
}


DEFINE_TEMPLATES(Binary_Classifier)
