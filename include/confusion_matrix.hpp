/*
 * confusion_matrix.hpp
 * 
 * Copyright 2020 Florian Thomas <flthomas@students.uni-mainz.de>
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
#include <cstdint>
#include <armadillo>

class Confusion_Matrix {
	
	public:
		
		Confusion_Matrix();
		Confusion_Matrix& operator+=(const Confusion_Matrix& rhs);
		Confusion_Matrix(const std::vector<bool>& test,
							const std::vector<bool>& truth);
							
		double get_FPR() const;
		double get_TPR() const;
		
		void print();
		
		static arma::Mat<double> ROC_curve(
                    const std::vector<Confusion_Matrix>& cm_matrices);
		
	private:
	
		uint64_t true_positives;
		uint64_t false_positives;
		uint64_t positive_labels;
		uint64_t negative_labels;
		
		void add_example(bool prediction, bool label);
};
