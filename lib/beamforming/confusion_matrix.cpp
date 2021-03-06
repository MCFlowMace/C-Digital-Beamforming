/*
 * confusion_matrix.cpp
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


#include "beamforming/confusion_matrix.hpp"

#include <stdexcept>
#include <iostream>

Confusion_Matrix::Confusion_Matrix(): 
true_positives{0}, 
false_positives{0}, 
positive_labels{0}, 
negative_labels{0}
{}

Confusion_Matrix::Confusion_Matrix(const std::vector<bool>& predictions,
									const std::vector<bool>& labels):
true_positives{0}, 
false_positives{0}, 
positive_labels{0}, 
negative_labels{0}
{
    if(predictions.size()!=labels.size())
        throw std::invalid_argument("predictions and labels must have same length!");

    for(unsigned int i=0; i<predictions.size(); ++i) {
		this->add_example(predictions[i], labels[i]);
    }
}

void Confusion_Matrix::add_example(bool prediction, bool label)
{
	if(label) {
		this->positive_labels++;

		if(prediction)
			this->true_positives++;
	} else {
		this->negative_labels++;

		if(prediction)
			this->false_positives++;
	}
}

double Confusion_Matrix::get_TPR() const
{
	//std::cerr << "tp: " << this->true_positives << " pl: " << this->positive_labels << std::endl;
	return (double) this->true_positives/this->positive_labels;
}

double Confusion_Matrix::get_FPR() const
{
	//std::cerr << "fp: " << this->false_positives << " nl: " << this->negative_labels << std::endl;
	return (double) this->false_positives/this->negative_labels;
}

double Confusion_Matrix::get_DTPR() const
{
	//std::cerr << "tp: " << this->true_positives << " pl: " << this->positive_labels << std::endl;
	
	double TP = this->true_positives;
	double P = this->positive_labels;
	
	//uncertainty
	double DTP = sqrt(this->true_positives);
	double DP = sqrt(this->positive_labels);
	
	double first = DTP/P;
	double second = DP*TP/(P*P);
	return (double) sqrt(first*first + second*second);
}

double Confusion_Matrix::get_DFPR() const
{
	//std::cerr << "fp: " << this->false_positives << " nl: " << this->negative_labels << std::endl;
	double FP = this->false_positives;
	double N = this->negative_labels;
	
	//uncertainty
	double DFP = sqrt(this->false_positives);
	double DN = sqrt(this->negative_labels);
	
	double first = DFP/N;
	double second = DN*FP/(N*N);
	return (double) sqrt(first*first + second*second);
}


void Confusion_Matrix::print()
{
	std::cout 	<< " tp: " << this->true_positives 
				<< " pl: " << this->positive_labels 
				<< " fp: " << this->false_positives
				<< " nl: " << this->negative_labels
				<< " TPR: " << this->get_TPR()
				<< " FPR: " << this->get_FPR() << std::endl;
}

Confusion_Matrix& Confusion_Matrix::operator +=(const Confusion_Matrix& rhs)
{
	this->true_positives += rhs.true_positives;
	this->false_positives += rhs.false_positives;
	this->positive_labels += rhs.positive_labels;
	this->negative_labels += rhs.negative_labels;
	
	return *this;
}

arma::Mat<double> Confusion_Matrix::ROC_curve(
                    const std::vector<Confusion_Matrix>& cm_matrices)
{
    arma::Mat<double> curve(4, cm_matrices.size());

    for(int i=0; i<cm_matrices.size(); ++i) {
        curve(0, i) = cm_matrices[i].get_FPR();
        curve(1, i) = cm_matrices[i].get_DFPR();
        curve(2, i) = cm_matrices[i].get_TPR();
        curve(3, i) = cm_matrices[i].get_DTPR();
    }

    return curve;
}

arma::Mat<double> Confusion_Matrix::to_arma(
                    const std::vector<Confusion_Matrix>& cm_matrices)
{
    arma::Mat<double> curve(4, cm_matrices.size());

    for(int i=0; i<cm_matrices.size(); ++i) {
        curve(0, i) = cm_matrices[i].false_positives;
        curve(1, i) = cm_matrices[i].negative_labels;
        curve(2, i) = cm_matrices[i].true_positives;
        curve(3, i) = cm_matrices[i].positive_labels;
    }

    return curve;
}
