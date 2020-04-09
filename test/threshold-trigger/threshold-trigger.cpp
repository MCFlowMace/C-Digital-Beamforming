/*
 * threshold-trigger.cpp
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
#include <vector>

#include "threshold_trigger.hpp"
#include "ROC_evaluator.hpp"

int main(int argc, char **argv)
{
    std::vector<float> data {0.51, 0.1, 0.3, 0.7, 0.8, 0.6, 0.4};
    std::vector<bool> label {true, false, false, true, true, false, true};

    float threshold=0.0f;

    while(threshold<1.0f) {


        Threshold_Trigger<float> trigger(threshold);

        ROC_Evaluator ev;

        ev.evaluate(trigger, data, label);

        std::cout << ev.get_FPR() << " " << ev.get_TPR() << std::endl;

        //~ std::cout << "FPR: " << ev.get_FPR() << std::endl;
        //~ std::cout << "TPR: " << ev.get_TPR() << std::endl;

        //~ for(auto val: ev.get_inference())
            //~ std::cout << val << " ";

        //~ std::cout << std::endl;

        //~ for(auto val: label)
            //~ std::cout << val << " ";

        //~ std::cout << std::endl;
        threshold +=0.1f;
    }
}

