/*
 * reconstruction.cpp
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
 
#include "reconstruction.hpp"
#include "reconstruction_gpu.hpp"
 
#ifdef PRECISION
	typedef double value_t;
#else
	typedef float value_t;
#endif
 
#ifdef USE_GPU
	typedef Reconstruction_GPU<value_t> Rec_Type;
#else
	typedef Reconstruction_CPU<value_t> Rec_Type;
#endif

#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdlib>
#include <memory>

#include "electron.hpp"
#include "antenna.hpp"
#include "antenna_array.hpp"
#include "simulation.hpp"
#include "hpc_helpers.hpp"

#include "threshold_trigger.hpp"
#include "binary_classifier.hpp"
#include "ROC_evaluator.hpp"

#include "confusion_matrix.hpp"

int main(int argc, char **argv)
{

    if(argc!=5) {
        std::cerr << "args: [grid_size] [N_samples] [snr] [n_packets] !" << std::endl;
        exit(0);
    }

    Simulation_Settings<value_t> settings;

    int grid_size = std::atoi(argv[1]);

    settings.n_events = 1;
    settings.w_min = 2*M_PI*24.6*1e9;
    settings.w_max = 2*M_PI*26.2*1e9;
    settings.R = 5.0f;
    settings.mean_event_lifetime = 1/(2*1e-4f);
    settings.trap_efficiency = 0.5f;
    settings.run_duration = 0.005f;

    //event observation and data generation
    settings.N = 30; //antennas
    settings.snr = std::stof(argv[3]); //0.05f;
    settings.sample_rate = 3.2*1e9;
    settings.w_mix = 2*M_PI*24.6*1e9;
    settings.n_samples = std::atoi(argv[2]); //for fourier transform

    int total_packets = std::atoi(argv[4]);

    settings.run_duration = total_packets*settings.n_samples/settings.sample_rate;
    value_t dt = 1/settings.sample_rate;

    Simulation<value_t> sim(settings);
    Antenna_Array<value_t> array(settings.N, settings.R, settings.snr, settings.w_mix, settings.sample_rate);
    arma::Col<value_t> frequency = Data_Packet<value_t>::get_frequency(settings.n_samples, dt);

    int n_packets = 100;

    Rec_Type rec(grid_size, n_packets, frequency, array);
    
    std::vector<std::unique_ptr<Binary_Classifier<float>>> triggers;
    std::vector<Confusion_Matrix> cm_matrices;

    float threshold=1e3f;
    while(threshold<3e6f) {
		
		//std::cerr << "collecting triggers " << threshold << std::endl;
        triggers.emplace_back(new Threshold_Trigger<float>(threshold));
        cm_matrices.push_back(Confusion_Matrix());
        threshold +=5.f;

    }

    int packet = 0;
    while(packet<total_packets) {
		value_t t_start = packet*settings.n_samples*dt;
		//fprintf(stderr,"t_start: %18.15f\n", t_start);
		TIMERSTART(SAMPLE)
		//std::vector<std::vector<Data_Packet<value_t>>> data = sim.observation(t_start,n_packets);
		std::vector<std::complex<value_t>> data = sim.observation_flat(t_start, n_packets);
		TIMERSTOP(SAMPLE)
		
		std::vector<bool> truth(n_packets);

		TIMERSTART(GET_TRUTH)
		for(int i=0; i<truth.size(); ++i) {
			truth[i]=false;
			int packet_0 = packet+i;
			for(int j=packet_0*settings.n_samples; j<(packet_0+1)*settings.n_samples; ++j) {
			   // std::cout << i << " " << j << " " << sim.w_mat.n_rows << " " << truth.size() << std::endl;
				if(sim.w_mat(j,0)!=0.0) {
					truth[i]=true;
					break;
				}
			}
		}
		TIMERSTOP(GET_TRUTH)
		
		packet+=n_packets;
		rec.run(data);
		
		std::vector<float> test_vals;

		TIMERSTART(FIND_MAX)
		for(int j=0; j<n_packets; ++j) {
			
			unsigned int index_max = rec.get_max_bin(j);
			
			//FOR DEBUGGING - I know I should start doing this properly ...
			//~ if(index_max>frequency.n_elem) {
				//~ exit(-1);
				//~ for(int i=0; i<settings.N; ++i) {
					//~ for(int k=0; k<frequency.n_elem; ++k)
						//~ std::cerr << data[(j*settings.N+i)*frequency.n_elem+k] << std::endl;
				//~ }
			//~ }
			
			float val_max = rec.get_max_val(index_max, j);
			test_vals.push_back(val_max);

			//std::cerr << "val_max: " << val_max << std::endl;
		}
		TIMERSTOP(FIND_MAX)
		
		//std::cerr << "trigger inference with " << triggers.size() << " triggers" << std::endl;
		TIMERSTART(TRIGGER_INFERENCE)
		for(int i=0; i<triggers.size(); ++i) {
			std::vector<bool> inference = triggers[i]->classify(test_vals);
			cm_matrices[i] += Confusion_Matrix(inference,truth);
		}
		TIMERSTOP(TRIGGER_INFERENCE)
	}
	
	arma::Mat<double> curve = Confusion_Matrix::ROC_curve(cm_matrices);
	
	std::cerr << "printing" << std::endl;
    curve.t().print();
/*

    ROC_Evaluator<float> ev;

	std::cerr << "evaluating ROC curve" << std::endl;
    arma::Mat<double> curve = ev.ROC_curve(triggers, test_vals, truth);

	std::cerr << "printing" << std::endl;
    curve.t().print();
*/
}

