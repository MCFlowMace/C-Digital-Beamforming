#!/bin/bash

../../bin/reconstruction 100 1000 0.1 -1 1 > beamforming_rec_test.out

python3 plot_result.py beamforming_rec_test.out 5

rm beamforming_rec_test.out
