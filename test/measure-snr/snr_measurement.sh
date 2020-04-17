#!/bin/bash

../../bin/measure-snr 50 100 1501 30 31 >snr_samples.dat 2>log1
../../bin/measure-snr 50 500 501 15 81 >snr_antennas.dat 2>log2

python3 plot_snr.py snr_samples.dat True
python3 plot_snr.py snr_antennas.dat False
