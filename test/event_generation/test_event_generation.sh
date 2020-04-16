#!/bin/bash

tmin=0.05
tmax=0.1
fmin=24.6
fmax=26.2
R=5.0

../../bin/generate_events $tmin $tmax $fmin $fmax $R

python plot_event.py event0_freq.dat $tmin $tmax $fmin $fmax event0_pos.dat $R

rm event0_freq.dat event0_pos.dat
