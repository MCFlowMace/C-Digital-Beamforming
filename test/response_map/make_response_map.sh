#!/bin/bash

../../bin/response_map 100 1000 10000 > response.out

python3 plot_result.py response.out 5

#rm response.out
