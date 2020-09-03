
import numpy
cimport numpy

from libcpp cimport bool

cdef extern from "reconstruction_test.hpp":
    void run_test_f(int, int, float, int, bool,	float, float, float, int, float*)

def py_run(int grid_size, int n_samples,float snr, int seed, bool weighted, 
                float e_r, float e_phi, float f0, int N) -> None:
					
    n_packets = 2 #bad hardcoded!!
    bins = n_samples//2 + 1
    cdef numpy.ndarray[ndim=1, dtype=numpy.float32_t] res = numpy.zeros(grid_size*grid_size*bins*n_packets, dtype=numpy.float32)

    run_test_f(grid_size, n_samples, snr, seed, weighted, e_r, e_phi, f0, N, &res[0])

    res_ = res.reshape(2, grid_size, grid_size, bins)
    
    ind = (numpy.isfinite(res_))^True
    res_[ind] = 0
    res_[res_==-1] = 0
    
    signal = res_[0]
    noise = res_[1]

    return signal, noise
    

