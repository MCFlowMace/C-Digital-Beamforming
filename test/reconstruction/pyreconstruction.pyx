
from libcpp cimport bool

cdef extern from "reconstruction_test.h":
    void run_test_f(int, int, float, int, bool,	float, float, float, int)

def py_run(int grid_size, int n_samples,float snr, int seed, bool weighted, 
				float e_r, float e_phi, float f0, int N) -> None:
    run_test_f(grid_size, n_samples, snr, seed, weighted, e_r, e_phi, f0, N)
