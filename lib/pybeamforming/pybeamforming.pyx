
import numpy
cimport numpy
import ctypes

from libcpp cimport bool
from libcpp.memory cimport make_shared, shared_ptr

cdef extern from "beamforming/simulation.hpp":
    cdef cppclass Simulation_Settingsf:
        Simulation_Settingsf() except +
		
        int n_events
        float w_min
        float w_max
        float R
        float run_duration
        float mean_event_lifetime
        float trap_efficiency
        long seed

        bool manual
        float e_r
        float e_phi
        float w0

        int N
        float snr
        float sample_rate
        float w_mix
        int n_samples
        
    cdef cppclass Simulation_Settingsd:
        Simulation_Settingsd() except +
		
        int n_events
        double w_min
        double w_max
        double R
        double run_duration
        double mean_event_lifetime
        double trap_efficiency
        long seed

        bool manual
        double e_r
        double e_phi
        double w0

        int N
        double snr
        double sample_rate
        double w_mix
        int n_samples
        
cdef extern from "beamforming/beamformer.hpp":
    cdef cppclass Beamformerf:
        Beamformerf(Simulation_Settingsf settings, int grid_size, 
                    int n_packets, bool weighted, bool full_frequency,
                    float r_grid) except +

        void get_next(float* dest)
        void get_result(float complex* src, float* dest)
        
    cdef cppclass Beamformerd:
        Beamformerf(Simulation_Settingsd settings, int grid_size, 
                    int n_packets, bool weighted, bool full_frequency,
                    double r_grid) except +

        void get_next(double* dest)
        void get_result(double complex* src, double* dest)

cdef class PySimulation_Settingsd:
    cdef Simulation_Settingsd c_settings

    def __cinit__(self):
        self.c_settings = Simulation_Settingsd()
        
    #here come the getters and setters FML ...
    
    @property
    def n_events(self):
        return self.c_settings.n_events

    @n_events.setter
    def n_events(self, val):
        self.c_settings.n_events = val

    @property
    def w_min(self):
        return self.c_settings.w_min

    @w_min.setter
    def w_min(self, val):
        self.c_settings.w_min = val

    @property
    def w_max(self):
        return self.c_settings.w_max

    @w_max.setter
    def w_max(self, val):
        self.c_settings.w_max = val

    @property
    def R(self):
        return self.c_settings.R

    @R.setter
    def R(self, val):
        self.c_settings.R = val

    @property
    def run_duration(self):
        return self.c_settings.run_duration

    @run_duration.setter
    def run_duration(self, val):
        self.c_settings.run_duration = val

    @property
    def mean_event_lifetime(self):
        return self.c_settings.mean_event_lifetime

    @mean_event_lifetime.setter
    def mean_event_lifetime(self, val):
        self.c_settings.mean_event_lifetime = val

    @property
    def trap_efficiency(self):
        return self.c_settings.trap_efficiency

    @trap_efficiency.setter
    def trap_efficiency(self, val):
        self.c_settings.trap_efficiency = val

    @property
    def seed(self):
        return self.c_settings.seed

    @seed.setter
    def seed(self, val):
        self.c_settings.seed = val

    @property
    def manual(self):
        return self.c_settings.manual

    @manual.setter
    def manual(self, val):
        self.c_settings.manual = val

    @property
    def e_r(self):
        return self.c_settings.e_r

    @e_r.setter
    def e_r(self, val):
        self.c_settings.e_r = val

    @property
    def e_phi(self):
        return self.c_settings.e_phi

    @e_phi.setter
    def e_phi(self, val):
        self.c_settings.e_phi = val

    @property
    def w0(self):
        return self.c_settings.w0

    @w0.setter
    def w0(self, val):
        self.c_settings.w0 = val

    @property
    def N(self):
        return self.c_settings.N

    @N.setter
    def N(self, val):
        self.c_settings.N = val

    @property
    def snr(self):
        return self.c_settings.snr

    @snr.setter
    def snr(self, val):
        self.c_settings.snr = val

    @property
    def sample_rate(self):
        return self.c_settings.sample_rate

    @sample_rate.setter
    def sample_rate(self, val):
        self.c_settings.sample_rate = val

    @property
    def w_mix(self):
        return self.c_settings.w_mix

    @w_mix.setter
    def w_mix(self, val):
        self.c_settings.w_mix = val

    @property
    def n_samples(self):
        return self.c_settings.n_samples

    @n_samples.setter
    def n_samples(self, val):
        self.c_settings.n_samples = val     
		
cdef class PySimulation_Settings:
    cdef Simulation_Settingsf c_settings

    def __cinit__(self):
        self.c_settings = Simulation_Settingsf()
        
    #here come the getters and setters FML ...
    
    @property
    def n_events(self):
        return self.c_settings.n_events

    @n_events.setter
    def n_events(self, val):
        self.c_settings.n_events = val

    @property
    def w_min(self):
        return self.c_settings.w_min

    @w_min.setter
    def w_min(self, val):
        self.c_settings.w_min = val

    @property
    def w_max(self):
        return self.c_settings.w_max

    @w_max.setter
    def w_max(self, val):
        self.c_settings.w_max = val

    @property
    def R(self):
        return self.c_settings.R

    @R.setter
    def R(self, val):
        self.c_settings.R = val

    @property
    def run_duration(self):
        return self.c_settings.run_duration

    @run_duration.setter
    def run_duration(self, val):
        self.c_settings.run_duration = val

    @property
    def mean_event_lifetime(self):
        return self.c_settings.mean_event_lifetime

    @mean_event_lifetime.setter
    def mean_event_lifetime(self, val):
        self.c_settings.mean_event_lifetime = val

    @property
    def trap_efficiency(self):
        return self.c_settings.trap_efficiency

    @trap_efficiency.setter
    def trap_efficiency(self, val):
        self.c_settings.trap_efficiency = val

    @property
    def seed(self):
        return self.c_settings.seed

    @seed.setter
    def seed(self, val):
        self.c_settings.seed = val

    @property
    def manual(self):
        return self.c_settings.manual

    @manual.setter
    def manual(self, val):
        self.c_settings.manual = val

    @property
    def e_r(self):
        return self.c_settings.e_r

    @e_r.setter
    def e_r(self, val):
        self.c_settings.e_r = val

    @property
    def e_phi(self):
        return self.c_settings.e_phi

    @e_phi.setter
    def e_phi(self, val):
        self.c_settings.e_phi = val

    @property
    def w0(self):
        return self.c_settings.w0

    @w0.setter
    def w0(self, val):
        self.c_settings.w0 = val

    @property
    def N(self):
        return self.c_settings.N

    @N.setter
    def N(self, val):
        self.c_settings.N = val

    @property
    def snr(self):
        return self.c_settings.snr

    @snr.setter
    def snr(self, val):
        self.c_settings.snr = val

    @property
    def sample_rate(self):
        return self.c_settings.sample_rate

    @sample_rate.setter
    def sample_rate(self, val):
        self.c_settings.sample_rate = val

    @property
    def w_mix(self):
        return self.c_settings.w_mix

    @w_mix.setter
    def w_mix(self, val):
        self.c_settings.w_mix = val

    @property
    def n_samples(self):
        return self.c_settings.n_samples

    @n_samples.setter
    def n_samples(self, val):
        self.c_settings.n_samples = val
        
cdef class PyBeamformer:
    cdef shared_ptr[Beamformerf] bf
    cdef public PySimulation_Settings settings
    cdef public int n_packets
    cdef public int grid_size

    def __cinit__(self, PySimulation_Settings settings, int grid_size,
                    int n_packets, bool weighted, bool full_frequency,
                    float r_grid):
        self.bf = make_shared[Beamformerf](settings.c_settings, grid_size,
                                            n_packets, weighted, 
                                            full_frequency, r_grid)
                                            
        self.settings = settings
        self.n_packets = n_packets
        self.grid_size = grid_size
                                            
                                            
    def get_next(self):

        n_samples = self.settings.n_samples
        bins = n_samples//2 + 1
        cdef numpy.ndarray[ndim=1, dtype=numpy.float32_t] res = numpy.zeros(
                        self.grid_size*self.grid_size*bins*self.n_packets, 
                        dtype=numpy.float32)

        self.bf.get().get_next(&res[0])

        res_ = res.reshape(self.n_packets, self.grid_size, self.grid_size, bins)
        res_ = numpy.moveaxis(res_,[1,2,3],[2,3,1])
        
        ind = (numpy.isfinite(res_))^True
        res_[ind] = 0
        res_[res_==-1] = 0

        return res_
        
    def get_result(self, data):

        n_samples = data.shape[-1]
        print(data.shape)
        bins = n_samples
        cdef numpy.ndarray[ndim=1, dtype=numpy.float32_t] res = numpy.zeros(
                        self.grid_size*self.grid_size*bins*self.n_packets, 
                        dtype=numpy.float32)

        cdef float complex[:] data_in = data.astype(numpy.complex64).flatten()
        self.bf.get().get_result(&data_in[0], &res[0])

        res_ = res.reshape(self.n_packets, self.grid_size, self.grid_size, bins)
        res_ = numpy.moveaxis(res_,[1,2,3],[2,3,1])
        
        ind = (numpy.isfinite(res_))^True
        
        res_[ind] = 0
        res_[res_==-1] = 0

        return res_
        
cdef class PyBeamformerd:
    cdef shared_ptr[Beamformerd] bf
    cdef public PySimulation_Settingsd settings
    cdef public int n_packets
    cdef public int grid_size

    def __cinit__(self, PySimulation_Settingsd settings, int grid_size,
                    int n_packets, bool weighted, bool full_frequency,
                    double r_grid):
        self.bf = make_shared[Beamformerd](settings.c_settings, grid_size,
                                            n_packets, weighted, 
                                            full_frequency, r_grid)
                                            
        self.settings = settings
        self.n_packets = n_packets
        self.grid_size = grid_size
                                            
                                            
    def get_next(self):

        n_samples = self.settings.n_samples
        bins = n_samples//2 + 1
        cdef numpy.ndarray[ndim=1, dtype=numpy.float64_t] res = numpy.zeros(
                        self.grid_size*self.grid_size*bins*self.n_packets, 
                        dtype=numpy.float64)

        self.bf.get().get_next(&res[0])

        res_ = res.reshape(self.n_packets, self.grid_size, self.grid_size, bins)
        res_ = numpy.moveaxis(res_,[1,2,3],[2,3,1])
        
        ind = (numpy.isfinite(res_))^True
        res_[ind] = 0
        res_[res_==-1] = 0

        return res_
        
    def get_result(self, data):

        n_samples = data.shape[-1]
        print(data.shape)
        bins = n_samples
        cdef numpy.ndarray[ndim=1, dtype=numpy.float64_t] res = numpy.zeros(
                        self.grid_size*self.grid_size*bins*self.n_packets, 
                        dtype=numpy.float64)

        cdef double complex[:] data_in = data.flatten()
        self.bf.get().get_result(&data_in[0], &res[0])

        res_ = res.reshape(self.n_packets, self.grid_size, self.grid_size, bins)
        res_ = numpy.moveaxis(res_,[1,2,3],[2,3,1])
        
        ind = (numpy.isfinite(res_))^True
        
        res_[ind] = 0
        res_[res_==-1] = 0

        return res_
        
