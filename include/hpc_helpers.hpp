#ifndef HPC_HELPERS_HPP
#define HPC_HELPERS_HPP

#include <iostream>
#include <cstdint>

static_assert(sizeof(unsigned long long int) == sizeof(uint64_t), "No 64bit datatype used!");
typedef unsigned long long int uint64_cu;

#ifndef __CUDACC__
    #include <chrono>
#endif

#ifndef __CUDACC__
    #define TIMERSTART(label)                                                  \
        std::chrono::time_point<std::chrono::system_clock> a##label, b##label; \
        a##label = std::chrono::system_clock::now();
#else
    #define TIMERSTART(label)                                                  \
        cudaEvent_t start##label, stop##label;                                 \
        float time##label;                                                     \
        cudaEventCreate(&start##label);                                        \
        cudaEventCreate(&stop##label);                                         \
        cudaEventRecord(start##label, 0);
#endif

#ifndef __CUDACC__
    #define TIMERSTOP(label)                                                   \
        b##label = std::chrono::system_clock::now();                           \
        std::chrono::duration<double> delta##label = b##label-a##label;        \
        std::cerr << "TIMING: " << delta##label.count() << " s "			   \
				  << #label << std::endl;
        //~ std::cerr << "# elapsed time ("<< #label <<"): "                       \
                  //~ << delta##label.count()  << "s" << std::endl;

#else
    #define TIMERSTOP(label)                                                   \
            cudaEventRecord(stop##label, 0);                                   \
            cudaEventSynchronize(stop##label);                                 \
            cudaEventElapsedTime(&time##label, start##label, stop##label);     \
            std::cerr << "TIMING: " << time##label << " ms "  				   \
                      << #label << std::endl;
	
	#define TIMERBW(memsize, label)											   \
			cudaEventRecord(stop##label, 0);                                   \
            cudaEventSynchronize(stop##label);                                 \
            cudaEventElapsedTime(&time##label, start##label, stop##label);     \
            std::cerr << "TIMING: " << time##label << " ms "  				   \
                      << #label << std::endl;								   \
			double bandwidth##label = (double) memsize / (double)1e9;		   \
			bandwidth##label = bandwidth##label / (time##label*1e-3);          \
			std::cerr << "BANDWIDTH: " << bandwidth##label << " GB/s "		   \
                      << #label << std::endl;	
#endif


#ifdef __CUDACC__
    #define CUERR {                                                            \
        cudaError_t err;                                                       \
        if ((err = cudaGetLastError()) != cudaSuccess) {                       \
            std::cerr << "CUDA error: " << cudaGetErrorString(err) << " : "    \
                      << __FILE__ << ", line " << __LINE__ << std::endl;       \
            exit(1);                                                           \
        }                                                                      \
    }
    
    #define FREE_CUDA(pointer) \
    do { \
		if(pointer) { \
			cudaFree(pointer); \
			CUERR \
			pointer=nullptr; \
		} \
    } while (false)
    
    #define FREE_HOST(pointer) \
    do { \
		if(pointer) { \
			cudaFreeHost(pointer); \
			CUERR \
			pointer=nullptr; \
		} \
    } while (false) 

    // transfer constants
    #define H2D (cudaMemcpyHostToDevice)
    #define D2H (cudaMemcpyDeviceToHost)
    #define H2H (cudaMemcpyHostToHost)
    #define D2D (cudaMemcpyDeviceToDevice)
#endif

// safe division
#define SDIV(x,y)(((x)+(y)-1)/(y))

#endif
