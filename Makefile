#parameters
EXEC=USE_GPU
ARCH=sm_70

#c++ compiler options
GCC=g++
GCCFLAGS= -O3 -std=c++17 -fopenmp -D$(EXEC) -fPIC

PYTHON=python

#-ffast-math

#NVCC options
NVCC=nvcc
NVCC_FLAGS= -arch=$(ARCH) -std=c++11 -I -dc --expt-extended-lambda --use_fast_math -D$(EXEC) -DARMA_ALLOW_FAKE_GCC -Xcompiler -fopenmp -Xcompiler -fPIC

#libpacker
LIBPACK=ar
LIBPACK_FLAGS = rsv

#-D_GLIBCXX_USE_CXX11_ABI=0

#all: release

SRCDIR = src
INCDIR = include
BUILDDIR = build
OBJDIR = $(BUILDDIR)/obj
BINDIR = $(BUILDDIR)/bin
TESTDIR = test
LIBDIR = lib/beamforming

INCLUDES = -I $(INCDIR)
#LIB_DIR = -L

TARGET = $(BUILDDIR)/lib/libbeamforming.a

#armadillo for FFT
LIBS = -larmadillo

SOURCES  := $(wildcard $(LIBDIR)/*.cpp)
OBJECTS  := $(SOURCES:$(LIBDIR)/%.cpp=$(OBJDIR)/%.o)

SOURCES_CU  := $(wildcard $(LIBDIR)/*.cu)
OBJECTS_CU  := $(SOURCES_CU:$(LIBDIR)/%.cu=$(OBJDIR)/%.cu.o)

REBUILDABLES = $(OBJECTS) $(OBJECTS_CU) $(TARGET)


#build the library
$(TARGET): $(OBJECTS) $(OBJECTS_CU)
	@mkdir -p $(BINDIR)
	@mkdir -p $(@D)
	@$(LIBPACK) $(LIBPACK_FLAGS) $@ $(OBJECTS) $(OBJECTS_CU)
	@echo "Library build!"

# Compile C++ source files to object files:
$(OBJECTS): $(OBJDIR)/%.o : $(LIBDIR)/%.cpp
	@mkdir -p $(@D)
	@$(GCC) $(INCLUDES) $(GCCFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

# Compile CUDA source files to object files:
$(OBJECTS_CU): $(OBJDIR)/%.cu.o : $(LIBDIR)/%.cu
	@mkdir -p $(@D)
	@$(NVCC) $(INCLUDES) $(NVCC_FLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

# Link c++ and CUDA compiled object files to target executable:
#$(BINDIR)/$(TARGET): $(OBJECTS) $(OBJECTS_CU)
#	@$(NVCC) $(NVCC_FLAGS) $(OBJECTS) $(OBJECTS_CU) $(LIB_DIR) $(LIBS) -Xcompiler \"$(RUNTIME_SEARCH_PATH)\" -o $@
#	@echo "Linking complete!"

# Tests

pybeamforming: $(TARGET)
	@$(PYTHON) lib/pybeamforming/setup.py build_ext --inplace
	@mv pybeamforming.*.so build/lib/
	@rm -rf build/temp*
	@rm -f lib/pybeamforming/pybeamforming.cpp
	@echo "Python module complete!"
	
reconstruction: $(TARGET)
	@$(NVCC) $(INCLUDES) $(NVCC_FLAGS) $(TESTDIR)/reconstruction/test_rec.cpp $(TARGET) $(LIBS) -o $(BINDIR)/test_rec
	@echo "Linking complete!"
	
response_map: $(TARGET)
	@$(NVCC) $(INCLUDES) $(NVCC_FLAGS) $(TESTDIR)/response_map/response_map.cpp $(TARGET) $(LIBS) -o $(BINDIR)/response_map
	@echo "Linking complete!"

threshold-trigger: $(TARGET)
	@$(NVCC) $(INCLUDES) $(NVCC_FLAGS) $(TESTDIR)/threshold-trigger/threshold-trigger.cpp $(TARGET) $(LIBS) -o $(BINDIR)/threshold-trigger
	@echo "Linking complete!"
	
roc-evaluation: $(TARGET)
	@$(NVCC) $(INCLUDES) $(NVCC_FLAGS) $(TESTDIR)/threshold-trigger/roc_evaluation.cpp $(TARGET) $(LIBS) -o $(BINDIR)/roc-evaluation
	@echo "Linking complete!"

event_generator: $(TARGET)
	@$(NVCC) $(INCLUDES) $(NVCC_FLAGS) $(TESTDIR)/event_generation/generate_events.cpp $(TARGET) $(LIBS) -o $(BINDIR)/generate_events
	@echo "Linking complete!"

e_loss: $(TARGET)
	@$(NVCC) $(INCLUDES) $(NVCC_FLAGS) $(TESTDIR)/event_generation/e_loss.cpp $(TARGET) $(LIBS) -o $(BINDIR)/e_loss
	@echo "Linking complete!"

measure-snr: $(TARGET)
	@$(NVCC) $(INCLUDES) $(NVCC_FLAGS) $(TESTDIR)/measure-snr/measure_snr.cpp $(TARGET) $(LIBS) -o $(BINDIR)/measure-snr
	@echo "Linking complete!"

clean :
	@echo " Cleaning...";
	@$(RM) -rf $(BUILDDIR) ;

.PHONY: clean
