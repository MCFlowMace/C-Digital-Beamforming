#parameters
EXEC=USE_GPU
ARCH=sm_70

#c++ compiler options
GCC=g++
GCCFLAGS= -O3 -std=c++17 -ffast-math -fopenmp -D$(EXEC)

#NVCC options
NVCC=nvcc
NVCC_FLAGS= -arch=$(ARCH) -std=c++11 -I -dc --expt-extended-lambda -D$(EXEC) -DARMA_ALLOW_FAKE_GCC -Xcompiler -fopenmp --use_fast_math

#-D_GLIBCXX_USE_CXX11_ABI=0

#all: release

SRCDIR = src
INCDIR = include
OBJDIR = build
BINDIR = bin
TESTDIR = test

TARGET = $(BINDIR)/Beamforming

INCLUDES = -I $(INCDIR)
LIB_DIR = -L lib

#armadillo for FFT
LIBS = -larmadillo

SOURCES  := $(wildcard $(SRCDIR)/*.cpp)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

SOURCES_CU  := $(wildcard $(SRCDIR)/*.cu)
OBJECTS_CU  := $(SOURCES_CU:$(SRCDIR)/%.cu=$(OBJDIR)/%.cu.o)

OBJ_TEST := $(filter-out $(OBJDIR)/main.o, $(OBJECTS))

REBUILDABLES = $(OBJECTS) $(OBJECTS_CU) $(BINDIR)/$(TARGET)

# Link c++ and CUDA compiled object files to target executable:
#@$(GCC) $(GCCFLAGS) $(OBJECTS) $(LIB_DIR) $(LIBS) -Xcompiler \"$(RUNTIME_SEARCH_PATH)\" -o $@

$(TARGET): $(OBJECTS) $(OBJECTS_CU)
	@$(NVCC) $(NVCC_FLAGS) $(OBJECTS) $(OBJECTS_CU) $(LIB_DIR) $(LIBS) -o $@
	@echo "Linking complete!"

# Compile C++ source files to object files:
$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	@$(GCC) $(INCLUDES) $(GCCFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

# Compile CUDA source files to object files:
$(OBJECTS_CU): $(OBJDIR)/%.cu.o : $(SRCDIR)/%.cu
	@mkdir -p $(@D)
	@$(NVCC) $(INCLUDES) $(NVCC_FLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

# Link c++ and CUDA compiled object files to target executable:
$(BINDIR)/$(TARGET): $(OBJECTS) $(OBJECTS_CU)
	@$(NVCC) $(NVCC_FLAGS) $(OBJECTS) $(OBJECTS_CU) $(LIB_DIR) $(LIBS) -Xcompiler \"$(RUNTIME_SEARCH_PATH)\" -o $@
	@echo "Linking complete!"

# Tests

reconstruction: $(OBJECTS) $(OBJECTS_CU)
	@$(NVCC) $(INCLUDES) $(NVCC_FLAGS) $(TESTDIR)/reconstruction/reconstruction.cpp $(OBJ_TEST) $(OBJECTS_CU) $(LIB_DIR) $(LIBS) -o $(BINDIR)/reconstruction
	@echo "Linking complete!"

threshold-trigger: $(OBJECTS) $(OBJECTS_CU)
	@$(NVCC) $(INCLUDES) $(NVCC_FLAGS) $(TESTDIR)/threshold-trigger/threshold-trigger.cpp $(OBJ_TEST) $(OBJECTS_CU) $(LIB_DIR) $(LIBS) -o $(BINDIR)/threshold-trigger
	@echo "Linking complete!"

event_generator: $(OBJECTS) $(OBJECTS_CU)
	@$(NVCC) $(INCLUDES) $(NVCC_FLAGS) $(TESTDIR)/event_generation/generate_events.cpp $(OBJ_TEST) $(OBJECTS_CU) $(LIB_DIR) $(LIBS) -o $(BINDIR)/generate_events
	@echo "Linking complete!"

e_loss: $(OBJECTS) $(OBJECTS_CU)
	@$(NVCC) $(INCLUDES) $(NVCC_FLAGS) $(TESTDIR)/event_generation/e_loss.cpp $(OBJ_TEST) $(OBJECTS_CU) $(LIB_DIR) $(LIBS) -o $(BINDIR)/e_loss
	@echo "Linking complete!"

measure-snr: $(OBJECTS) $(OBJECTS_CU)
	@$(NVCC) $(INCLUDES) $(NVCC_FLAGS) $(TESTDIR)/measure-snr/measure_snr.cpp $(OBJ_TEST) $(OBJECTS_CU) $(LIB_DIR) $(LIBS) -o $(BINDIR)/measure-snr
	@echo "Linking complete!"

clean :
	@echo " Cleaning...";
	@echo " $(RM) -r $(OBJDIR) $(TARGET)"; $(RM) -rf $(OBJDIR) $(BINDIR)/*
#rm -f $(REBUILDABLES)
#rm -rf $(BINDIR)
#echo Clean done

.PHONY: clean
