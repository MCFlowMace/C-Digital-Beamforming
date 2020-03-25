
#c++ compiler options
EXEC=SERIAL
GCC=g++
GCCFLAGS= -O3 -std=c++17 -ffast-math -fopenmp -D$(EXEC)

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

OBJ_TEST := $(filter-out $(OBJDIR)/main.o, $(OBJECTS))

REBUILDABLES = $(OBJECTS) $(BINDIR)/$(TARGET)

# Link c++ and CUDA compiled object files to target executable:
#@$(GCC) $(GCCFLAGS) $(OBJECTS) $(LIB_DIR) $(LIBS) -Xcompiler \"$(RUNTIME_SEARCH_PATH)\" -o $@

$(TARGET): $(OBJECTS)
	@$(GCC) $(GCCFLAGS) $(OBJECTS) $(LIB_DIR) $(LIBS) -o $@
	@echo "Linking complete!"

# Compile C++ source files to object files:
$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	@$(GCC) $(INCLUDES) $(GCCFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

# Tests

reconstruction: $(OBJECTS)
	@$(GCC) $(INCLUDES) $(GCCFLAGS) $(TESTDIR)/reconstruction/reconstruction.cpp $(OBJ_TEST) $(LIB_DIR) $(LIBS) -o $(BINDIR)/reconstruction
	@echo "Linking complete!"

measure-snr: $(OBJECTS)
	@$(GCC) $(INCLUDES) $(GCCFLAGS) $(TESTDIR)/measure-snr/measure-snr.cpp $(OBJ_TEST) $(LIB_DIR) $(LIBS) -o $(BINDIR)/measure-snr
	@echo "Linking complete!"

clean :
	@echo " Cleaning...";
	@echo " $(RM) -r $(OBJDIR) $(TARGET)"; $(RM) -rf $(OBJDIR) $(BINDIR)/*
#rm -f $(REBUILDABLES)
#rm -rf $(BINDIR)
#echo Clean done

.PHONY: clean
