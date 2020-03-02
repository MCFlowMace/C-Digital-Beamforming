
#c++ compiler options
EXEC=SERIAL
GCC=g++
GCCFLAGS= -O3 -std=c++17 -ffast-math -fopenmp -D$(EXEC)

#-D_GLIBCXX_USE_CXX11_ABI=0

all: release

TARGET = Beamforming

SRCDIR = src
INCDIR = include
OBJDIR = build/obj
BINDIR = build

INCLUDES = -I $(INCDIR)
#LIB_DIR = -L $(HDF5_LIB)
#LIBS = $(HDF5_LINK_LIBS)
LIB_DIR =

#armadillo for FFT
LIBS = -larmadillo
#RUNTIME_SEARCH_PATH = -Wl,-rpath,$(HDF5_LIB)

SOURCES  := $(wildcard $(SRCDIR)/*.cpp)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

REBUILDABLES = $(OBJECTS) $(BINDIR)/$(TARGET)

# Link c++ and CUDA compiled object files to target executable:
#@$(GCC) $(GCCFLAGS) $(OBJECTS) $(LIB_DIR) $(LIBS) -Xcompiler \"$(RUNTIME_SEARCH_PATH)\" -o $@

$(BINDIR)/$(TARGET): $(OBJECTS)
	@$(GCC) $(GCCFLAGS) $(OBJECTS) $(LIB_DIR) $(LIBS) -o $@
	@echo "Linking complete!"

# Compile C++ source files to object files:
$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	@$(GCC) $(INCLUDES) $(GCCFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"


release: $(BINDIR)/$(TARGET)

clean :
	rm -f $(REBUILDABLES)
	rm -rf $(BINDIR)
	echo Clean done
