# Compiler and flags
# PATH for required libraries
ARMADILLO_DIR := /opt/homebrew/opt/armadillo
ifneq ($(MKLROOT),)
    MKLROOT :=/opt/intel/oneapi/mkl/latest
endif
OPENBLAS_DIR := /opt/homebrew/Cellar/openblas/0.3.29
OPENMPI_DIR  := /opt/homebrew/opt/openmpi
LIBOMP_DIR   := /opt/homebrew/opt/libomp
NETCDF_DIR   := /opt/homebrew/opt/netcdf
NETCDFCXX_DIR:= /opt/homebrew/opt/netcdf-cxx
BIGWHAM_DIR  := ../HydroMech3D_arma/BigWham
BIGWHAM_BUILD_DIR := $(BIGWHAM_DIR)/build_openblas
BIGWHAM_INCLUDE := -I$(BIGWHAM_DIR)/src -I$(BIGWHAM_DIR)/il

CXX = mpicxx
USE_MKL=0

ifeq ($(USE_MKL), 1)
CXXFLAGS = -g -std=c++14 -Wall -Wextra -pedantic -O2 -m64 -DIL_MKL -DIL_BLAS

LDFLAGS  = -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -L$(LIBOMP_DIR)/lib -lomp -ldl -lpthread -lm -L$(ARMADILLO_DIR)/lib -larmadillo -L$(NETCDF_DIR)/lib -lnetcdf -L$(NETCDFCXX_DIR)/lib -lnetcdf-cxx4  

INCLUDES := -I./JSON -I$(MKLROOT)/include -I./ -I$(ARMADILLO_DIR)/include -I$(LIBOMP_DIR)/include -I$(NETCDF_DIR)/include -I$(NETCDFCXX_DIR)/include
else 
#Compiler flags
CXXFLAGS = -g -std=c++17 -Wall -Wextra -pedantic -O2 -DIL_OPENBLAS -Xpreprocessor -fopenmp -w 

#Linking flags
LDFLAGS  = -L$(LIBOMP_DIR)/lib -lomp -L$(OPENBLAS_DIR)/lib -lopenblas -L$(ARMADILLO_DIR)/lib -larmadillo -L$(NETCDF_DIR)/lib -lnetcdf -L$(NETCDFCXX_DIR)/lib -lnetcdf-cxx4  

#Include directories
INCLUDES := -I./JSON -I$(OPENBLAS_DIR)/include -I./ -I$(ARMADILLO_DIR)/include -I$(NETCDF_DIR)/include -I$(NETCDFCXX_DIR)/include -I$(LIBOMP_DIR)/include 
endif

# Define the path to hpro-config
#HPRO_CONFIG = ../HydroMech3D_arma/hlibpro-3.1/bin/hpro-config

# Capture compiler flags and linker flags
#INCLUDES += $(shell $(HPRO_CONFIG) --cflags)
#LDFLAGS += $(shell $(HPRO_CONFIG) --lflags)

# Add BigWham includes
INCLUDES += $(BIGWHAM_INCLUDE)
# Add BigWham library path and link flag
LDFLAGS += -L$(BIGWHAM_BUILD_DIR) -lBigWham -Wl,-rpath,/Users/zhenhuan/Documents/code/HydroMech3D_arma/BigWham/build_openblas

INCLUDES += -I/opt/homebrew/opt/lapack/include
LDFLAGS += -L/opt/homebrew/opt/lapack/lib -llapack

#INCLUDES += -I/opt/homebrew/opt/open-mpi/include
#LDFLAGS  += -L/opt/homebrew/opt/open-mpi/lib 
# Source files
SRC_DIR = src
SRC_FILES = $(SRC_DIR)/loadProgramArguments.cpp \
            $(SRC_DIR)/EQsolver.cpp \
            $(SRC_DIR)/ImportJsonInputData.cpp \
            $(SRC_DIR)/Mesh.cpp \
            $(SRC_DIR)/FullSpaceElasticity.cpp \
            $(SRC_DIR)/StressKernelsP0/StressKernelsDxP0.cpp \
            $(SRC_DIR)/StressKernelsP0/StressKernelsDyP0.cpp \
            $(SRC_DIR)/StressKernelsP0/StressKernelsDzP0.cpp \
            $(SRC_DIR)/Solution.cpp \
            $(SRC_DIR)/RK45.cpp \
            $(SRC_DIR)/Utils.cpp \
            $(SRC_DIR)/AssembleElasticityMatrix.cpp \
	    ${SRC_DIR}/Hmatrix_bigwham.cpp \
	    ${SRC_DIR}/ExportBackgroundStress.cpp \
            ${SRC_DIR}/ExportBackgroundStressDirect.cpp \
            main.cpp

# Object files
OBJ_FILES = $(SRC_FILES:.cpp=.o)

# Target executable
TARGET = 3DEQSim

# Rules
all: $(TARGET)

$(TARGET): $(OBJ_FILES)
	$(CXX) $(OBJ_FILES) $(LDFLAGS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJ_FILES) $(TARGET)

.PHONY: all clean
