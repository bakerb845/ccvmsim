# Define your compilers and your compiler flags - I recommend GNU
CC = gcc-5
#MPICC = mpicc
MPICC = gcc-5
F90 = gfortran-5
INCC = -I./include

CFLAGS = -O2 -Wall -fopenmp -std=gnu11 -g3

# Find your libraries
H5_INC = -I/home/bakerb3/Fortran/hdf5-1.8.16/include
H5_LIB = /home/bakerb3/Fortran/hdf5-1.8.16/lib/libhdf5.a
Z_LIB = /home/bakerb3/C/zlib-1.2.8/libz.a
INI_INC = -I/home/bakerb3/C/iniparser/src
INI_LIB = /home/bakerb3/C/iniparser/libiniparser.a
ISCL_INC = -I/home/bakerb3/C/libiscl/include
ISCL_LIB = /home/bakerb3/C/libiscl/libiscl.so
IPP_LIB = /opt/intel/ipp/lib/intel64/libipps.a \
	  /opt/intel/ipp/lib/intel64/libippcore.a
LAPACK_LIB = /usr/lib/liblapacke.so /usr/lib/liblapack.so /usr/lib/libblas.so
GSL_INC = -I/home/bakerb3/C/gsl-2.1/include
GSL_LIB = /home/bakerb3/C/gsl-2.1/lib/libgsl.a
GEO_LIB = /home/bakerb3/Fortran/GeographicLib-1.46/src/libGeographic.a
FFTW_LIB = /home/bakerb3/C/fftw-3.3.4/lib/libfftw3.a

OBJ = ./Obj

INCALL = $(INCC) $(H5_INC) $(INI_INC) $(ISCL_INC) $(GSL_INC)
LIBALL = $(ISCL_LIB) $(FFTW_LIB) $(GEO_LIB) $(LAPACK_LIB) \
	 $(GSL_LIB) $(INI_LIB) $(IPP_LIB) $(H5_LIB) $(Z_LIB) \
	 -lm -ldl -lgfortran -lstdc++

