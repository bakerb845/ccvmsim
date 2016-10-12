# Define your compilers and your compiler flags - I recommend GNU
CC = clang
F90 = gfortran
INCC = -I./include

CPRE_FLAGS = -DVISIT_WRITER
CFLAGS = -O2 -Weverything -fopenmp -std=gnu11 -g3 $(CPRE_FLAGS)

# Find your libraries
H5_INC = -I/home/bakerb3/C/hdf5-1.10.0-patch1/include
H5_LIB = /home/bakerb3/C/hdf5-1.10.0-patch1/lib/libhdf5.a
Z_LIB = /home/bakerb3/C/zlib-1.2.8/libz.a
INI_INC = -I/home/bakerb3/C/iniparser/src
INI_LIB = /home/bakerb3/C/iniparser/libiniparser.a
ISCL_INC = -I/home/bakerb3/C/libiscl/include
ISCL_LIB = /home/bakerb3/C/libiscl/libiscl_shared.so
IPP_LIB = /opt/intel/ipp/lib/intel64/libipps.a \
	  /opt/intel/ipp/lib/intel64/libippvm.a \
	  /opt/intel/ipp/lib/intel64/libippcore.a
LAPACK_LIB = /usr/lib/liblapacke.so /usr/lib/liblapack.so /usr/lib/libblas.so
MKLROOT = /opt/intel/mkl
LMKL = -Wl,--start-group \
	${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
	${MKLROOT}/lib/intel64/libmkl_core.a \
	${MKLROOT}/lib/intel64/libmkl_sequential.a \
	-Wl,--end-group -lpthread -lm -ldl
LAPACK_LIB = $(LMKL)
#GSL_INC = -I/home/bakerb3/C/gsl-2.1/include
#GSL_LIB = /home/bakerb3/C/gsl-2.1/lib/libgsl.a
GEO_LIB = /home/bakerb3/C/GeographicLib-1.46/src/libGeographic.so
FFTW_LIB = /home/bakerb3/C/fftw-3.3.5/lib/libfftw3.a
VISIT_LIB = /home/bakerb3/C/VisIt/libvisit.a
VISIT_INC = -I/home/bakerb3/C/VisIt 

OBJ = ./Obj

INCALL = $(INCC) $(H5_INC) $(INI_INC) $(ISCL_INC) $(GSL_INC) $(VISIT_INC)
LIBALL = $(ISCL_LIB) $(FFTW_LIB) $(GEO_LIB) $(LAPACK_LIB) \
	 $(GSL_LIB) $(INI_LIB) $(IPP_LIB) $(H5_LIB) $(Z_LIB) \
	 $(VISIT_LIB) -lm -ldl -lgfortran

