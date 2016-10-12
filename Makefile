include Makefile.inc

ifeq "$(wildcard $(OBJ) )" ""
-include $(shell mkdir $(OBJ)) $(wildcard $(OBJ)/*)
endif

#CVM_MESHER = cvm_mesher
CVM2H5 = cvm2h5

#EXECS = $(CVM_MESHER)
EXECS = $(CVM2H5)

OBJ_FILEIO = $(OBJ)/cvmio.o $(OBJ)/h5_cinter.o \
	     $(OBJ)/meshio.o $(OBJ)/readini.o $(OBJ)/unpack.o
OBJ_COMMON = $(OBJ)/density.o $(OBJ)/memory.o $(OBJ)/qualityFactor.o \
             $(OBJ)/extendBase.o $(OBJ)/utm_geo.o
OBJ_MESH_ROUTINES = $(OBJ)/element.o $(OBJ)/regmesh.o
OBJ_CVM_UTILS = $(OBJ)/layeredMesh.o \
	        $(OBJ)/qualityFactor.o $(OBJ)/template.o \
	        $(OBJ)/topo30.o $(OBJ_COMMON)

OBJ_CVM2H5 = $(OBJ)/cvm2h5.o $(OBJ_FILEIO) $(OBJ_MESH_ROUTINES) $(OBJ_COMMON)

OBJ_MESHER = $(OBJ)/cvm_mesher.o $(OBJ_FILEIO) $(OBJ_CVM_UTILS)

all: $(EXECS)

$(CVM_MESHER): $(OBJ_MESHER)
	$(CC) $(CFLAGS) -o $(CVM_MESHER) $(OBJ_MESHER) $(LIBALL)

$(CVM2H5): $(OBJ_CVM2H5)
	$(CC) $(CFLAGS) -o $(CVM2H5) $(OBJ_CVM2H5) $(LIBALL)

$(OBJ)/cvm2h5.o: cvm2h5.c
	$(CC) $(CFLAGS) $(INCALL) -c cvm2h5.c -o $(OBJ)/cvm2h5.o

$(OBJ)/cvmio.o: cvmio.c
	$(CC) $(CFLAGS) $(INCALL) -c cvmio.c -o $(OBJ)/cvmio.o

$(OBJ)/cvm_mesher.o: cvm_mesher.c
	$(CC) $(CFLAGS) $(INCALL) -c cvm_mesher.c -o $(OBJ)/cvm_mesher.o

$(OBJ)/density.o: density.c
	$(CC) $(CFLAGS) $(INCALL) -c density.c -o $(OBJ)/density.o

$(OBJ)/element.o: element.c
	$(CC) $(CFLAGS) $(INCALL) -c element.c -o $(OBJ)/element.o

$(OBJ)/extendBase.o: extendBase.c
	$(CC) $(CFLAGS) $(INCALL) -c extendBase.c -o $(OBJ)/extendBase.o

$(OBJ)/h5_cinter.o: h5_cinter.c
	$(CC) $(CFLAGS) $(INCALL) -c h5_cinter.c -o $(OBJ)/h5_cinter.o

$(OBJ)/layeredMesh.o: layeredMesh.c
	$(CC) $(CFLAGS) $(INCALL) -c layeredMesh.c -o $(OBJ)/layeredMesh.o

$(OBJ)/memory.o: memory.c
	$(CC) $(CFLAGS) $(INCALL) -c memory.c -o $(OBJ)/memory.o

$(OBJ)/meshio.o: meshio.c
	$(CC) $(CFLAGS) $(INCALL) -c meshio.c -o $(OBJ)/meshio.o

$(OBJ)/qualityFactor.o: qualityFactor.c
	$(CC) $(CFLAGS) $(INCALL) -c qualityFactor.c -o $(OBJ)/qualityFactor.o

$(OBJ)/readini.o: readini.c
	$(CC) $(CFLAGS) $(INCALL) -c readini.c -o $(OBJ)/readini.o

$(OBJ)/regmesh.o: regmesh.c
	$(CC) $(CFLAGS) $(INCALL) -c regmesh.c -o $(OBJ)/regmesh.o

$(OBJ)/template.o: template.c
	$(CC) $(CFLAGS) $(INCALL) -c template.c -o $(OBJ)/template.o

$(OBJ)/topo30.o: topo30.c
	$(CC) $(CFLAGS) $(INCALL) -c topo30.c -o $(OBJ)/topo30.o

$(OBJ)/unpack.o: unpack.c
	$(CC) $(CFLAGS) $(INCALL) -c unpack.c -o $(OBJ)/unpack.o

$(OBJ)/utm_geo.o: utm_geo.f90
	$(F90) $(FFLAGS) -c utm_geo.f90 -o $(OBJ)/utm_geo.o

clean:
	@$(RM) $(OBJ)/*.o $(EXECS) *.mod

