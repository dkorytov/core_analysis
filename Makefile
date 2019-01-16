#CC       = h5c++
CC       = mpicxx.mpich2
#-I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmp
GIOLIBDIR= /media/luna1/dkorytov/projects/hacc/trunk/datastar.cpu/frontend/lib
GIOLIBMPIDIR= /media/luna1/dkorytov/projects/hacc/trunk/datastar.cpu/mpi/lib
GIOINCDIR= /media/luna1/dkorytov/projects/hacc/trunk/genericio
hdf5_opts= /usr/lib/x86_64-linux-gnu/libhdf5_hl_cpp.a /usr/lib/x86_64-linux-gnu/libhdf5_cpp.a /usr/lib/x86_64-linux-gnu/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/libhdf5.a -Bsymbolic-functions -z relro -lpthread -lz -ldl -lstdc++ -lm -lgcc_s -lgcc -lc -lgcc_s -lgcc /usr/lib/gcc/x86_64-linux-gnu/4.8/../../../x86_64-linux-gnu/crtn.o
CFLAGS   =  -O3 -Isrc -I${GIOINCDIR} -Idtk -fPIC -fopenmp -g
#-L/usr/lib/openmpi/lib -lmpi_cxx -lmpi

#LDFLAGS  = -Lgenericio/genericio.build/libs -Ldtk/lib -lGenericIO -ldtk -lgsl -lgslcblas ${hdf5_opts} -ldtk -pthread
LDFLAGS  = -L${GIOLIBMPIDIR} -Ldtk/lib -lGenericIOMPI -ldtk -lgsl -lgslcblas ${hdf5_opts} -ldtk -pthread
SOURCES  = $(wildcard src/*.cpp)
HEADERS  = $(wildcard src/*.hpp)
OBJECTS  = $(SOURCES:src/%.cpp=obj/%.o)
EXE      = main
N2LIB    = lib/libn2merg.so

all: ${EXE} ${N2LIB}

${EXE}: ${OBJECTS} 
	${CC} ${CFLAGS} -o $@ $^ ${LDFLAGS}

${OBJECTS}: obj/%.o : src/%.cpp 
	${CC} ${CFLAGS} -c -o $@ $<

${N2LIB}: obj/n2_merger.o 
	${CC} -shared -o $@ $<

.PHONY:
clean:
	@rm -f obj/*.o
	@rm -f main

