#CC       = h5c++
CC       = mpicxx.mpich

#
GIOLIBDIR=/home/dkorytov/tools/genericio/frontend
GIOLIBMPIDIR=/home/dkorytov/tools/genericio/mpi
GIOINCDIR=/home/dkorytov/tools/genericio
hdf5_inc=/home/dkorytov/.local/include
#-I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmp
# GIOLIBDIR= /media/luna1/dkorytov/projects/hacc/trunk/datastar.cpu/frontend/lib
# GIOLIBMPIDIR= /media/luna1/dkorytov/projects/hacc/trunk/datastar.cpu/mpi/lib
# GIOINCDIR= /media/luna1/dkorytov/projects/hacc/trunk/genericio
hdf5_opts= /home/dkorytov/.local/lib/libhdf5_hl_cpp.a /home/dkorytov/.local/lib/libhdf5_cpp.a /home/dkorytov/.local/lib/libhdf5_hl.a /home/dkorytov/.local/lib/libhdf5.a -Bsymbolic-functions -z relro -lpthread -lz -ldl -lstdc++ -lm -lgcc_s -lgcc -lc -lgcc_s -lgcc
#/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../x86_64-linux-gnu/crtn.o
# ^ removed as a test
CFLAGS   =  -O3 -Isrc -I${GIOINCDIR} -Idtk -I${hdf5_inc} -fPIC -fopenmp -g
#-L/usr/lib/openmpi/lib -lmpi_cxx -lmpi

#LDFLAGS  = -Lgenericio/genericio.build/libs -Ldtk/lib -lGenericIO -ldtk -lgsl -lgslcblas ${hdf5_opts} -ldtk -pthread
LDFLAGS  = -L${GIOLIBMPIDIR} -Ldtk/lib -lGenericIO -ldtk -lgsl -lgslcblas ${hdf5_opts} -ldtk -pthread
SOURCES  = $(wildcard src/*.cpp)
HEADERS  = $(wildcard src/*.hpp)
OBJECTS  = $(SOURCES:src/%.cpp=obj/%.o)
EXE      = core_fit
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

