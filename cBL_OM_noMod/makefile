
#F77 = mpif90
#F77 = /Users/Eva/Programs/openmpi/bin/mpif90 -O0 -fdefault-real-8 -ffixed-line-length-none
#F77 = mpif90 -O0 -Wall -Wno-tabs -Wno-unused-dummy-argument -Wno-unused-variable -fdefault-real-8 -ffixed-line-length-none
F77 = mpif90 -O2 -fdefault-real-8 -ffixed-line-length-none
#F77 = mpif90 -O0 -g -fbounds-check -finit-local-zero -fdefault-real-8 -ffixed-line-length-none
#F77 = mpif90 -O0 -g -fbounds-check -Wall -Wno-unused-variable -fdefault-real-8 -ffixed-line-length-none
#F77 = mpif90 -O2 -r8 -132
#F77 = /usr/mpi/intel/openmpi-1.4.2/bin/mpif90  -r8 -132 -O5 #-fpe0 -traceback
#F77 = /opt/mvapich2/bin/mpif90 -r8 -132
#F77 = /opt/openmpi/bin/mpif90 -r8 -132
FC = $(F77) $(FLAGS)
#DBG = -g

LIBS = #-lmpi
RM = rm -f


PROGRAM = forced_real_VF
SRCS    = parpipe.f numerics.f fileio.f math.f mpistuff.f vfft.f
OBJS    = parpipe.o numerics.o fileio.o math.o mpistuff.o vfft.o

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(F77) $(DBG) $(FLAGS) -o $(PROGRAM) $(OBJS) $(LIBS)

parpipe.o: parpipe.f param.txt makefile
	$(F77) $(DBG) $(FLAGS) -c parpipe.f
numerics.o: numerics.f param.txt makefile
	$(F77) $(DBG) $(FLAGS) -c numerics.f
fileio.o: fileio.f param.txt makefile
	$(F77) $(DBG) $(FLAGS) -c fileio.f
math.o: math.f param.txt makefile
	$(F77) $(DBG) $(FLAGS) -c math.f
mpistuff.o: mpistuff.f param.txt makefile
	$(F77) $(DBG) $(FLAGS) -c mpistuff.f	
vfft.o: vfft.f makefile
	$(F77) $(DBG) $(FLAGS) -c vfft.f



changeCPU: changeNcpu.o
	$(F77) $(FLAGS) -o changeCPU changeNcpu.o $(LIBS)

changeCPU.o: changeCPU.f90 param.txt makefile
	$(F77) $(FLAGS) -c changeCPU.f90

clean:
	$(RM) a.out core *.mod *.o $(PROGRAM) $(OBJS)



