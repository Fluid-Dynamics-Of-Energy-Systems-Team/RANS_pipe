
#F77 = mpif90
#F77 = /Users/Eva/Programs/openmpi/bin/mpif90 -O0 -fdefault-real-8 -ffixed-line-length-none
#F77 = mpif90 -O0 -Wall -Wno-tabs -Wno-unused-dummy-argument -Wno-unused-variable -fdefault-real-8 -ffixed-line-length-none
#F77 = mpif90 -O2 -fdefault-real-8 -ffixed-line-length-none
#F77 = mpif90 -O0 -g -fbounds-check -finit-local-zero -fdefault-real-8 -ffixed-line-length-none
#F77 = mpif90 -O0 -g -fbounds-check -Wall -Wno-unused-variable -fdefault-real-8 -ffixed-line-length-none
F77 = mpif90 -O3 -fdefault-real-8 -ffixed-line-length-none
#F77 = mpif90 -O2 -r8 -132
#F77 = /usr/mpi/intel/openmpi-1.4.2/bin/mpif90  -r8 -132 -O5 #-fpe0 -traceback
#F77 = /opt/mvapich2/bin/mpif90 -r8 -132
#F77 = /opt/openmpi/bin/mpif90 -r8 -132
FC = $(F77) $(FLAGS)
DBG = -g

LIBS = -lmpi
RM = rm -f


PROGRAM = forced_real_VF
SRCS    = numerics.f90 mk.f90 vf.f90 sst.f90 sa.f90 turbmodels.f90 fileio.f math.f mpistuff.f vfft.f parpipe.f90
OBJS    = param.o common.o numerics.o mk.o vf.o sst.o sa.o turbmodels.o fileio.o math.o mpistuff.o vfft.o parpipe.o

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(F77) $(DBG) $(FLAGS) -o $(PROGRAM) $(OBJS) $(LIBS)

parpipe.o: parpipe.f90 param.f90 makefile
	$(F77) $(DBG) $(FLAGS) -c parpipe.f90
numerics.o: numerics.f90 param.f90 makefile
	$(F77) $(DBG) $(FLAGS) -c numerics.f90
mk.o: mk.f90 param.f90 makefile
	$(F77) $(DBG) $(FLAGS) -c mk.f90
vf.o: vf.f90 param.f90 makefile
	$(F77) $(DBG) $(FLAGS) -c vf.f90
sa.o: sa.f90 param.f90 makefile
	$(F77) $(DBG) $(FLAGS) -c sa.f90
sst.o: sst.f90 param.f90 makefile
	$(F77) $(DBG) $(FLAGS) -c sst.f90
turbmodels.o: param.f90 turbmodels.f90 makefile
	$(F77) $(DBG) $(FLAGS) -c turbmodels.f90
fileio.o: fileio.f90 param.f90 makefile
	$(F77) $(DBG) $(FLAGS) -c fileio.f90
math.o: math.f90 param.f90 makefile
	$(F77) $(DBG) $(FLAGS) -c math.f90
mpistuff.o: mpistuff.f90 param.f90 makefile
	$(F77) $(DBG) $(FLAGS) -c mpistuff.f90	
vfft.o: vfft.f makefile
	$(F77) $(DBG) $(FLAGS) -c vfft.f
param.o: param.f90
	$(F77) $(DBG) $(FLAGS) -c param.f90
common.o: common.f90
	$(F77) $(DBG) $(FLAGS) -c common.f90



changeCPU: changeNcpu.o
	$(F77) $(FLAGS) -o changeCPU changeNcpu.o $(LIBS)

changeCPU.o: changeCPU.f90 param.f90 makefile
	$(F77) $(FLAGS) -c changeCPU.f90

clean:
	$(RM) a.out core *.mod *.o $(PROGRAM) $(OBJS)



