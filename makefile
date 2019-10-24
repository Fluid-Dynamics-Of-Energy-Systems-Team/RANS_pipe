F77 = mpif90 -O5 -fdefault-real-8 -ffixed-line-length-none
FC = $(F77) $(FLAGS)

LIBS = -lmpi
RM = rm -f


PROGRAM = run 
OBJS    = param.o math.o tm.o eosmodels.o mesh.o  tm_ke.o tm_mk.o tm_vf.o tm_sa.o tm_sst.o   common.o numerics.o  fileio.o mpistuff.o vfft.o main.o

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(F77) $(DBG) $(FLAGS) -o $(PROGRAM) $(OBJS) $(LIBS)

main.o: main.f90 param.f90 makefile
	$(F77) $(DBG) $(FLAGS) -c main.f90
tm_vf.o: tm_vf.f90 makefile
	$(F77) $(DBG) $(FLAGS) -c tm_vf.f90
tm_ke.o: tm_ke.f90 makefile
	$(F77) $(DBG) $(FLAGS) -c tm_ke.f90
tm_mk.o: tm_mk.f90 makefile
	$(F77) $(DBG) $(FLAGS) -c tm_mk.f90
tm_sst.o: tm_sst.f90 makefile
	$(F77) $(DBG) $(FLAGS) -c tm_sst.f90
tm_sa.o: tm_sa.f90 makefile
	$(F77) $(DBG) $(FLAGS) -c tm_sa.f90
tm.o: tm.f90 makefile
	$(F77) $(DBG) $(FLAGS) -c tm.f90
mesh.o: mesh.f90 makefile
	$(F77) $(DBG) $(FLAGS) -c mesh.f90
eosmodels.o: eosmodels.f90 makefile
	$(F77) $(DBG) $(FLAGS) -c eosmodels.f90
numerics.o: numerics.f90 param.f90 makefile
	$(F77) $(DBG) $(FLAGS) -c numerics.f90
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



