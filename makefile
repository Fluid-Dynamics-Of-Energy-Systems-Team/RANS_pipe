F77 =mpif90
FLAGS =-fdefault-real-8 -Wall -Wno-unused-variable -Wno-unused-dummy-argument
FC = $(F77) $(FLAGS)
#DBG = -g
LIBS = -lmpi
RM = rm -f

PROGRAM = run 
OBJS = param.o common.o
OBJS += math.o
OBJS += mesh.o
OBJS += eos.o
OBJS += tm.o tm_ke.o tm_mk.o tm_sa.o tm_vf.o tm_sst.o tm_abe.o
OBJS += tdm.o tdm_vp.o tdm_ktet.o tdm_dwx.o tdm_nk.o 
OBJS += numerics.o
OBJS += fileio.o
OBJS += mpistuff.o vfft.o
OBJS += main.o
OBJS += blktri.o solver.o

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(F77) $(DBG) $(FLAGS) -o $(PROGRAM) $(OBJS) $(LIBS)
%.o: %.f90
	$(F77) $(DBG) $(FLAGS) -c $<
clean:
	$(RM) *.mod *.o $(PROGRAM) $(OBJS)



