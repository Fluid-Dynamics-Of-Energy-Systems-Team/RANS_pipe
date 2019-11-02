F77 =mpif90
FLAGS =-O5 -fdefault-real-8
FC = $(F77) $(FLAGS)
#DBG = -g
LIBS = -lmpi
RM = rm -f

PROGRAM = run 
OBJS = param.o common.o
OBJS += math.o
OBJS += eos.o
OBJS += mesh.o
OBJS += tm.o tm_ke.o tm_mk.o tm_sa.o tm_vf.o tm_sst.o
OBJS += tdm.o tdm_vp.o 
OBJS += numerics.o
OBJS += fileio.o
OBJS += mpistuff.o vfft.o
OBJS += main.o

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(F77) $(DBG) $(FLAGS) -o $(PROGRAM) $(OBJS) $(LIBS)
%.o: %.f90
	$(F77) $(DBG) $(FLAGS) -c $<
clean:
	$(RM) *.mod *.o $(PROGRAM) $(OBJS)



