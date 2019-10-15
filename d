mpif90 -O3 -fdefault-real-8 -ffixed-line-length-none -g  -c common.f90
mpif90 -O3 -fdefault-real-8 -ffixed-line-length-none -g  -o forced_real_VF  turbmodels_module.o ke_tm.o mk_tm.o vf_tm.o sa_tm.o sst_tm.o math_module.o eosmodels.o  common_module.o param.o common.o numerics.o mk.o vf.o sst.o sa.o turbmodels.o fileio.o math.o mpistuff.o vfft.o parpipe.o -lmpi
