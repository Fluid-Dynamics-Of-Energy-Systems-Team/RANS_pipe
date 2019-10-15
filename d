mpif90 -O3 -fdefault-real-8 -ffixed-line-length-none -g  -c vf_tm.f90
vf_tm.f90:64:32:

         ((u(i,k)+u(im,k))/(2.*Rp(i)))**2.) +  &
                                1
Error: Function ‘rp’ at (1) has no IMPLICIT type
makefile:32: recipe for target 'vf_tm.o' failed
make: *** [vf_tm.o] Error 1
