module mod_param
  implicit none
  public

  integer   i,k,imax,kmax,i1,k1,iwork,isave,px,nTab, kelem
  integer   EOSmode,periodic,kmaxper,k1Old
  integer   nstep,K_start_heat,select_init,systemSolve
  integer   turbmod,modifDiffTerm,modVF,profiling,isothermalBC,pressIsoThermal
  real*8    Re,Pr,Qwall,CFL,Tw,dTwall, LOD,x_start_heat
  real*8    ctheta,Fr_1, dtmax
  real*8    alphac,alphak,alphae,alphav2
  integer   Mt,Nx,Mx,Nt , turbdiffmod
  character(len=40) output_fname, filename2

  NAMELIST /input/ CFL, systemsolve, imax, K_start_heat, kelem, x_start_heat,      &
                 iwork, isave, nTab, nstep, LoD, periodic, Re, Pr, Qwall, &
                 isothermalBC, Tw, dTwall, pressIsoThermal, Fr_1, EOSmode, &
                 select_init, turbmod, ctheta, modVF, alphac,alphak,alphae, &
                 alphav2, modifDiffTerm, turbdiffmod, output_fname, dtmax, filename2

contains

subroutine read_parameters()
  character(len=30) inputfile
  if (command_argument_count().lt.1) then
    write(*,*) "No input file given"
  stop
  else
    call get_command_argument(1, inputfile)
  endif
  open(unit=10, file=inputfile)
  read (10, NML=input)
  close (10)
end subroutine read_parameters


end module mod_param
