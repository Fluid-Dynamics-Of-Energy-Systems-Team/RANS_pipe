module mod_param
  implicit none
  public

  integer   i,k,imax,kmax,i1,k1,iwork,isave,px,nTab, kelem
  integer   EOSmode,periodic,kmaxper,k1Old
  integer   nstep,K_start_heat,select_init,systemSolve
  integer   turbmod,modifDiffTerm,modVF,profiling,isothermalBC,pressIsoThermal
  real*8    Re,Pr,Qwall,CFL,PrT,dTwall, LOD,x_start_heat
  real*8    ctheta,Fr_1, dtmax,Tw_top, Tw_bot
  real*8    alphac,alphak,alphae,alphav2, Qsource, convCrit,alphakt,alphaet
  integer   Mt,Nx,Mx,Nt , turbdiffmod, restart, imax_old, kelem_old, bulkmod,noutput
  character(len=40) output_fname, output_fname_bl, read_fname,fluid_name
  character(len=500) table_loc

  NAMELIST /input/ CFL, systemsolve, imax, K_start_heat, kelem, x_start_heat,      &
                 iwork, isave, nTab, nstep, LoD, periodic, Re, Pr,PrT, Qwall, &
                 isothermalBC, Tw_top, Tw_bot, dTwall, pressIsoThermal, Fr_1, EOSmode, &
                 select_init, turbmod, ctheta, modVF, alphac,alphak,alphae, &
                 alphav2, modifDiffTerm, turbdiffmod, output_fname, dtmax, output_fname_bl, read_fname, &
                 imax_old, kelem_old, Qsource, bulkmod,table_loc,fluid_name, convCrit,alphakt,alphaet, noutput

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
