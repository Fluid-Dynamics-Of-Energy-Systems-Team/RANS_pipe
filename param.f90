module mod_param
  implicit none
  public

  integer   i,k,imax,kmax,i1,k1,iwork,isave,px,nTab,LoD
  integer   EOSmode,periodic,kmaxper,k1Old
  integer   nstep,K_start_heat,x_start_heat,select_init,systemSolve
  integer   turbmod,modifDiffTerm,modVF,profiling,isothermalBC,pressIsoThermal
  real*8    Re,Pr,Qwall,CFL,Tw,dTwall
  real*8    ctheta,Fr_1
  real*8    alphac,alphak,alphae,alphav2
  integer   Mt,Nx,Mx,Nt 

  parameter (CFL               = 100.0)   ! CFL number / time step
  parameter (systemSolve       = 1)       ! 1: pipe, 2: channel 3: boundary layer !**dpdz now defined in mkgrid should be 1 for channel/BL, 4 for the pipe**!!
  parameter (imax              = 96)      ! radial direction
  parameter (i1                = imax + 1)
  parameter (K_start_heat      = 3)       ! for isoflux:    after home many cells in x will you start heating
  parameter (x_start_heat      = 5)       ! for isothermal: lenght at which you start heating
  parameter (iwork             = 1)
  parameter (isave             = 1)
  parameter (nTab              = 2000)    ! 2000 for co2, 2500 for pH2
  parameter (nstep             = 100000)
  parameter (LoD               = 30)      ! pH2= 30
  parameter (periodic          = 1)       ! 1..periodic, 2..developing   
  parameter (Re                = 360)
  parameter (Pr                = 3.19457)
  parameter (Qwall             = 0.0)     ! Qwall=(qw*D/(k0*T0))
  parameter (isothermalBC      = 0)       ! isothermal wall: 1
  parameter (Tw                = 1.01)    ! scaled with the inflow temp (Tin=1)
  parameter (dTwall            = 0.05)
  parameter (pressIsoThermal   = 1)       ! constant pressure for BL cases (2 or 4)
  parameter (Fr_1              =  0.0)
  parameter (EOSmode           = 0)       ! 0..IG, 1..SCCO2/sPH2
  parameter (select_init       = 0)       ! 0..std initialization, 1..read inflow , 2..read restart file
  parameter (turbmod           = 2)       ! 0..laminar, 1..SA, 2..MK, 3..V2F, 4..SST
  parameter (ctheta            = 0.3)
  parameter (modVF             = 0)       ! 0..Original 1..LienKalitzin Time/Length scale VF model
  parameter (alphac            = 0.1)
  parameter (alphak            = 0.1)     ! for SA
  parameter (alphae            = 0.1) 
  parameter (alphav2           = 0.1) 
  parameter (modifDiffTerm     = 1)       ! 0..Standard, 1..invSLS, 2..Aupoix
     
end module mod_param
