      integer         i,k,imax,kmax,i1,k1,iwork,isave,px,rank,nTab,LoD
      integer         EOSmode,periodic,kmaxper,k1Old
      integer         nstep,K_start_heat,x_start_heat,select_init,systemSolve
      integer         turbmod,modifDiffTerm,modVF,profiling,isothermalBC,pressIsoThermal
      real*8          Re,Pr,Qwall,CFL,Tw,dTwall
      real*8          ctheta,Fr_1
      real*8    alphac,alphak,alphae,alphav2

      parameter (CFL               = 100.0)      ! CFL number / time step
     
      parameter (systemSolve       = 3)        ! 1: pipe, 2: channel 3: boundary layer ! dpdz now defined in mkgrid should be 1 for channel/BL, 4 for the pipe 

      parameter (px                = 4)        ! number of cores
      parameter (imax              = 96)       ! radial direction
      parameter (kmax              = 384/px)    ! 384/px)   ! axial direction
      parameter (kmaxper           = kmax*px/2)
      parameter (i1                = imax + 1)
      parameter (k1                = kmax + 1)

      parameter (k1old             = k1)    !previous kmax/px + 1  (to restart to a finer mesh in x)
      
      parameter (K_start_heat      = 3)    ! for isoflux:    after home many cells in x will you start heating
      parameter (x_start_heat      = 5)    ! for isothermal: lenght at which you start heating

      parameter (iwork             = 1)
      parameter (isave             = 1)

      parameter (nTab              = 2000)    ! 2000 for co2, 2500 for pH2
      parameter (nstep             = 100000)

      parameter (LoD               = 30)      ! pH2= 30
      parameter (periodic          = 2)               ! 1..periodic, 2..developing 
     
      parameter (Re                = 360)
      parameter (Pr                = 3.19457)

!     Q+=Qwall=(qw*D/(k0*T0))        
      parameter (Qwall             = 0.0)     !0.0) !2.4)
      parameter (isothermalBC      = 1)       ! isothermal wall: 1

      parameter (Tw                = 1.25)     ! scaled with the inflow temp (Tin=1)
      parameter (dTwall            = 0.05)
!     parameter (Tw                = 4)      ! case for Tw=100K
!     parameter (Tw                = 8)      ! case for Tw=200K  

      parameter (pressIsoThermal   = 2)       ! constant pressure for BL cases (2 or 4)
!     parameter (pressIsoThermal   = 2)       ! 2MPa: pH2_2MPa_table.dat
!     parameter (pressIsoThermal   = 4)       ! 4MPa: pH2_4MPa_table.dat 

!     1/Fr = Gr/(beta_o*T_o*Re**2.*Qwall)
      parameter (Fr_1              =  0.0) 
!      parameter (Fr_1              = -9.96)   ! case B
!      parameter (Fr_1              = -79.67)  ! case C
!      parameter (Fr_1              = -268.89) ! case D 
!      parameter (Fr_1              =  99.59)  ! downward
!      /(beta_o*T_o*Re**2.*Qwall))   

      parameter (EOSmode           = 1)             ! 0..IG, 1..SCCO2/sPH2

      parameter (select_init       = 1)            ! 0..std initialization, 1..read inflow , 2..read restart file
      parameter (turbmod           = 1)            ! 0..laminar, 1..SA, 2..MK, 3..V2F, 4..SST
      parameter (ctheta            = 0.3) 
      parameter (modVF             = 0)            ! 0..Original 1..LienKalitzin Time/Length scale VF model
      
      !MK: alphac=0.3, alphak=0.1 and alphae=0.1
      parameter (alphac            = 0.1)
      parameter (alphak            = 0.1) !0.05)         ! for SA
      parameter (alphae            = 0.1) !0.05) 
      parameter (alphav2           = 0.05) !0.03) 

      parameter (modifDiffTerm     = 0)              ! 0..Standard, 1..invSLS, 2..Aupoix
     

