
!>    ****************************************************
!!    Main file of the code
!!    This code simulates Super Critical Fluid in a
!!    heated pipe with constant heat flux
!!    ****************************************************
use mod_param
use mod_common
use mod_common2
use mod_math
use mod_eosmodels
use mod_tm
use mod_tm_mk
use mod_tm_sa
use mod_tm_sst
use mod_tm_vf

implicit none
include 'mpif.h'

integer :: rank,ierr,istart,noutput
real(8) :: bulk,stress,stime,time1
real(8) :: resC,resK,resE,resV2,resOm,resSA   
real(8) :: tempWall
real(8) :: start, finish
real(8), dimension(:), allocatable :: Win,ekmtin,dummy

call read_parameters()
call cpu_time(time1)
call mpi_init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)
call mpi_comm_size(MPI_COMM_WORLD,px,ierr)



resC=0;resV2=0;resK=0;resE=0;resOm=0;resSA=0;

dtmax = 1.e-3
dt = dtmax
istart = 1

!domain integers
kmax    = 200/px
kmaxper = kmax*px/2
k1      = kmax + 1
i1 = imax+1

allocate(Win(0:i1),ekmtin(0:i1))


!init memory
call initMem()
!numerical stuff
call init_transpose


!*********************!
!  Equation of State  !
!*********************!
if (EOSmode.eq.0) allocate(eos_model,    source=   IG_EOSModel(Re,Pr))
if (EOSmode.eq.1) allocate(eos_model,    source=Table_EOSModel(Re,Pr,2000,'tables/co2h_table.dat'))
if (EOSmode.eq.2) allocate(eos_model,    source=Table_EOSModel(Re,Pr,2499, 'tables/ph2_table.dat'))
call eos_model%init()

!*********************!
!   Turbulence model  !
!*********************!
if (turbmod.eq.0) allocate(turb_model,source=Laminar_TurbModel(i1, k1, imax, kmax,'lam'))
if (turbmod.eq.1) allocate(turb_model,source=     SA_TurbModel(i1, k1, imax, kmax,'SA'))
if (turbmod.eq.2) allocate(turb_model,source=init_MK_TurbModel(i1, k1, imax, kmax,'MK'))
if (turbmod.eq.3) allocate(turb_model,source=init_VF_TurbModel(i1, k1, imax, kmax,'VF'))
if (turbmod.eq.4) allocate(turb_model,source=    SST_TurbModel(i1, k1, imax, kmax,'SST'))
call turb_model%init()

!*********************!
!         Grid        !
!*********************!
call mkgrid(rank)

!*********************!
!    Initialization   !
!*********************!
if (select_init < 2) then
  call fkdat(rank)
  istart=1
else
  ! call loadRestart(istart,rank)
  istart = istart+1
endif
!periodic=1, turb flow generator,periodic=2, heated pipe
! if (periodic.ne.1) then
if (turbmod.eq.0) open(29,file ='pipe/Inflow_lam_00360.dat',form='unformatted')
if (turbmod.eq.1) open(29,file= 'pipe/Inflow_SA_00360.dat',form='unformatted')
if (turbmod.eq.2) open(29,file= 'pipe/Inflow_MK_00360.dat',form='unformatted')
if (turbmod.eq.3) open(29,file= 'pipe/Inflow_VF_00360.dat',form='unformatted')
if (turbmod.eq.4) open(29,file= 'pipe/Inflow_SST_00360.dat',form='unformatted')
  read(29) Win(:),dummy(:),dummy(:),dummy(:),dummy(:),dummy(:),ekmtin(:),dummy(:)
  close(29)
! endif
call turb_model%init_w_inflow(Re, systemsolve)
! cnew=0
! cold=0
! do k=0,k1
!   ! ekmt(:,k)=ekmtin(:)
!   wnew(:,k)=win(:)
! enddo
!write(*,*) cnew
! cnew = 0
call state_upd(cnew,rnew,ekm,ekh,temp,beta,istart,rank);

! stop
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! tempWall = 1.0
! if (isothermalBC.eq.1) then
!   tempWall = min(max(tempWall, (temp(imax,kmax)+temp(i1,kmax))*0.5),Tw)
!   call funcIsothermalEnthBC_upd(tempWall) ! calc. enthalpy at the wall (isothermal BC)
!   if (Qwall.ne.0) then
!     if (rank.eq.0) print '("Isothermal BC, Qwall should be 0  but it is ",f6.3,"... stopping")',Qwall
!     stop
!   else
!     if (rank.eq.0) print*,"*************SOLVING AN ISOTHERMAL WALL*************!"
!   endif
!   if (rank.eq.0) print '("temperature at the wall = ",f6.3," .")',tempWall
! endif
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ss

call bound_c(Tw,rank)
! cnew = 0
call turb_model%set_bc(ekm,rnew,walldist,centerBC,periodic,rank,px)

call state_upd(cnew,rnew,ekm,ekh,temp,beta,istart,rank);            !necessary to call it twice

rold = rnew

call calc_mu_eff(Unew,Wnew,rnew,ekm,ekmi,ekme,ekmt,ekmtin,rp,drp,dru,dz,walldist,rank) 
call bound_v(Unew,Wnew,Win,rank)
call chkdt(rank,istep)
call cpu_time(start)






!**************************
!       MAIN LOOP         !
!**************************

do istep=istart,nstep

  call calc_mu_eff(Unew,Wnew,rnew,ekm,ekmi,ekme,ekmt,ekmtin,rp,drp,dru,dz,walldist,rank) 
  call advanceC(resC,Unew,Wnew,Rnew,rank)
  call turb_model%advance_turb(uNew,wNew,rnew,ekm,ekmi,ekmk,ekmt,beta,temp, &
                     Ru,Rp,dru,drp,dz,walldist,alphak,alphae,alphav2, &
                     modifDiffTerm,rank,centerBC,periodic,resSA,resK, resV2)
  call bound_c(Tw,rank)
  if  (mod(istep,10) .eq. 0) then
    call output2d_upd2(rank,istep)
    call mpi_finalize(ierr)
stop
  endif
  call turb_model%set_bc(ekm,rnew,walldist,centerBC,periodic,rank,px)
  call state_upd(cnew,rnew,ekm,ekh,temp,beta,istep,rank)
  call advance(rank)
  call bound_m(dUdt,dWdt,wnew,rnew,Win,rank)
  call fillps(rank)
  call solvepois(p,Ru,Rp,dRu,dRp,dz,rank,centerBC)
  call correc(rank,1)
  call bound_v(Unew,Wnew,Win,rank)

  ! !ramping isothermal wall temperature
  ! if (mod(istep,2000).eq.0) then
  !   if (isothermalBC.eq.1 .AND.  tempWall.lt.Tw) then
  !     tempWall = min(tempWall+dTwall, Tw)
  !     call funcIsothermalEnthBC_upd(tempWall) ! calc. enthalpy at the wall (isothermal BC)
  !     if (rank.eq.0) print '("temperature at the wall ramped! = ",f6.3," .")',tempWall
  !   endif
  ! endif

  if  (mod(istep,10) .eq. 0)  call chkdiv(rank)
  call cmpinf(bulk,stress)
  call chkdt(rank,istep)
  if  (mod(istep,1000).eq.0)  call inflow_output_upd(rank)
  if  (mod(istep,1000).eq.0)  call output2d_upd2(rank,istep)
  noutput = 1
  if (rank.eq.0) then
    if (istep.eq.istart .or. mod(istep,noutput*20).eq.0) then
      write(6,'(A7,9A14)') 'istep'    ,'dt'      ,'bulk'   ,'stress' ,'cResid', &
                           'kineResid','epsResid','v2Resid','omResid','nuSAresid'
    endif
    if (istep.eq.istart .or. mod(istep,noutput).eq.0) then
      write(6,'(i7,9e14.5)') istep,dt,bulk,stress,resC,resK,resE,resV2,resOm,resSA
    endif
  end if
         
enddo
call cpu_time(finish)
print '("Time = ",f6.3," seconds.")',finish-start
call inflow_output_upd(rank)
call mpi_finalize(ierr)
stop
end




!>******************************************************************************************
!!      routine to estimate the effective viscosity
!!******************************************************************************************
subroutine calc_mu_eff(utmp,wtmp,rho,mu,mui,mue,mut,mutin,rp,drp,dru,dz,walldist,rank)
  use mod_param
  use mod_common2  
  implicit none
  real(8), dimension(0:i1,0:k1), intent(IN) :: utmp,wtmp,rho,mu,mui   
  real(8), dimension(0:i1),      intent(IN) :: mutin
  real(8), dimension(1:imax),    intent(IN) :: walldist
  real(8), dimension(0:i1),      intent(IN) :: Rp,dru,drp
  integer,                       intent(IN) :: rank
  real(8),                       intent(IN) :: dz
  real(8), dimension(0:i1,0:k1), intent(OUT):: mue,mut
  real(8), dimension(0:k1) :: tauw(0:k1)
  call turb_model%set_mut(utmp,wtmp,rho,mu,mui,walldist,rp,drp,dru,dz,mut)
  call turb_model%set_mut_bc(mut,periodic,px,rank)

  mue = mu + mut
end 


!>************************************************************************************
!!
!!     Performes time integration with second order
!!     Adams-Bashforth scheme, i.e
!!
!!
!!     n+1     n
!!     dU     U     - U                               n
!!     ---- = ------------ = 1.5*( -ADV + DIFF + Force)     -
!!     dt        dt                                   n-1
!!     0.5*( -ADV + DIFF + Force)
!!
!!     This scheme is weakly instabel for pure advection,
!!     and therefore a very small amount of physical diffusion
!!     is necessary.
!!     The timestep is limited (see routine chkdt)
!!
!!************************************************************************************
subroutine advanceC(resC,Utmp,Wtmp,Rtmp,rank)
  use mod_param
  use mod_common
  use mod_common2
  use mod_math
  implicit none
      
  real*8 dnew(0:i1,0:k1),tempArray(0:i1,0:k1),dimpl(0:i1,0:k1),tscl
  real*8 Utmp(0:i1,0:k1),Wtmp(0:i1,0:k1),Rtmp(0:i1,0:k1),sigmakSST(0:i1,0:k1)
  real*8 rho2(0:i1,0:k1), rho3(0:i1,0:k1), eknu(0:i1,0:k1),eknui(0:i1,0:k1),eknuk(0:i1,0:k1)
  real*8 cb3,Q,hbc
  integer rank,ierr
  real*8     a  (imax)
  real*8     b  (imax)
  real*8     c  (imax)
  real*8     rhs(imax)

  real*8 t1,t2,t3,t4,bc,scl
  real*8 term, adiffm, adiffp
  real*8 resC, sigmat
  sigmat = 0.9

  ! ------------------------------------------------------------------------
  !     ENTHAPLY ENTHAPLY ENTHAPLY ENTHAPLY
  ! ------------------------------------------------------------------------
  resC  = 0.0
  dnew=0.0; dimpl = 0.0;
  call advecc(dnew,dimpl,cnew,Utmp,Wtmp,Ru,Rp,dru,dz,i1,k1,rank,periodic,.true.)
  call diffc(dnew,cnew,ekh,ekhi,ekhk,ekmt,sigmat,Rtmp,Ru,Rp,dru,dz,rank,0)

  if (isothermalBC.eq.1) then
    !!!!!!!!! isothermal wall
    if (centerBC.eq.1) then
      do k=1,kmax
        do i=1,imax
          a(i) = -Ru(i-1)*(ekhi(i-1,k)+0.5*(ekmt(i,k)+ekmt(i-1,k))/sigmat)/(dRp(i-1)*Rp(i)*dru(i))/Rtmp(i,k)
          c(i) = -Ru(i  )*(ekhi(i  ,k)+0.5*(ekmt(i,k)+ekmt(i+1,k))/sigmat)/(dRp(i  )*Rp(i)*dru(i))/Rtmp(i,k)
          b(i) = (-a(i)-c(i) + dimpl(i,k) )/alphac        ! BUG
          rhs(i) = dnew(i,k) + (1-alphac)*b(i)*cnew(i,k)  ! BUG
        enddo

        i=1
        b(i)=b(i)+a(i)
         
        i=imax
        rhs(i) = dnew(i,k) - c(i)*cNew(i1,k) + (1-alphac)*b(i)*cNew(i,k)

        call matrixIdir(imax,a,b,c,rhs)
   
        do i=1,imax
          !resC = resC + (cnew(i,k) - rhs(i))**2.0
          resC = resC + ((cnew(i,k) - rhs(i))/(cnew(i,k)+1.0e-20))**2.0
          cnew(i,k) = max(rhs(i), 0.0)
        enddo
               
      enddo
            

    else
      if (rank.eq.0) print '("Isothermal boundary condition coded only for 1 wall.... stopping")'
      stop
    endif
     
  else
    !!!!!!!!! isoflux
    if (centerBC.eq.-1) then  ! wall both sides!!!!
      do k=1,kmax
        if (rank.eq.0.and.k.lt.K_start_heat) then
          Q=0.0
        else
          Q=Qwall
        endif

        do i=1,imax-1
          a(i) = -Ru(i-1)*(ekhi(i-1,k)+0.5*(ekmt(i,k)+ekmt(i-1,k))/sigmat)/(dRp(i-1)*Rp(i)*dru(i))/Rtmp(i,k)
          c(i) = -Ru(i  )*(ekhi(i  ,k)+0.5*(ekmt(i,k)+ekmt(i+1,k))/sigmat)/(dRp(i  )*Rp(i)*dru(i))/Rtmp(i,k)
          b(i) = (-a(i)-c(i) + dimpl(i,k) )/alphac        ! BUG
          rhs(i) = dnew(i,k) + (1-alphac)*b(i)*cnew(i,k)  ! BUG
        enddo
   
        i=1
        a(i)   = 0.0
        c(i)   = -Ru(i )*(ekhi(i ,k)+0.5*(ekmt(i,k)+ekmt(i+1,k))/sigmat)/(dRp(i )*Rp(i)*dru(i))/Rtmp(i,k)
        b(i)   =  (-a(i)-c(i) + dimpl(i,k) )/alphac
        rhs(i) = dnew(i,k) + Ru(i)*Q/(Re*Pr*Rtmp(i,k)*Rp(i)*dru(i)) + (1-alphac)*b(i)*cnew(i,k)

        i=imax
        a(i)   = -Ru(i-1)*(ekhi(i-1,k)+0.5*(ekmt(i,k)+ekmt(i-1,k))/sigmat)/(dRp(i-1)*Rp(i)*dru(i))/Rtmp(i,k)
        c(i)   =  0.0
        b(i)   =  (-a(i)-c(i) + dimpl(i,k) )/alphac
        rhs(i) = dnew(i,k) + Ru(i)*Q/(Re*Pr*Rtmp(i,k)*Rp(i)*dru(i)) + (1-alphac)*b(i)*cnew(i,k)

        call matrixIdir(imax,a,b,c,rhs)

        do i=1,imax
          !resC = resC + (cnew(i,k) - rhs(i))**2.0
          resC = resC + ((cnew(i,k) - rhs(i))/(cnew(i,k)+1.0e-20))**2.0
          cnew(i,k) = max(rhs(i), 0.0)
        enddo
      enddo
    else if (centerBC.eq.1) then
      do k=1,kmax
        if (rank.eq.0.and.k.lt.K_start_heat) then
          Q=0.0
        else
          Q=Qwall
        endif

        do i=1,imax-1
          a(i) = -Ru(i-1)*(ekhi(i-1,k)+0.5*(ekmt(i,k)+ekmt(i-1,k))/sigmat)/(dRp(i-1)*Rp(i)*dru(i))/Rtmp(i,k)
          c(i) = -Ru(i  )*(ekhi(i  ,k)+0.5*(ekmt(i,k)+ekmt(i+1,k))/sigmat)/(dRp(i  )*Rp(i)*dru(i))/Rtmp(i,k)
          b(i) = (-a(i)-c(i) + dimpl(i,k) )/alphac        ! BUG
          rhs(i) = dnew(i,k) + (1-alphac)*b(i)*cnew(i,k)  ! BUG
        enddo

        i=1
        b(i)=b(i)+a(i)
         
        i=imax
        a(i)   = -Ru(i-1)*(ekhi(i-1,k)+0.5*(ekmt(i,k)+ekmt(i-1,k))/sigmat)/(dRp(i-1)*Rp(i)*dru(i))/Rtmp(i,k)
        c(i)   =  0.0
        b(i)   =  (-a(i)-c(i) + dimpl(i,k) )/alphac
        rhs(i) = dnew(i,k) + Ru(i)*Q/(Re*Pr*Rtmp(i,k)*Rp(i)*dru(i)) + (1-alphac)*b(i)*cnew(i,k)
        call matrixIdir(imax,a,b,c,rhs)
   
        do i=1,imax
          !resC = resC + (cnew(i,k) - rhs(i))**2.0
          resC = resC + ((cnew(i,k) - rhs(i))/(cnew(i,k)+1.0e-20))**2.0
          cnew(i,k) = max(rhs(i), 0.0)
        enddo
      enddo
    endif
  endif


  if (periodic.eq.1) then
    cnew = 0.0; resC=0.0;
  endif



end


!>*************************************************************************************
!!
!!     Performes time integration with second order 
!!     Adams-Bashforth scheme, i.e
!!     
!!     
!!     n+1     n
!!     dU     U     - U                               n
!!     ---- = ------------ = 1.5*( -ADV + DIFF + Force)     -
!!     dt        dt                                   n-1
!!     0.5*( -ADV + DIFF + Force)
!!     
!!     This scheme is weakly instabel for pure advection,
!!     and therefore a very small amount of physical diffusion
!!     is necessary.
!!     The timestep is limited (see routine chkdt)
!!*************************************************************************************
subroutine advance(rank)
  use mod_param
  use mod_common
  use mod_math
  implicit none
      
  integer rank
  real*8, dimension(imax)   :: a, b, c, rhs
  real*8, dimension(imax-1) :: au, bu, cu, rhsu
  real*8 dnew(0:i1,0:k1)
  real*8 dif,alpha,rhoa,rhob,rhoc

  dif=0.0

  !>********************************************************************
  !!     CALCULATE advection, diffusion and Force in r-direction
  !!     at the old(n-1) and new(n) timelevels
  !!********************************************************************
  dnew = 0.0
  call advecu(dnew,Unew,Wnew,Rnew,Ru,Rp,dru,drp,dz,i1,k1) ! new
  call diffu (dnew,Unew,Wnew,ekme,Ru,Rp,dru,drp,dz,i1,k1,dif,numDomain) ! new

  if (centerBC == -1) then
    do k=1,kmax
      do i=1,imax-1
        au(i) = -dt*ekme(i  ,k)*Rp(i  )/(dRp(i)*Ru(i)*dru(i  ))
        cu(i) = -dt*ekme(i+1,k)*Rp(i+1)/(dRp(i)*Ru(i)*dru(i+1))
        bu(i) = -au(i)-cu(i)
        rhoa = 0.5*(rnew(i  ,k)+rnew(i-1,k))
        rhoc = 0.5*(rnew(i+1,k)+rnew(i+2,k))
        rhob = 0.5*(rnew(i+1,k)+rnew(i  ,k))
        au(i) = au(i)/rhoa
        bu(i) = bu(i)/rhob + 1.0
        cu(i) = cu(i)/rhoc
      enddo
   
      i = imax-1;    cu(i) = 0.0
      i = 1;
      au(i+1) = 0.0           ! BC wall
   
      do i=1,imax-1
        rhsu(i) = dt*dnew(i,k) + Unew(i,k)*(Rnew(i+1,k)+Rnew(i,k))*0.5
      enddo
   
      call matrixIdir(imax-1,au,bu,cu,rhsu)
      do i=1,imax-1
        dUdt(i,k)=rhsu(i)
      enddo
    enddo
  elseif (centerBC == 1) then
    do k=1,kmax
      do i=1,imax-1
        au(i) = -dt*ekme(i  ,k)*Rp(i  )/(dRp(i)*Ru(i)*dru(i  ))
        cu(i) = -dt*ekme(i+1,k)*Rp(i+1)/(dRp(i)*Ru(i)*dru(i+1))
        bu(i) = -au(i)-cu(i)
        rhoa = 0.5*(rnew(i  ,k)+rnew(i-1,k))
        rhoc = 0.5*(rnew(i+1,k)+rnew(i+2,k))
        rhob = 0.5*(rnew(i+1,k)+rnew(i  ,k))
        au(i) = au(i)/rhoa
        bu(i) = bu(i)/rhob + 1.0
        cu(i) = cu(i)/rhoc
      enddo

      i = imax-1;    cu(i) = 0.0
      i = 1;
      bu(i) = bu(i) + au(i)    ! BC at center: Neumann: cancel coeff a
   
      do i=1,imax-1
        rhsu(i) = dt*dnew(i,k) + Unew(i,k)*(Rnew(i+1,k)+Rnew(i,k))*0.5
      enddo
   
      call matrixIdir(imax-1,au,bu,cu,rhsu)
      do i=1,imax-1
        dUdt(i,k)=rhsu(i)
      enddo
    enddo
  endif

  !********************************************************************
  !     CALCULATE advection, diffusion and Force in z-direction
  !     at the old(n-1) and new(n) timelevels
  !********************************************************************

  dnew = 0.0
  call advecw(dnew,Unew,Wnew,Rnew,Ru,Rp,dru,dz,ekm,peclet)
  call diffw (dnew,Unew,Wnew,ekme,Ru,Rp,dru,drp,dz,i1,k1,dif,numDomain)  ! new

  if (periodic.eq.1) dnew = dnew + dpdz

  if (Qwall.ne.0) then
    dnew(1:imax,1:kmax) = dnew(1:imax,1:kmax) + 0.5*(Rnew(1:imax,1:kmax)+Rnew(1:imax,2:kmax+1))*Fr_1
  endif

  do k=1,kmax
    do i=1,imax
      a(i) = -0.25*dt*(ekme(i,k)+ekme(i,k+1)+ekme(i-1,k)+ekme(i-1,k+1))*Ru(i-1)/(dRp(i-1)*Rp(i)*dru(i))
      c(i) = -0.25*dt*(ekme(i,k)+ekme(i,k+1)+ekme(i+1,k)+ekme(i+1,k+1))*Ru(i  )/(dRp(i  )*Rp(i)*dru(i))
      b(i) = -a(i)-c(i)
      rhoa = 0.5*(rnew(i-1,k+1)+rnew(i-1,k))
      rhoc = 0.5*(rnew(i+1,k+1)+rnew(i+1,k))
      rhob = 0.5*(rnew(i  ,k+1)+rnew(i  ,k))
      a(i) = a(i)/rhoa
      b(i) = b(i)/rhob + 1.0
      c(i) = c(i)/rhoc
    enddo

    i=imax;    b(i) = b(i) - c(i)    ! BC at wall: zero vel: subtract one c
    i = 1;     b(i) = b(i) + centerBC*a(i)     ! BC at wall: zero vel: subtract one a

    do i=1,imax
      rhs(i) = dt*dnew(i,k) + Wnew(i,k)*(Rnew(i,k)+Rnew(i,k+1))*0.5
    enddo

    call matrixIdir(imax,a,b,c,rhs)
    do i=1,imax
      dWdt(i,k)=rhs(i)
    enddo
  enddo


end




!<*************************************************************************************
!!
!!  bound_c equation
!!
!!*************************************************************************************
subroutine bound_c(Twall,rank)
  use mod_param
  use mod_common
  use mod_common2
  implicit none
  include 'mpif.h'
  real(8), intent(IN) :: Twall
  integer, intent(IN) :: rank
  real*8 tmpShift(0:i1)
  real(8) enth_wall

  ! ISOTHERMAL
  if (isothermalBC.eq.1) then
    call eos_model%set_w_temp(Twall, "H", enth_wall)
    !pipe/bl
    if (centerBC.eq.1) then
      do k=0,k1
        if ((k+rank*kmax)*dz.lt.x_start_heat) then
          cnew(i1,k) = cnew(imax,k)
        else
          cnew(i1,k) = 2.0*enth_wall - cnew(imax,k)
        endif
      enddo
    !channel (TODO!!!)
    else
      if (rank.eq.0) print '("Isothermal boundary condition coded only for 1 wall.... stopping")'
      stop
    endif
  ! ISOFLUX
  else
    do k=0,k1
      if (rank.eq.0.and.k.lt.K_start_heat) then
        cnew(i1,:) = cnew(imax,:)
      else
        !pipe/bl
        call funcNewtonSolve_upd(cnew(i1,k), cnew(imax,k))
        !channel
        if (centerBC.eq.-1) call funcNewtonSolve_upd(cnew(0,k), cnew(1,k))
      endif
    enddo
  endif

  if (centerBC.eq.1) then     
    !center line BC
    cnew(0,:)    = cnew(1,:)
  endif
  call shiftf(cnew,   tmpShift,rank);       cnew(:,0) = tmpShift(:);
  call shiftb(cnew,   tmpShift,rank);      cnew(:,k1) = tmpShift(:);

  if (periodic.eq.1) return
  if (rank.eq.0) then
    cnew(:,0) = 0.0
  endif
  if (rank.eq.px-1) then
    cnew(:,k1) = 2.0*   cnew(:,kmax) -    cnew(:,kmax-1)
  endif  
end

!>*************************************************************************************
!!      bound_v(Ubound,Wbound,Win,rank)
!!
!!*************************************************************************************
subroutine bound_v(Ubound,Wbound,Win,rank)

  use mod_param
  use mod_common
  implicit none  
  include 'mpif.h'

  integer rank
  real*8 Ubound(0:i1,0:k1), Wbound(0:i1,0:k1)
  real*8 tmp(0:i1),Win(0:i1)

  ! channel
  if (centerBC.eq.-1) then 
    do k=0,k1
      Ubound(1,k)    =   0.0
      Ubound(0,k)    = - Ubound(2,k)
      Ubound(imax,k) =   0.0
      Ubound(i1,k)   = - Ubound(imax-1,k)
      Wbound(0,k)   = - Wbound(1,k)
      Wbound(i1,k)  = - Wbound(imax,k)
    enddo
  !pipe/BL
  else
    do k=0,k1
      Ubound(0,k)    =   Ubound(1,k)
      Ubound(imax,k) =   0.0
      Ubound(i1,k)   = - Ubound(imax-1,k)
      Wbound(0,k)   =   Wbound(1,k)
      Wbound(i1,k)  = - Wbound(imax,k)
    enddo
  endif

  call shiftf(Ubound,tmp,rank);     Ubound(:,0)  = tmp(:);
  call shiftf(Wbound,tmp,rank);     Wbound(:,0)  = tmp(:);
  call shiftb(Ubound,tmp,rank);     Ubound(:,k1) = tmp(:);
  call shiftb(Wbound,tmp,rank);     Wbound(:,k1) = tmp(:);

  !developing
  if (periodic.eq. 1) return
  if (rank.eq.0) then
    Ubound(:,0) = 0.0
    Wbound(:,0) = Win(:)
  endif
  if (rank.eq.px-1)then
    ubound(:,k1) = 2.*ubound(:,kmax)-ubound(:,kmax-1)
    wbound(:,k1) = 2.*wbound(:,kmax)-wbound(:,kmax-1)
  endif

end




!>*************************************************************************************
!!
!!     Subroutine bound sets the boundary conditions for all variables,
!!     except for the diffusion coefficients. These are set in submod.
!!     The common boundary conditions for the pressure are set in mkgrid.
!!
!!     Set boundary conditions for j=0 and j=j1. Because of the
!!     fact that the index j denotes the tangential direction,
!!     we have to set the values at j=0 equal to the values at
!!     j=jmax and the values at j=j1 equal to the values at j=1.
!!
!!*************************************************************************************
subroutine bound_m(Ubound,Wbound,W_out,Rbound,Win,rank)
  use mod_param
  use mod_common
  implicit none

      include 'mpif.h'
  character*5 inflow
  !
  real*8  W_out(0:i1,0:k1)
  integer rank,ierr
  real*8  y1,y2,y3,y4
  real*8  Ubound(0:i1,0:k1),Vbound(0:i1,0:k1), &
          Wbound(0:i1,0:k1),Rbound(0:i1,0:k1),Win(0:i1)
  real*8 Ub,flux,flux_tot,deltaW,rhob,wfunc,wr(1:imax)
  real*8 tmp(0:i1)
  integer ib,ie,kb,ke

  !channel
  if (centerBC.eq.-1) then ! channal bc
    do k=0,k1
      Ubound(1,k)    =   0.0
      Ubound(0,k)    = - Ubound(2,k)
      Ubound(imax,k) =   0.0
      Ubound(i1,k)   = - Ubound(imax-1,k)
      Wbound(0,k)   = - Wbound(1,k)
      Wbound(i1,k)  = - Wbound(imax,k)
    enddo
  !pipe/BL
  else
    do k=0,k1
      Ubound(0,k)    =   Ubound(1,k)
      Ubound(imax,k) =   0.0
      Ubound(i1,k)   = - Ubound(imax-1,k)

      Wbound(0,k)   =   Wbound(1,k)
      Wbound(i1,k)  = - Wbound(imax,k)
    enddo
  endif

  call shiftf(Ubound,tmp,rank);     Ubound(:,0)  = tmp(:);
  call shiftf(Wbound,tmp,rank);     Wbound(:,0)  = tmp(:);
  call shiftb(Ubound,tmp,rank);     Ubound(:,k1) = tmp(:);
  call shiftb(Wbound,tmp,rank);     Wbound(:,k1) = tmp(:);


  !developing
  if (periodic.eq.1) return

  if (rank.eq.0) then
    Rbound(:,0) = 1.0
    Ubound(:,0) = 0.0
    Wbound(:,0) = Win(:)
  endif

  wr = 0
  Ub = 0.
  flux = 0.0
  if (rank.eq.px-1) then
    Ub = 0.
    do i=1,imax
      wr(i) = W_out(i,kmax)
      Ub = max(Ub,2.0*Wbound(i,kmax)/(Rbound(i,kmax)+Rbound(i,k1)))
    enddo

    do i=0,i1
      Wbound(i,kmax) = 2.0*W_out(i,kmax-1) - W_out(i,kmax-2)
      Wbound(i,kmax) = Wbound(i,kmax)*0.5*(Rbound(i,kmax)+Rbound(i,k1))
    enddo

    Wbound(i1,kmax) = -Wbound(imax,kmax)
    Wbound(0,kmax)  = centerBC*Wbound(1,kmax)
  endif

  !     compute drho/dt*dvol
  do k=1,kmax
    do i=1,imax
      flux = flux - (rnew(i,k)-rold(i,k))/dt*Rp(i)*dru(i)*dz
    enddo
  enddo

  !     compute mf in
  if (rank.eq.0)then
    do i=1,imax
      flux = flux + Wbound(i,0)*dru(i)*rp(i)
    enddo
  endif

  if (rank.eq.px-1)then
    Ub = 0
    wfunc = 0
    do i=1,imax
      flux = flux - Wbound(i,kmax)*dru(i)*rp(i)
      wfunc = wfunc + wr(i)*dru(i)*rp(i) ! based on averaged outflow velocity
    enddo
  endif

  !      write(*,*) "delta flux: ", flux

  call mpi_allreduce(flux,flux_tot,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)


  if (rank.eq.px-1)then
    deltaW = (flux_tot)/wfunc
    !     write(*,*)'------------------------',deltaW
    do i=1,imax
      Wbound(i,kmax) = Wbound(i,kmax) + deltaW*wr(i) ! based on averaged outflow velocity
    enddo
    Wbound(i1,kmax) = -Wbound(imax,kmax)
    Wbound(0,kmax)  = centerBC*Wbound(1,kmax)

  !         flux = 0
  !         do i=1,imax
  !            flux = flux - Wbound(i,kmax)*dru(i)*rp(i)
  !         enddo
  !         write(*,*) "flux out: ", flux

  endif




  return
end






!>*************************************************************************************
!!
!!           fkdat(rank)
!!
!!*************************************************************************************
subroutine fkdat(rank)
  use mod_param
  use mod_common
  implicit none

  integer rank
  real*8 yplus,t1,t2,t3,in,chl,ran,Wvel,delta,gridSize
  real(8), dimension(0:i1) :: dummy
  character*5 inflow

  delta=0.5
  t1=3.5
  t2=5.
  t3=360.00
  Unew =0.
  Uold =0.
  Cnew =0.

  gridSize = y_fa(imax)
  rold =1.
  rnew =1.


  if (select_init.eq.1) then
    !initialized from inflow
    if (rank.eq.0)  write(*,*) 'Initializing flow with inflow = ', select_init

    do k=0,k1
      if (turbmod.eq.0) open(29,file= 'pipe/Inflow_lam_00360.dat',form='unformatted')
      if (turbmod.eq.1) open(29,file= 'pipe/Inflow_SA_00360.dat',form='unformatted')
      if (turbmod.eq.2) open(29,file= 'pipe/Inflow_MK_00360.dat',form='unformatted')
      if (turbmod.eq.3) open(29,file= 'pipe/Inflow_VF_00360.dat',form='unformatted')
      if (turbmod.eq.4) open(29,file= 'pipe/Inflow_SST_00360.dat',form='unformatted')
      read(29) Wnew(:,k),dummy,dummy,dummy,dummy,dummy,ekmt(:,k),dummy
      close(29)
    enddo
  else
    !initialized from scratch values
    if (rank.eq.0)  write(*,*) 'Initializing flow from scratch = ', select_init

    do i=1,imax
                            
      if (numDomain.eq.-1) then
        if (centerBC.eq.-1) then     ! channel
          Wnew(i,:)  = Re*dpdz*y_cv(i)*0.5*(gridSize-y_cv(i))
        elseif (centerBC.eq.1) then ! boundary layer
          Wnew(i,:)  = Re*dpdz*0.5*((gridSize*gridSize)-(y_cv(i)*y_cv(i)))
        endif
      else                           ! pipe
        Wnew(i,:)  = Re/6*3/2.*(1-(y_cv(i)/0.5)**2)
      endif
      ! TODO: BoundaryLayer: Analytical solution?

      ! knew(i,:)  = 0.1
      ! enew(i,:)  = 1.0
      ! omnew(i,:) = 0.001
      ! v2new(i,:) = 2./3.*knew(i,:)
      ! nuSAnew(i,:) = 0.001
    enddo
  endif


end
      


!>*************************************************************************************
!!           cmpinf(Bulk,Stress)
!!
!!*************************************************************************************
subroutine cmpinf(Bulk,Stress)
  use mod_param
  use mod_common
  implicit none
     
      include 'mpif.h'
  integer ierr
  real*8 waver(imax),waver2(imax)
     
  real*8  Bulk,Stress
     
  !     --- Initialization ---
   
  Waver = 0.0
     
  !     -------------------------------------------start i,j,k-loop
  do  i=1,imax
    do k=1,kmax
      Waver(i) = Waver(i) + wnew(i,k)
    enddo
  enddo
  call mpi_allreduce(waver,waver2,imax, &
    mpi_real8,mpi_sum,mpi_comm_world,ierr)
  !     -------------------------------------------end i,j,k-loop
  waver = waver2/(kmax*px)

     
  !     Stress =  Waver(imax) /(1.0-rp(imax))/Re

  !     *** use MIDPOINT INTEGRATION RULE ***
     
  Bulk = 0.0
  do i=1,imax
    Bulk = Bulk + 8.*Waver(i) *Rp(i)* dru(i)
  enddo
  Stress =  Waver(imax)/wallDist(imax)/Re
  return
end








