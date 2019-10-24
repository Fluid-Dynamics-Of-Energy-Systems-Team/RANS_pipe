!!    ****************************************************
!!    Main file of the code
!!    This code simulates Super Critical Fluid in a
!!    heated pipe/channel/bl with constant heat flux/isothermal wall
!!    ****************************************************

use mod_param
use mod_math
use mod_mesh
use mod_common
use mod_eos
use mod_tm
use sa_tm
use sst_tm
use mk_tm
use vf_tm
implicit none
include 'mpif.h'

integer ::  rank,ierr,istart,noutput
real(8) ::  bulk,stress,stime,time1
real(8) ::  resC,resK,resE,resV2,resOm,resSA   
real(8) ::  start, finish
resC=0;resV2=0;resK=0;resE=0;resOm=0;resSA=0;

!read parameters
call read_parameters()

!initalize mpi
call cpu_time(time1)
call mpi_init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)
call mpi_comm_size(MPI_COMM_WORLD,px,ierr)

!domain integers
kmax    = kelem/px
kmaxper = kmax*px/2
k1      = kmax + 1
k1old   = k1
i1 = imax+1
Mt=imax/px
Nx=kmax*px
Mx=kmax
Nt=imax

!***************************!
!      INITIALIZATION       !
!***************************!

!initialize variables
call initMem()

!initialize EOS
if (EOSmode.eq.0) allocate(eos_model,    source=IG_EOSModel(Re,Pr))
if (EOSmode.eq.1) allocate(eos_model,    source=Table_EOSModel(Re,Pr,2000, 'tables/co2h_table.dat'))
if (EOSmode.eq.2) allocate(eos_model,    source=Table_EOSModel(Re,Pr,2499, 'tables/ph2_table.dat'))
call eos_model%init()

!initialize turbomodel
if (turbmod.eq.0) allocate(turb_model,source=Laminar_TurbModel(i1, k1, imax, kmax,'lam'))
if (turbmod.eq.1) allocate(turb_model,source=     SA_TurbModel(i1, k1, imax, kmax,'SA'))
if (turbmod.eq.2) allocate(turb_model,source=init_MK_TurbModel(i1, k1, imax, kmax,'MK'))
if (turbmod.eq.3) allocate(turb_model,source=init_VF_TurbModel(i1, k1, imax, kmax,'VF'))
if (turbmod.eq.4) allocate(turb_model,source=    SST_TurbModel(i1, k1, imax, kmax,'SST'))
call turb_model%init()

!initialize numerical
call init_transpose

!initialize grid
call mkgrid(rank)


dtmax = 1.e-3
dt = dtmax
istart = 1

!initialize solution 
call initialize_solution(rank,wnew,unew,cnew,ekmt,win,ekmtin,i1,k1,y_fa,y_cv,dpdz,Re,systemsolve,select_init)
call state_upd(cnew,rnew,ekm,ekmi,ekmk,ekh,ekhi,ekhk,cp,cpi,cpk,temp,beta)
call bound_c(cnew, Tw, Qwall, dz, centerBC,rank)
call turb_model%set_bc(ekm,rnew,walldist,centerBC,periodic,rank,px)
call state_upd(cnew,rnew,ekm,ekmi,ekmk,ekh,ekhi,ekhk,cp,cpi,cpk,temp,beta) ! necessary to call it twice

rold = rnew
call calc_mu_eff(Unew,Wnew,rnew,ekm,ekmi,ekme,ekmt,ekmtin,rp,drp,dru,dz,walldist,rank) 
call bound_v(Unew,Wnew,Win,centerBC,rank)
call chkdt(rank,istep)
call cpu_time(start)


!***************************!
!        MAIN LOOP          !
!***************************!

do istep=istart,nstep

  call calc_mu_eff(Unew,Wnew,rnew,ekm,ekmi,ekme,ekmt,ekmtin,rp,drp,dru,dz,walldist,rank) 
  call advanceC(resC,Unew,Wnew,Rnew,rank)
  call turb_model%advance_turb(uNew,wNew,rnew,ekm,ekmi,ekmk,ekmt,beta,temp,           &
                               Ru,Rp,dru,drp,dz,walldist,alphak,alphae,alphav2,       &
                               modifDiffTerm,rank,centerBC,periodic,resSA,resK, resV2)
  call bound_c(Cnew, Tw, Qwall, dz, centerBC,rank)
  call turb_model%set_bc(ekm,rnew,walldist,centerBC,periodic,rank,px)
  call state_upd(cnew,rnew,ekm,ekmi,ekmk,ekh,ekhi,ekhk,cp,cpi,cpk,temp,beta);
  call advance(rank)
  call bound_m(dUdt,dWdt,wnew,rnew,Win,rank)
  call fillps(rank)
  call solvepois(p,Ru,Rp,dRu,dRp,dz,rank,centerBC)
  call correc(rank,1)
  call bound_v(Unew,Wnew,Win,centerBC,rank)

  if   (mod(istep,10) .eq. 0) call chkdiv(rank)
  call cmpinf(bulk,stress)
  call chkdt(rank,istep)
  if  ((mod(istep,1000).eq.0).and.(periodic .eq.1)) call inflow_output_upd(rank)
  if   (mod(istep,1000).eq.0)                       call output2d_upd2(rank,istep)

  noutput = 100
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
call output2d_upd2(rank,istep)
call mpi_finalize(ierr)
stop
end

!***********************************************************************************************************************************
!***********************************************************************************************************************************
!***********************************************************************************************************************************
!***********************************************************************************************************************************


!!******************************************************************************************
!!      routine to estimate the effective viscosity
!!******************************************************************************************
subroutine calc_mu_eff(utmp,wtmp,rho,mu,mui,mue,mut,mutin,rp,drp,dru,dz,walldist,rank)
  use mod_param
  use mod_tm

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
end subroutine calc_mu_eff

!!*************************************************************************************
!!  Apply the boundary conditions for the energy equation
!!*************************************************************************************
subroutine bound_c(c, Twall, Qwalll, dz, centerBC,rank)
  use mod_param
  use mod_eos
  implicit none
  include 'mpif.h'
  real(8),                       intent(IN) :: Twall, Qwalll,dz
  integer,                       intent(IN) :: centerBC, rank
  real(8), dimension(0:i1,0:k1), intent(OUT):: c
  real(8), dimension(0:i1) :: tmp
  real(8)                  :: enth_wall  
  
  !isothermal
  if (isothermalBC.eq.1) then
    call eos_model%set_w_temp(Twall, "H", enth_wall)
    do k=0,k1
      if ((k+rank*kmax)*dz.lt.x_start_heat) then
        c(i1,k) = c(imax,k)
      else
        if (centerBC.eq.-1) c(0,k) =  2.0*enth_wall - c(1,k)    !channel
                            c(i1,k) = 2.0*enth_wall - c(imax,k) !pipe/bl
      endif
    enddo

  !isoflux
  else
    do k=0,k1
      if (rank.eq.0.and.k.lt.K_start_heat) then
        c(i1,k) = c(imax,k)
      else
        if (centerBC.eq.-1) call funcNewtonSolve_upd(c(0,k), c(1,k))     !channel
                            call funcNewtonSolve_upd(c(i1,k), c(imax,k)) !pipe/bl
      endif
    enddo
  endif

  !pipe/bl
  if (centerBC.eq.1) then     
    c(0,:)    = c(1,:)
  endif

  call shiftf(c,tmp,rank);  c(:,0) = tmp(:);
  call shiftb(c,tmp,rank); c(:,k1) = tmp(:);

  !developing 
  if (periodic.eq.1) return
  if (rank.eq.0) then
    c(:,0) = 0.0
  endif
  if (rank.eq.px-1) then
    c(:,k1) = 2.0*   c(:,kmax) -    c(:,kmax-1)
  endif  
end subroutine bound_c

!!*************************************************************************************
!!  Apply the boundary conditions for the velocity
!!*************************************************************************************
subroutine bound_v(u,w,win,centerBC,rank)
  use mod_param
  implicit none  
  include 'mpif.h'
  integer,                       intent(IN) :: centerBC, rank
  real(8), dimension(0:i1),      intent(IN) :: win
  real(8), dimension(0:i1,0:k1), intent(OUT):: u, w
  real(8), dimension(0:i1)                  :: tmp

  ! channel
  if (centerBC.eq.-1) then 
    do k=0,k1
      u(1,k)    =   0.0
      u(0,k)    = - u(2,k)
      u(imax,k) =   0.0
      u(i1,k)   = - u(imax-1,k)
      w(0,k)    = - w(1,k)
      w(i1,k)   = - w(imax,k)
    enddo
  !pipe/BL
  else
    do k=0,k1
      u(0,k)    =   u(1,k)
      u(imax,k) =   0.0
      u(i1,k)   = - u(imax-1,k)
      w(0,k)    =   w(1,k)
      w(i1,k)   = - w(imax,k)
    enddo
  endif

  call shiftf(u,tmp,rank);     u(:,0)  = tmp(:);
  call shiftf(w,tmp,rank);     w(:,0)  = tmp(:);
  call shiftb(u,tmp,rank);     u(:,k1) = tmp(:);
  call shiftb(w,tmp,rank);     w(:,k1) = tmp(:);
  if (periodic.eq. 1) return

  !developing
  if (rank.eq.0) then
    u(:,0) = 0.0
    w(:,0) = win(:)
  endif
  if (rank.eq.px-1)then
    u(:,k1) = 2.*u(:,kmax)-u(:,kmax-1)
    w(:,k1) = 2.*w(:,kmax)-w(:,kmax-1)
  endif
end subroutine bound_v

!!*************************************************************************************
!!   Apply the boundary conditions for the velocity using the mass flux
!!*************************************************************************************
subroutine bound_m(Ubound,Wbound,W_out,Rbound,W_in,rank)
  use mod_param
  use mod_mesh
  use mod_common
  implicit none
  include 'mpif.h'
  integer,                       intent(IN) :: rank
  real(8), dimension(0:i1,0:k1), intent(IN) :: rbound
  real(8), dimension(0:i1),      intent(IN) :: W_in  
  real(8), dimension(0:i1,0:k1), intent(OUT):: ubound, wbound, w_out
  real(8), dimension(0:i1) :: tmp
  real(8), dimension(1:imax) :: wr
  integer :: ierr
  real(8) :: Ub,flux,flux_tot,deltaW,wfunc
  
  call bound_v(ubound,wbound,W_in,centerBC,rank)

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
    Wbound(0,kmax)  = centerBC*Wbound(1,kmax) !either symmetry or wall
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

  call mpi_allreduce(flux,flux_tot,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)

  if (rank.eq.px-1)then
    deltaW = (flux_tot)/wfunc
    do i=1,imax
      Wbound(i,kmax) = Wbound(i,kmax) + deltaW*wr(i) ! based on averaged outflow velocity
    enddo
    Wbound(i1,kmax) = -Wbound(imax,kmax)
    Wbound(0,kmax)  = centerBC*Wbound(1,kmax) !either symmetry or wall
    ! flux = 0
    ! do i=1,imax
    !    flux = flux - Wbound(i,kmax)*dru(i)*rp(i)
    ! enddo
    ! write(*,*) "flux out: ", flux
  endif
end subroutine bound_m

!!*************************************************************************************
!!           initialize the solution
!!*************************************************************************************

subroutine initialize_solution(rank, w, u,c, mut, win, mutin, i1,k1, y_fa, y_cv, dpdz,Re, systemsolve, select_init)
  use mod_tm
  implicit none
  integer,                        intent(IN) :: rank,systemsolve,i1,k1,select_init
  real(8),                        intent(IN) :: dpdz,Re
  real(8), dimension(0:i1),       intent(IN) :: y_cv,y_fa
  real(8), dimension(0:i1, 0:k1), intent(OUT):: w,u,c,mut
  real(8), dimension(0:i1),       intent(OUT):: win,mutin
  real(8), dimension(0:i1) :: dummy
  character(len=5)  :: Re_str
  character(len=7)  :: case
  integer           :: Re_int, i,k, imax
  real(8)           :: gridSize
    
  imax = i1-1
  ! inialize from inflow profile=
  if (select_init.eq.1) then
    if (rank.eq.0) write(*,*) 'Initializing flow with inflow'
    if (systemsolve .eq. 1) case = "pipe"
    if (systemsolve .eq. 2) case = "channel"
    if (systemsolve .eq. 3) case = "bl"
    Re_int = int(Re)
    write(Re_str,'(I5.5)') Re_int
    open(29,file =trim(case)//'/Inflow_'//trim(turb_model%name)//'_'//Re_str//'.dat',form='unformatted')
    read(29) win(:),dummy,dummy,dummy,dummy,dummy,mutin(:),dummy
    close(29)
    do k=0,k1
      w(:,k)  = win(:)
      mut(:,k)= mutin(:)
      c(:,k)= 0
    enddo
    call turb_model%init_w_inflow(Re, systemsolve)

  !initialize with laminar analytical solution
  else
    if (rank.eq.0) write(*,*) 'Initializing flow from scratch'
    gridSize = y_fa(imax)
    do i=1,imax         
      if (systemsolve.eq.1) w(i,:)  = Re/6*3/2.*(1-(y_cv(i)/0.5)**2)                      !pipe
      if (systemsolve.eq.2) w(i,:)  = Re*dpdz*y_cv(i)*0.5*(gridSize-y_cv(i))              !channel
      if (systemsolve.eq.3) w(i,:)  = Re*dpdz*0.5*((gridSize*gridSize)-(y_cv(i)*y_cv(i))) !bl
    enddo
  endif
end subroutine initialize_solution



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
  use mod_math
  use mod_mesh
  use mod_common
  implicit none
  real(8), dimension(0:i1,0:k1), intent(IN) :: Utmp, Wtmp, Rtmp
  integer,                       intent(IN) :: rank
  real(8),                       intent(OUT) :: resC
  real(8), dimension(0:i1,0:k1) :: dnew,dimpl
  real(8), dimension(imax)      :: a,b,c,rhs
  real(8)                       :: sigmat,Q
  
  sigmat = 0.9 !turbulent prandtl
  resC   = 0.0; dnew   = 0.0; dimpl = 0.0;

  call advecc(dnew,dimpl,cnew,Utmp,Wtmp,Ru,Rp,dru,dz,i1,k1,rank,periodic,.true.)
  call diffc(dnew,cnew,ekh,ekhi,ekhk,ekmt,sigmat,Rtmp,Ru,Rp,dru,dz,rank,0)

  !---------------------------------------   ISOTHERMAL
  if (isothermalBC.eq.1) then
    !pipe/bl
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
    !channel
    else
      do k=1,kmax
        do i=1,imax
          a(i) = -Ru(i-1)*(ekhi(i-1,k)+0.5*(ekmt(i,k)+ekmt(i-1,k))/sigmat)/(dRp(i-1)*Rp(i)*dru(i))/Rtmp(i,k)
          c(i) = -Ru(i  )*(ekhi(i  ,k)+0.5*(ekmt(i,k)+ekmt(i+1,k))/sigmat)/(dRp(i  )*Rp(i)*dru(i))/Rtmp(i,k)
          b(i) = (-a(i)-c(i) + dimpl(i,k) )/alphac        ! BUG
          rhs(i) = dnew(i,k) + (1-alphac)*b(i)*cnew(i,k)  ! BUG
        enddo

        i=1
        rhs(i) = dnew(i,k) - a(i)*cNew(0,k) + (1-alphac)*b(i)*cNew(i,k)
        
        i=imax
        rhs(i) = dnew(i,k) - c(i)*cNew(i1,k) + (1-alphac)*b(i)*cNew(i,k)

        call matrixIdir(imax,a,b,c,rhs)
   
        do i=1,imax
          !resC = resC + (cnew(i,k) - rhs(i))**2.0
          resC = resC + ((cnew(i,k) - rhs(i))/(cnew(i,k)+1.0e-20))**2.0
          cnew(i,k) = max(rhs(i), 0.0)
        enddo    
      enddo
    endif
  !---------------------------------------   ISOFLUX
  else
    !pipe/bl
    if (centerBC.eq.1) then
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
    !channel
    else
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
  use mod_math
  use mod_mesh
  use mod_common
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

  !channel
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
   
      i = imax-1; cu(i)   = 0.0
      i = 1;      au(i+1) = 0.0           ! BC wall
   
      do i=1,imax-1
        rhsu(i) = dt*dnew(i,k) + Unew(i,k)*(Rnew(i+1,k)+Rnew(i,k))*0.5
      enddo
   
      call matrixIdir(imax-1,au,bu,cu,rhsu)
      do i=1,imax-1
        dUdt(i,k)=rhsu(i)
      enddo
    enddo
  !bl/pipe
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

      i = imax-1; cu(i) = 0.0
      i = 1;      bu(i) = bu(i) + au(i)    ! BC at center: Neumann: cancel coeff a
      

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

    i=imax; b(i) = b(i) - c(i)              ! BC at wall: zero vel: subtract one c
    i=1;    b(i) = b(i) + centerBC*a(i)     ! BC at wall: zero vel: subtract one a

    do i=1,imax
      rhs(i) = dt*dnew(i,k) + Wnew(i,k)*(Rnew(i,k)+Rnew(i,k+1))*0.5
    enddo

    call matrixIdir(imax,a,b,c,rhs)
    do i=1,imax
      dWdt(i,k)=rhs(i)
    enddo
  enddo

end




!>*************************************************************************************
!!           cmpinf(Bulk,Stress)
!!
!!*************************************************************************************
subroutine cmpinf(Bulk,Stress)
  use mod_param
  use mod_mesh
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








