!!    ****************************************************
!!    Main file of the code
!!    This code simulates Super Critical Fluid in a
!!    heated pipe/channel/bl with constant heat flux/isothermal wall
!!    ****************************************************

use mod_param
use mod_math
use mod_common
use mod_eos
use mod_tm
use sa_tm
use sst_tm
use mk_tm
use vf_tm
use mod_tdm
use vp_tdm
use mod_mesh

implicit none
include 'mpif.h'

integer ::  rank,ierr,istart,noutput
real(8) ::  bulk,stress,stime,time1
real(8) ::  resC,resK,resE,resV2,resOm,resSA   
real(8) ::  start, finish
real(8) :: resU, resW


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

!initialize numerical
call init_transpose

!initialize EOS
if (EOSmode.eq.0) allocate(eos_model,    source=   IG_EOSModel(Re,Pr))
if (EOSmode.eq.1) allocate(eos_model,    source=Table_EOSModel(Re,Pr,2000, 'tables/co2h_table.dat'))
if (EOSmode.eq.2) allocate(eos_model,    source=Table_EOSModel(Re,Pr,2499, 'tables/ph2_table.dat' ))
call eos_model%init() 

!initialize turbulent viscosity model
if (turbmod.eq.0) allocate(turb_model,source=Laminar_TurbModel(i1, k1, imax, kmax,'lam'))
if (turbmod.eq.1) allocate(turb_model,source=     SA_TurbModel(i1, k1, imax, kmax,'SA' ))
if (turbmod.eq.2) allocate(turb_model,source=init_MK_TurbModel(i1, k1, imax, kmax,'MK' ))
if (turbmod.eq.3) allocate(turb_model,source=init_VF_TurbModel(i1, k1, imax, kmax,'VF' ))
if (turbmod.eq.4) allocate(turb_model,source=    SST_TurbModel(i1, k1, imax, kmax,'SST'))
call turb_model%init()


!initialize grid including bc
if (systemsolve .eq. 1) allocate(mesh, source=       Pipe_Mesh(i1,k1,imax,kmax))
if (systemsolve .eq. 2) allocate(mesh, source=    Channel_Mesh(i1,k1,imax,kmax))
if (systemsolve .eq. 3) allocate(mesh, source= SymChannel_Mesh(i1,k1,imax,kmax))
if (systemsolve .eq. 4) allocate(mesh, source=     BLayer_Mesh(i1,k1,imax,kmax))
call mesh%init(LoD, K_start_heat, x_start_heat, rank,px)
call mesh%discretize_streamwise2( LoD,rank, px)

!initialize turbulent diffusivity model
if (turbdiffmod.eq.0) allocate(turbdiff_model,source=        CPrt_TurbDiffModel(i1, k1, imax, kmax,'Pr', Pr))
if (turbdiffmod.eq.1) allocate(turbdiff_model,source=  Irrenfried_TurbDiffModel(i1, k1, imax, kmax,'IF'    ))
if (turbdiffmod.eq.2) allocate(turbdiff_model,source=        Tang_TurbDiffModel(i1, k1, imax, kmax,'Tang'  ))
if (turbdiffmod.eq.3) allocate(turbdiff_model,source=KaysCrawford_TurbDiffModel(i1, k1, imax, kmax,'KC'    ))
if (turbdiffmod.eq.4) allocate(turbdiff_model,source=        Kays_TurbDiffModel(i1, k1, imax, kmax,'Kays'  ))
if (turbdiffmod.eq.5) allocate(turbdiff_model,source=    init_Bae_TurbDiffModel(i1, k1, imax, kmax,'Bae', 70.,1.,mesh%walldist*Re))
call turbdiff_model%init()

! if (rank .eq. 0) then
! do i=0,i1
!   write(*,*) mesh%dru(i), dru(i)
! enddo
! endif

! call mpi_finalize(ierr)
! stop

dt = dtmax
istart = 1

!initialize solution 
call initialize_solution(rank,wnew,unew,cnew,ekmt,win,ekmtin,i1,k1,mesh%y_fa,mesh%y_cv,mesh%dpdz,Re,systemsolve,select_init)
! write(*,*) win
call bound_v(Unew,Wnew,Win,rank,istep)

  

call calc_prop(cnew,rnew,ekm,ekmi,ekmk,ekh,ekhi,ekhk,cp,cpi,cpk,temp,beta)
call bound_c(cnew, Tw, Qwall,rank)
call turb_model%set_bc(ekm,rnew,periodic,rank,px)
call calc_prop(cnew,rnew,ekm,ekmi,ekmk,ekh,ekhi,ekhk,cp,cpi,cpk,temp,beta) ! necessary to call it twice

rold = rnew
call calc_mu_eff(Unew,Wnew,rnew,ekm,ekmi,ekme,ekmt,rank) 
call bound_v(Unew,Wnew,Win,rank,istep)
call chkdt(rank,istep)
call cpu_time(start)

! call output2d_upd2(rank,istep)
! call mpi_finalize(ierr)
! stop

!***************************!
!        MAIN LOOP          !
!***************************!


do istep=istart,nstep
  call calc_mu_eff(Unew,Wnew,rnew,ekm,ekmi,ekme,ekmt,rank) 
  call turbdiff_model%set_alphat(unew,wnew,rnew,temp,ekm,ekmi,ekh,ekmt,alphat)
  call advanceC(resC,Unew,Wnew,Rnew,rank)


  ! if   (mod(istep,100).eq.0) then
  !   call output2d_upd(rank,istep)
  ! endif
  
  call turb_model%advance_turb(uNew,wNew,rnew,ekm,ekmi,ekmk,ekmt,beta,temp,           &
                              mesh%Ru,mesh%Rp,mesh%dru,mesh%drp,mesh%dz,mesh%walldist,alphak,alphae,alphav2,        &
                              modifDiffTerm,rank,mesh%centerBC,periodic,resSA,resK, resV2)
  call bound_c(Cnew, Tw, Qwall,rank)
  call turb_model%set_bc(ekm,rnew,periodic,rank,px)
  call calc_prop(cnew,rnew,ekm,ekmi,ekmk,ekh,ekhi,ekhk,cp,cpi,cpk,temp,beta);
  call advance(rank)


  call bound_m(dUdt,dWdt,wnew,rnew,Win,rank, istep)
  ! call output2d_upd2(rank,istep)
  ! call mpi_finalize(ierr)
  ! stop

  call fillps(rank)
  call solvepois_cr(p,0,rank,mesh%centerBC)
  ! call solvepois(p,Ru,Rp,dRu,dRp,dz,rank,centerBC)
  call correc(rank,1)
  call bound_v(Unew,Wnew,Win,rank, istep)
  if   (mod(istep,10) .eq. 0) call chkdiv(rank)

  call cmpinf(bulk,stress)
  call chkdt(rank,istep)
  if  ((mod(istep,500).eq.0).and.(periodic .eq.1)) call inflow_output_upd(rank);

  if   (mod(istep,500).eq.0) then
    call output2d_upd2(rank,istep) 
    if (systemSolve .eq. 4.) call write_output_bl(rank,istep) !extra output for the bl
  endif
  
  !write the screen output
  noutput = 100
  ! if (mod(istep,noutput) .eq. 0) then 
  !   call  calc_residual(unew, uold, wnew, wold, resU, resW)
  !   if ((resU .le. 1e-10) .and. (resW .le. 1e-10)) exit
  ! endif
      
  if (rank.eq.0) then
    if (istep.eq.istart .or. mod(istep,noutput*20).eq.0) then
      write(6,'(A7,9A14)') 'istep'    ,'dt'      ,'bulk'   ,'stress' ,'cResid', &
                           'kineResid','epsResid','v2Resid','omResid','nuSAresid'
    endif
    if (istep.eq.istart .or. mod(istep,noutput).eq.0) then
      write(6,'(i7,9e14.5)') istep,dt,bulk,stress,resC,resK,resE,resV2,resU,resW
    endif
  end if
         
enddo
call cpu_time(finish)
print '("Time = ",f6.3," seconds.")',finish-start
if (systemSolve .eq. 4.) call write_output_bl(rank,istep) !extra output for the bl
call output2d_upd2(rank,istep)
call mpi_finalize(ierr)
stop
end

!***********************************************************************************************************************************
!***********************************************************************************************************************************
!***********************************************************************************************************************************
!***********************************************************************************************************************************
subroutine get_max_2D(vector, sizei, sizek, max_global)
  implicit none
   include 'mpif.h'
  integer ierr,rank
  integer, intent(IN) :: sizei,sizek
  real(8), dimension(0:sizei,0:sizek), intent(IN) :: vector
  real(8), intent(OUT) ::  max_global 
  real(8) :: max_local
  integer :: i,k
  max_local = 1e-100
  do i=0,sizei
    do k=0,sizek
      if (vector(i,k) > max_local) then
        max_local = vector(i,k)
      endif
    enddo
  enddo
  call mpi_allreduce(max_local,max_global,1,mpi_real8,mpi_max,mpi_comm_world,ierr)
end subroutine get_max_2D

subroutine calc_residual(unew, uold, wnew, wold, resU, resW)
  use mod_param, only: k1, i1
  real(8), dimension(0:i1,0:k1), intent(IN) :: unew,uold,wNew,wold
  real(8), intent(OUT) :: resW, resU
  call get_max_2D(unew-uold, i1,k1,resU)
  call get_max_2D(wnew-wold, i1,k1,resW)
end subroutine  calc_residual
!!********************************************************************
!!     Calculates the thermodynamic properties
!!********************************************************************
subroutine calc_prop(enth,rho,mu,mui,muk,lam,lami,lamk,cp,cpi,cpk,tp,be)
  use mod_param
  use mod_eos
  implicit none
  real(8), dimension(0:i1, 0:k1), intent(OUT):: enth,rho,mu,mui,muk, &
                                                lam,lami,lamk,cp,cpk,cpi,tp,be
  real(8) :: enthface
  !centers
  do k=0,k1
    do i=0,i1
        call eos_model%set_w_enth(enth(i,k),"D", rho(i,k))
        call eos_model%set_w_enth(enth(i,k),"V", mu(i,k))
        call eos_model%set_w_enth(enth(i,k),"C", Cp(i,k))
        call eos_model%set_w_enth(enth(i,k),"L", lam(i,k))
        call eos_model%set_w_enth(enth(i,k),"T", tp(i,k))
        call eos_model%set_w_enth(enth(i,k),"B", be(i,k)) 
    enddo
  enddo
  !faces
  do k=0,kmax
    do i=0,imax
      enthface = 0.5*(enth(i,k)+enth(i+1,k))
      call eos_model%set_w_enth(enthface,"C", cpi(i,k))
      call eos_model%set_w_enth(enthface,"L", lami(i,k))
      call eos_model%set_w_enth(enthface,"V", mui(i,k)) 
      enthface = 0.5*(enth(i,k)+enth(i,k+1))
      call eos_model%set_w_enth(enthface,"C", cpk(i,k))
      call eos_model%set_w_enth(enthface,"L", lamk(i,k))
      call eos_model%set_w_enth(enthface,"V", muk(i,k)) 
    enddo
  enddo
  return
end subroutine calc_prop

!!******************************************************************************************
!!      routine to estimate the effective viscosity
!!******************************************************************************************
subroutine calc_mu_eff(utmp,wtmp,rho,mu,mui,mue,mut,rank)
  use mod_param, only : k1,i1,kmax,imax,periodic,px
  use mod_tm,    only : turb_model
  use mod_mesh,  only : mesh

  implicit none
  real(8), dimension(0:i1,0:k1), intent(IN) :: utmp,wtmp,rho,mu,mui   
  integer,                       intent(IN) :: rank
  real(8), dimension(0:i1,0:k1), intent(OUT):: mue,mut
  real(8), dimension(0:k1) :: tauw(0:k1)

  call turb_model%set_mut(utmp,wtmp,rho,mu,mui,mesh%walldist,mesh%rp,mesh%drp,mesh%dru,mesh%dz,mut)
  call turb_model%set_mut_bc(mut,periodic,px,rank)
  mue = mu + mut  
end subroutine calc_mu_eff

!!*************************************************************************************
!!  Apply the boundary conditions for the energy equation
!!*************************************************************************************
subroutine bound_c(c, Twall, Qwall,rank)
  use mod_param,   only : k1,i1,kmax,imax,i,k,isothermalBC,periodic,px
  use mod_eos,     only : eos_model
  use mod_mesh, only : mesh
  implicit none
  include 'mpif.h'
  real(8),                       intent(IN) :: Twall, Qwall
  integer,                       intent(IN) :: rank
  real(8), dimension(0:i1,0:k1), intent(OUT):: c
  real(8)                  :: enth_wall  
  real(8), dimension(0:i1) :: tmp,drp
  real(8), dimension(0:k1) :: bot_bcvalue1, top_bcvalue1

  top_bcvalue1 = mesh%top_bcvalue1
  bot_bcvalue1 = mesh%bot_bcvalue1
  drp = mesh%drp

  !isothermal
  if (isothermalBC.eq.1) then
    call eos_model%set_w_temp(Twall, "H", enth_wall)
    c(0,:) = (1-bot_bcvalue1(:))*(2.0*enth_wall - c(1,:))   +bot_bcvalue1(k)*c(1,:)    !pipe/bl
    c(i1,:)= (1-top_bcvalue1(:))*(2.0*enth_wall - c(imax,:))+top_bcvalue1(k)*c(imax,:) !pipe/bl

  !isoflux
  else
    do k=0,k1
      if (top_bcvalue1(k) .eq. 1) then
        c(i1,k) = c(imax,k)                                                !symmetry
      else
        call eos_model%set_enth_w_qwall(qwall,c(imax,k),drp(imax),c(i1,k)) !heatflux
      endif
      if (bot_bcvalue1(k) .eq. 1) then
        c(0,k) = c(1,k)                                                    !symmetry
      else
        call eos_model%set_enth_w_qwall(qwall,c(1,k),   drp(0),   c(0,k))  !heatflux
      endif
    enddo
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
subroutine bound_v(u,w,win,rank,step)
  use mod_param,   only : i1,k1,imax,kmax,periodic,px
  use mod_mesh, only : mesh
  implicit none  
  include 'mpif.h'
  
  integer,                       intent(IN) :: rank, step
  real(8), dimension(0:i1),      intent(IN) :: win
  real(8), dimension(0:i1,0:k1), intent(OUT):: u, w
  real(8), dimension(0:i1)                  :: tmp
  real(8) :: x, vfs, vfsm
  real(8), dimension(0:k1) :: dis

  u(0,:)    =  (1-mesh%ubot_bcvalue(:))*u(1,:) !wall and symmetry !pipe&chan: u=0, bl: du/dy=0
  u(imax,:) =   0.0                            !wall and symmetry
  u(i1,:)   = - u(imax-1,:)

  w(0,:)    = mesh%bot_bcnovalue(:)*w(1,:)    !wall (bot_bcnovalue=-1) or symmetry (bot_bcnovalue=1)
  w(i1,:)   = mesh%top_bcnovalue(:)*w(imax,:) !wall or symmetry
  
    
  call shiftf(u,tmp,rank);     u(:,0)  = tmp(:);
  call shiftb(u,tmp,rank);     u(:,k1) = tmp(:);
  call shiftf(w,tmp,rank);     w(:,0)  = tmp(:);
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
!!************************  *************************************************************
subroutine bound_m(Ubound,Wbound,W_out,Rbound,W_in,rank, step)
  use mod_param,   only : k,i,kmax,imax,k1,i1,px
  use mod_mesh, only : mesh
  use mod_common,  only : dt,rold,rnew
  implicit none
  include 'mpif.h'
  integer,                       intent(IN) :: rank, step
  real(8), dimension(0:i1,0:k1), intent(IN) :: rbound
  real(8), dimension(0:i1),      intent(IN) :: W_in  
  real(8), dimension(0:i1,0:k1), intent(OUT):: ubound, wbound, w_out
  real(8), dimension(0:i1) :: tmp,rp, dru
  real(8), dimension(0:k1) :: dzw
  real(8), dimension(1:imax) :: wr
  integer :: ierr
  real(8) :: Ub,flux,flux_tot,deltaW,wfunc
  
  dzw = mesh%dzw
  dru = mesh%dRu
  rp  = mesh%rp
  
  call bound_v(ubound,wbound,W_in,rank, step)
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
      Wbound(i,kmax) = 2.0*Wbound(i,kmax-1) - Wbound(i,kmax-2) !NOTE: CHANGE WITH SIMONE
      Wbound(i,kmax) = Wbound(i,kmax)*0.5*(Rbound(i,kmax)+Rbound(i,k1))
    enddo
  endif
  !compute drho/dt*dvol
  do k=1,kmax
    do i=1,imax
      flux = flux - (rnew(i,k)-rold(i,k))/dt*Rp(i)*dru(i)*dzw(k)      
    enddo
  enddo

  !     compute mf in
  if (rank.eq.0)then
    do i=1,imax
      flux = flux + Wbound(i,0)*dru(i)*rp(i)
    enddo
  endif

  !compute mf out to the bot (used for the bl)
  do k=1,kmax
     flux = flux + Ubound(0,k)*dzw(k)
  enddo

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
  endif

end subroutine bound_m

!!*************************************************************************************
!!           initialize the solution
!!*************************************************************************************

subroutine initialize_solution(rank, w, u,c, mut, win, mutin, i1,k1, y_fa, y_cv, dpdz, Re, systemsolve, select_init)
  use mod_tm
  use mod_param, only : imax_old, kelem_old,px
  implicit none
  include "mpif.h"
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
    if (systemsolve .eq. 3) case = "symchan"
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
  else if (select_init .eq. 2) then
    call interpolate_solution(imax_old+1, (kelem_old/px)+1, rank, px)
  else
    if (rank.eq.0) write(*,*) 'Initializing flow from scratch'
    gridSize = y_fa(imax)
    do i=0,i1!imax         
      if (systemsolve.eq.1) w(i,:)  = Re/6*3/2.*(1-(y_cv(i)/0.5)**2); 
      if (systemsolve.eq.2) w(i,:)  = Re*dpdz*y_cv(i)*0.5*(gridSize-y_cv(i))              !channel
      if (systemsolve.eq.3) w(i,:)  = Re*dpdz*0.5*((gridSize*gridSize)-(y_cv(i)*y_cv(i))) !bl
      if (systemsolve.eq.4) then
         w(i,:)  = 1.;u=0;  win=1.; mutin=0!bl
      endif
    enddo
    ! win=w(:,0);          !pipe
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
  use mod_param,   only : k1,i1,imax,kmax,alphac,k,i,periodic
  use mod_math,    only : matrixIdir
  use mod_mesh, only : mesh
  use mod_common,  only : cnew,ekhk,ekh,alphat,ekhi
  implicit none
  real(8), dimension(0:i1,0:k1), intent(IN) :: Utmp, Wtmp, Rtmp
  integer,                       intent(IN) :: rank
  real(8),                       intent(OUT):: resC
  real(8), dimension(0:i1,0:k1) :: dnew,dimpl
  real(8), dimension(imax)      :: a,b,c,rhs
  real(8)                       :: sigmat,Q,dz
  integer  :: ierr
  real(8), dimension(0:i1) :: rp, dru, ru, drp
  real(8), dimension(0:k1) :: dzw, top_bcvalue1, bot_bcvalue1

  dzw = mesh%dzw
  dru = mesh%dRu
  drp = mesh%drp
  rp  = mesh%rp
  ru  = mesh%ru
  top_bcvalue1 = mesh%top_bcvalue1
  bot_bcvalue1 = mesh%bot_bcvalue1

  resC   = 0.0; dnew   = 0.0; dimpl = 0.0;

  call advecc(dnew,dimpl,cnew,Utmp,Wtmp,Ru,Rp,dru,dz,i1,k1,rank,periodic,.true.)
  call diffc(dnew,cnew,ekh,ekhi,ekhk,alphat,1.,Rtmp,Ru,Rp,dru,dz,rank,0)

  do k=1,kmax
    do i=1,imax
      a(i) = -Ru(i-1)*(ekhi(i-1,k)+0.5*(alphat(i,k)+alphat(i-1,k)))/(dRp(i-1)*Rp(i)*dru(i))/Rtmp(i,k)
      c(i) = -Ru(i  )*(ekhi(i  ,k)+0.5*(alphat(i,k)+alphat(i+1,k)))/(dRp(i  )*Rp(i)*dru(i))/Rtmp(i,k)
      b(i) = (-a(i)-c(i) + dimpl(i,k) )        
      rhs(i) = dnew(i,k) + ((1-alphac)/alphac)*b(i)*cnew(i,k)  
    enddo

    i=1
    b(i)=b(i)+bot_bcvalue1(k)*a(i) !symmetry or nothing
    rhs(i) = dnew(i,k) - (1-bot_bcvalue1(k))*a(i)*cNew(i-1,k) + ((1-alphac)/alphac)*b(i)*cNew(i,k) !nothing or value

    i=imax
    b(i)=b(i)+top_bcvalue1(k)*c(i) !symmetry or nothing
    rhs(i) = dnew(i,k) - (1-top_bcvalue1(k))*c(i)*cNew(i+1,k) + ((1-alphac)/alphac)*b(i)*cNew(i,k) !nothing or value

    call matrixIdir(imax,a,b/alphac,c,rhs)

    do i=1,imax
      resC = resC + ((cnew(i,k) - rhs(i))/(cnew(i,k)+1.0e-20))**2.0
      cnew(i,k) = max(rhs(i), 0.0)
    enddo    
  enddo
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
  use mod_param,   only : k1,i1,kmax,imax,k,i,Fr_1,periodic,qwall
  use mod_math,    only : matrixIdir
  use mod_mesh, only : mesh
  use mod_common,  only : rnew,ekme,unew,wnew,dUdt,dWdt,ekm,dt,peclet
  implicit none
      
  integer rank
  real*8, dimension(imax)   :: a, b, c, rhs
  real*8, dimension(imax-1) :: au, bu, cu, rhsu
  real*8 dnew(0:i1,0:k1)
  real*8 dif,alpha,rhoa,rhob,rhoc
  real*8 :: x, dpdz, dz
  real(8), dimension(0:i1) :: rp, dru, ru, drp
  real(8), dimension(0:k1) :: dzw, top_bcnovalue, bot_bcnovalue, ubot_bcvalue

  dzw = mesh%dzw
  dru = mesh%dRu
  drp = mesh%drp
  rp  = mesh%rp
  ru  = mesh%ru
  dpdz= mesh%dpdz
  top_bcnovalue = mesh%top_bcnovalue
  bot_bcnovalue = mesh%bot_bcnovalue
  ubot_bcvalue  = mesh%ubot_bcvalue

  dif=0.0

  !>********************************************************************
  !!     CALCULATE advection, diffusion and Force in r-direction
  !!     at the old(n-1) and new(n) timelevels
  !!********************************************************************
  dnew = 0.0
  call advecu(dnew,Unew,Wnew,Rnew,Ru,Rp,dru,drp,dz,i1,k1) ! new
  call diffu (dnew,Unew,Wnew,ekme,Ru,Rp,dru,drp,dz,i1,k1,dif,mesh%numDomain) ! new

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
      rhsu(i) = dt*dnew(i,k) + Unew(i,k)*(Rnew(i+1,k)+Rnew(i,k))*0.5
    enddo
 
    i = imax-1; cu(i)   = 0.0           ! BC wall and symmetry
    i=1;   bu(i) = bu(i) + (1.-ubot_bcvalue(k))*au(i); au(i) = (1.-ubot_bcvalue(k))*au(i) !symmetry with 0 or with the derivative
    !ubot_bcvalue=1: au(i)=0 and bu(i)=bu(i), ubot_bcvalue=0: au(i)=au(i), bu(i)=bu(i)+au(i)
   
    call matrixIdir(imax-1,au,bu,cu,rhsu)
    do i=1,imax-1
      dUdt(i,k)=rhsu(i)
    enddo
  enddo

  !********************************************************************
  !     CALCULATE advection, diffusion and Force in z-direction
  !     at the old(n-1) and new(n) timelevels
  !********************************************************************

  dnew = 0.0
  call advecw(dnew,Unew,Wnew,Rnew,Ru,Rp,dru,dz,ekm,peclet)
  call diffw (dnew,Unew,Wnew,ekme,Ru,Rp,dru,drp,dz,i1,k1,dif,mesh%numDomain)  ! new

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

    i=1;    b(i) = b(i) + bot_bcnovalue(k)*a(i)   !BC wall or symmetry
    i=imax; b(i) = b(i) + top_bcnovalue(k)*c(i)   !BC wall or symmetry

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
  use mod_common
  use mod_mesh, only : mesh
  implicit none   
  include 'mpif.h'
  integer ierr
  real*8 waver(imax),waver2(imax)
  real(8), dimension(0:i1)   :: rp, dru 
  real(8), dimension(1:imax) :: wallDist
  real*8  Bulk,Stress
  
  rp       = mesh%rp
  dru      = mesh%dru
  walldist = mesh%wallDist

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








