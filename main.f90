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
use abe_tm
use vf_tm
use mod_tdm
use vp_tdm
use mod_mesh
use dwx_tdm
use nk_tdm

implicit none
include 'mpif.h'

integer ::  rank,ierr,istart,noutput
real(8) ::  bulk,stress,stime,time1,hbulk,tbulk,massflow, Re_bulk,vbulk
real(8) ::  resC,resTV1,resTV2,resTV3,resTD1,resTD2,totResC
real(8) ::  start, finish
real(8) :: resU, resW
real(8) :: fr_1_goal


resC=0;resTV1=0;resTV2=0;resTV3=0;resTD1=0;resTD2=0;totResC=0!;resKt=0;resEpst=0

!read parameters
call read_parameters()
!Fr_1_goal = Fr_1
!Fr_1 = 0.
! write(*,*) trim(test)
! call mpi_finalize(ierr)
! stop

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

!initialize grid including bc
if (systemsolve .eq. 1) allocate(mesh, source=       init_PipeMesh("pipe"))
if (systemsolve .eq. 2) allocate(mesh, source=    init_ChannelMesh("channel"))
if (systemsolve .eq. 3) allocate(mesh, source= init_SymChannelMesh("symchan"))
if (systemsolve .eq. 4) allocate(mesh, source=     init_BLayerMesh("bl"))
call mesh%init(LoD, K_start_heat, x_start_heat, rank,px)
call mesh%discretize_streamwise2( LoD,rank, px)

!initialize EOS
if (EOSmode.eq.0) allocate(eos_model,    source=init_ConstProp_EOSModel(Re,Pr))
if (EOSmode.eq.1) allocate(eos_model,    source=       init_IG_EOSModel(Re,Pr))
if (EOSmode.eq.2) allocate(eos_model,    source=    init_Table_EOSModel(Re,Pr,nTab, table_loc,fluid_name))
! if (EOSmode.eq.3) allocate(eos_model,    source=    init_Table_EOSModel(Re,Pr,2499, 'tables/ph2_table.dat' ,'h2-'))
! if (EOSmode.eq.4) allocate(eos_model,    source=    init_Table_EOSModel(Re,Pr,2000, "tables/h2o_table.dat",'h20' ))

call eos_model%init() 

!initialize turbulent viscosity model
if (turbmod.eq.0) allocate(turb_model,source= Laminar_TurbModel('lam'))
if (turbmod.eq.1) allocate(turb_model,source=      SA_TurbModel('SA' ))
if (turbmod.eq.2) allocate(turb_model,source= init_MK_TurbModel('MK' ))
if (turbmod.eq.3) allocate(turb_model,source= init_VF_TurbModel('VF' ))
if (turbmod.eq.4) allocate(turb_model,source=init_SST_TurbModel('SST'))
if (turbmod.eq.5) allocate(turb_model,source=init_Abe_TurbModel('Abe'))
call turb_model%init()

!initialize turbulent diffusivity model
if (turbdiffmod.eq.0) allocate(turbdiff_model,source=   init_CPrt_TurbDiffModel('cPr', PrT))
if (turbdiffmod.eq.1) allocate(turbdiff_model,source=  Irrenfried_TurbDiffModel('IF'    ))
if (turbdiffmod.eq.2) allocate(turbdiff_model,source=        Tang_TurbDiffModel('Tang'  ))
if (turbdiffmod.eq.3) allocate(turbdiff_model,source=KaysCrawford_TurbDiffModel('KC'    ))
if (turbdiffmod.eq.4) allocate(turbdiff_model,source=        Kays_TurbDiffModel('Kays'  ))
if (turbdiffmod.eq.5) allocate(turbdiff_model,source=    init_Bae_TurbDiffModel('Bae', 70.,20.))
if (turbdiffmod.eq.6) allocate(turbdiff_model,source=    init_DWX_TurbDiffModel('DWX'))
if (turbdiffmod.eq.7) allocate(turbdiff_model,source=     init_NK_TurbDiffModel('NK'))
if (((turbdiffmod.eq.6).or.(turbdiffmod.eq.7)).and.((turbmod.eq.1).or.(turbmod.eq.4))) then
  if (rank .eq. 0)  write(*,*) "Combination of eddy viscosity model and turbulent diffusivity model not valid", turbmod, turbdiffmod
  call mpi_finalize(ierr)
  stop
endif
call turbdiff_model%init()

dt = dtmax
istart = 1

!initialize solution 
call initialize_solution(rank,wnew,unew,cnew,ekmt,alphat,win,Re,systemsolve,select_init)

!#call debug(rank)
!#call mpi_finalize(ierr)
!#stop

call bound_v(Unew,Wnew,Win,Wnew,rank,istep)
call calc_prop(cnew,rnew,ekm,ekmi,ekmk,ekh,ekhi,ekhk,cp,cpi,cpk,temp,beta)

!call set_qwall(Fr_1, Fr_1_goal, -1.)
call bound_c(cnew, Tw_top,Tw_bot, Qwall,rank)
call turb_model%set_bc(ekm,rnew,periodic,rank,px)
call calc_prop(cnew,rnew,ekm,ekmi,ekmk,ekh,ekhi,ekhk,cp,cpi,cpk,temp,beta) ! necessary to call it twice

rold = rnew
call calc_mu_eff(Unew,Wnew,rnew,ekm,ekmi,ekme,ekmt,rank) 
call calc_ekh_eff(Unew,Wnew,rnew,temp,ekm,ekmi,ekh,ekhi,ekhk,ekmt,ekhe,alphat,rank)

! call turbdiff_model%set_alphat(unew,wnew,rnew,temp,ekm,ekmi,ekh,ekmt,alphat)
call bound_v(Unew,Wnew,Win,Wnew,rank,istep)
call chkdt(rank,istep)
call cpu_time(start)


  
!***************************!
!        MAIN LOOP          !
!***************************!


do istep=istart,nstep
  
  !calc the 
  call calc_mu_eff(Unew,Wnew,rnew,ekm,ekmi,ekme,ekmt,rank) 

  call calc_ekh_eff(Unew,Wnew,rnew,temp,ekm,ekmi,ekh,ekhi,ekhk,ekmt,ekhe,alphat,rank)
                
  !scalar equations
  call advanceC(resC,Unew,Wnew,Rnew,rank)
  call turb_model%advance_turb(uNew,wNew,rnew,ekm,ekmi,ekmk,ekmt,beta,temp,           &
                              alphak,alphae,alphav2,        &
                              modifDiffTerm,rank,periodic,resTV1,resTV2,resTV3)


  call turbdiff_model%advance_turbdiff(unew,wnew,cnew,temp,rnew,ekm,ekh,ekhi,ekhk,alphat, &
                                      alphakt,alphaet,rank,periodic,resTD1,resTD2)
  !apply bc
  call bound_c(Cnew, Tw_top, Tw_bot, Qwall,rank)
  call turb_model%set_bc(ekm,rnew,periodic,rank,px)
  call turbdiff_model%set_bc(ekh,rnew,periodic,rank,px)


  call calc_prop(cnew,rnew,ekm,ekmi,ekmk,ekh,ekhi,ekhk,cp,cpi,cpk,temp,beta);

  if (bulkmod .eq. 1) call set_dpdz_wbulk(wnew,rank)
  call advance(rank)

  call bound_m(dUdt,dWdt,wnew,rnew,Win,rank, istep)
  
  call fillps(rank)
  ! call solvepois_cr(p,0,rank,mesh%centerBC)
  call solvepois(p,rank,mesh%centerBC)
  call correc(rank,1)
  call bound_v(Unew,Wnew,Win,Wnew,rank, istep)
 ! if (mod(istep,10000).eq.0) call set_qwall(Fr_1,Fr_1_goal,-1.)
  if   (mod(istep,10) .eq. 0) call chkdiv(rank)

  call cmpinf(bulk,stress)
  call chkdt(rank,istep)
  if  ((mod(istep,100).eq.0).and.(periodic .eq.1)) call inflow_output_upd(rank);

  if   (mod(istep,100).eq.0) then
    call output2d_upd2(rank,istep) 
    !call MPI_ALLREDUCE (resC,totResC, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    !if ((periodic .eq.2) .and.(totResC .lt. convCrit )) exit
    if (systemSolve .eq. 4.) call write_output_bl(rank,istep) !extra output for the bl
  endif
  
  !write the screen output
  noutput = 100

!   if (mod(istep,noutput).eq.0) then
!     call debug(rank)
!   endif

  if (rank.eq.0) then
    if (istep.eq.istart .or. mod(istep,noutput*20).eq.0) then
      write(6,'(A7,11A13)') 'istep'  ,'dt' ,'bulk','tbulk','Re_b'   ,'stress' ,'cResid', &
                           'resTV1','resTV2','resTV3', &
                           'resTD1', 'resTD2'
    endif
    if (istep.eq.istart .or. mod(istep,noutput).eq.0) then
      call calc_avg_quantities(wnew, rnew, cnew,massflow, hbulk, tbulk,vbulk, Re_bulk,kmax)
      write(6,'(i7,11e13.4)') istep,dt,bulk,vbulk,Re_bulk,stress,resC,resTV1,resTV2,resTV3,resTD1,resTD2
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

subroutine set_qwall(qwall,qwall_goal,dq)
  implicit none 
  real(8),  intent(IN) :: qwall_goal
  real(8),  intent(OUT) :: qwall
  real(8) :: dq
  qwall = min(qwall+dq,qwall_goal)
end subroutine
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
  use mod_param, only : k1,i1,kmax,imax,k,i
  use mod_eos, only : eos_model
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
subroutine calc_ekh_eff(utmp,wtmp,rho,temp,mu,mui,lam_cp,lam_cpi,lam_cpk,mut,ekhe,alphat,rank)
  use mod_param, only : k1,i1,kmax,imax,periodic,px,i,k
  use mod_tdm,    only : turbdiff_model
  use mod_common, only : ekhek, ekhei
  implicit none
  real(8), dimension(0:i1,0:k1), intent(IN) :: utmp,wtmp,rho,temp,mu,mui,lam_cp,lam_cpi,lam_cpk,mut
  integer,                       intent(IN) :: rank
  real(8), dimension(0:i1,0:k1), intent(OUT):: ekhe,alphat
  real(8), dimension(0:k1) :: tauw(0:k1)
  !calc alphat
  
  call turbdiff_model%set_alphat(utmp,wtmp,rho,temp,mu,mui,lam_cp,mut,alphat)
  call turbdiff_model%set_alphat_bc(alphat,periodic,px,rank)

  do k=0,k1
    do i=0,i1
      ekhe(i,k) = lam_cp(i,k) + alphat(i,k)*rho(i,k)
    enddo
  enddo
  do k=0,kmax
    do i=0,imax
      ekhek(i,k) = lam_cpk(i,k) + 0.25*(alphat(i,k)+alphat(i,k+1))*(rho(i,k)+rho(i,k+1))
      ekhei(i,k) = lam_cpi(i,k) + 0.25*(alphat(i,k)+alphat(i+1,k))*(rho(i,k)+rho(i+1,k))
    enddo
  enddo
end subroutine calc_ekh_eff

!!******************************************************************************************
!!      routine to estimate the effective viscosity
!!******************************************************************************************
subroutine calc_mu_eff(utmp,wtmp,rho,mu,mui,mue,mut,rank)
  use mod_param, only : k1,i1,kmax,imax,periodic,px
  use mod_tm,    only : turb_model
  implicit none
  real(8), dimension(0:i1,0:k1), intent(IN) :: utmp,wtmp,rho,mu,mui   
  integer,                       intent(IN) :: rank
  real(8), dimension(0:i1,0:k1), intent(OUT):: mue,mut
  real(8), dimension(0:k1) :: tauw(0:k1)

  call turb_model%set_mut(utmp,wtmp,rho,mu,mui,mut)
  call turb_model%set_mut_bc(mut,periodic,px,rank)
  mue = mu + mut
end subroutine calc_mu_eff


!!*************************************************************************************
!!  Apply the boundary conditions for the energy equation
!!*************************************************************************************
subroutine bound_c(c, Twall_top, Twall_bot, Qwall,rank)
  use mod_param,only : k1,i1,kmax,imax,i,k,isothermalBC,periodic,px
  use mod_eos,  only : eos_model
  use mod_mesh, only : top_bcvalue1,bot_bcvalue1,drp
  use mod_common, only: wnew
  implicit none
  include 'mpif.h'
  real(8),                       intent(IN) ::  Twall_top, Twall_bot, Qwall
  integer,                       intent(IN) :: rank
  real(8), dimension(0:i1,0:k1), intent(OUT):: c
  real(8)                  :: enth_wall  
  real(8), dimension(0:i1) :: tmp
  
  !isothermal
  if (isothermalBC.eq.1) then
    call eos_model%set_w_temp(Twall_bot, "H", enth_wall)
    c(0,:) = (1-bot_bcvalue1(:))*(2.0*enth_wall - c(1,:))   +bot_bcvalue1(k)*c(1,:)    !pipe/bl
    call eos_model%set_w_temp(Twall_top, "H", enth_wall)
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
 !extrapolation outlet
 if (rank.eq.px-1) then
   c(:,k1) = 2.0*   c(:,kmax) -    c(:,kmax-1)
 endif 
  ! call bound_conv(c,wnew,rank,0) 

end subroutine bound_c

subroutine set_dpdz_wbulk(w,rank)

  use mod_param, only : i1,kmax,px,k1,i,imax, systemsolve
  use mod_mesh, only : dRu, mesh,rp
  implicit none
  include "mpif.h"
  real(8), dimension(0:i1,0:k1), intent(IN) :: w
  real(8) bulk,bulk_tot, pi
  integer :: ierror,rank

  bulk = 0.
  ! if (rank .eq. 0) then
  pi = 4.*atan(1.)

  if (systemsolve .eq. 2) then
      do i=1,imax
        bulk = bulk + dru(i)*w(i,kmax)
      enddo
    ! endif
      
    call mpi_allreduce(bulk,  bulk_tot,1, MPI_REAL8,MPI_SUM,MPI_COMM_WORLD, ierror)
    bulk = (bulk_tot/2.)/px

  else
    if (systemsolve .eq. 1) then
      ! 2*pi;

      do i=1,imax
        bulk = bulk + 2*pi*rp(i)*dru(i)*w(i,kmax)
      enddo
      call mpi_allreduce(bulk,  bulk_tot,1, MPI_REAL8,MPI_SUM,MPI_COMM_WORLD, ierror)
      bulk = (bulk_tot/(pi/4.0))/px
    !
    endif
  endif

  mesh%dpdz     = (0.99*mesh%dpdz + (1.-bulk));

end subroutine

!!*************************************************************************************
!!  Apply the boundary conditions for the velocity
!!*************************************************************************************
subroutine bound_v(u,w,win,wout,rank,step)
  use mod_param,only : i1,k1,imax,kmax,periodic,px
  use mod_mesh, only : ubot_bcvalue,top_bcnovalue,bot_bcnovalue
  implicit none  
  include 'mpif.h'
  
  integer,                       intent(IN) :: rank, step
  real(8), dimension(0:i1),      intent(IN) :: win
  real(8), dimension(0:i1,0:k1), intent(OUT):: u, w, wout
  real(8), dimension(0:i1)                  :: tmp
  real(8)                                   :: bulk

  u(0,:)    =  (1-ubot_bcvalue(:))*u(1,:) !wall and symmetry !pipe&chan: u=0, bl: du/dy=0
  u(imax,:) =   0.0                            !wall and symmetry
  u(i1,:)   = - u(imax-1,:)

  w(0,:)    = bot_bcnovalue(:)*w(1,:)    !wall (bot_bcnovalue=-1) or symmetry (bot_bcnovalue=1)
  w(i1,:)   = top_bcnovalue(:)*w(imax,:) !wall or symmetry
  
    
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
  !if (rank.eq.px-1)then
  !  !extrapolation
  !  u(:,k1) = 2.*u(:,kmax)-u(:,kmax-1)
  !  w(:,kmax) = 2.*w(:,kmax-1)-w(:,kmax-2)
  !endif
  !convective outlet
  call bound_conv(u,wout,rank,0)
  call bound_conv(w,wout,rank,1)
end subroutine bound_v


subroutine bound_conv(phi,w,rank,flagw)
  use mod_param, only : i1,k1,imax,kmax,periodic,px
  use mod_common,only : dt
  use mod_mesh,  only : dru, dzp,dzw
  implicit none  
  include 'mpif.h'
  integer                       :: rank,i,flagw
  real(8), dimension(0:i1,0:k1) :: phi,w
  real(8)                       :: bulk

  if (rank.eq.px-1)then
    bulk=0
    do i=1,imax
      bulk = bulk + dru(i)*w(i,kmax)/2.0
    enddo
    if (flagw .eq. 0) then
      phi(:,k1)  =phi(:,k1)  -bulk*dt*(phi(:,k1)  -phi(:,kmax)  )/(dzp(kmax))
    else 
      phi(:,kmax)=phi(:,kmax)-bulk*dt*(phi(:,kmax)-phi(:,kmax-1))/(dzw(kmax-1))
    endif
  endif
end subroutine 

!!*************************************************************************************
!!   Apply the boundary conditions for the velocity using the mass flux
!!************************  *************************************************************
subroutine bound_m(Ubound,Wbound,W_out,Rbound,W_in,rank, step)
  use mod_param, only : k,i,kmax,imax,k1,i1,px, periodic
  use mod_mesh,  only : dzw,dru,rp
  use mod_common,only : dt,rold,rnew
  implicit none
  include 'mpif.h'
  integer,                       intent(IN) :: rank, step
  real(8), dimension(0:i1,0:k1), intent(IN) :: rbound
  real(8), dimension(0:i1),      intent(IN) :: W_in  
  real(8), dimension(0:i1,0:k1), intent(OUT):: ubound, wbound, w_out
  real(8), dimension(0:i1) :: tmp
  real(8), dimension(1:imax) :: wr
  integer :: ierr
  real(8) :: Ub,flux,flux_tot,deltaW,wfunc
  
  call bound_v(ubound,wbound,W_in,W_out,rank, step)

  if (periodic .eq. 1) return
  
  wr = 0
  Ub = 0.
  flux = 0.0
  if (rank.eq.px-1) then
    Ub = 0.
    do i=1,imax
      wr(i) = W_out(i,kmax)
    enddo
!    do i=0,i1
!      Wbound(i,kmax) = 2.0*Wbound(i,kmax-1) - Wbound(i,kmax-2) !NOTE: CHANGE WITH SIMONE
!      Wbound(i,kmax) = Wbound(i,kmax)*0.5*(Rbound(i,kmax)+Rbound(i,k1))
!    enddo
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

subroutine initialize_solution(rank, w, u,c, mut,alphat, &
                               win, Re, systemsolve, select_init)
  use mod_tm,    only : turb_model
  use mod_tdm,   only : turbdiff_model
  use mod_eos,   only : eos_model
  use mod_param, only : k1,i1,imax,imax_old, kelem_old,px, i,k,Pr,Tw_bot,Tw_top, modifDiffTerm
  use mod_mesh,  only : y_fa,y_cv, mesh
  implicit none
  integer,                        intent(IN) :: rank,systemsolve,select_init
  real(8),                        intent(IN) :: Re
  real(8), dimension(0:i1, 0:k1), intent(OUT):: w,u,c,mut,alphat
  real(8), dimension(0:i1),       intent(OUT):: win
  real(8), dimension(0:i1) :: nuSAin,pkin,kin,epsin,omin,mutin,v2in, &
                              Prtin,alphatin,ktin,epstin,pktin
  character(len=5)  :: Re_str
  character(len=1)  :: mod_str
  integer           :: Re_int
  real(8)           :: gridSize
  character(len=100) :: fname

  ! inialize from inflow profile
  if (select_init.eq.1) then
    if (rank.eq.0) write(*,*) 'Initializing flow with inflow'
    Re_int = int(eos_model%Re)
    write(Re_str,'(I5.5)') Re_int
    write(mod_str,'(I1.1)') modifDiffTerm
    
    fname = 'Inflow_'//trim(turb_model%name)     &
          //'_'//trim(turbdiff_model%name) &
          //'_'//'mod'//mod_str            &
          //"_"//trim(eos_model%name)      &
          //"_"//Re_str
    Re_int = int(Re)
    open(29,file=trim(mesh%name)//'/'//trim(fname)//'.dat',form='unformatted')

    read(29) win(:),kin,epsin,v2in,omin,nuSAin,mutin,pkin, alphatin,prtin, ktin, epstin,pktin
    close(29)
    do k=0,k1
      w(:,k)  = win(:)
      mut(:,k)= mutin(:)
      c(:,k)= 0.
      alphat(:,k) = alphatin(:)
    enddo
    
    call turb_model%init_w_inflow(nuSAin,pkin,kin,epsin,omin,mutin,v2in)
    call turbdiff_model%init_w_inflow(Prtin,alphatin,ktin,epstin,pktin)

  !initialize with laminar analytical solution
  else if (select_init .eq. 2) then
    call interpolate_solution(imax_old+1, (kelem_old/px)+1, rank, px)
  else
    if (rank.eq.0) write(*,*) 'Initializing flow from scratch'
    gridSize = y_fa(imax)
    do i=0,i1!imax         
      if (systemsolve.eq.1) w(i,:)  = Re*mesh%dpdz/6*3/2.*(1-(y_cv(i)/0.5)**2); c(i,:) = 0.0; 
      if (systemsolve.eq.2) then
        w(i,:)  = Re*mesh%dpdz*y_cv(i)*0.5*(gridSize-y_cv(i))              !channel
        c(i,:)= 0.!Tw_bot -Tw_bot* (y_cv(i)/2.0)!(Tw_bot+Tw_top)*0.5!.05!1. !- (y_cv(i)/2.0)
        ! if (rank .eq. 1) 
        ! write(*,*) c(i,1), i
        ! c(i,:)  = Pr*w(i,:)
      endif
      if (systemsolve.eq.3) w(i,:)  = Re*mesh%dpdz*0.5*((gridSize*gridSize)-(y_cv(i)*y_cv(i))) !bl
      if (systemsolve.eq.4) then
         w(i,:)  = 1.; u=0.;  win=1.; mutin=1.!bl 
      endif
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
  use mod_param, only : k1,i1,imax,kmax,alphac,k,i,periodic,Qsource,Re,Pr
  use mod_math,  only : matrixIdir
  use mod_mesh,  only : dzw,dru,drp,rp,ru,top_bcvalue1,bot_bcvalue1
  use mod_common,only : cnew,ekh,ekhi,ekhk,alphat,ekhe,ekhek,ekhei
  implicit none
  real(8), dimension(0:i1,0:k1), intent(IN) :: Utmp, Wtmp, Rtmp
  integer,                       intent(IN) :: rank
  real(8),                       intent(OUT):: resC
  real(8), dimension(0:i1,0:k1) :: dnew,dimpl
  real(8), dimension(imax)      :: a,b,c,rhs
  integer  :: ierr
  
  resC   = 0.0; dnew   = 0.0; dimpl = 0.0;
  ! call calc_avg_quantities(wnew, rnew, cnew,massflow, hbulk, tbulk,vbulk, Re_bulk,kmax)

  call advecc(dnew,dimpl,cnew,Utmp,Wtmp,rank,periodic,.true.)
  ! call diffc(dnew,cnew,ekh,ekhi,ekhk,alphat,1.,Rtmp,0)
  call diffc(dnew,cnew,ekhe,ekhei,ekhek,alphat,1.,Rtmp,10)


  do k=1,kmax
    do i=1,imax
      a(i) = -Ru(i-1)*ekhei(i-1,k)/(dRp(i-1)*Rp(i)*dru(i))/Rtmp(i,k)
      c(i) = -Ru(i  )*ekhei(i  ,k)/(dRp(i  )*Rp(i)*dru(i))/Rtmp(i,k)

      ! a(i) = -Ru(i-1)*(ekhi(i-1,k)+0.5*(alphat(i,k)+alphat(i-1,k)))/(dRp(i-1)*Rp(i)*dru(i))/Rtmp(i,k)
      ! c(i) = -Ru(i  )*(ekhi(i  ,k)+0.5*(alphat(i,k)+alphat(i+1,k)))/(dRp(i  )*Rp(i)*dru(i))/Rtmp(i,k)

      ! a(i) = -Ru(i-1)*(0.5*(ekhe(i,k)+ekhe(i-1,k)))/(dRp(i-1)*Rp(i)*dru(i))/Rtmp(i,k)
      ! c(i) = -Ru(i  )*(0.5*(ekhe(i,k)+ekhe(i+1,k)))/(dRp(i  )*Rp(i)*dru(i))/Rtmp(i,k)

      b(i) = (-a(i)-c(i) + dimpl(i,k) )        
      rhs(i) = dnew(i,k) + ((1-alphac)/alphac)*b(i)*cnew(i,k) + Qsource/(Re*Pr)
    enddo

    i=1
    b(i)=b(i)+bot_bcvalue1(k)*a(i) !symmetry or nothing
    rhs(i) = dnew(i,k) - (1-bot_bcvalue1(k))*a(i)*cNew(i-1,k) + ((1-alphac)/alphac)*b(i)*cNew(i,k) + Qsource/(Re*Pr) !nothing or value

    i=imax
    b(i)=b(i)+top_bcvalue1(k)*c(i) !symmetry or nothing
    rhs(i) = dnew(i,k) - (1-top_bcvalue1(k))*c(i)*cNew(i+1,k) + ((1-alphac)/alphac)*b(i)*cNew(i,k) + Qsource/(Re*Pr) !nothing or value

    call matrixIdir(imax,a,b/alphac,c,rhs)

    do i=1,imax
      resC = resC + ((cnew(i,k) - rhs(i))/(cnew(i,k)+1.0e-20))**2.0
      cnew(i,k) = max(rhs(i), 0.0)
    enddo    
  enddo
  ! if (periodic.eq.1) then
  !   cnew = 0.0; resC=0.0;
  ! endif

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
  use mod_param, only : k1,i1,kmax,imax,k,i,Fr_1,periodic,qwall
  use mod_math,  only : matrixIdir
  use mod_mesh,  only : mesh,dzw,dru,drp,rp,ru,top_bcnovalue,bot_bcnovalue,ubot_bcvalue
  use mod_common,only : rnew,ekme,unew,wnew,dUdt,dWdt,ekm,dt,peclet
  implicit none
      
  integer rank
  real(8), dimension(imax)   :: a, b, c, rhs
  real(8), dimension(imax-1) :: au, bu, cu, rhsu
  real(8) dnew(0:i1,0:k1)
  real(8) dif,alpha,rhoa,rhob,rhoc
  dif=0.0

  !>********************************************************************
  !!     CALCULATE advection, diffusion and Force in r-direction
  !!     at the old(n-1) and new(n) timelevels
  !!********************************************************************
  dnew = 0.0
  call advecu(dnew,Unew,Wnew,Rnew) ! new
  call diffu (dnew,Unew,Wnew,ekme,dif,mesh%numDomain) ! new

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
  call advecw(dnew,Unew,Wnew,Rnew,ekm,peclet)
  call diffw (dnew,Unew,Wnew,ekme,dif,mesh%numDomain)  ! new

  if (periodic.eq.1) dnew = dnew + mesh%dpdz

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

subroutine debug(rank)
  use mod_param, only :k1,i1,kmax,imax
  use mod_common, only : wnew,cnew,rnew,p,temp
  implicit none
  integer rank

  call write_2D_vector(wnew,i1,k1,rank,'w')


end subroutine


!>*************************************************************************************
!!           cmpinf(Bulk,Stress)
!!
!!*************************************************************************************
subroutine cmpinf(Bulk,Stress)
  use mod_param, only : kmax,imax,k,i,px,Re
  use mod_common
  use mod_mesh, only : rp,dru,walldist
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








