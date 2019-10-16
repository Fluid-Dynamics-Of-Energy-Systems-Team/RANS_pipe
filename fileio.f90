
!********************************************************************
!     write inflow file
!********************************************************************
subroutine Inflow_output(rank,istap) !tubstress,rank,istap)
  use mod_param
  use mod_common
  implicit none
  
  integer rank,istap
  character*5 inflow

  if (rank.eq.px/2) then
    if (turbmod.eq.0) open(29,file='0/Inflow', form='unformatted')
    if (turbmod.eq.1) open(29,file='SA/Inflow',form='unformatted')
    if (turbmod.eq.2) open(29,file='MK/Inflow',form='unformatted')
    if (turbmod.eq.3) open(29,file='VF/Inflow',form='unformatted')
    if (turbmod.eq.4) open(29,file='OM/Inflow',form='unformatted')
    write(29) Wnew(:,kmax/2),knew(:,kmax/2),enew(:,kmax/2),v2new(:,kmax/2),omNEW(:,kmax/2), &
      nuSAnew(:,kmax/2),ekmt(:,kmax/2),Pk(:,kmax/2)
    close(29)
  endif

end

subroutine inflow_output_upd(rank,istap) !tubstress,rank,istap)
  use mod_param
  use mod_common
  use mod_common2
  implicit none
  
  integer rank,istap
  real(8), dimension(0:i1) :: p_nuSA,p_k,p_eps,p_om,p_v2,p_Pk,p_bF1
  real(8), dimension(1:imax) :: p_bF2
  
  character*5 inflow

  if (rank.eq.px/2) then
    if (systemsolve.eq.1) open(29,file='pipe/Inflow.dat', form='unformatted')
    if (systemsolve.eq.2) open(29,file='channel/Inflow.dat', form='unformatted')
    if (systemsolve.eq.3) open(29,file='bl/Inflow.dat', form='unformatted')
  
    k = kmax/2
    call turb_model%get_profile(p_nuSA,p_k,p_eps,p_om,p_v2,p_Pk,p_bF1,p_bF2,k)

    write(29) Wnew(:,kmax/2),p_k,p_eps,p_v2,p_om,p_nuSA,ekmt(:,kmax/2),p_pk
    close(29)

    if (systemsolve.eq.1) open(29,file='pipe/Inflow.csv')
    if (systemsolve.eq.2) open(29,file='channel/Inflow.csv')
    if (systemsolve.eq.3) open(29,file='bl/Inflow.csv')

    write(29, '(15a20)' ) 'y','u','w','h','T', &
                          'rho','k','eps','v2','om', &
                          'nuSA','mut','Pk','bF1','bF2'
    do i=1,imax
      write(29, '(E20.12,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12)') &
                y_cv(i),unew(i,k),Wnew(i,k),cnew(i,k),temp(i,k), &
                rnew(i,k),p_k(i),p_eps(i),p_v2(i),p_om(i), &
                p_nuSA(i),ekmt(i,k),p_Pk(i),p_bf1(i),p_bf2(i)
    enddo
    close(29)
  endif

  

end


!***************************************************************************************
!     read istep from the restart file to identify the corresponding inflow profile
!***************************************************************************************
subroutine loadRestart(startStep,rank)
  use mod_param
  use mod_common
  implicit none

  character*5 cha
  integer startStep,ierr,rank

  real*8     UNEWR(0:i1,0:k1old)
  real*8     UOLDR(0:i1,0:k1old)
  real*8     WNEWR(0:i1,0:k1old)
  real*8     WOLDR(0:i1,0:k1old)
  real*8     CNEWR(0:i1,0:k1old)
  real*8     KNEWR(0:i1,0:k1old)
  real*8     ENEWR(0:i1,0:k1old)
  real*8    V2NEWR(0:i1,0:k1old)
  real*8  NUSANEWR(0:i1,0:k1old)
  real*8    omNEWR(0:i1,0:k1old)
  real*8    PKR(0:i1,0:k1old)

  write(cha,'(I5.5)')rank
  if (turbmod.eq.0) open(19,file= '0/Restart/start_stop.'//cha,form='unformatted')
  if (turbmod.eq.1) open(19,file='SA/Restart/start_stop.'//cha,form='unformatted')
  if (turbmod.eq.2) open(19,file='MK/Restart/start_stop.'//cha,form='unformatted')
  if (turbmod.eq.3) open(19,file='VF/Restart/start_stop.'//cha,form='unformatted')
  if (turbmod.eq.4) open(19,file='OM/Restart/start_stop.'//cha,form='unformatted')
  read(19) startStep
  read(19) UNEWR,UOLDR,WNEWR,WOLDR,CNEWR,KNEWR,ENEWR,V2NEWR,omNEWR,nuSANEWR, PKR
  close(19)

  do i=0,i1
    do k=0,k1
      UNEW(i,k) =    UNEWR(i,k1old*k/k1)
      UOLD(i,k) =    UOLDR(i,k1old*k/k1)
      WNEW(i,k) =    WNEWR(i,k1old*k/k1)
      WOLD(i,k) =    WOLDR(i,k1old*k/k1)
      CNEW(i,k) =    CNEWR(i,k1old*k/k1)
      KNEW(i,k) =    KNEWR(i,k1old*k/k1)
      ENEW(i,k) =    ENEWR(i,k1old*k/k1)
      V2NEW(i,k) =   V2NEWR(i,k1old*k/k1)
      nuSAnew(i,k) = nuSAnewR(i,k1old*k/k1)
      omNew(i,k) =   omNewR(i,k1old*k/k1)
      Pk(i,k) =      PKR(i,k1old*k/k1)
    enddo
  enddo


end

!***************************************************************************************
!     store istep in the restart file to identify the corresponding inflow profile
!***************************************************************************************
subroutine saveRestart(rank)
  use mod_param
  use mod_common
  implicit none

  character*5 cha
  integer rank
  write(cha,'(I5.5)')rank
  if (turbmod.eq.0) open(19,file= '0/Restart/start_stop.'//cha,form='unformatted')
  if (turbmod.eq.1) open(19,file='SA/Restart/start_stop.'//cha,form='unformatted')
  if (turbmod.eq.2) open(19,file='MK/Restart/start_stop.'//cha,form='unformatted')
  if (turbmod.eq.3) open(19,file='VF/Restart/start_stop.'//cha,form='unformatted')
  if (turbmod.eq.4) open(19,file='OM/Restart/start_stop.'//cha,form='unformatted')
  write(19) istep
  write(19) UNEW,UOLD,WNEW,WOLD,CNEW,KNEW,ENEW,V2NEW,omNEW,nuSANEW,Pk
  close(19)
end




!***************************************************************************************
!   writes the inflow profile
!***************************************************************************************
subroutine outputProfile(rank)
  use mod_param
  use mod_common
  implicit none
  
  integer rank
  ! radius= 1, velocity=3, temperature=5, density=6, turb. kine=7,
  ! epsilon=8, v2=9, nuSA=10, mut=11
  if (rank == 0) then
    open(19,file='profile')
    k = 1
    do i=1,imax
      write(19,'(14E18.6)') y_cv(i),  unew(i,k),Wnew(i,k), cnew(i,k), temp(i,k), &
                            rnew(i,k),knew(i,k),enew(i,k),v2new(i,k),omNew(i,k), &
                            nuSAnew(i,k),ekmt(i,k),bf1(i,k), bf2(i,k)
    enddo
    close(19)
  endif

end

!***************************************************************************************
!   writes the output file in 2D gnuplot style
!***************************************************************************************
subroutine outputX_h_upd(rank,istap)
  use mod_param
  use mod_common
  use mod_common2
  implicit none

  include 'mpif.h'
  character*5 cha
  integer rank,ierr,istap
  real*8 massflow(kmax),enthflow(kmax),enth_b(kmax),Twall(kmax),Tbulk(kmax),Qp(kmax), &
    tmp(kmax),massfl,enth,laminter,tempinter,cpinter,Gb(kmax),rhob(kmax),ekmb(kmax),Ub(kmax), &
    addedHeat,addedHeatTot,subheat,subheatTot,w_c

  twall    = 0.0
  Qp       = 0.0
  massflow = 0.0
  enthflow = 0.0
  w_c      = 0.0
  Gb       = 0.0

  do k=1,kmax
    do i=1,imax
      massfl = 0.5*rnew(i,k)*(Wnew(i,k)+Wnew(i,k-1))*rp(i)*dru(i)
      massflow(k) = massflow(k) + massfl
      enthflow(k) = enthflow(k) + massfl*Cnew(i,k)
    enddo
  enddo
  
  enth_b=enthflow/massflow
  do k=1,kmax
    w_c=(Cnew(i1,k)+Cnew(imax,k))/2.
    call eos_model%set_w_enth(enth_b(k), 'D', rhob(k))
    call eos_model%set_w_enth(enth_b(k), 'V', ekmb(k))
    call eos_model%set_w_enth(enth_b(k), 'V', ekmb(k))
    call eos_model%set_w_enth(w_c,       "T", Twall(k))
    call eos_model%set_w_enth(enth_b(k), "T", Tbulk(k))
    Gb(k) =massflow(k)/(4.*atan(1.)*Ru(imax)**2.)
  enddo

  ! bulk stuff
  write(cha,'(I5.5)')rank
  if (turbmod.eq.0) open(9,file='0/profX.'//cha)
  if (turbmod.eq.1) open(9,file='SA/profX.'//cha)
  if (turbmod.eq.2) open(9,file='MK/profX.'//cha)
  if (turbmod.eq.3) open(9,file='VF/profX.'//cha)
  if (turbmod.eq.4) open(9,file='OM/profX.'//cha)
  do k=1,kmax
    write(9,'(8E18.6)') (k+rank*kmax)*dz, massflow(k),enthflow(k),enth_b(k), &
                        Twall(k),Tbulk(k),Gb(k)/Gb(1),Gb(k)/Gb(1)/rhob(k)
  enddo
  close(9)
  return
end

!***************************************************************************************
!   writes the output file in 2D tecplot style
!***************************************************************************************
subroutine output2d_upd(rank,istap)
  use mod_param
  use mod_common
  use mod_common2
  implicit none
  include 'mpif.h'
  character*5 cha
  real*8 pecletx,peclety,pecletz
  real*8 massflow(kmax),enthflow(kmax),enth_b(kmax),Twall(kmax),Tbulk(kmax), &
         massfl,w_c,ndt
  integer rank,istap,jstart
  write(cha,'(I5.5)')rank


  twall    = 0.0
  massflow = 0.0
  enthflow = 0.0
  w_c      = 0.0

  do k=1,kmax
    do i=1,imax
      massfl = 0.5*rnew(i,k)*(Wnew(i,k)+Wnew(i,k-1))*rp(i)*dru(i)
      massflow(k) = massflow(k) + massfl
      enthflow(k) = enthflow(k) + massfl*Cnew(i,k)
    enddo
  enddo

  enth_b=enthflow/massflow
  do k=1,kmax
    w_c=(Cnew(i1,k)+Cnew(imax,k))/2.
    call eos_model%set_w_enth(w_c,      "T", Twall(k))
    call eos_model%set_w_enth(enth_b(k),"T", Tbulk(k))
  enddo

  if (turbmod.eq.0) open(15,file='0/tecp.'//cha)
  if (turbmod.eq.1) open(15,file='SA/tecp.'//cha)
  if (turbmod.eq.2) open(15,file='MK/tecp.'//cha)
  if (turbmod.eq.3) open(15,file='VF/tecp.'//cha)
  if (turbmod.eq.4) open(15,file='OM/tecp.'//cha)

  if (rank.eq.0) then
    write(15,*) 'VARIABLES ="X","Y","U","W","P","C","T","k","eps", "v2","omega","nuSA","yplus","RHO","Pe","mu","mut"'
    write(15,*) 'ZONE I=  ', imax+2,' J=  ',(kmax+2)*px,' F=POINT '
  endif

  do k=0,k1
    do i=0,i1
      write(15,'(17ES24.10E3)')  (k+rank*kmax)*dz, y_cv(i), &
        0.5*(unew(max(i-1,0),k)+unew(i,k)), &
        Wnew(i,k), p(min(max(i,1),imax),min(max(k,1),kmax)), cnew(i,k), temp(i,k), &
        knew(i,k),enew(i,k),v2new(i,k), &
        omNew(i,k),nuSAnew(i,k),yp(i,k),rnew(i,k),peclet(i,k),ekm(i,k),ekmt(i,k)
    enddo
  enddo
  
  close(15)

end