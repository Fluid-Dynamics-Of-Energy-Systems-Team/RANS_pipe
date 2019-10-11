
!!*********************************************************************
!!     Calculate enthalpy at the wall boundary condition for isothermal
!!********************************************************************
subroutine funcIsothermalEnthBC_upd(Twall_bc)
  use mod_param
  use mod_common
  use mod_common2
  implicit none
  real*8 Twall_bc
  call eos_model%set_w_temp(Twall_bc,"H",enth_wall)
end


!!********************************************************************
!!     Newton solver for wall boundary condition
!!********************************************************************
subroutine funcNewtonSolve_upd(enth_i1, enth_imax)
  use mod_param
  implicit none
  real*8 enth_i1, enth_imax
  if (EOSmode.eq.0) call funcNewtonSolveIG(enth_i1, enth_imax)
  if (EOSmode.eq.1) call funcNewtonSolveRG_upd(enth_i1, enth_imax)
end

subroutine funcNewtonSolveIG(enth_i1, enth_imax)
use mod_param
use mod_common
implicit none

real*8 enth_i1, enth_imax, ekh_imax

  ekh_imax = 1./(Re*Pr)
  enth_i1 = enth_imax + dRp(imax)*Qwall/(ekh_imax*Re*Pr) ! new
  
end


!!********************************************************************
!!     Newton solver for wall boundary condition RG
!!********************************************************************
subroutine funcNewtonSolveRG_upd(enth_i1, enth_imax)
  implicit none
  integer nIterNewton,success
  real*8 enth_i1, enth_imax, fxValue, fxValue1, ekh_imax
  success = 1
  fxValue = 1000.0
  nIterNewton = 0

  if (enth_imax.gt.2)    enth_imax =  2.0
  if (enth_imax.lt.-0.1) enth_imax = -0.1

  enth_i1 = enth_imax

  do while (abs(fxValue).gt.1.0e-10)
    call funcNewtonBC_upd(enth_i1,        enth_imax, fxValue)
    call funcNewtonBC_upd(enth_i1+1.0e-8, enth_imax, fxValue1)

    enth_i1 = enth_i1 - fxValue/((fxValue1-fxValue)/1.0e-8)

    if (nIterNewton.gt.200) then
      fxValue = 0.0
      success = 0
    endif

    nIterNewton = nIterNewton + 1
    ! write (*,*) 'newton iter: ', nIterNewton,enth_i1,enth_imax,fxValue
  enddo

  if (success.eq.0) then
    write (*,*) 'newton didnt converge, enthimax= ',enth_imax,', ', nIterNewton, ', ', enth_i1
    stop
  endif
end

!>********************************************************************
!!     function for wall boundary condition used by Newton solver
!!********************************************************************
subroutine funcNewtonBC_upd(enth, enthIMAX, fxValue)
  use mod_param
  use mod_common
  use mod_common2
  implicit none
  real*8 enth,lamOcpinter,enthIMAX,fxValue
  call eos_model%set_w_enth(0.5*(enth+enthIMAX), 'L', lamocpinter)
  lamocpinter = lamocpinter*(eos_model%Re*eos_model%Pr)
  fxValue = enth - enthIMAX - dRp(imax)*Qwall/lamocpinter 
end

!>********************************************************************
!!     PIG equation of state
!!********************************************************************
subroutine state_upd(enth,rho,mu,lam,tp,be,istap,rank)
  use mod_param
  use mod_common2
  use mod_common
  implicit none
  integer istap, rank
  real(8) enth(0:i1,0:k1),rho(0:i1,0:k1),mu (0:i1,0:k1), &
          lam(0:i1,0:k1),tp(0:i1,0:k1),be(0:i1,0:k1)
  real    enthface

  if(isothermalBC.eq.1) then
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
  else 
    do k=0,k1
      call funcNewtonSolveRG_upd(enth(i1,k), enth(imax,k))
      if (rank.eq.0.and.k.lt.K_start_heat) enth(i1,k)=enth(imax,k)
      do i=0,i1
          call eos_model%set_w_enth(enth(i,k),"D", rho(i,k))
          call eos_model%set_w_enth(enth(i,k),"V", mu(i,k))
          call eos_model%set_w_enth(enth(i,k),"C", Cp(i,k))
          call eos_model%set_w_enth(enth(i,k),"L", lam(i,k))
          call eos_model%set_w_enth(enth(i,k),"T", tp(i,k))
          call eos_model%set_w_enth(enth(i,k),"B", be(i,k)) 
      enddo
    enddo
  endif

  do k=0,kmax
    do i=0,imax
      enthface = 0.5*(enth(i,k)+enth(i+1,k))
      call eos_model%set_w_enth(enthface,"C", Cpi(i,k))
      call eos_model%set_w_enth(enthface,"L", ekhi(i,k))
      call eos_model%set_w_enth(enthface,"V", ekmi(i,k)) 
      enthface = 0.5*(enth(i,k)+enth(i,k+1))
      call eos_model%set_w_enth(enthface,"C", Cpk(i,k))
      call eos_model%set_w_enth(enthface,"L", ekhk(i,k))
      call eos_model%set_w_enth(enthface,"V", ekmk(i,k)) 
    enddo
  enddo
  return

end subroutine state_upd




!>********************************************************************
!!     poisson solver
!!********************************************************************
subroutine fillps(rank)
  use mod_param
  use mod_common
  implicit none

  include 'mpif.h'

  integer ierr,rank

  real*8 sumps,sumps_tot
  !
  !     *** Fill the right hand for the poisson solver. ***
  !
  sumps = 0.0

  do  k=1,kmax
    do i=1,imax
      p(i,k)  = ( &
        ( Ru(i)*dUdt(i,k) - Ru(i-1)*dUdt(i-1,k) )/( Rp(i)*dru(i)) &
        + (     dWdt(i,k) -         dWdt(i,k-1) )/( dz         ) )/dt &
        +      (rnew(i,k)-rold(i,k))/(dt*dt)

      qcrit(i,k) = p(i,k)*dt

      sumps = sumps + p(i,k)*dru(i)*dz
    enddo
  enddo


  call mpi_allreduce(sumps,sumps_tot,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)

end








!>********************************************************************
!!     correc
!!********************************************************************
subroutine correc(rank,setold)
  use mod_param
  use mod_common
  implicit none
      
  integer rank,setold
  real*8 pplus_w(imax)

  do k=1,kmax
    do i=1,imax-1
      dUdt(i,k)=dUdt(i,k)-dt*(p(i+1,k)-p(i,k))/dRp(i) !(Rp(i+1)-Rp(i))
    enddo
  enddo

  do k=1,kmax-1
    do i=1,imax
      dWdt(i,k)=dWdt(i,k)-dt*(p(i,k+1)-p(i,k))/dz
    enddo
  enddo

  !     if periodic is true, don't overwrite pplus on px-1
  call pshiftb_w(p,pplus_w,rank)
  if ((rank.eq. px-1).and.(periodic.ne.1)) then
    do i=1,imax
      pplus_w(i)=p(i,kmax)
    enddo
  endif
      
  do i=1,imax
    dWdt(i,kmax)=dWdt(i,kmax)-dt*(pplus_w(i)-p(i,kmax))/dz
  enddo
  if (setold.eq.0)then
    do k=0,kmax
      do i=0,imax
        dudt(i,k)=dUdt(i,k)/(0.5*(rnew(i,k)+rnew(i+1,k))) !dudt(i,k)=dUdt(i,k)/(0.5*(drdt(i,k)+drdt(i+1,k)))
        dwdt(i,k)=dWdt(i,k)/(0.5*(rnew(i,k)+rnew(i,k+1))) !dwdt(i,k)=dWdt(i,k)/(0.5*(drdt(i,k)+drdt(i,k+1)))
      enddo
    enddo
  endif
      
  if (setold.eq.1)then
    Uold = Unew
    Wold = Wnew

    do k=0,kmax
      do i=0,imax
        Unew(i,k)=dUdt(i,k)/(0.5*(rnew(i,k)+rnew(i+1,k))) !Unew(i,k)=dUdt(i,k)/(0.5*(drdt(i,k)+drdt(i+1,k)))
        Wnew(i,k)=dWdt(i,k)/(0.5*(rnew(i,k)+rnew(i,k+1))) !Wnew(i,k)=dWdt(i,k)/(0.5*(drdt(i,k)+drdt(i,k+1)))
      enddo
    enddo

  endif


  rold = rnew


end




!>********************************************************************
!!     chkdt
!!********************************************************************
subroutine chkdt(rank,istap)
  use mod_param
  use mod_common
  implicit none
      
      include 'mpif.h'
  integer rank,ierr,istap
  real*8  tmp,Courant,dtmp,tmp1,tmp2,tmp3,dr2,dz2,kcoeff

  dt = dtmax

  do k=1,kmax
    do i=1,imax
      tmp = ( abs(Unew(i,k)) /  dRp(i) ) +   &  ! ( Rp(i+1)-Rp(i) ) ) + new
        ( abs(Wnew(i,k)) /         dz        )
      tmp = CFL/tmp
      dt  = min(dt, tmp)
    enddo
  enddo

  dtmp = dt

  call mpi_allreduce(dtmp,dt,1,mpi_real8,mpi_min,mpi_comm_world,ierr)

end




!>********************************************************************
!!     chkdiv
!!********************************************************************
subroutine chkdiv(rank)
  use mod_param
  use mod_common
  implicit none
      
      include 'mpif.h'
  integer rank,ierr,ll
  real*8   div,divmax,divbar,divmax_tot,divbar_tot,rhoip,rhoim,rhokp,rhokm
  divbar = 0.0
  divmax = 0.0

  do k=1,kmax
    do i=1,imax
      rhoip = 0.5*(rNew(i,k)+rNew(i+1,k))
      rhoim = 0.5*(rNew(i,k)+rNew(i-1,k))

      rhokp = 0.5*(rNew(i,k)+rNew(i,k+1))
      rhokm = 0.5*(rNew(i,k)+rNew(i,k-1))

      div = &
        (Ru(i)*Unew(i,k)*rhoip-Ru(i-1)*Unew(i-1,k)*rhoim)*dz/Rp(i)+ &
        (Wnew(i,k)*rhokp-Wnew(i,k-1)*rhokm)*dru(i)+ &
        (rNew(i,k)-rold(i,k))/dt*dru(i)*dz

      !     if (abs(div).gt.10e-6) write(6,*) i,k+kmax*rank,div

      divbar = divbar+div
      div    = abs(div)
      divmax = max(divmax,div)
    enddo
  enddo

  call mpi_allreduce(divbar,divbar_tot,1,mpi_real8,mpi_max,mpi_comm_world,ierr)
  call mpi_allreduce(divmax,divmax_tot,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)


!     if (rank.eq.0) write(6,100) divbar_tot,divmax_tot
! 100  format('Mass_loss/gain:Tot= ',e15.4,' Max= ',e15.4)

end






!>********************************************************************
!!     mkgrid
!!********************************************************************
subroutine mkgrid(rank)
  use mod_param
  use mod_common
  implicit none
  
  integer rank
  real*8  delta(imax),Yplus,X,rmax,drr,Rei
  real*8  pii,y,y1,y2,fA,fB,fC,fact,gridSize

  !     const = 0.25
  !      const = 0.65
  !     const = 0.50
  !******************************************************************
  pi    = 4.0*atan(1.0)
  Rei   = 1.0/Re
  dz    = 1.0*LoD/(kmax*px)
  ru(0) = 0

  fA = 0.12
  fB = 2.4
  if (systemSolve.eq.1) then
    numDomain = 1
    centerBC  = 1
    gridSize = 0.5
    dpdz =4.0
    if (rank.eq.0) print*,"************* SOLVING A PIPE FLOW *************!"
    write(*,*) centerBC
  elseif (systemSolve.eq.2) then
    numDomain = -1
    centerBC  = -1
    gridSize  = 2.0
    fA = 0.5
    fB = 4.6
    dpdz= 1.0
    if (rank.eq.0) print*,"************* SOLVING A CHANNEL FLOW *************!"
  elseif (systemSolve.eq.3) then
    numDomain = -1
    centerBC  = 1
    gridSize = 1.0
    dpdz= 1.0
    if (rank.eq.0) print*,"************* SOLVING A BOUNDARY LAYER FLOW *************!"
  else
    if (rank.eq.0) print '("systemSolve is ",i7," when it should be either 1 (pipe), 2(channel) or 3(BL)")',systemSolve
    stop
  endif

      

  do i = 1,imax
    fact = (i-0.)/(imax-0.)
    ru(i) = (1.-tanh(fB*(fA-fact))/tanh(fA*fB))
    ru(i) = ru(i)/(1.-tanh(fB*(fA-1.))/tanh(fA*fB))
    delta(i)=0.5*(ru(i)-ru(i-1))
  enddo

  do i=0,imax
    ru(i)=ru(i)/ru(imax)*gridSize
  enddo

  do i = 1 , imax
    Rp(i)  = (Ru(i)+Ru(i-1))/2.0
    dru(i) = (Ru(i)-Ru(i-1))
  enddo

  dru(i1) = dru(imax)
  Ru(i1) = Ru(imax) + dru(i1)
  Rp(i1) = Ru(imax) + dru(i1)/2.0
  dru(0)  = dru(1)
  Rp(0)  = Ru(0) - dru(0)/2.0

  do i = 0,imax
    drp(i) = Rp(i+1) - Rp(i)
  enddo

  if (centerBC.eq.-1) then ! two walls
    do i = 1,imax
      if (rp(i).le.1) then
        wallDist(i) = rp(i)
      else
        wallDist(i) = gridSize-rp(i)
      endif
    enddo
  else ! one wall
    do i = 1,imax
      wallDist(i) = gridSize - rp(i)
    enddo
  endif
  do i=0,i1
    y_cv(i)=rp(i)
    y_fa(i)=ru(i)
  enddo

  if (numDomain.eq.-1) then
    do i=0,i1
      ru(i)=1.0
      rp(i)=1.0
    enddo
  endif

  if (rank.eq.0) then
    open(11,file = 'grid.txt')
    write(11,*) Re, imax
    do i=1,imax
      Yplus = wallDist(i)*Re
      write(11,'(i5,4F12.6)') i,yplus,y_fa(i),y_cv(i),delta(max(1,i))
    enddo
    close(11)
  endif
  return
end



!>********************************************************************
!!     readOldGrid
!!********************************************************************
subroutine readOldGrid(rank)
  implicit none
  integer i,i1,imax,tmp,rank
  real*8 Re,Ru(0:200),Rp(0:200),yplus

  if (rank.eq.0) then
    open(11,file = 'grid.txt')
    read(11,*) Re, i1
    do i=0,i1
      Yplus = (0.5-Rp(i))*Re
      read(11,'(i5,4F12.6)') tmp,yplus,Ru(i),Rp(i)
    enddo
    close(11)
  endif
  stop
end


!>********************************************************************
!!     diffusion term in the z-direction, set as a source term...
!!********************************************************************
subroutine diffc(putout,putin,ek,eki,ekk,ekmt,sigma,rho,Ru,Rp,dru,dz,rank1,diffVersion)
  use mod_param
  implicit none
  integer   km,kp,rank1,diffVersion
  real*8     putout(0:i1,0:k1),putin(0:i1,0:k1), &
    rho(0:i1,0:k1),ek(0:i1,0:k1),eki(0:i1,0:k1),ekk(0:i1,0:k1), &
    ekmt(0:i1,0:k1),dru(0:i1),dz,Ru(0:i1),Rp(0:i1),sigma

  if (diffVersion == 1) then       ! Inverse SLS
    do k=1,kmax
      kp=k+1
      km=k-1
      do i=1,imax
        putout(i,k) = putout(i,k) + 1.0/rho(i,k)/sqrt(rho(i,k))*( &
   ( (ekk(i,k ) + 0.5*(ekmt(i,k)+ekmt(i,kp))/sigma)/sqrt(0.5*(rho(i,k)+rho(i,kp)))*(rho(i,kp)*putin(i,kp)-rho(i,k )*putin(i,k )) &
    -(ekk(i,km) + 0.5*(ekmt(i,k)+ekmt(i,km))/sigma)/sqrt(0.5*(rho(i,k)+rho(i,km)))*(rho(i,k )*putin(i,k )-rho(i,km)*putin(i,km)) &
          )/(dz*dz)   )
      enddo
    enddo
  elseif (diffVersion == 2) then   ! Aupoix
    do k=1,k1-1
      kp=k+1
      km=k-1
      do i=1,i1-1
        putout(i,k) = putout(i,k) + 1.0/rho(i,k)*( &
      ( (ekk(i,k ) + 0.5*(ekmt(i,k)+ekmt(i,kp))/sigma)/(0.5*(rho(i,k)+rho(i,kp)))*(rho(i,kp)*putin(i,kp)-rho(i,k )*putin(i,k )) &
       -(ekk(i,km) + 0.5*(ekmt(i,k)+ekmt(i,km))/sigma)/(0.5*(rho(i,k)+rho(i,km)))*(rho(i,k )*putin(i,k )-rho(i,km)*putin(i,km)) &
          )/(dz*dz)   )
      enddo
    enddo
  else                               ! Standard
    do k=1,k1-1
      kp=k+1
      km=k-1
      do i=1,i1-1
        putout(i,k) = putout(i,k) + 1.0/rho(i,k)*( &
          ( (ekk(i,k ) + 0.5*(ekmt(i,k)+ekmt(i,kp))/sigma)*(putin(i,kp)-putin(i,k )) &
          -(ekk(i,km) + 0.5*(ekmt(i,k)+ekmt(i,km))/sigma)*(putin(i,k )-putin(i,km))   )/(dz*dz)   )
      enddo
    enddo
  endif
end




     
!>*****************************************************************
!!     
!!     diffu calculates the diffusion of u-velocity, which is
!!     the velocity in the radial direction.
!!     
!!
!!     In formula:  (4 terms)
!!
!!     1  d                  1  d                     d
!!     - -- (r Sigma(r r)) + - ---- (Sigma(phi r)) + -- (Sigma(z r)) -
!!     r dr                  r dphi                  dz
!!
!!
!!     1
!!     - Sigma(phi phi)
!!     r
!!
!!     r   : direction  ---> explicit (subroutine diffu)
!!     phi : direction  ---> implicit (subroutine predic)
!!     z   : direction  ---> explicit (subroutine diffu)
!!     
!!     on input :
!!     
!!     putout            : advection part
!!     Uvel,Vvel,Wvel    : contain velocities at n-1
!!     ekm               : diffusion coefficients (for velocity) in
!!     center points of the grid cells
!!     dr,dphi,dz        : grid spacing in r, phi and z-direction
!!     i1,j1,k1          : parameters for array-dimensions
!!     ib,ie,jb,je,kb,ke : range of gridpoints for which the
!!     diffusion part has to be calculated
!!     Ru,Rp             : radial positions of the U-velocity
!!     component and the pressure location
!!     respectively
!!     
!!     on output :
!!     
!!     putout            : advection and diffusion part
!!     other parameters  : all unchanged
!!     
!!*****************************************************************
subroutine diffu (putout,Uvel,Wvel,ekme,Ru,Rp,dru,drp,dz,i1,k1,dif,numDom)
  implicit none

  integer  i,k,im,ip,km,kp,i1,k1,numDom
  real*8     putout(0:i1,0:k1),Uvel(0:i1,0:k1),Wvel(0:i1,0:k1), &
    ekme(0:i1,0:k1),dru(0:i1),drp(0:i1),dz,Ru(0:i1),Rp(0:i1), &
    epop,epom,dzi,divUim,divUip,divUi,dif

  if (numDom == -1) then

    dzi =1./dz
    do k=1,k1-1
      kp=k+1
      km=k-1
      do i=1,i1-1
        ip=i+1
        im=i-1

        epop = 0.25*(ekme(i,k)+ekme(ip,k) + ekme(ip,kp) + ekme(i,kp))
        epom = 0.25*(ekme(i,k)+ekme(ip,k) + ekme(ip,km) + ekme(i,km))

        putout(i,k) = putout(i,k) + &
          2.0*( Rp(ip)*ekme(ip,k)*(dif*(Uvel(ip,k)-Uvel(i ,k))/dru(ip)) &
          -Rp(i )*ekme(i ,k)*(dif*(Uvel(i ,k)-Uvel(im,k))/dru(i )) &
          )/(Ru(i)*drp(i)) &
          + &
          ( epop * ( (Uvel(i,kp)-Uvel(i,k))*dzi &
          + (Wvel(ip,k)-Wvel(i,k))/drp(i)) &
          - &
          epom * ( (Uvel(i,k)  -Uvel(i,km))*dzi &
          + (Wvel(ip,km)-Wvel(i,km))/drp(i)) &
          )*dzi
      enddo
    enddo

  elseif (numDom == 1) then
      
    dzi =1./dz
    do k=1,k1-1
      kp=k+1
      km=k-1
      do i=1,i1-1
        ip=i+1
        im=i-1

        epop = 0.25*(ekme(i,k)+ekme(ip,k) + ekme(ip,kp) + ekme(i,kp))
        epom = 0.25*(ekme(i,k)+ekme(ip,k) + ekme(ip,km) + ekme(i,km))

        divUim = (Ru(i)*Uvel(i,k) - Ru(im)*Uvel(im,k))/(Rp(i)*dru(i)) &
          + (        Wvel(i,k) -        Wvel(i,km))/dz

        divUip = (Ru(ip)*Uvel(ip,k)-Ru(i)*Uvel(i ,k ))/(Rp(ip)*dru(ip)) &
          + (         Wvel(ip,k)-      Wvel(ip,km))/dz

        divUi = ( Rp(ip)*(Uvel(ip,k)+Uvel(i,k)) - Rp(i)*(Uvel(i,k)+Uvel(im,k)) &
          )/(2.*Ru(i)*drp(i)) &!(Rp(ip)-Rp(i))) new
          + ((Wvel(ip,k)+Wvel(i,k))-(Wvel(ip,km)+Wvel(i,km)))/(2.*dz)

        putout(i,k) = putout(i,k) + &
          2.0*( Rp(ip)*ekme(ip,k)*(dif*(Uvel(ip,k)-Uvel(i ,k))/dru(ip) -1./3.*divUip) &
          -Rp(i )*ekme(i ,k)*(dif*(Uvel(i ,k)-Uvel(im,k))/dru(i ) -1./3.*divUim) &
          )/(Ru(i)*drp(i)) &
          + &
          ( epop * ( (Uvel(i,kp)-Uvel(i,k))*dzi &
          + (Wvel(ip,k)-Wvel(i,k))/drp(i)) &
          - &
          epom * ( (Uvel(i,k)  -Uvel(i,km))*dzi &
          + (Wvel(ip,km)-Wvel(i,km))/drp(i)) &
          )*dzi &
          - &
          (ekme(i,k)+ekme(ip,k))/Ru(i)*(Uvel(i,k)/Ru(i)-1./3.*divUi)
      enddo
    enddo
  endif



  return
end




!>*****************************************************************
!!     
!!     diffw calculates the diffusion of w-velocity, which is
!!     the velocity in the axial direction.
!!     
!!     
!!     In formula:  (3 terms)
!!     
!!     1  d                  1  d                     d
!!     - -- (r Sigma(r z)) + - ---- (Sigma(phi z)) + -- (Sigma(z z))
!!     r dr                  r dphi                  dz
!!     
!!     r   : direction  ---> explicit (subroutine diffw)
!!     phi : direction  ---> implicit (subroutine predic)
!!     z   : direction  ---> explicit (subroutine diffw)
!!     
!!     on input :
!!     
!!     putout            : advection part
!!     Uvel,Vvel,Wvel    : contain velocities at n-1
!!     ekm               : diffusion coefficients (for velocity) in
!!     center points of the grid cells
!!     dr,dphi,dz        : grid spacing in r, phi and z-direction
!!     i1,j1,k1          : parameters for array-dimensions
!!     ib,ie,jb,je,kb,ke : range of gridpoints for which the
!!     diffusion part has to be calculated
!!     Ru,Rp             : radial positions of the U-velocity
!!     component and the pressure location
!!     respectively
!!     
!!     on output :
!!     
!!     putout            : advection and diffusion part
!!     other parameters  : all unchanged
!!     
!!*****************************************************************
subroutine diffw(putout,Uvel,Wvel,ekme,Ru,Rp,dru,drp,dz,i1,k1,dif,numDom)
  implicit none
     

  integer  i,k,im,ip,km,kp,i1,k1,numDom
  real*8     putout(0:i1,0:k1),Uvel(0:i1,0:k1),Wvel(0:i1,0:k1),  &
    ekme(0:i1,0:k1),dru(0:i1),drp(0:i1),dz,Ru(0:i1),Rp(0:i1), &
    epop,emop,divUkm,divUkp,dif

  if (numDom == -1) then
    do k=1,k1-1
      kp=k+1
      km=k-1
      do i=1,i1-1
        ip=i+1
        im=i-1

        epop = 0.25*(ekme(i,k)+ekme(i,kp) + ekme(ip,k) + ekme(ip,kp) )
        emop = 0.25*(ekme(i,k)+ekme(i,kp) + ekme(im,k) + ekme(im,kp) )

        putout(i,k) = putout(i,k) + &
          (Ru(i )*epop*(  (Uvel(i,kp) - Uvel(i,k))/dz &
          +dif*(Wvel(ip,k) - Wvel(i,k))/drp(i)) &!(Rp(ip)-Rp(i))) new
          - &
          Ru(im)*emop*(  (Uvel(im,kp) - Uvel(im,k))/dz &
          +dif*(Wvel(i,k)   - Wvel(im,k))/drp(im)) &!(Rp(i)-Rp(im))) new
          )/(Rp(i)*dru(i)) &
          + &
          (2.*ekme(i,kp)*((Wvel(i,kp)-Wvel(i,k ))/dz )- &
          2.*ekme(i,k )*((Wvel(i,k )-Wvel(i,km))/dz ))/dz
      enddo
    enddo

  elseif (numDom == 1) then
    do k=1,k1-1
      kp=k+1
      km=k-1
      do i=1,i1-1
        ip=i+1
        im=i-1

        epop = 0.25*(ekme(i,k)+ekme(i,kp) + ekme(ip,k) + ekme(ip,kp) )
        emop = 0.25*(ekme(i,k)+ekme(i,kp) + ekme(im,k) + ekme(im,kp) )

        divUkm = (Ru(i)*Uvel(i,k) - Ru(im)*Uvel(im,k))/(Rp(i)*dru(i)) &
          + (      Wvel(i,k) -        Wvel(i,km))/dz

        divUkp = (Ru(i)*Uvel(i,kp)- Ru(im)*Uvel(im, kp))/(Rp(i)*dru(i)) &
          + (      Wvel(i,kp)-        Wvel(i,  k ))/dz

        putout(i,k) = putout(i,k) + &
          (Ru(i )*epop*(  (Uvel(i,kp) - Uvel(i,k))/dz &
          +dif*(Wvel(ip,k) - Wvel(i,k))/drp(i)) &!(Rp(ip)-Rp(i))) new
          - &
          Ru(im)*emop*(  (Uvel(im,kp) - Uvel(im,k))/dz &
          +dif*(Wvel(i,k)   - Wvel(im,k))/drp(im)) &!(Rp(i)-Rp(im))) new
          )/(Rp(i)*dru(i)) &
          + &
          (2.*ekme(i,kp)*((Wvel(i,kp)-Wvel(i,k ))/dz - 1./3.*divUkp)- &
          2.*ekme(i,k )*((Wvel(i,k )-Wvel(i,km))/dz - 1./3.*divUkm))/dz
      enddo
    enddo
  endif
  return
end




    
!>********************************************************************
!!     
!!     advecc calculates the advection for a scalar variable, which is
!!     situated in the center point of the grid cell.
!!     
!!     In formula:
!!     
!!     1 d(ruC)     1 d(vC)     d(wC)
!!     - (  - ------  +  - -----  +  -----  )
!!     r   dr       r  dphi      dz
!!     
!!     on input :
!!     
!!     putout            : "empty" (initialised to zero)
!!     putin             : variable for which the advection has
!!     to be calculated
!!     U,V,W             : contain velocities at former timestep
!!     putinn            : contains subgrid energy at oldest timestep
!!     dr,dphi,dz        : grid spacing in r, phi and z-direction
!!     i1,j1,k1          : parameters for array-dimensions
!!     ib,ie,jb,je,kb,ke : range of gridpoints for which the
!!     advection has to be calculated
!!     Ru,Rp             : radial positions of the U-velocity
!!     component and the pressure location
!!     respectively
!!     
!!     on output :
!!     
!!     putout            : advection part
!!     other parameters  : all unchanged
!!     
!!********************************************************************
subroutine advecc(putout,dimpl,putin,U,W,Ru,Rp,dru,dz,i1,k1,rank,periodic,flagImpl)

  implicit none

  integer  i,k,im,ip,km,kp,i1,k1,ib,ie,kb,ke,rank,periodic
  real*8 putout(0:i1,0:k1),putin(0:i1,0:k1),dimpl(0:i1,0:k1)
  real*8 U(0:i1,0:k1),W(0:i1,0:k1),dru(0:i1),dz,Ru(0:i1),Rp(0:i1)
  real*8 eps,r1,rk1,re1,phi1,phik1,phie1,r2,rk2,re2,phi2,phik2,phie2,r3,rk3,re3,phi3,phik3,phie3,fak
  real*8 rhoip,rhoim,rhokp,rhokm
  real*8 dcubf(0:i1),dcwbf(0:i1)
  real*8 dcubb(0:i1),dcwbb(0:i1)
  real*8 dcu(0:i1,0:k1),dcw(0:i1,0:k1)
  real*8 cu(0:i1,0:k1),cw(0:i1,0:k1)
  logical flagImpl
  !
  cu = 0.0
  cw = 0.0

  dcu = 0.0
  dcw = 0.0

  ib = 1
  ie = i1-1

  kb = 1
  ke = k1-1
  !
  !     compute delta C and distribute cpu boundaries
  do k=0,ke
    do i=0,ie
      dcu(i,k) = putin(i+1,k)-putin(i,k)
      dcw(i,k) = putin(i,k+1)-putin(i,k)
    enddo
  enddo

  !     if (periodic.eq.1) then

  call shiftf(dcu,dcubf,rank)
  call shiftf(dcw,dcwbf,rank)

  do i=ib,ie
    dcu(i,0) = dcubf(i)
    dcw(i,0) = dcwbf(i)
  enddo

  !     endif

  if ((periodic.ne.1).and.(rank.eq.0)) then
    dcu(:,0) = 0.0
    dcw(:,0) = 0.0
  endif

  eps=1.0e-16
  fak=1.0/3.0

  !     calculate face value for C

  do k=kb,ke
    kp=k+1
    km=k-1
    do i=ib,ie
      ip=i+1
      im=i-1
      if (U(i,k).ge.(0.0)) then
        r1=(dcu(i,k)+eps)/(dcu(im,k)+eps)
        phi1=max(0.,min(2.*r1,min(fak*(1.+2.*r1),2.)))
        Cu(i,k)=putin(i,k)+0.5*phi1*(dcu(im,k))
      else
        r1=(dcu(i,k)+eps)/(dcu(ip,k)+eps)
        phi1=max(0.,min(2.*r1,min(fak*(1.+2.*r1),2.)))
        Cu(i,k)=putin(ip,k)+0.5*phi1*(-dcu(ip,k))
      endif
      if (W(i,k).ge.(0.0)) then
        r3=(dcw(i,k)+eps)/(dcw(i,km)+eps)
        phi3=max(0.,min(2.*r3,min(fak*(1.+2.*r3),2.)))
        Cw(i,k)=putin(i,k) + 0.5*phi3*(dcw(i,km))
      else
        r3=(dcw(i,k)+eps)/(dcw(i,kp)+eps)
        phi3=max(0.,min(2.*r3,min(fak*(1.+2.*r3),2.)))
        Cw(i,k)=putin(i,kp)+0.5*phi3*(-dcw(i,kp))
      endif
    enddo
  enddo

  !     if (periodic.eq.1) then

  call shiftf(cu,dcubf,rank)
  call shiftf(cw,dcwbf,rank)


  do i=ib,ie
    cu(i,0)  = dcubf(i)
    cw(i,0)  = dcwbf(i)
  enddo

  !     endif


  if ((periodic.ne.1).and.(rank.eq.0)) then
    do i=0,ie
      cu(i,0) = putin(i,0)
      cw(i,0) = putin(i,0)
    enddo
  endif

  !     adv = u dc/dz = d(u c)/dz -c du/dz
  do k=kb,ke
    km=k-1
    do i=ib,ie
      im=i-1
      putout(i,k) =     - (Ru(i)*U(i,k)*cu(i,k) - Ru(im)*U(im,k)*cu(im,k))/(Rp(i)*dru(i))  &
        - (      W(i,k)*cw(i,k) -        W(i,km)*cw(i,km))/(dz) &
        + putin(i,k)*(   (Ru(i)*U(i,k) - Ru(i-1)*U(i-1,k ))/(Rp(i)*dru(i)) &
        + (      W(i,k) -         W(i  ,km))/dz )
    enddo
  enddo

  !     Implicit part
  if (flagImpl .eqv. .true.) then
    do k=kb,ke
      km=k-1
      do i=ib,ie
        im=i-1
        putout(i,k) = putout(i,k) + W(i,k)*putin(i,k)/dz
        dimpl(i,k) = dimpl(i,k)  + W(i,k)/dz
      enddo
    enddo
  endif

end


!>********************************************************************
!!     
!!     advecu calculates the advection of the u-velocity, which is
!!     the velocity in the radial direction.
!!     
!!     In formula:
!!     
!!          1 d(ruu)     1 d(uv)     d(uw)     vv
!!     - (  - ------  +  - -----  +  -----  -  --  )
!!          r   dr       r  dphi      dz        r
!!     
!!     on input :
!!     
!!     putout            : "empty" (initialised to zero)
!!     Uvel,Vvel,Wvel    : contain velocities at former timestep
!!     Utmp              : contains velocity at oldest timestep
!!     dr,dphi,dz        : grid spacing in r, phi and z-direction
!!     i1,j1,k1          : parameters for array-dimensions
!!     ib,ie,jb,je,kb,ke : range of gridpoints for which the
!!     advection has to be calculated
!!     Ru,Rp             : radial positions of the U-velocity
!!     component and the pressure location
!!     respectively
!!     
!!     on output :
!!     
!!     putout            : advection part
!!     other parameters  : all unchanged
!!     
!!********************************************************************
subroutine advecu(putout,Uvel,Wvel,RHO,Ru,Rp,dru,drp,dz,i1,k1)
  implicit none


  integer  i,k,im,ip,km,kp,i1,k1,ib,ie,kb,ke
  real*8     putout(0:i1,0:k1),Uvel(0:i1,0:k1),Wvel(0:i1,0:k1), &
    dru(0:i1),drp(0:i1),dz,Ru(0:i1),Rp(0:i1)
  real*8 rho(0:i1,0:k1)
  real*8 rhoip,rhoim,rhokp,rhokm

  !     if (adv.eq.1) Uin=Uvel
  ib = 1
  ie = i1-1
  kb = 1
  ke = k1-1

  do k=kb,ke
    kp=k+1
    km=k-1
    do  i=ib,ie
      ip=i+1
      im=i-1

      rhokp=0.25*(rho(i,k)+rho(i,kp)+rho(ip,k)+rho(ip,kp))
      rhokm=0.25*(rho(i,k)+rho(i,km)+rho(ip,k)+rho(ip,km))

      putout(i,k) = 0.0

      putout(i,k) = - 0.25 * (  &
        (Rp(ip)*(Uvel(i,k)+Uvel(ip,k))*(Uvel(i,k)+Uvel(ip,k)) &
        *rho(ip,k) - &
        Rp(i )*(Uvel(im,k)+Uvel(i,k))*(Uvel(i,k)+Uvel(im,k)) &
        *rho(i ,k)  ) &
        / ( Ru(i) * drp(i))  &   ! ( Rp(ip)-Rp(i) ) ) new
        + &
        ( (Wvel(i,k) +Wvel(ip,k) )*(Uvel(i,k)+Uvel(i,kp))*rhokp - &
        (Wvel(i,km)+Wvel(ip,km))*(Uvel(i,k)+Uvel(i,km))*rhokm  ) &
        / ( dz ) &
        )
    enddo
  enddo
  return
end




!>********************************************************************
!!     
!!     advecw calculates the advection of the w-velocity, which is
!!     the velocity in the axial direction.
!!     
!!     In formula:
!!     
!!     1 d(ruw)     1 d(wv)     d(ww)
!!     - (  - ------  +  - -----  +  -----  )
!!     r   dr       r  dphi      dz
!!     
!!     on input :
!!     
!!     putout            : "empty" (initialised to zero)
!!     Uvel,Vvel,Wvel    : contain velocities at former timestep
!!     Wtmp              : contains velocity at oldest timestep
!!     dr,dphi,dz        : grid spacing in r, phi and z-direction
!!     i1,j1,k1          : parameters for array-dimensions
!!     ib,ie,jb,je,kb,ke : range of gridpoints for which the
!!     advection has to be calculated
!!     Ru,Rp             : radial positions of the U-velocity
!!     component and the pressure location
!!     respectively
!!     
!!     on output :
!!     
!!     putout            : advection part
!!     other parameters  : all unchanged
!!     
!!********************************************************************
subroutine advecw(putout,Uvel,Wvel,RHO,Ru,Rp,dru,dz,ekm,peclet_z)
  use mod_param
  implicit none

  integer   im,ip,km,kp,ib,ie,kb,ke
  real*8     putout(0:i1,0:k1),Uvel(0:i1,0:k1),Wvel(0:i1,0:k1), &
    dru(0:i1),dz,Ru(0:i1),Rp(0:i1)
  real*8 rho(0:i1,0:k1),ekm(0:i1,0:k1),peclet_z(0:i1,0:k1)
  real*8 rhoip,rhoim,advcecw_w

  ib = 1
  ie = i1-1
  kb = 1
  ke = k1-1


  do k=kb,ke
    kp=k+1
    km=k-1
    do i=ib,ie
      ip=i+1
      im=i-1

      rhoip=0.25*(rho(i,k)+rho(i,kp)+rho(ip,k)+rho(ip,kp))
      rhoim=0.25*(rho(i,k)+rho(i,kp)+rho(im,k)+rho(im,kp))

      peclet_z(i,k)=Wvel(i,k)*rho(i,k)*dz/(ekm(i,k))
      if (peclet_z(i,k).gt.2.)then
        advcecw_w= 2.0*(rho(i,k)+rho(i,k))*Wvel(i,k)*(Wvel(i,k)-Wvel(i,km))/ dz+ &
          2.0*Wvel(i,k)*((rho(i,k)+rho(i,k))*Wvel(i,k)-(rho(i,km)+rho(i,km))*Wvel(i,km))/ dz
      else
        advcecw_w= ((Wvel(i,k) +Wvel(i,kp) )*(Wvel(i,k)+Wvel(i,kp))*rho(i,kp) - &
          (Wvel(i,km)+Wvel(i,k)  )*(Wvel(i,k)+Wvel(i,km))*rho(i,k )  ) / dz
      endif

      putout(i,k) = 0.0
      putout(i,k) = - 0.25 * (  &
        (Ru(i )*(Uvel(i,k) +Uvel(i,kp) )*(Wvel(i,k)+Wvel(ip,k)) &
        *rhoip - &
        Ru(im)*(Uvel(im,k)+Uvel(im,kp))*(Wvel(i,k)+Wvel(im,k)) &
        *rhoim ) &
        / ( Rp(i) * dru(i) ) +advcecw_w)
    enddo
  enddo
  return
end

