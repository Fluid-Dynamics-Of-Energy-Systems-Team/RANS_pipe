
!>    ****************************************************
!!    Main file of the code
!!    This code simulates Super Critical Fluid in a
!!    heated pipe with constant heat flux
!!    ****************************************************

use mod_param
use mod_common
use mod_eosmodels
use mod_common2
implicit none

include      'mpif.h'             !> mpi stuff

integer      rank,ierr,istart,noutput
integer      iTabFoundHi,iTabFoundLo
real*8       bulk,stress,stime,time1,time2,timer,time3,dif,adv
real*8       newTemp,newRho,newMu,newLam,newCP,newEnth,newTemp2,enth_i1,enth_imax,fxvalue,str,str_tot
real*8       resC,resK,resE,resV2,resOm,resSA   ! Residuals for energy, kine, eps, v2, omega, nuSA
real(8), dimension(0:i1) :: Win,kin,ein,ekmtin,v2in,omIn,nuSAin
real*8       tempWall
real :: start, finish

call cpu_time(time1)
call mpi_init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)
call mpi_comm_size(MPI_COMM_WORLD,px,ierr)

kmax    = 384/px
kmaxper = kmax*px/2
k1      = kmax + 1
k1old   = k1
Mt=imax/px
Nx=kmax*px
Mx=kmax
Nt=imax

call initMem()


!initialize EOS
if (EOSmode.eq.0) then
  allocate(eos_model,    source=IG_EOSModel(Re,Pr))
else
  allocate(eos_model,    source=Table_EOSModel(Re,Pr,2000, 'co2h_table.dat'))
endif
call eos_model%init()






call init_transpose
call mkgrid(rank)

dtmax = 1.e-3
dt = dtmax
       
istart = 1

if (select_init < 2) then
  call fkdat(rank)
  istart=1
else
  call loadRestart(istart,rank)
  istart = istart+1
endif

! periodic = 1, turb flow generator
! periodic = 2, heated pipe
if (periodic.ne.1) then
  if (turbmod.eq.0) open(29,file =  '0/Inflow',form='unformatted')
  if (turbmod.eq.1) open(29,file = 'SA/Inflow',form='unformatted')
  if (turbmod.eq.2) open(29,file = 'MK/Inflow',form='unformatted')
  if (turbmod.eq.3) open(29,file = 'VF/Inflow',form='unformatted')
  if (turbmod.eq.4) open(29,file = 'OM/Inflow',form='unformatted')
  read(29) Win(:),kin(:),ein(:),v2in(:),omIn(:),nuSAin(:),ekmtin(:),Pk(:,0)
  close(29)
endif

call state_upd(cnew,rnew,ekm,ekh,temp,beta,istart,rank);

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
tempWall = 1.0
if (isothermalBC.eq.1) then
  tempWall = min(max(tempWall, (temp(imax,kmax)+temp(i1,kmax))*0.5),Tw)
  call funcIsothermalEnthBC_upd(tempWall) ! calc. enthalpy at the wall (isothermal BC)
  if (Qwall.ne.0) then
    if (rank.eq.0) print '("Isothermal BC, Qwall should be 0  but it is ",f6.3,"... stopping")',Qwall
    stop
  else
    if (rank.eq.0) print*,"*************SOLVING AN ISOTHERMAL WALL*************!"
  endif
  if (rank.eq.0) print '("temperature at the wall = ",f6.3," .")',tempWall
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ss
      
call bound_h_upd(kin,ein,v2in,omIn,nuSAin,rank)
call state_upd(cnew,rnew,ekm,ekh,temp,beta,istart,rank);! necessary to call it twice
rold = rnew


call turbprop(Unew,Wnew,ekme,ekmt,ekmtin,rank,istep)
call bound_v(Unew,Wnew,Win,rank)
      


call chkdt(rank,istep)

call cpu_time(start)
! simulation loop
do istep=istart,nstep


  ! calculating turbulent viscosity
  call turbprop(Unew,Wnew,ekme,ekmt,ekmtin,rank,istep)
         
  if (turbmod.eq.3)  then
    call fillhm(rank)
    call SOLVEhelm(fv2,Ru,Rp,dRu,dRp,dz,rank,Lh,centerBC)
  endif
  
  call advanceScalar(resC,resK,resE,resV2,resOm,resSA,Unew,Wnew,Rnew,fv2,rank)
         
  call bound_h_upd(kin,ein,v2in,omIn,nuSAin,rank)
  call state_upd(cnew,rnew,ekm,ekh,temp,beta,istep,rank)
  
  call advance(rank)
  call bound_m(dUdt,dWdt,wnew,rnew,Win,rank)
  call fillps(rank)

  call SOLVEpois(p,Ru,Rp,dRu,dRp,dz,rank,centerBC)
  call correc(rank,1)
  call bound_v(Unew,Wnew,Win,rank)

  !ramping isothermal wall temperature
  if (mod(istep,2000).eq.0) then
    if (isothermalBC.eq.1 .AND.  tempWall.lt.Tw) then
      tempWall = min(tempWall+dTwall, Tw)
      call funcIsothermalEnthBC_upd(tempWall) ! calc. enthalpy at the wall (isothermal BC)
      if (rank.eq.0) print '("temperature at the wall ramped! = ",f6.3," .")',tempWall
    endif
  endif

  if (mod(istep,10) .eq. 0)      call chkdiv(rank)

  call cmpinf(bulk,stress)
  call chkdt(rank,istep)

  if (mod(istep,1000 ).eq.0)     call outputProfile(rank)
  if (mod(istep,1000 ).eq.0)     call outputX_h_upd(rank,istep)
  if (mod(istep,1000 ).eq.0)     call output2d_upd(rank,istep)
  if (mod(istep,5000 ).eq.0)     call saveRestart(rank)

  if ((periodic.eq.1) .and. mod(istep,5000).eq.0) then
    call Inflow_output(rank,istep)
  endif

  noutput = 100
  if (rank.eq.0) then
    if (istep.eq.istart .or. mod(istep,noutput*20).eq.0) then
      write(6,'(A7,9A14)') 'istep','dt','bulk', &
        'stress','cResid','kineResid','epsResid','v2Resid','omResid','nuSAresid'
    endif
    if (istep.eq.istart .or. mod(istep,noutput).eq.0) then
      write(6,'(i7,9e14.5)') istep,dt,bulk,stress,resC,resK,resE,resV2,resOm,resSA
    endif
  end if
         
enddo
call cpu_time(finish)
print '("Time = ",f6.3," seconds.")',finish-start
call outputProfile(rank)
call output2d_upd(rank,istep)
call mpi_finalize(ierr)
stop
end



!>******************************************************************************************
!!      turbprop routine to estimate the eddy viscosity
!!******************************************************************************************
subroutine turbprop(U,W,ekmetmp,ekmttmp,ekmtin,rank,step)

  use mod_param
  use mod_common
  implicit none
      
  integer  rank,im,ip,km,kp,step
  real*8   tauwLoc, tauw(0:k1)
  real*8, dimension(0:i1,0:k1) :: U,W,ekmetmp,ekmttmp
  real*8, dimension(0:i1) :: ekmtb,ekmtf,ekmtin


  if (turbmod.eq.0) then
    do k=1,kmax
      tauw(k) = ekmi(imax,k)*0.5*(W(imax,k-1)+W(imax,k))/wallDist(imax)
      do i=1,imax
        ekmttmp(i,k) = 0.
      enddo
    enddo
  elseif (turbmod.eq.1) then
    call calculate_mut_SA(U,W,ekmetmp,ekmttmp,ekmtin,step)
  elseif (turbmod.eq.2) then
    call calculate_mut_MK(U,W,ekmetmp,ekmttmp,ekmtin,step)
  elseif (turbmod.eq.3) then
    call calculate_mut_VF(U,W,ekmetmp,ekmttmp,ekmtin,step)
  elseif (turbmod.eq.4) then
    call calculate_mut_SST(U,W,ekmetmp,ekmttmp,ekmtin,step)

    ! Boundary condition for bF1 of sst model
    bF1(i1,:) =  bF1(imax,:)
    bF1(0,:)  =  bF1(1,:)

    call shiftf(bF1,ekmtf,rank)
    call shiftb(bF1,ekmtb,rank)
    bF1(:,0)  = ekmtf(:)
    bF1(:,k1) = ekmtb(:)

    if ((periodic.ne.1).and.(rank.eq.0)) then
      bF1(:,0) = bF1(:,1)  ! ATTENTION
    endif

    if ((periodic.ne.1).and.(rank.eq.px-1)) then
      bF1(:,k1) = 2.*bF1(:,kmax)-bF1(:,kmax-1)
    endif

  endif

  sigmat = 0.9

  ekmttmp(i1,:) = -ekmttmp(imax,:)
  ekmttmp(0,:)  =  ekmttmp(1,:)

  call shiftf(ekmttmp,ekmtf,rank)
  call shiftb(ekmttmp,ekmtb,rank)
  ekmttmp(:,0)  = ekmtf(:)
  ekmttmp(:,k1) = ekmtb(:)

  if ((periodic.ne.1).and.(rank.eq.0)) then
    ekmttmp(:,0) = ekmtin(:)
  endif

  if ((periodic.ne.1).and.(rank.eq.px-1)) then
    ekmttmp(:,k1) = 2.*ekmt(:,kmax)-ekmt(:,kmax-1)
  endif

  ekmetmp = ekm + ekmttmp
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
subroutine advanceScalar(resC,resK,resE,resV2,resOm,resSA,Utmp,Wtmp,Rtmp,ftmp,rank)
  use mod_param
  use mod_common
  implicit none
      
  real*8 dnew(0:i1,0:k1),tempArray(0:i1,0:k1),dimpl(0:i1,0:k1),tscl
  real*8 Utmp(0:i1,0:k1),Wtmp(0:i1,0:k1),Rtmp(0:i1,0:k1),ftmp(imax,kmax),sigmakSST(0:i1,0:k1)
  real*8 rho2(0:i1,0:k1), rho3(0:i1,0:k1), eknu(0:i1,0:k1),eknui(0:i1,0:k1),eknuk(0:i1,0:k1)
  real*8 cb3,Q,hbc
  integer rank,ierr
  real*8     a  (imax)
  real*8     b  (imax)
  real*8     c  (imax)
  real*8     rhs(imax)

  real*8 t1,t2,t3,t4,bc,scl
  real*8 term, adiffm, adiffp
  real*8 resC, resK, resE, resV2, resOm, resSA


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

  resK = 0.0; resE = 0.0;  resOm = 0.0; resSA = 0.0;  resV2 = 0.0;
  if (turbmod.eq.1) then
    call advanceScalar_SA(resSA,Utmp,Wtmp,Rtmp,rank)
  elseif (turbmod.eq.2) then
    call advanceScalar_MK(resK,resE,Utmp,Wtmp,Rtmp,ftmp,rank)
  elseif (turbmod.eq.3) then
    call advanceScalar_VF(resK,resE,resV2,Utmp,Wtmp,Rtmp,ftmp,rank)
  elseif (turbmod.eq.4) then
    call advanceScalar_SST(resK,resOm,Utmp,Wtmp,Rtmp,rank)
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
!!  bound_h equation
!!
!!*************************************************************************************
subroutine bound_h_upd(kin,ein,v2in,omin,nuSAin,rank)
  use mod_param
  use mod_common
  implicit none
     
      
      include 'mpif.h'
  integer rank
  real*8 unin,flux,Ub,BCvalue(0:k1)
  real*8 Sk(0:k1),kin(0:i1),ein(0:i1),v2in(0:i1),omin(0:i1),nuSAin(0:i1)
  real*8 tmpShift(0:i1)


  !     Radial boundary condition for enthalpy c
  if (isothermalBC.eq.1) then
    !!!!!!!!!!!! isothermal
    if (centerBC.eq.1) then
      do k=0,k1
        if ((k+rank*kmax)*dz.lt.x_start_heat) then
          cnew(i1,k) = cnew(imax,k)
        else
          cnew(i1,k) = 2.0*enth_wall - cnew(imax,k)
        endif
      enddo

    else
      if (rank.eq.0) print '("Isothermal boundary condition coded only for 1 wall.... stopping")'
      stop
    endif
     
  else
    !!!!!!!!!!!! isoflux
    do k=0,k1
      if (rank.eq.0.and.k.lt.K_start_heat) then
        cnew(i1,:) = cnew(imax,:)
      else
        call funcNewtonSolve_upd(cnew(i1,k), cnew(imax,k))
        if (centerBC.eq.-1) call funcNewtonSolve_upd(cnew(0,k), cnew(1,k))
      endif

    enddo

  endif

  !     Radial boundary condition
  if (turbmod.eq.1) then
    ! SA
    nuSAnew(i1,:) = -nuSAnew(imax,:)

    if (centerBC.eq.-1) then
      nuSAnew(0,:) = -nuSAnew(1,:)
    endif

  elseif (turbmod.eq.2) then
    ! MK
    knew(i1,:) = -knew(imax,:)
    BCvalue(:) = 2.0*ekm(imax,:)/rNew(imax,:)*knew(imax,:)/wallDist(imax)**2
    enew(i1,:) = 2.0*BCvalue(:) - enew(imax,:)

    if (centerBC.eq.-1) then
      knew(0,:)  = -knew(1,:)
      BCvalue(:) = 2.0*ekm(1,:)/rNew(1,:)*knew(1,:)/wallDist(1)**2
      enew(0,:)  = 2.0*BCvalue(:) - enew(1,:)
    endif

  elseif (turbmod.eq.3) then
    ! VF
    knew(i1,:)  = -knew(imax,:)
    v2new(i1,:) = -v2new(imax,:)
    BCvalue(:)  = 2.0*ekm(imax,:)/rNew(imax,:)*knew(imax,:)/wallDist(imax)**2
    enew(i1,:)  = 2.0*BCvalue(:) - enew(imax,:)

    if (centerBC.eq.-1) then
      knew(0,:)  = -knew(1,:)
      v2new(0,:) = -v2new(1,:)
      BCvalue(:) = 2.0*ekm(1,:)/rNew(1,:)*knew(1,:)/wallDist(1)**2
      enew(0,:)  = 2.0*BCvalue(:) - enew(1,:)
    endif

  elseif (turbmod.eq.4) then
    ! SST
    knew(i1,:)  = -knew(imax,:)
    BCvalue(:)  = 60.0/0.075*ekm(imax,:)/rNew(imax,:)/wallDist(imax)**2
    omNew(i1,:) = 2.0*BCvalue(:) - omNew(imax,:)

    if (centerBC.eq.-1) then
      knew(0,:)  = -knew(1,:)
      BCvalue(:) = 60.0/0.075*ekm(1,:)/rNew(1,:)/wallDist(1)**2
      omNew(0,:) = 2.0*BCvalue(:) - omNew(1,:)
    endif

  endif

  if (centerBC.eq.1) then
    !     center line BC
    cnew(0,:)    = cnew(1,:)
    knew(0,:)    = knew(1,:)
    enew(0,:)    = enew(1,:)
    v2new(0,:)   = v2new(1,:)
    omNew(0,:)   = omNew(1,:)
    nuSAnew(0,:) = nuSAnew(1,:)
  endif
  !     ghost cells between CPU
  call shiftf(cnew,   tmpShift,rank);       cnew(:,0) = tmpShift(:);
  call shiftf(knew,   tmpShift,rank);       knew(:,0) = tmpShift(:);
  call shiftf(enew,   tmpShift,rank);       enew(:,0) = tmpShift(:);
  call shiftf(v2new,  tmpShift,rank);      v2new(:,0) = tmpShift(:);
  call shiftf(omNew,  tmpShift,rank);      omNew(:,0) = tmpShift(:);
  call shiftf(nuSAnew,tmpShift,rank);    nuSAnew(:,0) = tmpShift(:);

  call shiftb(cnew,   tmpShift,rank);      cnew(:,k1) = tmpShift(:);
  call shiftb(knew,   tmpShift,rank);      knew(:,k1) = tmpShift(:);
  call shiftb(enew,   tmpShift,rank);      enew(:,k1) = tmpShift(:);
  call shiftb(v2new,  tmpShift,rank);     v2new(:,k1) = tmpShift(:);
  call shiftb(omNew,  tmpShift,rank);     omNew(:,k1) = tmpShift(:);
  call shiftb(nuSAnew,tmpShift,rank);   nuSAnew(:,k1) = tmpShift(:);
      
  if (periodic.eq.1) return


  !     set inlet and outlet BC for developing flow
  if (rank.eq.0) then
    cnew(:,0) = 0.0
    knew(:,0) = kin(:)
    enew(:,0) = ein(:)
    v2new(:,0) = v2in(:)
    omNew(:,0) = omin(:)
    nuSAnew(:,0) = nuSAin(:)
  endif

  if (rank.eq.px-1) then
    cnew(:,k1) = 2.0*   cnew(:,kmax) -    cnew(:,kmax-1)
    knew(:,k1) = 2.0*   knew(:,kmax) -    knew(:,kmax-1)
    enew(:,k1) = 2.0*   enew(:,kmax) -    enew(:,kmax-1)
    v2new(:,k1) = 2.0*  v2new(:,kmax) -   v2new(:,kmax-1)
    omNew(:,k1) = 2.0*  omNew(:,kmax) -   omNew(:,kmax-1)
    nuSAnew(:,k1) = 2.0*nuSAnew(:,kmax) - nuSAnew(:,kmax-1)
  endif

end

!>*************************************************************************************
!!      bound_v(Ubound,Wbound,Win,rank)
!!
!!
!!*************************************************************************************
subroutine bound_v(Ubound,Wbound,Win,rank)

  use mod_param
  use mod_common
  implicit none
    
      include 'mpif.h'
  character*5 inflow
  integer rank,ierr,tabkhi,tabklo
  real*8  y1,y2,y3,y4
  real*8  Ubound(0:i1,0:k1), Wbound(0:i1,0:k1),flux,Ub,Win(0:i1)
  real*8 Rbb(0:i1)
  real*8 ubb(0:i1)
  real*8 wbb(0:i1)
  real*8 Rbf(0:i1)
  real*8 ubf(0:i1)
  real*8 wbf(0:i1)
  integer ib,ie,kb,ke



  !     Radial Boundary condition
  if (centerBC.eq.-1) then ! channal bc
    do k=0,k1
      Ubound(1,k)    =   0.0
      Ubound(0,k)    = - Ubound(2,k)
      Ubound(imax,k) =   0.0
      Ubound(i1,k)   = - Ubound(imax-1,k)

      Wbound(0,k)   = - Wbound(1,k)
      Wbound(i1,k)  = - Wbound(imax,k)
    enddo
  else
    do k=0,k1
      Ubound(0,k)    =   Ubound(1,k)
      Ubound(imax,k) =   0.0
      Ubound(i1,k)   = - Ubound(imax-1,k)

      Wbound(0,k)   =   Wbound(1,k)
      Wbound(i1,k)  = - Wbound(imax,k)
    enddo
  endif

  call shiftf(Ubound,ubf,rank);     Ubound(:,0)  = Ubf(:);
  call shiftf(Wbound,wbf,rank);     Wbound(:,0)  = Wbf(:);
  call shiftb(Ubound,ubb,rank);     Ubound(:,k1) = Ubb(:);
  call shiftb(Wbound,wbb,rank);     Wbound(:,k1) = Wbb(:);


  !     if periodic is true, no need to overwrite k=0 proofile on rank=0 and to
  !     apply advective outflow BC
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
  real*8      W_out(0:i1,0:k1)
  integer rank,ierr
  real*8  y1,y2,y3,y4
  real*8  Ubound(0:i1,0:k1),Vbound(0:i1,0:k1), &
    Wbound(0:i1,0:k1),Rbound(0:i1,0:k1),Win(0:i1)
  real*8 Ub,flux,flux_tot,deltaW,rhob,wfunc,wr(1:imax)
  real*8 Rbb(0:i1)
  real*8 ubb(0:i1)
  real*8 vbb(0:i1)
  real*8 wbb(0:i1)
  real*8 Rbf(0:i1)
  real*8 ubf(0:i1)
  real*8 vbf(0:i1)
  real*8 wbf(0:i1)
  integer   ib,ie,kb,ke

  if (centerBC.eq.-1) then ! channal bc
    do k=0,k1
      Ubound(1,k)    =   0.0
      Ubound(0,k)    = - Ubound(2,k)
      Ubound(imax,k) =   0.0
      Ubound(i1,k)   = - Ubound(imax-1,k)

      Wbound(0,k)   = - Wbound(1,k)
      Wbound(i1,k)  = - Wbound(imax,k)
    enddo
  else
    do k=0,k1
      Ubound(0,k)    =   Ubound(1,k)
      Ubound(imax,k) =   0.0
      Ubound(i1,k)   = - Ubound(imax-1,k)

      Wbound(0,k)   =   Wbound(1,k)
      Wbound(i1,k)  = - Wbound(imax,k)
    enddo
  endif


  call shiftf(Ubound,ubf,rank)
  call shiftf(Wbound,wbf,rank)
  call shiftb(Ubound,ubb,rank)
  call shiftb(Wbound,wbb,rank)
  !!!!!!!!!!!!!!!!!!!!!!!

  do i=0,i1
    Ubound(i,k1) = Ubb(i)
    Wbound(i,k1) = Wbb(i)
    Ubound(i,0)  = Ubf(i)
    Wbound(i,0)  = Wbf(i)
  enddo




  !     if periodic is true, no need to overwrite k=0 proofile on rank=0 and to
  !     apply advective outflow BC
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
      if (turbmod.eq.0) open(29,file=  '0/Inflow',form='unformatted')
      if (turbmod.eq.1) open(29,file= 'SA/Inflow',form='unformatted')
      if (turbmod.eq.2) open(29,file= 'MK/Inflow',form='unformatted')
      if (turbmod.eq.3) open(29,file= 'VF/Inflow',form='unformatted')
      if (turbmod.eq.4) open(29,file= 'OM/Inflow',form='unformatted')
      read(29) Wnew(:,k),knew(:,k),enew(:,k),v2new(:,k),omNew(:,k),nuSAnew(:,k),ekmt(:,k),Pk(:,k)
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

      knew(i,:)  = 0.1
      enew(i,:)  = 1.0
      omnew(i,:) = 0.001
      v2new(i,:) = 2./3.*knew(i,:)
      nuSAnew(i,:) = 0.001
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








