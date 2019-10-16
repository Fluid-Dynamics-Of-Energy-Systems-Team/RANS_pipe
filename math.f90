SUBROUTINE SOLVEpois(rhs,Ru,Rp,dRu,dRp,dz,rank,centerBC)
  !
  !     FAST POISSON SOLVER IN CYLINDRICAL COORDINATES
  !     BENDIKS JAN BOERSMA
  !     LABORATORIUM VOOR AERO EN HYDRODYNAMICA
  !     ROTTERDAMSEWEG 145
  !     2628 AL DELFT
  !     email ::::   b.j.boersma@wbmt.tudelft.nl
  !
 
  use mod_param
  implicit none
        
  real*8      RHS(IMAX,KMAX),Ru(0:IMAX+1),Rp(0:IMAX+1)
  real*8      dz,dzi,pi,d(IMAX,kmax),bbb,z,dru(0:IMAX+1),drp(0:IMAX+1)
  real*8      a(imax),b(imax),c(imax)
  real*8      zrt(kmax*px)
  real*8      vfftk(imax/px,kmax*px)
  real*8      wk(4*px*kmax+15),bb(imax),rtmp(imax/px,kmax*px)
  integer     rank, centerBC

  !     generate tridiagonal systems

  pi = 4.*atan(1.)
  if (periodic.eq.1) then
    call vrffti(kmax*px,wk)
  else
    call vcosqi(kmax*px,wk)
  endif
  call t2fp(rhs,rtmp,rank)
  !     K --> direction
  do i=1,imax/px
    do k=1,kmax*px
      vfftk(i,k)=rtmp(i,k)
    enddo
  enddo
  if (periodic.eq.1) then
    call vrfftf(imax/px,kmax*px,vfftk,rtmp,imax/px,wk) ! RENE: commented Sept_21_2011
  else
    call vcosqb(imax/px,kmax*px,vfftk,rtmp,imax/px,wk)
  endif
  do i=1,imax/px
    do k=1,kmax*px
      rtmp(i,k)=vfftk(i,k)
    enddo
  enddo
  call t2np(rtmp,rhs,rank)


  do i=1,imax
    a(i)= Ru(I-1)/(dRp(I-1)*Rp(I)*dRu(I))
    b(i)=-(Ru(I)/(dRp(I))+Ru(I-1)/dRp(I-1))/ &
      (Rp(I)*dRu(I))
    c(i)= Ru(I) /(dRp(I)*Rp(I)*dRu(I))
  enddo

  if (centerBC.eq.-1) then
    b(1)    = b(1)+a(1)
  else
    b(1)    =-Ru(1) /(dRp(1)*Rp(1)*dRu(1))
  endif
  b(imax) = b(imax)+c(imax)

  c(imax)=0.
  a(1)=0.


  !     Generate Eigenvalues

  !     K --> direction      (zrt)
  dzi = 1./dz

  !     RENE: added Sept21_2011
  !     RENE: commented Sept21_2011
  if (periodic.eq.1) then
    zrt(1)=0.
    do k=3,kmax*px,2
      zrt(k-1)=-4.*dzi*dzi*(sin(float((k-1))*pi/(2.*kmax*px)))**2
      zrt(k)=zrt(k-1)
    enddo
    zrt(kmax*px)=-4.*dzi*dzi
  else
    do k=1,kmax*px
      zrt(k)=-4.*dzi*dzi*(sin(float((k-1))*pi/(2.*kmax*px)))**2
    enddo
  endif



  !     solve tridiagonal systems with Gaussian elemination
  do k=1,kmax
    bbb        = b(1)+zrt(k+rank*kmax)
    z          = 1./bbb
    d(1,k)   = c(1)*z
    rhs(1,k) = rhs(1,k)*z
  enddo
  do k=1,kmax
    do i=2,imax-1
      bb(i)      = b(i)+zrt(k+rank*kmax)
      z          = 1./(bb(i)-a(i)*d(i-1,k))
      d(i,k)   = c(i)*z
      rhs(i,k) = (rhs(i,k)-a(i)*rhs(i-1,k))*z
    enddo
  enddo
  do k=1,kmax
    z = b(imax)+zrt(k+rank*kmax)-a(imax)*d(imax-1,k)
    if (z .ne. 0.) then
      rhs(imax,k) = (rhs(imax,k)-a(imax)*rhs(imax-1,k))/z
    else
      rhs(imax,k) = 0.
    endif
  enddo
  do k=1,kmax
    do  i=imax-1,1,-1
      rhs(i,k) = rhs(i,k)-d(i,k)*rhs(i+1,k)
    enddo
  enddo

  call t2fp(rhs,rtmp,rank)
  !     BACKWARD FFT ---> K direction
  do i=1,imax/px
    do k=1,kmax*px
      vfftk(i,k)=rtmp(i,k)
    enddo
  enddo
  if (periodic.eq.1) then
    call vrfftb(imax/px,kmax*px,vfftk,rtmp,imax/px,wk) ! RENE: commented Sept_21_2011
  else
    call vcosqf(imax/px,kmax*px,vfftk,rtmp,imax/px,wk)
  endif

  do i=1,imax/px
    do k=1,kmax*px
      rtmp(i,k)=vfftk(i,k)
    enddo
  enddo

  call t2np(rtmp,rhs,rank)
  return
end



SUBROUTINE SOLVEhelm(rhs,Ru,Rp,dRu,dRp,dz,rank,hterm,centerBC)
  !
  !     FAST POISSON SOLVER IN CYLINDRICAL COORDINATES
  !     BENDIKS JAN BOERSMA
  !     LABORATORIUM VOOR AERO EN HYDRODYNAMICA
  !     ROTTERDAMSEWEG 145
  !     2628 AL DELFT
  !     email ::::   b.j.boersma@wbmt.tudelft.nl
  !
  use mod_param
  implicit none
  real*8      RHS(IMAX,KMAX),Ru(0:IMAX+1),Rp(0:IMAX+1),hterm(IMAX,KMAX)
  real*8      dz,dzi,pi,d(IMAX,kmax),bbb,z,dru(0:IMAX+1),drp(0:IMAX+1)
  real*8      a(imax),b(imax),c(imax)
  real*8      zrt(kmax*px)
  real*8      vfftk(imax/px,kmax*px)
  real*8      wk(4*px*kmax+15),bb(imax),rtmp(imax/px,kmax*px)
  integer     rank, centerBC
  !     generate tridiagonal systems

  pi = 4.*atan(1.)
  if (periodic.eq.1) then
    call vrffti(kmax*px,wk)
  else
    call vcosqi(kmax*px,wk)
  endif
  call t2fp(rhs,rtmp,rank)
  !     K --> direction
  do i=1,imax/px
    do k=1,kmax*px
      vfftk(i,k)=rtmp(i,k)
    enddo
  enddo
  if (periodic.eq.1) then
    call vrfftf(imax/px,kmax*px,vfftk,rtmp,imax/px,wk) ! RENE: commented Sept_21_2011
  else
    call vcosqb(imax/px,kmax*px,vfftk,rtmp,imax/px,wk)
  endif
  do i=1,imax/px
    do k=1,kmax*px
      rtmp(i,k)=vfftk(i,k)
    enddo
  enddo
  call t2np(rtmp,rhs,rank)

  do i=1,imax
    a(i)= Ru(I-1)/(dRp(I-1)*Rp(I)*dRu(I))     ! new
    b(i)=-(Ru(I)/(dRp(I))+Ru(I-1)/dRp(I-1))/  &! new
      (Rp(I)*dRu(I))
    c(i)= Ru(I) /(dRp(I)*Rp(I)*dRu(I))        ! new
  enddo
  if (centerBC.eq.-1) then
    b(1)    = b(1)-a(1)
  else
    b(1)    =-Ru(1) /(dRp(1)*Rp(1)*dRu(1))
  endif
  b(imax) = b(imax)-c(imax)
  c(imax)=0.
  a(1)=0.


  !     Generate Eigenvalues

  !     K --> direction      (zrt)
  dzi = 1./dz

  !     RENE: added Sept21_2011
  ! RENE: commented Sept21_2011
  if (periodic.eq.1) then
    zrt(1)=0.
    do k=3,kmax*px,2
      zrt(k-1)=-4.*dzi*dzi*(sin(float((k-1))*pi/(2.*kmax*px)))**2
      zrt(k)=zrt(k-1)
    enddo
    zrt(kmax*px)=-4.*dzi*dzi
  else
    do k=1,kmax*px
      zrt(k)=-4.*dzi*dzi*(sin(float((k-1))*pi/(2.*kmax*px)))**2
    enddo
  endif

  !     solve tridiagonal systems with Gaussian elemination

  do k=1,kmax
    bbb        = b(1)+zrt(k+rank*kmax)-1./hterm(1,k)**2.
    z          = 1./bbb
    d(1,k)   = c(1)*z
    rhs(1,k) = rhs(1,k)*z
  enddo
  do k=1,kmax
    do i=2,imax-1
      bb(i)      = b(i)+zrt(k+rank*kmax)-1./hterm(i,k)**2.
      z          = 1./(bb(i)-a(i)*d(i-1,k))
      d(i,k)   = c(i)*z
      rhs(i,k) = (rhs(i,k)-a(i)*rhs(i-1,k))*z
    enddo
  enddo
  do k=1,kmax
    z = b(imax)+zrt(k+rank*kmax)-a(imax)*d(imax-1,k)-1./hterm(imax,k)**2.
    if (z .ne. 0.) then
      rhs(imax,k) = (rhs(imax,k)-a(imax)*rhs(imax-1,k))/z
    else
      rhs(imax,k) = 0.
    endif
  enddo
  do k=1,kmax
    do  i=imax-1,1,-1
      rhs(i,k) = rhs(i,k)-d(i,k)*rhs(i+1,k)
    enddo
  enddo


  call t2fp(rhs,rtmp,rank)
  !     BACKWARD FFT ---> K direction
  do i=1,imax/px
    do k=1,kmax*px
      vfftk(i,k)=rtmp(i,k)
    enddo
  enddo
  if (periodic.eq.1) then
    call vrfftb(imax/px,kmax*px,vfftk,rtmp,imax/px,wk) ! RENE: commented Sept_21_2011
  else
    call vcosqf(imax/px,kmax*px,vfftk,rtmp,imax/px,wk)
  endif
  do i=1,imax/px
    do k=1,kmax*px
      rtmp(i,k)=vfftk(i,k)
    enddo
  enddo

  call t2np(rtmp,rhs,rank)
  return
end

subroutine matrixIdir(imax,A,B,C,RHS)
  implicit none
  integer imax,i
  real*8  A(imax),B(imax),C(imax),RHS(imax)
  do  i=1,imax-1
    B(i+1)    =  B  (i+1)  - (A(i+1)/B(i))  * C  (i)
    RHS(i+1)  =  RHS(i+1)  - (A(i+1)/B(i))  * RHS(i)
  enddo

  RHS(imax) = RHS(imax)/B(imax)
  do i=imax-1,1,-1
    RHS(i) = ( RHS(i) - C(i)*RHS(i+1))/ B(i)
  enddo
end




