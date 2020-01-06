subroutine make_vector(matrix_in, vector_out, imax, kmax)
  implicit none
  integer, intent(IN) :: imax, kmax
  real(8), dimension(1:imax,1:kmax), intent(IN) :: matrix_in
  real(8), dimension(1:imax*kmax), intent(OUT) :: vector_out
  integer :: index, i, k

  index = 1
  do k=1,kmax
    do i=1,imax
      vector_out(index) = matrix_in(i,k)
      index = index+1
    enddo
  enddo
end subroutine make_vector

subroutine make_matrix(vector_in, matrix_out, imax, kmax)
  implicit none
  integer, intent(IN) :: imax, kmax
    real(8), dimension(1:imax*kmax), intent(IN) :: vector_in
  real(8), dimension(1:imax,1:kmax), intent(OUT) :: matrix_out
  integer :: index,i,k
  k=1
  i=1
  do index=1,imax*kmax
    if (i > imax) then
      i = 1
      k = k+1
    endif
    matrix_out(i,k) = vector_in(index)
    i = i +1
  enddo
end subroutine make_matrix


SUBROUTINE SOLVEpois(rhs,rank,centerBC)
  !
  !     FAST POISSON SOLVER IN CYLINDRICAL COORDINATES
  !     BENDIKS JAN BOERSMA
  !     LABORATORIUM VOOR AERO EN HYDRODYNAMICA
  !     ROTTERDAMSEWEG 145
  !     2628 AL DELFT
  !     email ::::   b.j.boersma@wbmt.tudelft.nl
  !
 
  use mod_param
  use mod_mesh, only : dru,ru,rp,dz,drp
  implicit none
        
  real*8      RHS(IMAX,KMAX)
  real*8      dzi,pi,d(IMAX,kmax),bbb,z
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

!*********************************************************************
!*********************************************************************
!  This routine solves a Poisson equation on a staggered
!  grid in spherical coordinates.
!
!  METHOD:  FOURIER TRANSFORMATION IN THE PHI DIRECTION
!           CYCLIC REDUCTION IN THE RADIAL AND CIRCUMFE-
!           RENTIAL DIRECTION.
!
! AUTHOR:   Bendiks Jan Boersma
! DATE  :   11 October 1996
!
! USES  :   VFFTPACK and FISHPACK  (netlib)
!
!
!*********************************************************************
!
!
!   P(imax,jmax,kmax)  On entry: The right hand side of the poisson
!                                equation
!                      On exit : The solution of the poisson equation
!
!   Ru(0:i1),Rp(0:i1)  Arrays containing the grid locations
!                      of the Uvelocity point and the pressure
!
!   dtheta,dphi        Gridspacings
!
!
!********************************************************************
!********************************************************************
!********************************************************************

subroutine solvepois_cr(rhs,ini,rank,centerBC)
  use mod_param, only : k1,i1,kmax,imax,k,i,px,periodic
  use mod_mesh, only : dzp,dzw,dru,drp,ru,rp
  implicit none
  include 'mpif.h'  
  real*8      rhs(imax,kmax)
  integer     centerBC
  real(8), dimension(kmax*px*imax*8):: work
  real(8), dimension(1:imax)        :: am,bm,cm
  real(8), dimension(1:kmax)        :: an,bn,cn
  real(8), dimension(1:kmax*px)     :: an_t,bn_t,cn_t
  real(8), dimension(imax*kmax)     :: pvec
  real(8), dimension(imax*kmax*px)  :: pvec_t
  real(8), dimension(imax,kmax*px)  :: y
  integer ierr,ini,rank, ier
  

  !wall normal-direction
  do i=1,imax
    am(i)= Ru(I-1)/(dRp(I-1)*Rp(I)*dRu(I))
    bm(i)=-(Ru(I)/(dRp(I))+Ru(I-1)/dRp(I-1))/ &
      (Rp(I)*dRu(I))  ! twall    = 0.0
  ! massflow = 0.0
  ! enthflow = 0.0
  ! w_c      = 0.0

  ! do k=1,kmax
  !   do i=1,imax
  !     massfl = 0.5*rnew(i,k)*(Wnew(i,k)+Wnew(i,k-1))*rp(i)*dru(i)
  !     massflow(k) = massflow(k) + massfl
  !     enthflow(k) = enthflow(k) + massfl*Cnew(i,k)
  !   enddo
  ! enddo

  ! enth_b=enthflow/massflow
  ! do k=1,kmax
  !   w_c=(Cnew(i1,k)+Cnew(imax,k))/2.
  !   call eos_model%set_w_enth(w_c,      "T", Twall(k))
  !   call eos_model%set_w_enth(enth_b(k),"T", Tbulk(k))
  ! enddo

    cm(i)= Ru(I) /(dRp(I)*Rp(I)*dRu(I))
  enddo
  !apply bc for wall normal-direction
  if (centerBC.eq.-1) then
    bm(1)    = bm(1)+am(1)
  else
    bm(1)    =-Ru(1) /(dRp(1)*Rp(1)*dRu(1))
  endif
  am(1)=0.
  bm(imax) = bm(imax)+cm(imax)
  cm(imax)=0.
   
  !streamwise
  do k=1,kmax
    an(k)=  1.0/(dzp(k-1)*dzw(k))
    cn(k)=  1.0/(dzp(k)  *dzw(k))
    bn(k)= -( an(k) + cn(k) )
  enddo
  !gather all the coefficients to 1 coordinates
  call MPI_ALLGATHER(an,   kmax,      MPI_REAL8, an_t,  kmax, MPI_REAL8, MPI_COMM_WORLD, ierr)  
  call MPI_ALLGATHER(bn,   kmax,      MPI_REAL8, bn_t,  kmax, MPI_REAL8, MPI_COMM_WORLD, ierr) 
  call MPI_ALLGATHER(cn,   kmax,      MPI_REAL8, cn_t,  kmax, MPI_REAL8, MPI_COMM_WORLD, ierr)  
  !apply bc streamwise
  bn_t(1)=bn_t(1)-cn_t(1)      
  an_t(1)=0.
  if (periodic.eq.1) then
    bn_t(kmax*px)=bn_t(kmax*px)-cn_t(kmax*px) !periodic
  else 
    bn_t(kmax*px)=bn_t(kmax*px)+cn_t(kmax*px) !for developing 
  endif
  cn_t(kmax*px)=0.

  call make_vector(rhs,pvec,imax,kmax)
  call MPI_ALLGATHER(pvec, imax*kmax, MPI_REAL8, pvec_t, imax*kmax, MPI_REAL8, MPI_COMM_WORLD, ierr)
  call make_matrix(pvec_t,y,imax,kmax*px)

  !Call cyclic reduction algorithm     
  if (ini.eq.0) call blktri(0,1,kmax*px,an_t,bn_t,cn_t,1,imax,am,bm,cm,imax,y,ier,work)
  call blktri(1,1,kmax*px,an_t,bn_t,cn_t,1,imax,am,bm,cm,imax,y,ier,work)

  if (ier .ne. 0) then
     write(6,*) 'There is something wrong with the solution of the Poisson equation!'
     write(6,*) 'Result are not reliable!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     stop
  endif

  do k=1+rank*kmax,(kmax + rank*kmax)
     do i=1,imax
        rhs(i,k-rank*kmax)=y(i,k)
     enddo
  enddo  

end subroutine  solvepois_cr



