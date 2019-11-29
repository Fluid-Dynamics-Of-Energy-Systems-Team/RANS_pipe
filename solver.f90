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

subroutine solver(p,Ru,RP,thetav,thetap,dphi,ini,rank)
  use mod_param, only : kmax, imax, i1, k1
  implicit none
  include 'mpif.h'
  real(8), dimension(0:i1,0:k1) :: p
  real(8), dimension(kmax) :: an,bn,cn
  real(8), dimension(imax) :: am,bm,cm
  real(8), dimension(imax,kmax) :: y  
  real(8), dimension(3*imax*kmax) :: work
  
  !real*8 p(jmax,kmax,Mx),p2(imax,jmax,kmax/Px)
  real*8 Ru(0:i1),Rp(0:i1),dphi
  real*8 thetav(-1:k1),thetap(0:k1)
!  real*8 am(imax),bm(imax),cm(imax),y(imax,jmax), &
!         work(3*imax*jmax,kmax/Px)

!  ,pi,tpp(jmax,kmax)
  integer rank

  integer ier,it,ini,i,j,k,ipro

  do i=1,imax
    cm(i) = 1.0/((Rp(i+1)-Rp(i))*(Ru(i)-Ru(i-1)))
    am(i) = 1.0/((Rp(i)-Rp(i-1))*(Ru(i)-Ru(i-1)))
    bm(I) =-(am(i)+cm(i))
  enddo
      
  do j=1,kmax
    cn(j)=  1.0/((thetap(j+1)-thetap(j))*(thetav(j)-thetav(j-1) ) )
    an(j)=  1.0/((thetap(j)-thetap(j-1))*(thetav(j)-thetav(j-1) ) )
    bn(j)= -( an(j) + cn(j) )
  enddo

! 
! Set the boundary conditions
! 
      
  bm(1)=bm(1)+am(1)      
  am(1)=0.
  bm(imax)=bm(imax)+cm(imax)    
  cm(imax)=0.
  
  bn(1)    = bn(1)-an(1)
  an(1)    = 0.
  bn(kmax) = bn(kmax)-cn(kmax)
  cn(kmax) = 0.


  do j=1,kmax
    do i=1,imax
      y(i,j)=p(i,j)
    enddo
  enddo 
! 
!     Call cyclic reduction algorithm
!     
  if (ini.eq.0) call blktri(0,1,kmax,an,bn,cn,1,imax,am,bm,cm,imax,y,ier,work(1))
  call blktri(1,1,kmax,an,bn,cn,1,imax,am,bm,cm,imax,y,ier,work(1))
  if (ier .ne. 0) then
     write(6,*) 'There is something wrong with the solution of the Poisson equation!'
     write(6,*) 'Result are not reliable!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     stop
  endif
      
      
  do j=1,kmax
     do i=1,imax
        p(i,j)=y(i,j)
     enddo
  enddo          

end subroutine  solver



