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

subroutine solvepois_cr(rhs,ini,Ru,Rp,dRu,dRp,dz,rank,centerBC)
  use mod_param, only : kmax, imax, i1, k1, px, periodic
  use module_mesh, only : mesh
  implicit none
  include 'mpif.h'
  
  real*8      RHS(IMAX,KMAX),Ru(0:IMAX+1),Rp(0:IMAX+1)
  real*8      dz,dru(0:IMAX+1),drp(0:IMAX+1)
  integer     centerBC

  ! real(8), dimension(1:imax,1:kmax) :: rhs
  real(8), dimension(kmax*px*imax*8):: work
  real(8), dimension(1:imax)        :: am,bm,cm
  real(8), dimension(1:kmax)        :: an,bn,cn
  real(8), dimension(1:kmax*px)     :: an_t,bn_t,cn_t
  real(8), dimension(imax*kmax)     :: pvec
  real(8), dimension(imax*kmax*px)  :: pvec_t
  real(8), dimension(imax,kmax*px)  :: y
  real(8), dimension(0:k1) :: dzw, dzp
  
  integer ierr,ini,i,j,rank, ier
  character*5 cha
 
  dzp = mesh%dzp
  dzw = mesh%dzw
! write(*,*) px
  
  ! create the coefficients
  ! do i=1,imax
  !   cm(i) = 1.0/((Rp(i+1)-Rp(i))*(Ru(i)-Ru(i-1)))
  !   am(i) = 1.0/((Rp(i)-Rp(i-1))*(Ru(i)-Ru(i-1)))
  !   bm(I) =-(am(i)+cm(i))
  ! enddo
  
  !wall normal-direction
  do i=1,imax
    am(i)= Ru(I-1)/(dRp(I-1)*Rp(I)*dRu(I))
    bm(i)=-(Ru(I)/(dRp(I))+Ru(I-1)/dRp(I-1))/ &
      (Rp(I)*dRu(I))
    cm(i)= Ru(I) /(dRp(I)*Rp(I)*dRu(I))
  enddo
  if (centerBC.eq.-1) then
    bm(1)    = bm(1)+am(1)
  else
    bm(1)    =-Ru(1) /(dRp(1)*Rp(1)*dRu(1))
  endif
  am(1)=0.
  bm(imax) = bm(imax)+cm(imax)
  cm(imax)=0.

    
  !streamwise-direction  
  ! do j=1,kmax
  !   cn(j)=  1.0/((dz)*(dz ) )
  !   an(j)=  1.0/((dz)*(dz ) )
  !   bn(j)= -( an(j) + cn(j) )
  ! enddo
  
  do j=1,kmax
    an(j)=  1.0/(dzp(j-1)*dzw(j))
    cn(j)=  1.0/(dzp(j)  *dzw(j))
    bn(j)= -( an(j) + cn(j) )
  enddo


  !gather all the coefficients to 1 coordinates
  call MPI_ALLGATHER(an,   kmax,      MPI_REAL8, an_t,  kmax, MPI_REAL8, MPI_COMM_WORLD, ierr)  
  call MPI_ALLGATHER(bn,   kmax,      MPI_REAL8, bn_t,  kmax, MPI_REAL8, MPI_COMM_WORLD, ierr) 
  call MPI_ALLGATHER(cn,   kmax,      MPI_REAL8, cn_t,  kmax, MPI_REAL8, MPI_COMM_WORLD, ierr)  

  !apply bc
  bn_t(1)=bn_t(1)-cn_t(1)      
  an_t(1)=0.
  if (periodic.eq.1) then
    bn_t(kmax*px)=bn_t(kmax*px)-cn_t(kmax*px)   !for developing 
  else 
    bn_t(kmax*px)=bn_t(kmax*px)+cn_t(kmax*px)
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

  do j=1+rank*kmax,(kmax + rank*kmax)
     do i=1,imax
        rhs(i,j-rank*kmax)=y(i,j)
     enddo
  enddo  

end subroutine  solvepois_cr



