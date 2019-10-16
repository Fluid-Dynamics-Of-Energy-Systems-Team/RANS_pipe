subroutine shiftb(UT,UP,rank)
  use mod_param
  use mod_common
  implicit none
  include 'mpif.h'
  integer rank,ileng,rankb,rankf,ierr
  integer itag,status(MPI_STATUS_SIZE),l
  real*8 ut(0:i1,0:k1)
  real*8 up(0:i1),UTMP(0:i1)
  parameter (ileng=i1+1 )
  itag = 11
  do i=0,i1
    utmp(i) =UT(i,1)
  enddo
  rankf=rank+1
  rankb=rank-1
  if(rank.eq.px-1)rankf=0
  if(rank.eq.   0)rankb=px-1
  call mpi_sendrecv(utmp ,ileng,MPI_REAL8,rankb,itag, &
    up,ileng,MPI_REAL8,rankf,itag, &
    MPI_COMM_WORLD,status,ierr)
end

subroutine shiftf(UT,UP,rank)
  use mod_param
  use mod_common
  implicit none
  include 'mpif.h'
  integer rank,ileng,rankb,rankf,ierr
  integer  itag,status(MPI_STATUS_SIZE),l
  real*8 UT(0:i1,0:k1),UP(0:i1),UTMP(0:i1)
  parameter (ileng= i1+1)
  itag = 10
  do i=0,i1
    UTMP(i) =UT(i,kmax)
  enddo
  rankf=rank+1
  rankb=rank-1
  if(rank.eq.px-1)rankf=0
  if(rank.eq.   0)rankb=px-1
  call mpi_sendrecv(utmp ,ileng,MPI_REAL8,rankf,itag, &
    up,ileng,MPI_REAL8,rankb,itag, &
    MPI_COMM_WORLD,status,ierr)
end

subroutine pshiftb_w(UT,UP,rank)
  use mod_param
  use mod_common
  implicit none
  include 'mpif.h'
  integer rank,ileng,rankb,rankf,ierr
  integer itag,status(MPI_STATUS_SIZE),l
  real*8 ut(imax,kmax)
  real*8 up(imax),UTMP(imax)
  parameter (ileng= imax )
  itag = 11
  do i=1,imax
    utmp(i) =UT(i,1)
  enddo
  rankf=rank+1
  rankb=rank-1
  if(rank.eq.px-1)rankf=0
  if(rank.eq.   0)rankb=px-1
  call mpi_sendrecv(utmp ,ileng,MPI_REAL8,rankb,itag, &
    up ,ileng,MPI_REAL8,rankf,itag, &
    MPI_COMM_WORLD,status,ierr)
end

SUBROUTINE init_transpose
  use mod_param
  use mod_common
  do k = 1,Mt
    do i = 1,Nx
      Xii(i) = MOD(i-1,Mx)+1
      Xkk(i,k) = INT((i-1)/Mx)*Mt+k
    end do
  end do
  do i = 1,Mx
    do k = 1,Nt
      Tkk(k) = MOD(k-1,Mt) + 1
      Tii(k,i) = INT((k-1)/Mt)*Mx + i
    end do
  end do
  return
end

SUBROUTINE t2fp(Q,Qt,rank)
  use mod_param
  use mod_common
  include 'mpif.h'
  integer N
  integer ii,kk,l
  real*8  Q(Nt,Mx),Qt(Mt,Nx)
  do i = 1,Mx
    do k = 1,Nt
      W1t(i,k) = Q(k,i)
    end do
  end do
  call MPI_ALLTOALL(W1t,Mt*Mx,MPI_REAL8,W2t,Mt*Mx, &
    MPI_REAL8,MPI_COMM_WORLD,ierr)
  do k = 1,Mt
    do i = 1,Nx
      Qt(k,i) = W2t(Xii(i),Xkk(i,k))
    end do
  end do
  return
end

SUBROUTINE t2np(Q,Qt,rank)
  use mod_param
  use mod_common
  include 'mpif.h'
  integer N
  integer ii,kk,l
  real*8 Q(Mt,Nx),Qt(Nt,Mx)
  do k = 1,Mt
    do i = 1,Nx
      W1(k,i) = Q(k,i)
    end do
  end do
  call MPI_ALLTOALL(W1,Mt*Mx,MPI_REAL8,W2,Mt*Mx, &
    MPI_REAL8,MPI_COMM_WORLD,ierr)
  do i = 1,Mx
    do k = 1,Nt
      Qt(k,i) = W2(Tkk(k),Tii(k,i))
    end do
  end do
  return
end
