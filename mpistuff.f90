

      subroutine shiftb(UT,UP,rank)
      implicit none
      include 'param.f90'
      include 'common.f90'
      include 'mpif.h'
      integer ileng,rankb,rankf,ierr
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
      implicit none
      include 'param.f90'
      include 'common.f90'
      include 'mpif.h'
      integer ileng,rankb,rankf
      integer  itag,status(MPI_STATUS_SIZE),l,ierr
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
      implicit none
      include 'param.f90'
      include 'common.f90'
      include 'mpif.h'
      integer ileng,rankb,rankf,ierr
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
      include 'param.f90'
      parameter (mt=imax/px,nx=kmax*px,mx=kmax,NT=imax)
      integer Xii(Nx),Xkk(Nx,Mt)
      common /XPOSE/ Xii,Xkk
      integer Tkk(Nt),Tii(Nt,Mx)
      common /TPOSE/ Tkk,Tii

!     ...  Locals

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
      include 'param.f90'
      include 'mpif.h'
      parameter (mt=imax/px,nx=kmax*px,mx=kmax,NT=imax)
      integer N
      integer ii,kk,l
      integer Xii(Nx),Xkk(Nx,Mt)
      common /XPOSE/ Xii,Xkk
      real*8  Q(Nt,Mx),Qt(Mt,Nx)

!     ...  Locals

!     ...  Globals

      real*8  W1(Mx,Nt),W2(Mx,Nt)
      common  W1,W2


      do i = 1,Mx
         do k = 1,Nt
               W1(i,k) = Q(k,i)
            end do
         end do

      call MPI_ALLTOALL(W1,Mt*Mx,MPI_REAL8,W2,Mt*Mx, &
           MPI_REAL8,MPI_COMM_WORLD,ierr)

      do k = 1,Mt
            do i = 1,Nx
               Qt(k,i) = W2(Xii(i),Xkk(i,k))
            end do
      end do


      return
      end


      SUBROUTINE t2fpVector(Q,Qt,rank)
      include 'param.f90'
      include 'mpif.h'
      parameter (mt=imax/px,nx=kmax*px,mx=kmax,NT=imax)
      integer n
      integer ii,kk,l
      integer Xii(Nx),Xkk(Nx,Mt)
      common /XPOSE/ Xii,Xkk
      real*8  Q(Nt,Mx,3),Qt(Mt,Nx,3)

!     ...  Locals

!     ...  Globals

      real*8  W1(3,Mx,Nt),W2(3,Mx,Nt)
      common  W1,W2


      do i = 1,Mx
         do k = 1,Nt
              do n=1,3
               W1(n,i,k) = Q(k,i,n)
            enddo
         enddo
      enddo

      call MPI_ALLTOALL(W1,Mt*Mx*3,MPI_REAL8,W2,Mt*Mx*3, &
           MPI_REAL8,MPI_COMM_WORLD,ierr)

      do k = 1,Mt
            do i = 1,Nx
               do n=1,3
               Qt(k,i,n) = W2(n,Xii(i),Xkk(i,k))
            enddo
         enddo
      enddo


      return
      end
      SUBROUTINE t2fpTensor(Q,Qt,rank)
      include 'param.f90'
      include 'mpif.h'
      parameter (mt=imax/px,nx=kmax*px,mx=kmax,NT=imax)
      integer n,s
      integer ii,kk,l
      integer Xii(Nx),Xkk(Nx,Mt)
      common /XPOSE/ Xii,Xkk
      real*8  Q(Nt,Mx,3,3),Qt(Mt,Nx,3,3)

!     ...  Locals

!     ...  Globals

      real*8  W1(3,3,Mx,Nt),W2(3,3,Mx,Nt)
      common  W1,W2


      do i = 1,Mx
         do k = 1,Nt
              do n=1,3
                do s=1,3
               W1(n,s,i,k) = Q(k,i,n,s)
            enddo
         enddo
      enddo
      enddo

      call MPI_ALLTOALL(W1,Mt*Mx*9,MPI_REAL8,W2,Mt*Mx*9, &
           MPI_REAL8,MPI_COMM_WORLD,ierr)

      do k = 1,Mt
            do i = 1,Nx
              do n=1,3
                do s=1,3
               Qt(k,i,n,s) = W2(n,s,Xii(i),Xkk(i,k))
            enddo
         enddo
      enddo
      enddo


      return
      end
!     ... Performs data transform across the nodes.  The data
!     starts out (i,j,k) node distributed in theta and is transformed
!     to (j,k,i) node distributed in x.

      SUBROUTINE t2np(Q,Qt,rank)
      include 'param.f90'
      include 'mpif.h'
      parameter (mt=imax/px,nx=kmax*px,mx=kmax,NT=imax)
      integer N
      integer ii,kk,l
      integer Tkk(Nt),Tii(Nt,Mx)
      common /TPOSE/ Tkk,Tii
      real*8 Q(Mt,Nx),Qt(Nt,Mx)
      real*8  W1(Mt,Nx),W2(Mt,Nx)
      common  W1,W2


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
