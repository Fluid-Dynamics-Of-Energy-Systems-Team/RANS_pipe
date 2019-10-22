
!*****************************************************************!
!                                                                 !
!         Outputs the profile halfway along the domain            !
!                                                                 !
!*****************************************************************!

subroutine inflow_output_upd(rank,istap)
  use mod_param
  use mod_common
  use mod_common2
  implicit none
  
  integer rank,istap
  real(8), dimension(0:i1)   :: p_nuSA,p_k,p_eps,p_om,p_v2,p_Pk,p_bF1,p_yp
  real(8), dimension(1:imax) :: p_bF2
  character(len=50) :: fname
  character(len=5)  :: Re_str
  integer           :: Re_int

  Re_int = int(eos_model%Re)
  write(Re_str,'(I5.5)') Re_int

  fname = 'Inflow_'//trim(turb_model%name)//'_'//Re_str
  if (rank.eq.px/2) then
    k = kmax/2
    call turb_model%get_profile(p_nuSA,p_k,p_eps,p_om,p_v2,p_Pk,p_bF1,p_bF2,p_yp,k)
    !binary
    if (systemsolve.eq.1) open(29,file='pipe/'   //trim(fname)//'.dat',form='unformatted')
    if (systemsolve.eq.2) open(29,file='channel/'//trim(fname)//'.dat',form='unformatted')
    if (systemsolve.eq.3) open(29,file='bl/'     //trim(fname)//'.dat',form='unformatted')
    write(29) Wnew(:,kmax/2),p_k,p_eps,p_v2,p_om,p_nuSA,ekmt(:,kmax/2),p_pk
    close(29)
    !fixed width file
    if (systemsolve.eq.1) open(29,file='pipe/'   //trim(fname)//'.csv')
    if (systemsolve.eq.2) open(29,file='channel/'//trim(fname)//'.csv')
    if (systemsolve.eq.3) open(29,file='bl/'     //trim(fname)//'.csv')
    write(29, '(16a20)' ) 'y'   ,'u'  ,'w'  ,'h'  ,'T',  &
                          'rho' ,'k'  ,'eps','v2' ,'om', &
                          'nuSA','mut','Pk' ,'bF1','bF2', 'yp'
    do i=1,imax
      write(29, '(16E20.12)')                                           &
                y_cv  (i)  , unew(i,k),Wnew (i,k),cnew (i,k),temp (i,k), &
                rnew  (i,k),p_k (i)  ,p_eps(i)  ,p_v2 (i)  ,p_om (i),   &
                p_nuSA(i),  ekmt(i,k),p_Pk (i)  ,p_bf1(i)  ,p_bf2(i), p_yp(i)
    enddo
    close(29)
  endif
end subroutine inflow_output_upd

!*****************************************************************!
!                                                                 !
!       Writes the variables in a fixed width formatted file      !
!                                                                 !
!*****************************************************************!

subroutine write_mpiio_formatted(filename, x, y, u,w, rho,T,p,mu, mut, yp, &
                                 k, eps, v2, om,nuSA, i1, k1,rank,px)
  implicit none 
  include "mpif.h"
  character(*),                  intent(IN) :: filename
  real(8), dimension(0:i1,0:k1), intent(IN) :: x,y,u,w,rho,T,p,mu,mut,yp, &
                                               k,eps,v2,om,nuSA
  integer,                       intent(IN) :: i1,k1,rank,px
  integer nvar,i,j,index,k_max,k_min,size,fh,ierr
  integer(kind=MPI_OFFSET_KIND) disp 
  character(len=301), dimension(:), allocatable :: lines, lines2
  character(len=300) :: test
  character(len=301) :: line
  nvar = 15
  index=1

  !first core write from 0 to k1-1
  if (rank .eq. 0) then
    k_min = 0
    k_max = k1-1
    allocate(lines(1:(i1+1)*k1+1)) !+1 for header
    disp = 0
    size = ((i1+1)*(k1)+1)*(nvar*20+1)
    write(test,'(15(A20))') 'x','y','u','w','rho','T','p','mu','mut','yp','k','eps','v2','om','nuSA'
    write(line, '(A)') test // NEW_LINE("A")
    lines(index) = line
    index = index+1
  !last core write from 1 to k1
  else if (rank .eq. px-1) then
    k_min = 1
    k_max = k1
    allocate(lines(1:(i1+1)*k1))
    size =           (i1+1)*(k1)*(nvar*20+1)
   disp = ((i1+1)*(k1)+1)*(nvar*20+1) + (rank-1)*(i1+1)*(k1-1)*(nvar*20+1)
  !other core write from 1 to k1-1
  else
    k_min = 1
    k_max = k1-1
    allocate(lines(1:(i1+1)*(k1-1)))
    size =           (i1+1)*(k1-1)*(nvar*20+1)
    disp = ((i1+1)*(k1)+1)*(nvar*20+1) + (rank-1)*(i1+1)*(k1-1)*(nvar*20+1)
  endif
  do i = 0,i1
    do j = k_min,k_max
      write(test,'(15(E20.12))') x(i,j), y(i,j),u(i,j),w(i,j),rho(i,j),T(i,j),p(i,j),mu(i,j),mut(i,j),yp(i,j), &
                                k(i,j),eps(i,j),v2(i,j),om(i,j),nuSA(i,j)
      write(line, '(A)') test // NEW_LINE("A")
      lines(index) = line
      index=index+1
    enddo
  enddo
  call MPI_FILE_OPEN(MPI_COMM_WORLD, filename,MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL, fh, ierr) 
  call MPI_FILE_SET_VIEW(fh, disp, MPI_CHAR, MPI_CHAR, 'native', MPI_INFO_NULL, ierr) 
  call MPI_FILE_WRITE(fh, lines, size, MPI_CHAR,MPI_STATUS_IGNORE, ierr) 
  call MPI_FILE_CLOSE(fh, ierr) 
end subroutine write_mpiio_formatted

!*****************************************************************!
!                                                                 !
!       Outputs the variables in a fixed width formatted file     !
!                                                                 !
!*****************************************************************!

subroutine output2d_upd2(rank,istap)
  use mod_param
  use mod_common
  use mod_common2
  implicit none
  include 'mpif.h'
  character*5 cha
  real*8 pecletx,peclety,pecletz
  real*8 massflow(kmax),enthflow(kmax),enth_b(kmax),Twall(kmax),Tbulk(kmax), &
         massfl,w_c,ndt
  real(8), dimension(0:i1,0:k1) ::  xvec, yvec,nuSA_sol,k_sol,eps_sol,om_sol,v2_sol,yp_sol
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

  do k=0,k1
    do i=0,i1
      xvec(i,k)=(k+rank*kmax)*dz
      yvec(i,k) = y_cv(i)
    enddo
  enddo
  call turb_model%get_sol(nuSA_sol,k_sol,eps_sol,om_sol,v2_sol,yp_sol)
  call write_mpiio_formatted("test", xvec, yvec, unew,wnew, rnew,cnew,cp,ekm, ekmt,yp_sol,     &
                                 k_sol, eps_sol, v2_sol, om_sol,nuSA_sol, i1, k1,rank,px)
end subroutine output2d_upd2
