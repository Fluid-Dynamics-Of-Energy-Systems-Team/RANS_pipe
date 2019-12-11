subroutine write_output_bl(rank, istap)
  use mod_common, only : wnew, ekmi, unew
  use mod_param, only : k1, kmax, K_start_heat, output_fname_bl,px,i1, LoD
  use mod_mesh, only : mesh
  use mod_tdm, only : turbdiff_model
  implicit none
  include "mpif.h"
  integer, intent(IN) :: rank, istap
  real(8), dimension(0:k1) :: mom_th, dis_th, bl_th, wstress, sfriction,x
  integer(kind=MPI_OFFSET_KIND) disp 
  character(len=141), dimension(:), allocatable :: lines, lines2
  character(len=140) :: test
  character(len=141) :: line
  integer :: index, nvar, k, fh,ierr, k_max, k_min,size
  real(8), dimension(0:k1) :: zw
  ! integer K_start_heat
  zw = mesh%zw

  nvar = 7
  index=1
  call postprocess_bl(wnew, ekmi, 1., 1., mom_th, dis_th, bl_th, wstress,sfriction)

  do k=0,k1
      !@x(k)=(k+rank*kmax)*dz - K_start_heat*dz !!
      x(k)=zw(k)-mesh%start
  enddo
 
 !first core write from 0 to k1-1
  if (rank .eq. 0) then
    k_min = 0
    k_max = k1-1
    allocate(lines(1:k1+1)) !+1 for header
    disp = 0
    size = (k1+1)*(nvar*20+1)
    write(test,'(7(A20))') 'x','tau','Cf','blth','momth','disth','vinf'
    write(line, '(A)') test // NEW_LINE("A")
    lines(index) = line
    index = index+1
  !last core write from 1 to k1
  else if (rank .eq. px-1) then
    k_min = 1
    k_max = k1
    allocate(lines(1:k1))
    size =           (k1)*(nvar*20+1)
   disp = (k1+1)*(nvar*20+1) + (rank-1)*(k1-1)*(nvar*20+1)
  !other core write from 1 to k1-1
  else
    k_min = 1
    k_max = k1-1
    allocate(lines(1:(k1-1)))
    size =           (k1-1)*(nvar*20+1)
    disp = (k1+1)*(nvar*20+1) + (rank-1)*(k1-1)*(nvar*20+1)
  endif

  do k = k_min,k_max
      write(test,'(7(E20.10e3))') x(k), wstress(k),sfriction(k),bl_th(k),mom_th(k),  &
                                   dis_th(k), -unew(0,k)
      write(line, '(A)') test // NEW_LINE("A")
      lines(index) = line
      index=index+1
  enddo

  call MPI_FILE_OPEN(MPI_COMM_WORLD, output_fname_bl,MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL, fh, ierr) 
  call MPI_FILE_SET_VIEW(fh, disp, MPI_CHAR, MPI_CHAR, 'native', MPI_INFO_NULL, ierr) 
  call MPI_FILE_WRITE(fh, lines, size, MPI_CHAR,MPI_STATUS_IGNORE, ierr) 
  call MPI_FILE_CLOSE(fh, ierr) 

  call write_vector(turbdiff_model%Pkt, i1, k1, rank)

end subroutine


subroutine inflow_output_upd(rank,istap)
  use mod_param,   only : kmax,i1,imax,px,k,i,systemsolve
  use mod_tm,      only : turb_model
  use mod_tdm,     only : turbdiff_model
  use mod_eos,     only : eos_model
  use mod_common,  only : wnew, unew, cnew, temp, rnew, ekmt, alphat
  use mod_mesh, only : mesh
  implicit none
  integer, intent(IN) :: rank,istap
  real(8), dimension(0:i1)   :: p_nuSA,p_k,p_eps,p_om,p_v2,p_Pk,p_bF1,p_yp, &
                                p_prt,p_kt,p_epst,p_pkt
  real(8), dimension(1:imax) :: p_bF2
  character(len=100) :: fname
  character(len=5)  :: Re_str
  integer           :: Re_int

  Re_int = int(eos_model%Re)
  write(Re_str,'(I5.5)') Re_int

  fname = 'Inflow_'//trim(turb_model%name)//'_'//trim(turbdiff_model%name)//'_'//Re_str

  if (rank.eq.px/2) then
    k = kmax/2
    call turb_model%get_profile(p_nuSA,p_k,p_eps,p_om,p_v2,p_Pk,p_bF1,p_bF2,p_yp,k)
    call turbdiff_model%get_profile(p_prt,p_kt,p_epst,p_pkt,k)
    
    !binary
    if (systemsolve.eq.1) open(29,file='pipe/'   //trim(fname)//'.dat',form='unformatted')
    if (systemsolve.eq.2) open(29,file='channel/'//trim(fname)//'.dat',form='unformatted')
    if (systemsolve.eq.3) open(29,file='symchan/'//trim(fname)//'.dat',form='unformatted')
    write(29) Wnew(:,kmax/2),p_k,p_eps,p_v2,p_om,p_nuSA,ekmt(:,kmax/2),p_pk, &
              alphat(:,kmax/2), p_prt, p_kt, p_epst, p_pkt
    close(29)
    !fixed width file
    if (systemsolve.eq.1) open(29,file='pipe/'   //trim(fname)//'.csv')
    if (systemsolve.eq.2) open(29,file='channel/'//trim(fname)//'.csv')
    if (systemsolve.eq.3) open(29,file='symchan/'//trim(fname)//'.csv')
    write(29, '(21a20)' ) 'y'   ,'u'  ,'w'  ,'h'  ,'T',  &
                          'rho' ,'k'  ,'eps','v2' ,'om', &
                          'nuSA','mut','Pk' ,'bF1','bF2', 'yp', &
                          'alphat','prt','kt','epst','pkt'
    do i=1,imax
      write(29, '(21E20.12)')                                            &
                mesh%y_cv(i),unew(i,k),Wnew (i,k),cnew (i,k),temp (i,k),  &
                rnew  (i,k) ,p_k (i)  ,p_eps(i)  ,p_v2 (i)  ,p_om (i),    &
                p_nuSA(i)   ,ekmt(i,k),p_Pk (i)  ,p_bf1(i)  ,p_bf2(i), p_yp(i), &
                alphat(i,k), p_prt(i), p_kt(i), p_epst(i), p_pkt(i)
    enddo
    close(29)
  endif

end

subroutine read_mpiio_formatted(filename, x, y, u,w, rho,T,p,mu, mut, yp, &
                                 k, eps, v2, om,nuSA,alphat, i1, k1,rank,px)
  use mod_param, only : LoD
  implicit none 
  include "mpif.h"
  character(*),                  intent(IN) :: filename
  real(8), dimension(0:i1,0:k1), intent(OUT) :: x,y,u,w,rho,T,p,mu,mut,yp, &
                                               k,eps,v2,om,nuSA, &
                                               alphat
  integer,                       intent(IN) :: i1,k1,rank,px
  integer nvar,i,j,index,k_max,k_min,size,fh,ierr
  integer(kind=MPI_OFFSET_KIND) disp 
  character(len=321), dimension(:), allocatable :: lines, lines2
  character(len=320) :: test
  character(len=321) :: line
  real(8) :: x_v, y_v
  nvar =16
  index=1

  !first core write from 0 to kmax
  if (rank .eq. 0) then
    k_min = 0
    k_max = k1-1
    allocate(lines(1:(i1+1)*k1+1)) !+1 for header
    disp = 0
    size = ((i1+1)*(k1)+1)*(nvar*20+1)
  !last core write from 1 to k1
  else if (rank .eq. px-1) then
    k_min = 1
    k_max = k1
    allocate(lines(1:(i1+1)*k1))
    size =           (i1+1)*(k1)*(nvar*20+1)
   disp = ((i1+1)*(k1)+1)*(nvar*20+1) + (rank-1)*(i1+1)*(k1-1)*(nvar*20+1)
  !other core write from 1 to kmax
  else
    k_min = 1
    k_max = k1-1
    allocate(lines(1:(i1+1)*(k1-1)))
    size =           (i1+1)*(k1-1)*(nvar*20+1)
    disp = ((i1+1)*(k1)+1)*(nvar*20+1) + (rank-1)*(i1+1)*(k1-1)*(nvar*20+1)
  endif

  call MPI_FILE_OPEN(MPI_COMM_WORLD, filename,MPI_MODE_RDONLY,MPI_INFO_NULL, fh, ierr) 
  call MPI_FILE_SET_VIEW(fh, disp, MPI_CHAR, MPI_CHAR, 'native', MPI_INFO_NULL, ierr) 
  call MPI_FILE_READ(fh, lines, size, MPI_CHAR,MPI_STATUS_IGNORE, ierr) 
  call MPI_FILE_CLOSE(fh, ierr)
  
  if (rank .eq. 0) index=2 !skip header
  do j = k_min,k_max
    do i = 0,i1
      read(lines(index)(1  :20) ,*) x   (i,j)
      read(lines(index)(21 :40) ,*) y   (i,j)
      read(lines(index)(41 :60) ,*) u   (i,j)
      read(lines(index)(61 :80) ,*) w   (i,j)
      read(lines(index)(81 :100),*) rho (i,j)
      read(lines(index)(101:120),*) T   (i,j)
      read(lines(index)(121:140),*) p   (i,j)
      read(lines(index)(141:160),*) mu  (i,j)
      read(lines(index)(161:180),*) mut (i,j)
      read(lines(index)(181:200),*) yp  (i,j)
      read(lines(index)(201:220),*) k   (i,j)
      read(lines(index)(221:240),*) eps (i,j)
      read(lines(index)(241:260),*) v2  (i,j)
      read(lines(index)(261:280),*) om  (i,j)
      read(lines(index)(281:300),*) nuSA(i,j)      
      read(lines(index)(301:320),*) alphat(i,j)      
      index=index+1
    enddo
  enddo
end subroutine read_mpiio_formatted

subroutine write_mpiio_formatted(filename, x, y, u,w, rho,T,p,mu, mut, yp, &
                                 k, eps, v2, om,nuSA,alphat,Prt,kt,epst, i1, k1,rank,px)
  use mod_param, only : LoD
  implicit none 
  include "mpif.h"
  character(*),                  intent(IN) :: filename
  real(8), dimension(0:i1,0:k1), intent(IN) :: x,y,u,w,rho,T,p,mu,mut,yp, &
                                               k,eps,v2,om,nuSA, &
                                               alphat,Prt,kt,epst
  integer,                       intent(IN) :: i1,k1,rank,px
  integer nvar,i,j,index,k_max,k_min,size,fh,ierr
  integer(kind=MPI_OFFSET_KIND) disp 
  character(len=381), dimension(:), allocatable :: lines, lines2
  character(len=380) :: test
  character(len=381) :: line
  real(8) :: x_v, y_v
  nvar = 19
  index=1

  !first core write from 0 to kmax
  if (rank .eq. 0) then
    k_min = 0
    k_max = k1-1
    allocate(lines(1:(i1+1)*k1+1)) !+1 for header
    disp = 0
    size = ((i1+1)*(k1)+1)*(nvar*20+1)
    write(test,'(19(A20))') 'x','y','u','w','rho','T','p','mu','mut','yp','k','eps','v2','om','nuSA','alphat',"Prt","kt","epst" !write the header
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
  !other core write from 1 to kmax
  else
    k_min = 1
    k_max = k1-1
    allocate(lines(1:(i1+1)*(k1-1)))
    size =           (i1+1)*(k1-1)*(nvar*20+1)
    disp = ((i1+1)*(k1)+1)*(nvar*20+1) + (rank-1)*(i1+1)*(k1-1)*(nvar*20+1)
  endif

  do j = k_min,k_max      
    do i = 0,i1
      y_v = min(1., max(0.,y(i,j))) !take the zero in case negative, take the one in case bigger than 1
      x_v = min(LoD,max(0.,x(i,j))) !takes the LOD or zero in case the x is negative
      write(test,'(19(E20.10e3))') x_v, y_v,u(i,j),w(i,j),rho(i,j),T(i,j), p(i,j), mu(i,j), mut(i,j),  &
                                   yp(i,j),k(i,j), eps(i,j), v2(i,j), om(i,j), nuSA(i,j), alphat(i,j), &
                                   Prt(i,j), kt(i,j), epst(i,j)
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

subroutine set_scalar_to_coords(vector, i1, k1, output)
  implicit None
  integer, intent(IN) :: i1, k1
  real(8), dimension(0:i1,0:k1), intent(IN) :: vector
  real(8), dimension(0:i1,0:k1), intent(OUT) :: output
  real(8), dimension(0:i1) :: tmp
  integer :: i,k

  do k=0,k1
    if (k .eq. 0) then
      tmp = (vector(:,0)+ vector(:,1))/2.
    else if (k .eq. k1) then
      tmp = (vector(:,k1-1)+ vector(:,k1))/2.
    else
      tmp = vector(:,k)
    endif
    do i=0,i1
      if (i .eq. 0) then
        output(i,k) = (tmp(i)+tmp(i+1))/2.
      else if (i .eq. i1) then
        output(i,k) = (tmp(i1-1)+tmp(i1))/2.
      else
        output(i,k) = tmp(i)
      endif
    enddo
  enddo
end subroutine set_scalar_to_coords

subroutine set_wvector_to_coords(vector, i1, k1, output)
implicit None
  integer, intent(IN) :: i1, k1
  real(8), dimension(0:i1,0:k1), intent(IN) :: vector
  real(8), dimension(0:i1,0:k1), intent(OUT) :: output
  real(8), dimension(0:i1) :: tmp
  integer :: i,k
  do k=0,k1
    if (k .eq. 0) then
      tmp = vector(:,0) 
    else if (k .eq. k1) then
      tmp = vector(:,k1-1)
    else
      tmp = (vector(:,k)+vector(:,k-1))/2.
    endif
    do i=0,i1
      if (i .eq. 0) then
        output(i,k) = (tmp(i)+tmp(i+1))/2.
      else if (i .eq. i1) then
        output(i,k) = (tmp(i1-1)+tmp(i1))/2.
      else
        output(i,k) = tmp(i)
      endif
    enddo
  enddo
end subroutine set_wvector_to_coords

subroutine set_uvector_to_coords(vector, i1,k1,output)
implicit None
  integer, intent(IN) :: i1, k1
  real(8), dimension(0:i1,0:k1), intent(IN) :: vector
  real(8), dimension(0:i1,0:k1), intent(OUT) :: output
  real(8), dimension(0:i1) :: tmp
  integer :: i,k
  do k=0,k1
    if (k .eq. 0) then
      tmp = (vector(:,0) + vector(:,1))/2.0
    else if (k .eq. k1) then
      tmp = (vector(:,k1-1) + vector(:,k1))/2.0
    else
      tmp = vector(:,k)
    endif
    do i=0,i1
      if (i .eq. 0) then
        output(i,k) = tmp(i)
      else if (i .eq. i1) then
        output(i,k) = tmp(i1-1)
      else
        output(i,k) = (tmp(i)+tmp(i-1))/2.0
      endif
    enddo
  enddo
end subroutine set_uvector_to_coords




!***************************************************************************************
!   writes the output file in 2D tecplot style
!***************************************************************************************


subroutine output2d_upd2(rank,istap)
  use mod_param
  use mod_common
  use mod_tm
  use mod_eos, only : eos_model
  use mod_tdm, only : turbdiff_model
  use mod_mesh , only : mesh
  implicit none
  include 'mpif.h'
  real*8 pecletx,peclety,pecletz
  real*8 massflow(kmax),enthflow(kmax),enth_b(kmax),Twall(kmax),Tbulk(kmax), &
         massfl,w_c,ndt
  real(8), dimension(0:i1,0:k1) ::  xvec, yvec,nuSA_sol,k_sol,eps_sol,om_sol,v2_sol,yp_sol, &
                                    kt_sol, epst_sol, Prt_sol
  integer rank,istap,jstart
  real(8), dimension(0:i1,0:k1) :: x_plt, y_plt, u_plt,w_plt, rho_plt, c_plt, p_plt, mu_plt, mut_plt, &
                                   k_plt, eps_plt, v2_plt, om_plt, nuSA_plt,alphat_plt, &
                                   prt_plt, kt_plt, epst_plt
  real(8), dimension(0:i1) :: dru, rp
  dru = mesh%dru
  rp = mesh%rp
  
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
      ! xvec(i,k)=(k+rank*kmax)*dz-(1./2.)*dz 
      xvec(i,k) = mesh%zp(k)
      yvec(i,k) = mesh%y_cv(i)
    enddo
  enddo

  call turb_model%get_sol(nuSA_sol,k_sol,eps_sol,om_sol,v2_sol,yp_sol)
  call turbdiff_model%get_sol(prt_sol,epst_sol,kt_sol)

  call set_uvector_to_coords(unew,i1,k1,u_plt)
  call set_wvector_to_coords(wnew,i1,k1,w_plt)
  call set_scalar_to_coords (rnew,i1,k1,rho_plt)
  call set_scalar_to_coords (cnew,i1,k1,c_plt)
  call set_scalar_to_coords (p,i1,k1,p_plt)
  call set_scalar_to_coords (ekm,i1,k1,mu_plt)
  call set_scalar_to_coords (ekmt,i1,k1,mut_plt)
  call set_scalar_to_coords (k_sol,i1,k1,k_plt)
  call set_scalar_to_coords (eps_sol,i1,k1,eps_plt)
  call set_scalar_to_coords (v2_sol,i1,k1,v2_plt)
  call set_scalar_to_coords (om_sol,i1,k1,om_plt)
  call set_scalar_to_coords (nuSA_sol,i1,k1,nuSA_plt)
  call set_scalar_to_coords (alphat,i1,k1,alphat_plt)
  call set_scalar_to_coords (prt_sol,i1,k1,prt_plt)
  call set_scalar_to_coords (kt_sol,i1,k1,kt_plt)
  call set_scalar_to_coords (epst_sol,i1,k1,epst_plt)

  ! call write_vector(turbdiff_model%Pkt,i1,k1,rank)
  
  call write_mpiio_formatted(trim(output_fname), xvec, yvec, u_plt,w_plt, rho_plt,c_plt,p_plt,mu_plt, mut_plt,yp_sol,     &
                                 k_plt, eps_plt, v2_plt, om_plt,nuSA_plt,alphat_plt, prt_plt, kt_plt, epst_plt,           &
                                 i1, k1,rank,px)
end

subroutine output2d_upd(rank,istap)
  use mod_param
  use mod_mesh, only : mesh
  use mod_common
  ! use mod_common2
  implicit none
  include 'mpif.h'
  character*5 cha
  real*8 pecletx,peclety,pecletz
  real*8 massflow(kmax),enthflow(kmax),enth_b(kmax),Twall(kmax),Tbulk(kmax), &
         massfl,w_c,ndt
  integer rank,istap,jstart
  write(cha,'(I5.5)')rank


  ! twall    = 0.0
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

  open(15,file=cha)
  
  ! if (rank.eq.0) then
  !   write(15,*) 'VARIABLES ="X","Y","U","W","P","C","T","k","eps", "v2","omega","nuSA","yplus","RHO","Pe","mu","mut"'
  !   write(15,*) 'ZONE I=  ', imax+2,' J=  ',(kmax+2)*px,' F=POINT '
  ! endif

  do k=0,k1
    do i=0,i1
      write(15,'(17ES24.10E3)')  mesh%zp(k), mesh%y_cv(i), &
        unew(i,k), &
        Wnew(i,k), p(min(max(i,1),imax),min(max(k,1),kmax)), cnew(i,k), cnew(i,k), &
        cnew(i,k),cnew(i,k),cnew(i,k), &
        cnew(i,k),cnew(i,k),cnew(i,k),rnew(i,k),peclet(i,k),ekm(i,k),ekmt(i,k)
    enddo
  enddo
  
  close(15)

end

subroutine write_vector(vector, i1, k1, rank)
  implicit none
  real(8), dimension(0:i1,0:k1), intent(IN) :: vector
  integer, intent(IN) :: i1, k1, rank
  integer :: i,k
  character*4 cha
  write(cha,'(I4.4)')rank
  open(15,file=cha)

  do i=0,i1
    do k=0,k1
      write(15,'(3ES24.10E3)') i*1., k*1., vector(i,k)
    enddo
  enddo
  close(15)


end subroutine