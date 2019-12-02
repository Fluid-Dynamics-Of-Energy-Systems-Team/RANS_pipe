!!********************************************************************
!!     poisson solver
!!********************************************************************
subroutine fillps(rank)
  use mod_param, only : i, k, kmax, imax,k1
  use mod_common,only : dudt,dwdt,rnew,rold,p,qcrit,dt
  use mod_mesh,  only : dz, dru, rp, ru
  use module_mesh, only : mesh
  implicit none
  include 'mpif.h'
  integer ierr,rank
  real*8 sumps,sumps_tot
  real*8,dimension(0:k1) ::  dzw, dzp
  dzw = mesh%dzw
  dzp = mesh%dzp


  !
  !     *** Fill the right hand for the poisson solver. ***
  !
  sumps = 0.0

  do  k=1,kmax
    do i=1,imax
      ! p(i,k)  = (                                                         &
      !             (Ru(i)*dUdt(i,k) - Ru(i-1)*dUdt(i-1,k))/( Rp(i)*dru(i)) &
      !             +     (dWdt(i,k) -         dWdt(i,k-1))/( dz         )  &
      !           )/dt                                                      &
      !         + (rnew(i,k)-rold(i,k))/(dt*dt)

      p(i,k)  = (                                                         &
                  (Ru(i)*dUdt(i,k) - Ru(i-1)*dUdt(i-1,k))/( Rp(i)*dru(i)) &
                  +     (dWdt(i,k) -         dWdt(i,k-1))/( dzw(k)     )  &
                )/dt                                                      &
              + (rnew(i,k)-rold(i,k))/(dt*dt)

      qcrit(i,k) = p(i,k)*dt

      ! sumps = sumps + p(i,k)*dru(i)*dz
      sumps = sumps + p(i,k)*dru(i)*dzw(k)

    enddo
  enddo


  call mpi_allreduce(sumps,sumps_tot,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)

end

subroutine calc_momentum_thickness(w, w_fs, dru, i1,imax, mom_thickness)
  implicit none
  real(8), dimension(0:i1),   intent(IN) :: w
  real(8), dimension(1:imax), intent(IN) :: dru
  real(8),                    intent(IN) :: w_fs
  integer,                    intent(IN) :: i1, imax
  real(8),                    intent(OUT):: mom_thickness
  integer :: i
  mom_thickness = 0.
  do i=1,imax
    mom_thickness = mom_thickness + (w(i)/w_fs)*(1.-(w(i)/w_fs))*dru(i)
  enddo
end subroutine

subroutine calc_displacement_thickness(w, w_fs, dru, i1, imax, dis_thickness)
  implicit none
  integer,                    intent(IN) :: i1, imax
  real(8), dimension(0:i1),   intent(IN) :: w
  real(8), dimension(1:imax), intent(IN) :: dru
  real(8),                    intent(IN) :: w_fs
  real(8),                    intent(OUT):: dis_thickness
  integer :: i
  dis_thickness = 0.
  do i=1,imax
    dis_thickness = dis_thickness + (1.-(w(i)/w_fs))*dru(i)
  enddo
end subroutine

subroutine calc_bl_thickness(w, w_fs, y_vec, i1, imax, bl_thickness)
  use mod_math, only : spline, splint
  implicit none
  integer,                  intent(IN) :: i1, imax
  real(8), dimension(0:i1), intent(IN) :: w
  real(8), dimension(0:i1), intent(IN) :: y_vec
  real(8),                  intent(IN) :: w_fs
  real(8),                  intent(OUT):: bl_thickness
  real(8), dimension(:),allocatable :: w_int, y_int, y2_int
  integer :: i, elem
  integer :: tabkhi,tabklo = 0 
  real(8), dimension(0:i1) :: w_inv, y_inv

  !inverse beause bl is at top
  do i=0, i1
    w_inv(i) = w(i1-i)
    y_inv(i) = 1.-y_vec(i1-i)
  enddo
  do i=1,i1
     if (w_inv(i) .ge. 0.99*w_fs) then
      elem = i
      exit
     endif
 enddo
 
 bl_thickness = y_inv(elem-1)+(0.99*w_fs-w_inv(elem-1))* &
               ((y_inv(elem)-y_inv(elem-1))/(w_inv(elem) -w_inv(elem-1)))
   
end subroutine

subroutine calc_shear_stress(w, mui,drp, i1,imax,stress)
  implicit none
  integer,                  intent(IN) :: i1,imax
  real(8), dimension(0:i1), intent(IN) :: w,mui, drp
  real(8),                  intent(OUT):: stress
  stress = mui(imax)*((w(imax)-w(i1))/drp(imax))
end subroutine

subroutine calc_skin_friction(w, mui,drp,rho_fs, w_fs,i1,imax, sfriction)
  implicit none
  integer,                  intent(IN) :: i1,imax
  real(8), dimension(0:i1), intent(IN) :: w,mui, drp
  real(8),                  intent(IN) :: rho_fs, w_fs
  real(8),                  intent(OUT):: sfriction
  real(8)                              :: stress
  call calc_shear_stress(w,mui,drp,i1,imax, stress)
  sfriction = stress/(0.5*rho_fs*w_fs*w_fs)
end subroutine

subroutine postprocess_bl(w, mui, rho_fs, w_fs, &
                          mom_th, dis_th, bl_th, stress,sfriction)
  use mod_param, only : k1,i1,imax
  use mod_mesh, only : top_bcnovalue, y_cv, drp, dru
  implicit none
  real(8), dimension(0:i1,0:k1), intent(IN) :: w, mui
  real(8), intent(IN) :: rho_fs, w_fs
  real(8), dimension(0:k1), intent(OUT) :: mom_th, dis_th, bl_th, stress, sfriction
  integer :: k

  do k=0,k1
    if (top_bcnovalue(k) .eq. -1) then
      call calc_momentum_thickness(w(:,k),w_fs,dru, i1, imax, mom_th(k))
      call calc_displacement_thickness(w(:,k), w_fs,dru, i1, imax, dis_th(k))
      call calc_bl_thickness(w(:,k), w_fs, y_cv, i1, imax, bl_th(k))
      call calc_shear_stress(w(:,k), mui(:,k),drp, i1,imax, stress(k))
      call calc_skin_friction(w(:,k), mui(:,k),drp, rho_fs, w_fs,i1,imax, sfriction(k))
    else
      mom_th(k)=0; dis_th(k)=0; bl_th(k)=0; stress(k)=0; sfriction(k)=0;
    endif
  enddo
end subroutine

  subroutine interpolate_vector(x_old, y_old, x_new, y_new, i1_old, k1_old, i1, k1, vector, output,rank)
  use mod_math, only : spline, splint
  implicit none
  integer, intent(IN) :: i1_old, i1, k1_old, k1, rank
  real(8), dimension(0:i1_old, 0:k1_old), intent(IN) :: vector, x_old, y_old
  real(8), dimension(0:i1, 0:k1), intent(IN) :: output, x_new, y_new
  real(8), dimension(0:i1_old) :: vectory2  
  real(8), dimension(0:k1_old) :: vectorx2  
  real(8), dimension(0:i1) :: tmpvecy
  real(8), dimension(0:i1,0:k1_old) :: tmpvec
  integer :: k,i

  integer :: tabkhi,tabklo = 0

  ! interpolate on the new y coordinates same x coords
  do k=0,k1_old
    call spline(y_old(:,k),vector(:,k),i1_old+1,vectory2)
    do i=0,i1
      call splint(y_old(:,k),vector(:,k),vectory2,i1_old+1,y_new(i,1),tmpvec(i,k),tabkhi,tabklo)
    enddo
  enddo

 !interpolate on hte new x coordinates with the already interpolated y coords
  do i=0,i1
    call spline(x_old(1,:),tmpvec(i,:),k1_old+1,vectorx2)
      do k=0,k1
      call splint(x_old(1,:),tmpvec(i,:),vectorx2,k1_old+1,x_new(i,k),output(i,k),tabkhi,tabklo)
    enddo
  enddo

end subroutine interpolate_vector

subroutine interpolate_solution(i1_old, k1_old, rank, px)
  use mod_common, only : wnew, unew, rnew, ekm, ekmt, cnew, win, ekmtin, uin
  use mod_mesh, only : y_cv, y_fa, dz
  use mod_param, only : read_fname,i1,k1, kmax, periodic
  use module_mesh, only : mesh

  implicit none
  include "mpif.h"
  integer,                  intent(IN) :: rank, px, i1_old, k1_old
  real(8),dimension(0:i1_old,0:k1_old) :: xold,yold,uold,wold,dummy,rold, ekmold, ekmtold,cold, pold
  real(8),dimension(0:i1,0:k1) :: xw, x, yu, y
  real(8), dimension(0:i1_old) :: tmp
  integer ::  i,k
  integer :: ierror
  real*8,dimension(0:k1) ::  zw, zp
  zw = mesh%zw
  zp = mesh%zp

  call read_mpiio_formatted(trim(read_fname), xold, yold, uold,wold,rold,cold,pold,ekmold, ekmtold,dummy,     &
                                 dummy, dummy, dummy, dummy,dummy,dummy, i1_old, k1_old,rank,px)

  !interpolate on the y values of the new grid
  do k=0,k1
    do i=0,i1
      ! x (i,k)= (k+rank*kmax)*dz-(1./2.)*dz  
      x(i,k) = zp(i)
      y (i,k)= y_cv(i)
      ! xw(i,k)= (k+rank*kmax)*dz !
      xw(i,k) = zw(i)
      yu(i,k)= y_fa(i)
    enddo
  enddo
    
  !commuicate the values that are on the previous or next core
  call shiftf1(yold,tmp,rank,i1_old,k1_old);     yold(:,0)      = tmp(:);
  call shiftf1(xold,tmp,rank,i1_old,k1_old);     xold(:,0)      = tmp(:);
  call shiftb1(yold,tmp,rank,i1_old,k1_old);     yold(:,k1_old) = tmp(:);
  call shiftb1(xold,tmp,rank,i1_old,k1_old);     xold(:,k1_old) = tmp(:);

  if ((rank .ne. 0) .or. (periodic .eq. 1)) then 
    write(*,*) "started",rank
    call shiftf1(wold,tmp,rank,i1_old,k1_old);     wold(:,0)      = tmp(:);
    write(*,*) "finished",rank

    call MPI_Barrier(MPI_COMM_WORLD);
    write(*,*) "started2",rank
    call shiftf1(uold,tmp,rank,i1_old,k1_old);     uold(:,0)      = tmp(:);
    write(*,*) "finished2",rank
    
  endif  


  if ((rank .ne. px-1) .or. (periodic.eq.1)) then
    call shiftb1(uold,tmp,rank,i1_old,k1_old);     uold(:,k1_old) = tmp(:);
    call shiftb1(wold,tmp,rank,i1_old,k1_old);     wold(:,k1_old) = tmp(:);
  endif


  !interpolate on the new grid
  call interpolate_vector(xold,yold,x, yu,i1_old, k1_old,i1,k1,uold,unew,rank)
  call interpolate_vector(xold,yold,xw,y, i1_old, k1_old,i1,k1,wold,wnew,rank)  
  call interpolate_vector(xold,yold,x, y, i1_old, k1_old,i1,k1,rold,rnew,rank)
  call interpolate_vector(xold,yold,x, y, i1_old, k1_old,i1,k1,cold,cnew,rank)
  call interpolate_vector(xold,yold,x, y, i1_old, k1_old,i1,k1,ekmtold,ekmt,rank)
  
  if (periodic .eq. 1) return
  if (rank .eq. 0) then
    win(:) = wnew(:,0)
    ekmtin(:) = ekmt(:,0)
    uin(:) = unew(:,0)
  endif

  call MPI_Bcast(win,   i1+1, MPI_REAL8,0,MPI_COMM_WORLD, ierror)
  call MPI_Bcast(ekmtin,i1+1, MPI_REAL8,0,MPI_COMM_WORLD, ierror)
  call MPI_Bcast(uin,   i1+1, MPI_REAL8,0,MPI_COMM_WORLD, ierror)
  

end subroutine

!!********************************************************************
!!     correc
!!********************************************************************
subroutine correc(rank,setold)
  use mod_param
  use mod_common
  use mod_mesh
  use module_mesh, only : mesh
  implicit none
      
  integer rank,setold
  real*8 pplus_w(imax)
  real*8,dimension(0:k1) ::  dzw, dzp
  dzw = mesh%dzw
  dzp = mesh%dzp

  do k=1,kmax
    do i=1,imax-1
      dUdt(i,k)=dUdt(i,k)-dt*(p(i+1,k)-p(i,k))/dRp(i) !(Rp(i+1)-Rp(i))
    enddo
  enddo

  do k=1,kmax-1
    do i=1,imax
      ! dWdt(i,k)=dWdt(i,k)-dt*(p(i,k+1)-p(i,k))/dz
      dWdt(i,k)=dWdt(i,k)-dt*(p(i,k+1)-p(i,k))/dzp(k)

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
    ! dWdt(i,kmax)=dWdt(i,kmax)-dt*(pplus_w(i)-p(i,kmax))/dz
    dWdt(i,kmax)=dWdt(i,kmax)-dt*(pplus_w(i)-p(i,kmax))/dzp(kmax)

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

!!********************************************************************
!!     chkdt
!!********************************************************************
subroutine chkdt(rank,istap)
  use mod_param
  use mod_common
  use mod_mesh
  use module_mesh, only : mesh
  implicit none
  include 'mpif.h'
  integer rank,ierr,istap
  real*8  tmp,Courant,dtmp,tmp1,tmp2,tmp3,dr2,dz2,kcoeff
  real*8,dimension(0:k1) ::  dzw, dzp
  dzw = mesh%dzw
  dzp = mesh%dzp

  dt = dtmax

  do k=1,kmax
    do i=1,imax
      ! tmp = ( abs(Unew(i,k)) /  dRp(i) ) &
      !      +( abs(Wnew(i,k)) /  dz     )
      tmp = ( abs(Unew(i,k)) /  dRp(i) ) &
           +( abs(Wnew(i,k)) /  dzp(k) )
      
      tmp = CFL/tmp
      dt  = min(dt, tmp)
    enddo
  enddo

  dtmp = dt

  call mpi_allreduce(dtmp,dt,1,mpi_real8,mpi_min,mpi_comm_world,ierr)

end

!!********************************************************************
!!     chkdiv
!!********************************************************************
subroutine chkdiv(rank)
  use mod_param, only : kmax, imax, i, k,k1
  use mod_common,only : rnew, unew, wnew, dt, rold
  use mod_mesh,  only : dz, dru, ru, rp
  use module_mesh, only : mesh
  implicit none    
  include 'mpif.h'
  integer rank,ierr,ll
  real*8   div,divmax,divbar,divmax_tot,divbar_tot,rhoip,rhoim,rhokp,rhokm
  real*8,dimension(0:k1) ::  dzw, dzp
  dzw = mesh%dzw
  dzp = mesh%dzp

  divbar = 0.0
  divmax = 0.0

  do k=1,kmax
    do i=1,imax
      rhoip = 0.5*(rNew(i,k)+rNew(i+1,k))
      rhoim = 0.5*(rNew(i,k)+rNew(i-1,k))

      rhokp = 0.5*(rNew(i,k)+rNew(i,k+1))
      rhokm = 0.5*(rNew(i,k)+rNew(i,k-1))

      ! div = &
      !   (1./Rp(i))*(Ru(i)*Unew(i,k)*rhoip-Ru(i-1)*Unew(i-1,k)*rhoim)*dz      +     &
      !              (      Wnew(i,k)*rhokp-        Wnew(i,k-1)*rhokm)*dru(i)  +     &
      !              (      rNew(i,k)      -        rold(i,k)        )/(dt*dru(i)*dz)

      div = &
        (1./Rp(i))*(Ru(i)*Unew(i,k)*rhoip-Ru(i-1)*Unew(i-1,k)*rhoim)*dzw(k)  +     &
                   (      Wnew(i,k)*rhokp-        Wnew(i,k-1)*rhokm)*dru(i)  +     &
                   (      rNew(i,k)      -        rold(i,k)        )/(dt*dru(i)*dzw(k))


        ! if (abs(div).gt.10e-10) write(6,*) i,k+kmax*rank,div

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

!!********************************************************************
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


!!********************************************************************
!!     diffusion term in the z-direction, set as a source term...
!!********************************************************************
subroutine diffc(putout,putin,ek,eki,ekk,ekmt,sigma,rho,Ru,Rp,dru,dz,rank1,diffVersion)
  use mod_param
  use module_mesh, only : mesh
  implicit none
  integer   km,kp,rank1,diffVersion
  real*8     putout(0:i1,0:k1),putin(0:i1,0:k1), &
    rho(0:i1,0:k1),ek(0:i1,0:k1),eki(0:i1,0:k1),ekk(0:i1,0:k1), &
    ekmt(0:i1,0:k1),dru(0:i1),dz,Ru(0:i1),Rp(0:i1),sigma
  real(8), dimension(0:k1) :: dzw, dzp

  dzw = mesh%dzw
  dzp = mesh%dzp

  if (diffVersion == 1) then       ! Inverse SLS
    do k=1,kmax
      kp=k+1
      km=k-1
      do i=1,imax
        ! putout(i,k) = putout(i,k) + 1.0/rho(i,k)/sqrt(rho(i,k))*( &
        !   (  (ekk(i,k ) + 0.5*(ekmt(i,k)+ekmt(i,kp))/sigma)/sqrt(0.5*(rho(i,k)+rho(i,kp))) &
        !       *(rho(i,kp)*putin(i,kp)-rho(i,k )*putin(i,k )) &
        !     -(ekk(i,km) + 0.5*(ekmt(i,k)+ekmt(i,km))/sigma)/sqrt(0.5*(rho(i,k)+rho(i,km))) &
        !       *(rho(i,k )*putin(i,k )-rho(i,km)*putin(i,km)) &
        !   )/(dz*dz))

        putout(i,k) = putout(i,k) + 1.0/rho(i,k)/sqrt(rho(i,k))*( &
          (  (ekk(i,k ) + 0.5*(ekmt(i,k)+ekmt(i,kp))/sigma)/sqrt(0.5*(rho(i,k)+rho(i,kp))) &
              *(rho(i,kp)*putin(i,kp)-rho(i,k )*putin(i,k ))/dzp(k) &
            -(ekk(i,km) + 0.5*(ekmt(i,k)+ekmt(i,km))/sigma)/sqrt(0.5*(rho(i,k)+rho(i,km))) &
              *(rho(i,k )*putin(i,k )-rho(i,km)*putin(i,km))/dzp(km) &
          )/dzw(km))
      enddo
    enddo
  elseif (diffVersion == 2) then   ! Aupoix
    do k=1,k1-1
      kp=k+1
      km=k-1
      do i=1,i1-1
        ! putout(i,k) = putout(i,k) + 1.0/rho(i,k)*( &
        !   (  (ekk(i,k ) + 0.5*(ekmt(i,k)+ekmt(i,kp))/sigma)/(0.5*(rho(i,k)+rho(i,kp))) &
        !       *(rho(i,kp)*putin(i,kp)-rho(i,k )*putin(i,k )) &
        !     -(ekk(i,km) + 0.5*(ekmt(i,k)+ekmt(i,km))/sigma)/(0.5*(rho(i,k)+rho(i,km))) &
        !       *(rho(i,k )*putin(i,k )-rho(i,km)*putin(i,km)) &
        !   )/(dz*dz))

        putout(i,k) = putout(i,k) + 1.0/rho(i,k)*( &
          (  (ekk(i,k ) + 0.5*(ekmt(i,k)+ekmt(i,kp))/sigma)/(0.5*(rho(i,k)+rho(i,kp))) &
              *(rho(i,kp)*putin(i,kp)-rho(i,k )*putin(i,k ))/dzp(k)  &
            -(ekk(i,km) + 0.5*(ekmt(i,k)+ekmt(i,km))/sigma)/(0.5*(rho(i,k)+rho(i,km))) &
              *(rho(i,k )*putin(i,k )-rho(i,km)*putin(i,km))/dzp(km) &
          )/dzw(km))

      enddo
    enddo
  else                             ! Standard
    do k=1,k1-1
      kp=k+1
      km=k-1
      do i=1,i1-1
        ! putout(i,k) = putout(i,k) + 1.0/rho(i,k)*( &
        ! (  (ekk(i,k ) + 0.5*(ekmt(i,k)+ekmt(i,kp))/sigma)*(putin(i,kp)-putin(i,k )) &
        !   -(ekk(i,km) + 0.5*(ekmt(i,k)+ekmt(i,km))/sigma)*(putin(i,k )-putin(i,km)) &
        ! )/(dz*dz))

        putout(i,k) = putout(i,k) + 1.0/rho(i,k)*( &
        (  (ekk(i,k ) + 0.5*(ekmt(i,k)+ekmt(i,kp))/sigma)*(putin(i,kp)-putin(i,k ))/dzp(k) &
          -(ekk(i,km) + 0.5*(ekmt(i,k)+ekmt(i,km))/sigma)*(putin(i,k )-putin(i,km))/dzp(km) &
        )/dzw(km))

      enddo
    enddo
  endif
end

     
!!*****************************************************************
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
  use module_mesh, only : mesh
  implicit none
  integer  i,k,im,ip,km,kp,i1,k1,numDom
  real*8     putout(0:i1,0:k1),Uvel(0:i1,0:k1),Wvel(0:i1,0:k1), &
    ekme(0:i1,0:k1),dru(0:i1),drp(0:i1),dz,Ru(0:i1),Rp(0:i1), &
    epop,epom,dzi,divUim,divUip,divUi,dif
  real*8,dimension(0:k1) ::  dzw, dzp
  dzw = mesh%dzw
  dzp = mesh%dzp

  !pipe
  if (numDom == -1) then

    dzi =1./dz
    do k=1,k1-1
      kp=k+1
      km=k-1
      do i=1,i1-1
        ip=i+1
        im=i-1

        epop = 0.25*(ekme(i,k)+ekme(ip,k) + ekme(ip,kp) + ekme(i,kp)) !right top
        epom = 0.25*(ekme(i,k)+ekme(ip,k) + ekme(ip,km) + ekme(i,km)) !left  top

        putout(i,k) = putout(i,k) + &
          2.0*( Rp(ip)*ekme(ip,k)*(dif*(Uvel(ip,k)-Uvel(i ,k))/dru(ip)) &
               -Rp(i )*ekme(i ,k)*(dif*(Uvel(i ,k)-Uvel(im,k))/dru(i )) &
              )/(Ru(i)*drp(i)) &
          +                                              &
          ! ( epop * ( (Uvel(i,kp )-Uvel(i,k ))/dz         &
          !           +(Wvel(ip,k )-Wvel(i,k ))/drp(i))    &
          ! - epom * ( (Uvel(i,k  )-Uvel(i,km))/dz         &
          !           +(Wvel(ip,km)-Wvel(i,km))/drp(i))    &  
          ! )/dz
          ( epop * ( (Uvel(i,kp )-Uvel(i,k ))/dzp(k)     &
                    +(Wvel(ip,k )-Wvel(i,k ))/drp(i))    &
          - epom * ( (Uvel(i,k  )-Uvel(i,km))/dzp(km)    &
                    +(Wvel(ip,km)-Wvel(i,km))/drp(i))    &  
          )/dzw(k)
      enddo
    enddo

  !channel
  elseif (numDom == 1) then
      
    do k=1,k1-1
      kp=k+1
      km=k-1
      do i=1,i1-1
        ip=i+1
        im=i-1

        epop = 0.25*(ekme(i,k)+ekme(ip,k) + ekme(ip,kp) + ekme(i,kp)) !right top
        epom = 0.25*(ekme(i,k)+ekme(ip,k) + ekme(ip,km) + ekme(i,km)) !left  top

        ! divUim = (Ru(i )*Uvel(i ,k)              &
        !          -Ru(im)*Uvel(im,k)              &
        !          )/(Rp(i )*dru(i ))              &
        !        + ( Wvel(i ,k) - Wvel(i ,km))/dz

        divUim = (Ru(i )*Uvel(i ,k)              &
                 -Ru(im)*Uvel(im,k)              &
                 )/(Rp(i )*dru(i ))              &
               + ( Wvel(i ,k) - Wvel(i ,km))/dzw(k)

        ! divUip = (Ru(ip)*Uvel(ip,k)              &
        !          -Ru(i )*Uvel(i ,k)              &
        !          )/(Rp(ip)*dru(ip))              &
        !        + ( Wvel(ip,k) - Wvel(ip,km))/dz

        divUip = (Ru(ip)*Uvel(ip,k)              &
                 -Ru(i )*Uvel(i ,k)              &
                 )/(Rp(ip)*dru(ip))              &
               + ( Wvel(ip,k) - Wvel(ip,km))/dzw(k)


        ! divUi =  (Rp(ip)*(Uvel(ip,k)+Uvel(i ,k)) &
        !          -Rp(i )*(Uvel(i ,k)+Uvel(im,k)) &
        !          )/(2.*Ru(i)*drp(i))             &
        !        + ((Wvel(ip,k)+Wvel(i,k))-(Wvel(ip,km)+Wvel(i,km)))/(2.*dz)

        divUi =  (Rp(ip)*(Uvel(ip,k)+Uvel(i ,k)) &
                 -Rp(i )*(Uvel(i ,k)+Uvel(im,k)) &
                 )/(2.*Ru(i)*drp(i))             &
               + ((Wvel(ip,k)+Wvel(i,k))-(Wvel(ip,km)+Wvel(i,km)))/(2.*dzw(k))

        
        putout(i,k) = putout(i,k) + &
          2.0*( Rp(ip)*ekme(ip,k)*(dif*(Uvel(ip,k)-Uvel(i ,k))/dru(ip) -1./3.*divUip) &
               -Rp(i )*ekme(i ,k)*(dif*(Uvel(i ,k)-Uvel(im,k))/dru(i ) -1./3.*divUim) &
              )/(Ru(i)*drp(i)) &
          +(                                            &
          !   epop * ( (Uvel(i ,kp)-Uvel(i,k ))/dz        &
          !           +(Wvel(ip,k )-Wvel(i,k ))/drp(i))   &
          ! - epom * ( (Uvel(i ,k )-Uvel(i,km))/dz        &
          !           +(Wvel(ip,km)-Wvel(i,km))/drp(i))   &
          !  )/dz                                         &
            epop * ( (Uvel(i ,kp)-Uvel(i,k ))/dzp(k)    &
                    +(Wvel(ip,k )-Wvel(i,k ))/drp(i))   &
          - epom * ( (Uvel(i ,k )-Uvel(i,km))/dzp(km)   &
                    +(Wvel(ip,km)-Wvel(i,km))/drp(i))   &
           )/dzw(k)                                     &
          - (ekme(i,k)+ekme(ip,k))/Ru(i)*(Uvel(i,k)/Ru(i)-1./3.*divUi)
      enddo
    enddo
  endif



  return
end


!!*****************************************************************
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
  use module_mesh, only : mesh
  implicit none
     

  integer  i,k,im,ip,km,kp,i1,k1,numDom
  real*8     putout(0:i1,0:k1),Uvel(0:i1,0:k1),Wvel(0:i1,0:k1),  &
    ekme(0:i1,0:k1),dru(0:i1),drp(0:i1),dz,Ru(0:i1),Rp(0:i1), &
    epop,emop,divUkm,divUkp,dif
  real*8,dimension(0:k1) ::  dzw, dzp
  dzw = mesh%dzw
  dzp = mesh%dzp

  if (numDom == -1) then
    do k=1,k1-1
      kp=k+1
      km=k-1
      do i=1,i1-1
        ip=i+1
        im=i-1

        epop = 0.25*(ekme(i,k)+ekme(i,kp) + ekme(ip,k) + ekme(ip,kp) ) !right top
        emop = 0.25*(ekme(i,k)+ekme(i,kp) + ekme(im,k) + ekme(im,kp) ) !right bot

        ! putout(i,k) = putout(i,k) + &
        !   (Ru(i )*epop*((Uvel(i ,kp)-Uvel(i ,k))/dz      &
        !            +dif*(Wvel(ip,k) -Wvel(i ,k))/drp(i ))&
        !   -                                              &
        !    Ru(im)*emop*((Uvel(im,kp)-Uvel(im,k))/dz      &
        !            +dif*(Wvel(i ,k )-Wvel(im,k))/drp(im))&
        !   )/(Rp(i)*dru(i))                               &
        !   +                                              &
        !   (2.*ekme(i,kp)*((Wvel(i,kp)-Wvel(i,k ))/dz)    &
        !   -2.*ekme(i,k )*((Wvel(i,k )-Wvel(i,km))/dz)    &
        !   )/dz
        putout(i,k) = putout(i,k) + &
          (Ru(i )*epop*((Uvel(i ,kp)-Uvel(i ,k))/dzp(k)  &
                   +dif*(Wvel(ip,k) -Wvel(i ,k))/drp(i ))&
          -                                              &
           Ru(im)*emop*((Uvel(im,kp)-Uvel(im,k))/dzp(k) &
                   +dif*(Wvel(i ,k )-Wvel(im,k))/drp(im)) &
          )/(Rp(i)*dru(i))                                &
          +                                               &
          (2.*ekme(i,kp)*((Wvel(i,kp)-Wvel(i,k ))/dzw(kp)) &
          -2.*ekme(i,k )*((Wvel(i,k )-Wvel(i,km))/dzw(k))&
          )/dzp(k)
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

        ! divUkm = (Ru(i)*Uvel(i,k )-Ru(im)*Uvel(im,k ))/(Rp(i)*dru(i)) &
        !        + (Wvel(i,k )      -       Wvel(i,km))/dz

        divUkm = (Ru(i)*Uvel(i,k )-Ru(im)*Uvel(im,k ))/(Rp(i)*dru(i)) &
               + (Wvel(i,k )      -       Wvel(i,km))/dzw(k)


        ! divUkp = (Ru(i)*Uvel(i,kp)-Ru(im)*Uvel(im,kp))/(Rp(i)*dru(i)) &
        !        + (Wvel(i,kp)      -       Wvel(i ,k ))/dz

        divUkp = (Ru(i)*Uvel(i,kp)-Ru(im)*Uvel(im,kp))/(Rp(i)*dru(i)) &
               + (Wvel(i,kp)      -       Wvel(i ,k ))/dzw(kp)


        ! putout(i,k) = putout(i,k) + &
        !   (                                                      &
        !   Ru(i )*epop*(  (Uvel(i ,kp)-Uvel(i ,k))/dz             &
        !             +dif*(Wvel(ip,k )-Wvel(i ,k))/drp(i))        &!(Rp(ip)-Rp(i))) new
        !   -                                                      &
        !   Ru(im)*emop*(  (Uvel(im,kp)-Uvel(im,k))/dz             &
        !   +       dif*(   Wvel(i ,k )-Wvel(im,k))/drp(im))       &!(Rp(i)-Rp(im))) x_new
        !   )/(Rp(i)*dru(i))                                       &
        !   +                                                      &
        !   (2.*ekme(i,kp)*((Wvel(i,kp)-Wvel(i,k ))/dz - 1./3.*divUkp) &
        !   -2.*ekme(i,k )*((Wvel(i,k )-Wvel(i,km))/dz - 1./3.*divUkm) &
        !   )/dz

        putout(i,k) = putout(i,k) + &
          (                                                      &
          Ru(i )*epop*(  (Uvel(i ,kp)-Uvel(i ,k))/dzp(k)         &
                    +dif*(Wvel(ip,k )-Wvel(i ,k))/drp(i))        &!(Rp(ip)-Rp(i))) new
          -                                                      &
          Ru(im)*emop*(  (Uvel(im,kp)-Uvel(im,k))/dzp(k)         &
          +       dif*(   Wvel(i ,k )-Wvel(im,k))/drp(im))       &!(Rp(i)-Rp(im))) x_new
          )/(Rp(i)*dru(i))                                       &
          +                                                      &
          (2.*ekme(i,kp)*((Wvel(i,kp)-Wvel(i,k ))/dzw(kp) - 1./3.*divUkp) &
          -2.*ekme(i,k )*((Wvel(i,k )-Wvel(i,km))/dzw(k) - 1./3.*divUkm) &
          )/dzp(k)
      enddo
    enddo
  endif
  return
end


!!********************************************************************
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
  use module_mesh, only : mesh
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
  real(8), dimension(0:k1) :: dzw, dzp
  logical flagImpl
  !
  dzw = mesh%dzw
  dzp = mesh%dzp
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
      ! putout(i,k) =    - (Ru(i)*U(i,k)*cu(i,k)-Ru(im)*U(im,k)*cu(im,k))/(Rp(i)*dru(i))  &
      !                  - (      W(i,k)*cw(i,k)-       W(i,km)*cw(i,km))/dz              &
      !       + putin(i,k)*(                                                              &
      !                    (Ru(i)*U(i,k)        -Ru(im)*U(im,k)         )/(Rp(i)*dru(i))  &
      !                  + (      W(i,k)        -       W(i,km)         )/dz              &
      !       )


      putout(i,k) =    - (Ru(i)*U(i,k)*cu(i,k)-Ru(im)*U(im,k)*cu(im,k))/(Rp(i)*dru(i))  &
                       - (      W(i,k)*cw(i,k)-       W(i,km)*cw(i,km))/dzw(k)          &
            + putin(i,k)*(                                                              &
                         (Ru(i)*U(i,k)        -Ru(im)*U(im,k)         )/(Rp(i)*dru(i))  &
                       + (      W(i,k)        -       W(i,km)         )/dzw(k)         &
      )
    enddo
  enddo

  !     Implicit part
  if (flagImpl .eqv. .true.) then
    do k=kb,ke
      km=k-1
      do i=ib,ie
        im=i-1
        ! putout(i,k)= putout(i,k)  + W(i,k)*putin(i,k)/dz
        ! dimpl(i,k) =  dimpl(i,k)  + W(i,k)/dz
        putout(i,k)= putout(i,k)  + W(i,k)*putin(i,k)/dzw(k)
        dimpl(i,k) =  dimpl(i,k)  + W(i,k)/dzw(k)
      enddo
    enddo
  endif

end


!!********************************************************************
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
  use module_mesh, only : mesh
  implicit none
  integer  i,k,im,ip,km,kp,i1,k1,ib,ie,kb,ke
  real*8     putout(0:i1,0:k1),Uvel(0:i1,0:k1),Wvel(0:i1,0:k1), &
    dru(0:i1),drp(0:i1),dz,Ru(0:i1),Rp(0:i1)
  real*8 rho(0:i1,0:k1)
  real*8 rhoip,rhoim,rhokp,rhokm
  real*8,dimension(0:k1) ::  dzw, dzp
  dzw = mesh%dzw
  dzp = mesh%dzp

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

      rhokp=0.25*(rho(i,k)+rho(i,kp)+rho(ip,k)+rho(ip,kp)) !right top
      rhokm=0.25*(rho(i,k)+rho(i,km)+rho(ip,k)+rho(ip,km)) !left top

      putout(i,k) = 0.0
      ! putout(i,k) = - 0.25 * (                                            &
      !   (Rp(ip)*(Uvel(i, k)+Uvel(ip,k))*(Uvel(i,k)+Uvel(ip,k))*rho(ip,k)  &
      !   -Rp(i )*(Uvel(im,k)+Uvel(i ,k))*(Uvel(i,k)+Uvel(im,k))*rho(i ,k)  &
      !   )/(Ru(i)*drp(i))                                                  &   
      !   + &
      !   ((Wvel(i,k )+Wvel(ip,k ))*(Uvel(i,k)+Uvel(i,kp))*rhokp            &
      !   -(Wvel(i,km)+Wvel(ip,km))*(Uvel(i,k)+Uvel(i,km))*rhokm            &
      !   )/(dz)                                                            &
      !   )
      putout(i,k) = - 0.25 * (                                            &
        (Rp(ip)*(Uvel(i, k)+Uvel(ip,k))*(Uvel(i,k)+Uvel(ip,k))*rho(ip,k)  &
        -Rp(i )*(Uvel(im,k)+Uvel(i ,k))*(Uvel(i,k)+Uvel(im,k))*rho(i ,k)  &
        )/(Ru(i)*drp(i))                                                  &   
        + &
        ((Wvel(i,k )+Wvel(ip,k ))*(Uvel(i,k)+Uvel(i,kp))*rhokp            &
        -(Wvel(i,km)+Wvel(ip,km))*(Uvel(i,k)+Uvel(i,km))*rhokm            &
        )/(dzw(k))                                                        &
        )
      
    enddo
  enddo
end




!!********************************************************************
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
  use module_mesh, only : mesh
  implicit none

  integer   im,ip,km,kp,ib,ie,kb,ke
  real*8    putout(0:i1,0:k1),Uvel(0:i1,0:k1),Wvel(0:i1,0:k1), &
            dru(0:i1),dz,Ru(0:i1),Rp(0:i1)
  real*8    rho(0:i1,0:k1),ekm(0:i1,0:k1),peclet_z(0:i1,0:k1)
  real*8    rhoip,rhoim,advcecw_w
  real*8,dimension(0:k1) ::  dzw, dzp
  dzw = mesh%dzw
  dzp = mesh%dzp



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

      rhoip=0.25*(rho(i,k)+rho(i,kp)+rho(ip,k)+rho(ip,kp)) !right top
      rhoim=0.25*(rho(i,k)+rho(i,kp)+rho(im,k)+rho(im,kp)) !right bottom

      ! peclet_z(i,k)= Wvel(i,k)*rho(i,k)*dz/(ekm(i,k))
      peclet_z(i,k)= Wvel(i,k)*rho(i,k)*dzw(k)/(ekm(i,k))
      
      ! if (peclet_z(i,k).gt.2.)then
      !   advcecw_w= 2.0*(rho(i,k)+rho(i,k))*Wvel(i,k)*(Wvel(i,k)-Wvel(i,km))/dz & 
      !              +                                                           & 
      !              2.0*Wvel(i,k)*(                                             &
      !               (rho(i,k )+rho(i,k ))*Wvel(i,k )                           &
      !              -(rho(i,km)+rho(i,km))*Wvel(i,km)                           &
      !              )/dz
      ! else
      !   advcecw_w= ((Wvel(i,k )+Wvel(i,kp))*(Wvel(i,k)+Wvel(i,kp))*rho(i,kp)   &
      !              -(Wvel(i,km)+Wvel(i,k ))*(Wvel(i,k)+Wvel(i,km))*rho(i,k )   &
      !              )/dz
      ! endif
      if (peclet_z(i,k).gt.2.)then
        advcecw_w= 2.0*(rho(i,k)+rho(i,k))*Wvel(i,k)*(Wvel(i,k)-Wvel(i,km))/dzw(k)& 
                   +                                                              & 
                   2.0*Wvel(i,k)*(                                                &
                    (rho(i,k )+rho(i,k ))*Wvel(i,k )                              &
                   -(rho(i,km)+rho(i,km))*Wvel(i,km)                              &
                   )/dzw(k)
      else
        advcecw_w= ((Wvel(i,k )+Wvel(i,kp))*(Wvel(i,k)+Wvel(i,kp))*rho(i,kp)   &
                   -(Wvel(i,km)+Wvel(i,k ))*(Wvel(i,k)+Wvel(i,km))*rho(i,k )   &
                   )/dzp(k)
      endif

      putout(i,k) = 0.0
      putout(i,k) = - 0.25 * (                                          &
        (Ru(i )*(Uvel(i ,k)+Uvel(i ,kp))*(Wvel(i,k)+Wvel(ip,k))*rhoip   &
        -Ru(im)*(Uvel(im,k)+Uvel(im,kp))*(Wvel(i,k)+Wvel(im,k))*rhoim   &
        )/(Rp(i)*dru(i))                                                &
        +advcecw_w)
    enddo
  enddo
end

