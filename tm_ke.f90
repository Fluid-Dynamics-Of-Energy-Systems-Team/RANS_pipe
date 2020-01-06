module ke_tm
  use mod_tm
  implicit none

!****************************************************************************************

  !************************!
  !         KE class       !
  !************************!
  
  type,abstract,extends(TurbModel), public :: KE_TurbModel
  real(8), dimension(:,:), allocatable :: Gk,f1,f2,fmu,Lh,fv2, feps
  real(8), dimension(:),   allocatable :: epsin, kin, v2in
  real(8) :: sigmak,sigmae,cmu
  contains
    procedure(set_mut_KE), deferred :: set_mut
    procedure(set_bc_KE), deferred :: set_bc
    procedure(set_constants), deferred :: set_constants
    procedure :: init_w_inflow => init_w_inflow_KE
    procedure :: advance_turb => advance_KE
    procedure :: production_KE
    procedure :: init => init_KE
    procedure :: rhs_k_KE
    procedure :: rhs_eps_KE
    procedure :: diffusion_eps_KE
    procedure :: solve_eps_KE
    procedure :: solve_k_KE
    procedure :: init_mem_KE
    procedure :: init_sol => init_sol_KE
    procedure :: get_profile => get_profile_KE
    procedure :: get_sol => get_sol_KE

  end type KE_TurbModel

  interface
    subroutine set_constants(this)
      import :: KE_TurbModel
      class(KE_TurbModel) :: this
    end subroutine set_constants
    subroutine set_mut_KE(this,u,w,rho,mu,mui,mut)
      use mod_param, only: k1,i1
      import :: KE_TurbModel
      class(KE_TurbModel) :: this
      real(8), dimension(0:i1,0:k1),intent(IN) :: u,w,rho,mu,mui
      real(8), dimension(0:i1,0:k1),intent(OUT):: mut
    end subroutine set_mut_KE
    subroutine set_bc_KE(this,mu,rho,periodic,rank,px)
      use mod_param, only: k1,i1
      import :: KE_TurbModel
      class(KE_TurbModel) :: this
      real(8),dimension(0:i1,0:k1),intent(IN) :: rho,mu
      integer,                               intent(IN) :: periodic, rank, px
    end subroutine set_bc_KE
  end interface

contains
!****************************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************************

  !************************!
  !       KE routines      !
  !************************!

subroutine init_KE(this)
  class(KE_TurbModel) :: this
  call this%init_mem_KE()
  call this%set_constants()
  call this%init_sol()
end subroutine init_KE

subroutine init_sol_KE(this)
  use mod_param, only : i1,i
  class(KE_TurbModel) :: this
  do i=0,i1
    this%k(i,:)  =0.1
    this%eps(i,:)=1.0
    this%v2(i,:) =2./3.*this%k(i,:)    
    this%Pk(i,:) = 0 
    this%Gk(i,:) = 0
    this%kin(i) = 0.1
    this%epsin(i) = 1.0
    this%v2in(i) =2./3.*this%kin(i)
    this%Pkin(i) = 0 
  enddo
end subroutine init_sol_KE

subroutine init_w_inflow_KE(this,nuSAin,pkin,kin,epsin,omin,mutin,v2in)
  use mod_param, only : i1,k1,k
  class(KE_TurbModel) :: this
  real(8), dimension(0:i1), intent(IN) :: nuSAin,pkin,kin,epsin,omin,mutin,v2in
  this%epsin = epsin
  this%kin   = kin
  this%mutin = mutin
  do k=0,k1
    this%eps(:,k) = this%epsin(:)
    this%k(:,k)   = this%kin(:)
  enddo
end subroutine init_w_inflow_KE

subroutine init_mem_KE(this)
  use mod_param, only : kmax,imax,k1,i1
  class(KE_TurbModel) :: this
  allocate(this%eps(0:i1,0:k1),this%k (0:i1,0:k1), this%Gk (0:i1,0:k1),this%Pk(0:i1,0:k1), &
           this%f1 (0:i1,0:k1),this%f2(0:i1,0:k1), this%fmu(0:i1,0:k1),this%Tt(0:i1,0:k1), &
           this%v2 (0:i1,0:k1),this%yp(0:i1,0:k1), this%fv2(imax,kmax),this%Lh(imax,kmax), &
           this%feps(0:i1,0:k1))
  allocate(this%mutin(0:i1),this%Pkin (0:i1),this%epsin(0:i1),this%kin(0:i1),this%v2in(0:i1))
end subroutine init_mem_KE

subroutine production_KE(this,u,w,temp,rho,mu,mut,beta)
  use mod_param, only : i1,k1
  implicit none
  class(KE_TurbModel) :: this
  real(8), dimension(0:i1,0:k1), intent(IN) :: u,w,temp,rho,mu,mut,beta
  real(8), dimension(0:i1,0:k1) :: div
end subroutine production_KE

subroutine get_sol_KE(this,nuSA,k,eps,om,v2,yp)
  use mod_param, only : i1,k1  
  class(KE_TurbModel) :: this
  real(8),dimension(0:i1,0:k1), intent(OUT):: nuSA,k,eps,om,v2,yp
  nuSA=0
  k   =this%k    
  eps =this%eps
  v2  =this%v2
  om  =0
  yp  = this%yp
end subroutine get_sol_KE

subroutine get_profile_KE(this,p_nuSA,p_k,p_eps,p_om,p_v2,p_Pk,p_bF1,p_bF2,yp,k)
  use mod_param, only : i1,imax
  class(KE_TurbModel) :: this
  integer,                               intent(IN) :: k
  real(8),dimension(0:i1),          intent(OUT):: p_nuSA,p_k,p_eps,p_om,p_v2,p_Pk, p_bF1,yp
  real(8),dimension(1:imax),        intent(OUT):: p_bF2

  p_nuSA(:)=0
  p_k(:)   =this%k(:,k)
  p_eps(:) =this%eps(:,k)
  p_v2(:)  =this%v2(:,k)
  p_om(:)  =0
  p_Pk(:)  =this%Pk(:,k)
  p_bF1(:) =0
  p_bF2(:) =0
  yp(:)    =this%yp(:,k)
end subroutine get_profile_KE


subroutine advance_KE(this,u,w,rho,mu,mui,muk,mut,beta,temp, &
                      alpha1,alpha2,alpha3,                  &
                      modification,rank,periodic,   &
                      residual1, residual2, residual3)
  use mod_param, only : i1,k1,imax,kmax
  class(KE_TurbModel) :: this
  real(8), dimension(0:i1,0:k1),intent(IN) :: u,w,rho,mu,mui,muk,mut,beta,temp
  real(8),                                intent(IN) :: alpha1,alpha2, alpha3
  integer,                                intent(IN) :: modification,rank,periodic
  real(8),                                intent(OUT):: residual1,residual2, residual3
  real(8), dimension(0:i1,0:k1) :: rho_mod

  !1, our modification, 2, Aupoix modification
  if ((modification == 1) .or. (modification == 2)) then
    rho_mod = rho
  else
    rho_mod = 1.0
  endif

  call this%production_KE(u,w,temp,rho,mu,mut,beta)
  call this%solve_eps_KE(residual2,u,w,rho,mu,mui,muk,mut,rho_mod, &
                         alpha2,modification,rank,periodic)
  call this%solve_k_KE(residual1,u,w,rho,mu,mui,muk,mut,rho_mod, &
                       alpha1,modification,rank,periodic)
end

subroutine solve_k_KE(this,resK,u,w,rho,mu,mui,muk,mut,rho_mod, &
                       alphak,modification,rank,periodic)
  use mod_param, only : k1,i1,imax,kmax,i,k
  use mod_mesh, only : Ru,Rp,dru,drp,top_bcnovalue,bot_bcnovalue
  use mod_math, only : matrixIdir
  implicit none
  class(KE_TurbModel) :: this
  real(8),dimension(0:i1,0:k1), intent(IN) :: u, w, rho,mu,mui,muk,mut,rho_mod
  real(8),                      intent(IN) :: alphak
  integer,                      intent(IN) :: modification,rank,periodic
  real(8),                      intent(OUT):: resK
  real(8), dimension(0:i1,0:k1)      :: dnew,dimpl
  real(8), dimension(imax)           :: a,b,c,rhs
  real(8)                            :: dz
  
  
  resK  = 0.0;  dnew  = 0.0; dimpl = 0.0;

  call advecc(dnew,dimpl,this%k,u,w,rank,periodic,.true.)
  call this%rhs_k_KE(dnew,dimpl,rho) 
  call diffc(dnew,this%k,mu,mui,muk,mut,this%sigmak,rho,modification)

  do k=1,kmax
    do i=1,imax
      if ((modification == 0) .or. (modification == 1)) then
        a(i) = (mui(i-1,k)+0.5*(mut(i,k)+mut(i-1,k))/this%sigmak)/((0.5*(rho_mod(i-1,k)+rho_mod(i,k)))**0.5)
        a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/rho(i,k)/(rho_mod(i,k)**0.5)
        c(i) = (mui(i  ,k)+0.5*(mut(i,k)+mut(i+1,k))/this%sigmak)/((0.5*(rho_mod(i+1,k)+rho_mod(i,k)))**0.5)
        c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/rho(i,k)/(rho_mod(i,k)**0.5)
      else if (modification == 2) then
        a(i) = (mui(i-1,k)+0.5*(mut(i,k)+mut(i-1,k))/this%sigmak)/(0.5*(rho_mod(i-1,k)+rho_mod(i,k)))
        a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/rho(i,k)
        c(i) = (mui(i  ,k)+0.5*(mut(i,k)+mut(i+1,k))/this%sigmak)/(0.5*(rho_mod(i+1,k)+rho_mod(i,k)))
        c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/rho(i,k)
      endif
      b(i) = (rho_mod(i,k)*(-a(i)-c(i)) + dimpl(i,k))
      a(i) = a(i)*rho_mod(i-1,k)
      c(i) = c(i)*rho_mod(i+1,k)
      rhs(i) = dnew(i,k) + ((1-alphak)/alphak)*b(i)*this%k(i,k)
    enddo

     i=1
     b(i) = b(i) + bot_bcnovalue(k)*a(i)                      !symmetry = -1 ; wall = 1 
     rhs(i) = dnew(i,k) + ((1-alphak)/alphak)*b(i)*this%k(i,k)
       
     i=imax
     b(i) = b(i) + top_bcnovalue(k)*c(i)
     rhs(i) = dnew(i,k) + ((1-alphak)/alphak)*b(i)*this%k(i,k)

     call matrixIdir(imax,a,b/alphak,c,rhs)

    do i=1,imax
      resK = resK + ((this%k(i,k) - rhs(i))/(this%k(i,k)+1.0e-20))**2.0
      this%k(i,k) = max(rhs(i), 1.0e-8)
    enddo
  enddo
end

subroutine solve_eps_KE(this,resE,u,w,rho,mu,mui,muk,mut,rho_mod, &
                        alphae,modification,rank,periodic)
  use mod_param, only :k1,i1,imax,kmax,i,k
  use mod_math
  use mod_mesh, only : mesh,Ru,Rp,dru,drp,bot_bcvalue,top_bcvalue
  implicit none
  class(KE_TurbModel) :: this
  real(8),dimension(0:i1,0:k1), intent(IN) :: u, w, rho,mu,mui,muk,mut,rho_mod
  real(8),                                intent(IN) :: alphae
  integer,                                intent(IN) :: modification,rank,periodic
  real(8),                                intent(OUT):: resE
  real(8), dimension(0:i1,0:k1) :: dnew,dimpl
  real(8), dimension(imax)           :: a,b,c,rhs
  real(8) dz
  
  resE  = 0.0; dnew  = 0.0; dimpl = 0.0;

  call advecc(dnew,dimpl,this%eps,u,w,rank,periodic,.true.)
  call this%rhs_eps_KE(dnew,dimpl,rho)  
  call this%diffusion_eps_KE(dnew,this%eps,muk,mut,this%sigmae,rho,modification)

  do k=1,kmax
    do i=1,imax
      a(i) = (mui(i-1,k)+0.5*(mut(i,k)+mut(i-1,k))/this%sigmae)/sqrt(0.5*(rho_mod(i-1,k)+rho_mod(i,k)))
      a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/rho(i,k)/rho_mod(i,k)
      c(i) = (mui(i  ,k)+0.5*(mut(i,k)+mut(i+1,k))/this%sigmae)/sqrt(0.5*(rho_mod(i+1,k)+rho_mod(i,k)))
      c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/rho(i,k)/rho_mod(i,k)
      b(i) = ((-a(i)-c(i))*(rho_mod(i,k)**1.5) + dimpl(i,k)  )  
      a(i) = a(i)*(rho_mod(i-1,k)**1.5)
      c(i) = c(i)*(rho_mod(i+1,k)**1.5)
      rhs(i) = dnew(i,k) + ((1-alphae)/alphae)*b(i)*this%eps(i,k)
    enddo

    i=1
    b(i) = b(i)+bot_bcvalue(k)*a(i)
    rhs(i) = dnew(i,k) - (1.-bot_bcvalue(k))*a(i)*this%eps(i-1,k) + ((1-alphae)/alphae)*b(i)*this%eps(i,k)   !wall with value

    i=imax
    b(i) = b(i)+top_bcvalue(k)*c(i)
    rhs(i) = dnew(i,k) - (1.-top_bcvalue(k))*c(i)*this%eps(i+1,k) + ((1-alphae)/alphae)*b(i)*this%eps(i,k)   !wall with value
      
    call matrixIdir(imax,a,b/alphae,c,rhs)
  
    do i=1,imax
      resE = resE + ((this%eps(i,k) - rhs(i))/(this%eps(i,k)+1.0e-20))**2.0
      this%eps(i,k) = max(rhs(i), 1.0e-8)
    enddo
  enddo
end subroutine solve_eps_KE

subroutine rhs_k_KE(this,putout,dimpl,rho)
  use mod_param, only : k1,i1,kmax,imax,k,i
  implicit none
  class(KE_TurbModel) :: this
  real(8), dimension(0:i1,0:k1), intent(IN) :: rho
  real(8), dimension(0:i1,0:k1), intent(OUT):: putout,dimpl
  
  !k equation  
  do k=1,kmax
    do i=1,imax
      putout(i,k) = putout(i,k)+(this%Pk(i,k)+this%Gk(i,k))/rho(i,k)
      dimpl(i,k)  = dimpl(i,k) + this%eps(i,k)/this%k(i,k) ! note, rho*epsilon/(rho*k), set implicit and divided by density
    enddo
  enddo
end subroutine rhs_k_KE

subroutine rhs_eps_KE(this,putout,dimpl,rho)
  use mod_param, only :k1,i1,imax,kmax,k,i
  implicit none
  class(KE_TurbModel) :: this
  real(8), dimension(0:i1,0:k1), intent(IN) :: rho
  real(8), dimension(0:i1,0:k1), intent(OUT):: putout,dimpl
  
  do k=1,kmax
    do i=1,imax
      putout(i,k) = putout(i,k) +(this%ce1*this%f1(i,k)*this%Pk(i,k)/this%Tt(i,k) &
                                + this%ce1*this%f1(i,k)*this%Gk(i,k)/this%Tt(i,k) )/rho(i,k)
      dimpl(i,k)  = dimpl(i,k)  + this%ce2*this%f2(i,k)/this%Tt(i,k)   ! note, ce2*f2*rho*epsilon/T/(rho*epsilon), set implicit and divided by density
    enddo
  enddo
end subroutine rhs_eps_KE

subroutine diffusion_eps_KE(this,putout,putin,muk,mut,sigma,rho,modification)
  use mod_param, only : k1,i1,kmax,imax,k,i
  use mod_mesh,  only : dzw,dzp
  implicit none
  class(KE_TurbModel) :: this
  real(8), dimension(0:i1, 0:k1), intent(IN)  :: putin, muk, mut, rho
  real(8),                        intent(IN)  :: sigma
  integer,                        intent(IN)  :: modification
  real(8), dimension(0:i1, 0:k1), intent(OUT) :: putout
  integer :: kp,km
  real(8) :: difcp,difcm
  
  if ((modification == 1) .or. (modification == 2)) then       ! Inverse SLS  and Aupoix
    do k=1,kmax
      kp=k+1
      km=k-1
      do i=1,imax
        difcp = (muk(i,k ) + 0.5*(mut(i,k)+mut(i,kp))/sigma)/sqrt(0.5*(rho(i,k)+rho(i,kp)))
        difcm = (muk(i,km) + 0.5*(mut(i,k)+mut(i,km))/sigma)/sqrt(0.5*(rho(i,k)+rho(i,km)))
        putout(i,k) = putout(i,k) + 1.0/rho(i,k)/rho(i,k)*( &
          (     difcp * ((rho(i,kp)**1.5)*putin(i,kp)-(rho(i,k )**1.5)*putin(i,k ))/dzp(k) &
               -difcm * ((rho(i,k )**1.5)*putin(i,k )-(rho(i,km)**1.5)*putin(i,km))/dzp(km)&
          )/dzw(k)   )
      enddo
    enddo
  else                               ! Standard
    do k=1,kmax
      kp=k+1
      km=k-1
      do i=1,imax
        putout(i,k) = putout(i,k) + 1.0/rho(i,k)*( &
          ( (muk(i,k ) + 0.5*(mut(i,k)+mut(i,kp))/sigma)*(putin(i,kp)-putin(i,k ))/dzp(k) &
          - (muk(i,km) + 0.5*(mut(i,k)+mut(i,km))/sigma)*(putin(i,k )-putin(i,km))/dzp(km)&
          )/dzw(k)   )
      enddo
    enddo
  endif
end subroutine diffusion_eps_KE


end module