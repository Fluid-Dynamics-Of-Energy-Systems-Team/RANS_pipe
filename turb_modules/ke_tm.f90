module ke_tm
  use mod_turbmodels
  implicit none

!****************************************************************************************

  !************************!
  !         KE class       !
  !************************!
  
  type,abstract,extends(TurbModel), public :: KE_TurbModel
  real(8), dimension(:,:), allocatable :: Gk,Tt,f1,f2,fmu,Lh,fv2
  real(8) :: sigmak,sigmae,ce1,ce2,cmu
  contains

    procedure(set_mut_KE), deferred :: set_mut
    procedure(advance_KE), deferred :: advance_turb
    procedure(set_bc_KE), deferred :: set_bc
    procedure(set_constants), deferred :: set_constants
    procedure :: init => init_KE
    procedure :: rhs_k_KE
    procedure :: rhs_eps_KE
    procedure :: diffusion_eps_KE
    procedure :: solve_eps_KE
    procedure :: solve_k_KE
    procedure :: init_mem_KE
    procedure :: init_sol => init_sol_KE

  end type KE_TurbModel

  interface
    subroutine set_constants(this)
      import :: KE_TurbModel
      class(KE_TurbModel) :: this
    end subroutine set_constants
    subroutine set_mut_KE(this,u,w,rho,mu,mui,walldist,Rp,dRp,dru,dz,mut)
      import :: KE_TurbModel
      class(KE_TurbModel) :: this
      real(8), dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,mu,mui
      real(8), dimension(1:this%imax),        intent(IN) :: walldist
      real(8), dimension(0:this%i1),          intent(IN) :: Rp,dRp,dru
      real(8),                                intent(IN) :: dz
      real(8), dimension(0:this%i1,0:this%k1),intent(OUT):: mut
    end subroutine set_mut_KE
    subroutine advance_KE(this,u,w,rho,mu,mui,muk,mut,beta,temp, &
                             Ru,Rp,dru,drp,dz,walldist,          &
                             alpha1,alpha2,alpha3,               &
                             modification,rank,centerBC,periodic,&
                             residual1,residual2,residual3)
      import :: KE_TurbModel
      class(KE_TurbModel) :: this
      real(8), dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,mu,mui,muk,mut,beta,temp
      real(8), dimension(0:this%i1),          intent(IN) :: Ru,Rp,dru,drp
      real(8), dimension(1:this%i1),          intent(IN) :: walldist
      real(8),                                intent(IN) :: dz,alpha1,alpha2,alpha3
      integer,                                intent(IN) :: modification,rank,centerBC,periodic
      real(8),                                intent(OUT):: residual1,residual2,residual3
    end subroutine advance_KE
    subroutine set_bc_KE(this,mu,rho,walldist,centerBC,periodic,rank,px)
      import :: KE_TurbModel
      class(KE_TurbModel) :: this
      real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: rho,mu
      real(8),dimension(1:this%imax),        intent(IN) :: walldist
      integer,                               intent(IN) :: centerBC,periodic, rank, px
    end subroutine set_bc_KE
  end interface

contains

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
  class(KE_TurbModel) :: this
  integer :: i
  do i=1,this%imax
    this%k(i,:)  =0.1
    this%eps(i,:)=1.0
    this%v2(i,:) =2./3.*this%k(i,:)
  enddo
end subroutine init_sol_KE


subroutine init_mem_KE(this)
  class(KE_TurbModel) :: this
  allocate(this%eps(0:this%i1,0:this%k1),this%k (0:this%i1,0:this%k1), &
           this%Gk (0:this%i1,0:this%k1),this%Pk(0:this%i1,0:this%k1), &
           this%f1 (0:this%i1,0:this%k1),this%f2(0:this%i1,0:this%k1), &
           this%fmu(0:this%i1,0:this%k1),this%Tt(0:this%i1,0:this%k1), &
           this%v2 (0:this%i1,0:this%k1),                              &
           this%fv2(this%imax,this%kmax),this%Lh(this%imax,this%kmax))
end subroutine init_mem_KE

subroutine solve_k_KE(this,resK,u,w,rho,mu,mui,muk,mut,rho_mod, &
                       Ru,Rp,dru,drp,dz, &
                       alphak,modification,rank,centerBC,periodic)
  implicit none
  class(KE_TurbModel) :: this
  real(8),dimension(0:this%i1,0:this%k1), intent(IN) :: u, w, rho,mu,mui,muk,mut,rho_mod
  real(8),dimension(0:this%i1),           intent(IN) :: dru,ru,rp,drp
  real(8),                                intent(IN) :: dz, alphak
  integer,                                intent(IN) :: modification,rank,centerBC,periodic
  real(8),                                intent(OUT):: resK
  real(8), dimension(0:this%i1,0:this%k1) :: dnew,dimpl
  real(8), dimension(this%imax)           :: a,b,c,rhs
  integer                                 :: i,k

  resK  = 0.0;  dnew  = 0.0; dimpl = 0.0;

  call advecc(dnew,dimpl,this%k,u,w,Ru,Rp,dru,dz,this%i1,this%k1,rank,periodic,.true.)
  call this%rhs_k_KE(dnew,dimpl,rho) 
  call diffc(dnew,this%k,mu,mui,muk,mut,this%sigmak,rho,Ru,Rp,dru,dz,rank,modification)

  do k=1,this%kmax
    do i=1,this%imax
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

      b(i) = (rho_mod(i,k)*(-a(i)-c(i)) + dimpl(i,k))/alphak

      a(i) = a(i)*rho_mod(i-1,k)
      c(i) = c(i)*rho_mod(i+1,k)

      rhs(i) = dnew(i,k) + (1-alphak)*b(i)*this%k(i,k)
    enddo

    i=1
    b(i) = b(i) + centerBC*a(i)
           
    i=this%imax
    b(i) = b(i) - (c(i) /alphak)
    rhs(i) = dnew(i,k) + (1-alphak)*b(i)*this%k(i,k)

    call matrixIdir(this%imax,a,b,c,rhs)

    do i=1,this%imax
      resK = resK + ((this%k(i,k) - rhs(i))/(this%k(i,k)+1.0e-20))**2.0
      this%k(i,k) = max(rhs(i), 1.0e-8)
    enddo
  enddo
end


subroutine solve_eps_KE(this,resE,u,w,rho,mu,mui,muk,mut,rho_mod, &
                       Ru,Rp,dru,drp,dz, &
                       alphae,modification,rank,centerBC,periodic)
  implicit none
  class(KE_TurbModel) :: this
  real(8),dimension(0:this%i1,0:this%k1), intent(IN) :: u, w, rho,mu,mui,muk,mut,rho_mod
  real(8),dimension(0:this%i1),           intent(IN) :: dru,ru,rp,drp
  real(8),                                intent(IN) :: dz, alphae
  integer,                                intent(IN) :: modification,rank,centerBC,periodic
  real(8),                                intent(OUT):: resE
  real(8), dimension(0:this%i1,0:this%k1) :: dnew,dimpl
  real(8), dimension(this%imax)           :: a,b,c,rhs
  integer                                 :: i,k
  
  resE  = 0.0; dnew  = 0.0; dimpl = 0.0;

  call advecc(dnew,dimpl,this%eps,u,w,Ru,Rp,dru,dz,this%i1,this%k1,rank,periodic,.true.)
  call this%rhs_eps_KE(dnew,dimpl,rho)  
  call this%diffusion_eps_KE(dnew,this%eps,muk,mut,this%sigmae,rho,dz,modification)

  do k=1,this%kmax
    do i=1,this%imax
                  
      a(i) = (mui(i-1,k)+0.5*(mut(i,k)+mut(i-1,k))/this%sigmae)/sqrt(0.5*(rho_mod(i-1,k)+rho_mod(i,k)))
      a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/rho(i,k)/rho_mod(i,k)
  
      c(i) = (mui(i  ,k)+0.5*(mut(i,k)+mut(i+1,k))/this%sigmae)/sqrt(0.5*(rho_mod(i+1,k)+rho_mod(i,k)))
      c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/rho(i,k)/rho_mod(i,k)
  
      b(i) = ((-a(i)-c(i))*(rho_mod(i,k)**1.5) + dimpl(i,k)  )/alphae
  
      a(i) = a(i)*(rho_mod(i-1,k)**1.5)
      c(i) = c(i)*(rho_mod(i+1,k)**1.5)
  
      rhs(i) = dnew(i,k) + (1-alphae)*b(i)*this%eps(i,k)
    enddo
    i=1
    if (centerBC.eq.-1) then
      rhs(i) = dnew(i,k) - a(i)*this%eps(i-1,k) + (1-alphae)*b(i)*this%eps(i,k)
    else
      b(i)=b(i)+a(i)
    endif

    i=this%imax
    rhs(i) = dnew(i,k) - c(i)*this%eps(i+1,k) + (1-alphae)*b(i)*this%eps(i,k)
  
    call matrixIdir(this%imax,a,b,c,rhs)
  
    do i=1,this%imax
      resE = resE + ((this%eps(i,k) - rhs(i))/(this%eps(i,k)+1.0e-20))**2.0
      this%eps(i,k) = max(rhs(i), 1.0e-8)
  
    enddo
  enddo
end subroutine solve_eps_KE

subroutine rhs_k_KE(this,putout,dimpl,rho)
  implicit none
  class(KE_TurbModel) :: this
  real(8), dimension(0:this%i1,0:this%k1), intent(IN) :: rho
  real(8), dimension(0:this%i1,0:this%k1), intent(OUT):: putout,dimpl
  integer ib,ie,kb,ke,i,k
  
  ib = 1
  ie = this%i1-1

  kb = 1
  ke = this%k1-1

  do k=kb,ke
    do i=ib,ie
      !k equation
      putout(i,k) = putout(i,k) + ( this%Pk(i,k) + this%Gk(i,k) )/rho(i,k)
      dimpl(i,k)  = dimpl(i,k) + this%eps(i,k)/this%k(i,k) ! note, rho*epsilon/(rho*k), set implicit and divided by density
    enddo
  enddo
end subroutine rhs_k_KE

subroutine rhs_eps_KE(this,putout,dimpl,rho)
  implicit none
  class(KE_TurbModel) :: this
  real(8), dimension(0:this%i1,0:this%k1), intent(IN) :: rho
  real(8), dimension(0:this%i1,0:this%k1), intent(OUT):: putout,dimpl
  integer ib,ie,kb,ke,i,k 
 
  ib = 1
  ie = this%i1-1

  kb = 1
  ke = this%k1-1

  do k=kb,ke
    do i=ib,ie
      !epsilon equation
      putout(i,k) = putout(i,k) +(this%ce1*this%f1(i,k)*this%Pk(i,k)/this%Tt(i,k) &
                  + this%ce1*this%f1(i,k)*this%Gk(i,k)/this%Tt(i,k) )/rho(i,k)
      dimpl(i,k)  = dimpl(i,k)  + this%ce2*this%f2(i,k)/this%Tt(i,k)   ! note, ce2*f2*rho*epsilon/T/(rho*epsilon), set implicit and divided by density
    enddo
  enddo
end subroutine rhs_eps_KE

subroutine diffusion_eps_KE(this,putout,putin,muk,mut,sigma,rho,dz,modification)
  implicit none
  class(KE_TurbModel) :: this
  real(8), dimension(0:this%i1, 0:this%k1), intent(IN) :: putin, muk, mut, rho
  real(8),                                  intent(IN) :: sigma, dz
  integer,                                  intent(IN) :: modification
  real(8), dimension(0:this%i1, 0:this%k1), intent(OUT):: putout
  integer i,k,kp,km
  real(8) difcp,difcm

  if ((modification == 1) .or. (modification == 2)) then       ! Inverse SLS  and Aupoix
    do k=1,this%k1-1
      kp=k+1
      km=k-1
      do i=1,this%i1-1
        difcp = (muk(i,k ) + 0.5*(mut(i,k)+mut(i,kp))/sigma)/sqrt(0.5*(rho(i,k)+rho(i,kp)))
        difcm = (muk(i,km) + 0.5*(mut(i,k)+mut(i,km))/sigma)/sqrt(0.5*(rho(i,k)+rho(i,km)))
        putout(i,k) = putout(i,k) + 1.0/rho(i,k)/rho(i,k)*( &
          (     difcp * ((rho(i,kp)**1.5)*putin(i,kp)-(rho(i,k )**1.5)*putin(i,k )) &
          -difcm * ((rho(i,k )**1.5)*putin(i,k )-(rho(i,km)**1.5)*putin(i,km))  )/(dz*dz)   )
      enddo
    enddo
  else                               ! Standard
    do k=1,this%k1-1
      kp=k+1
      km=k-1
      do i=1,this%i1-1
        putout(i,k) = putout(i,k) + 1.0/rho(i,k)*( &
          ( (muk(i,k ) + 0.5*(mut(i,k)+mut(i,kp))/sigma)*(putin(i,kp)-putin(i,k )) &
          - (muk(i,km) + 0.5*(mut(i,k)+mut(i,km))/sigma)*(putin(i,k )-putin(i,km)))/(dz*dz)   )
      enddo
    enddo
  endif
end subroutine diffusion_eps_KE


end module