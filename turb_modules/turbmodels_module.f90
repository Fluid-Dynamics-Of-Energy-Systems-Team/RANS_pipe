module mod_turbmodels
  implicit none


  !************************!
  !     Abstract class     !
  !************************!
  
  type, abstract, public :: TurbModel
  integer i1,k1,imax,kmax
  real(8), dimension(:,:), allocatable :: nuSA,Pk,om,k,bF1,bF2,eps,v2

  contains
    procedure(init_tm), deferred :: init
    procedure(set_mut_tm), deferred :: set_mut
    procedure(advance_turb_tm), deferred :: advance_turb
    procedure(set_bc_tm), deferred :: set_bc
    procedure(init_sol_tm), deferred :: init_sol
    procedure :: set_mut_bc

  end type TurbModel

  interface
    subroutine init_tm(this)
      import :: TurbModel
      class(TurbModel) :: this
    end subroutine init_tm
    subroutine init_sol_tm(this)
      import :: TurbModel
      class(TurbModel) :: this
    end subroutine init_sol_tm
    
    subroutine set_mut_tm(this,u,w,rho,mu,mui,walldist,Rp,dRp,dru,dz,mut)
      import :: TurbModel
      class(TurbModel) :: this
      real(8), dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,mu,mui
      real(8), dimension(1:this%imax),        intent(IN) :: walldist
      real(8), dimension(0:this%i1),          intent(IN) :: Rp,dRp,dru
      real(8),                                intent(IN) :: dz
      real(8), dimension(0:this%i1,0:this%k1),intent(OUT):: mut
    end subroutine set_mut_tm
    subroutine advance_turb_tm(this,u,w,rho,mu,mui,muk,mut,beta,temp,&
                             Ru,Rp,dru,drp,dz,walldist,              &
                             alpha1,alpha2,alpha3,                   &
                             modification,rank,centerBC,periodic,    &
                             residual1, residual2, residual3)
      import :: TurbModel
      class(TurbModel) :: this
      real(8), dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,mu,mui,muk,mut,beta,temp
      real(8), dimension(0:this%i1),          intent(IN) :: Ru,Rp,dru,drp
      real(8), dimension(1:this%i1),          intent(IN) :: walldist
      real(8),                                intent(IN) :: dz,alpha1,alpha2,alpha3
      integer,                                intent(IN) :: modification,rank,centerBC,periodic
      real(8),                                intent(OUT):: residual1,residual2,residual3
    end subroutine advance_turb_tm
    subroutine set_bc_tm(this,mu,rho,walldist,centerBC,periodic,rank,px)
      import :: TurbModel
      class(TurbModel) :: this
      real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: rho,mu
      real(8),dimension(1:this%imax),        intent(IN) :: walldist
      integer,                               intent(IN) :: centerBC,periodic, rank, px
    end subroutine set_bc_tm

  end interface


  !************************!
  !      Laminar class     !
  !************************!
  
  type, extends(TurbModel), public :: Laminar_TurbModel
  contains
    procedure :: init => init_laminar
    procedure :: init_mem => init_mem_laminar
    procedure :: set_mut  => set_mut_laminar
    procedure :: advance_turb => advance_laminar
    procedure :: init_sol => init_sol_laminar
    procedure :: set_bc => set_bc_laminar
  end type Laminar_TurbModel



contains

subroutine set_mut_bc(this,mut,periodic, px, rank)
  class(TurbModel) :: this
  integer,                                 intent(IN) :: periodic,px,rank
  real(8), dimension(0:this%i1,0:this%k1), intent(OUT) :: mut
  real(8), dimension(0:this%i1) :: tmp
  
  mut(this%i1,:) = -mut(this%imax,:)
  mut(0,:)  =  mut(1,:)

  call shiftf(mut,tmp,rank); mut(:,0)  = tmp(:);
  call shiftb(mut,tmp,rank); mut(:,this%k1) = tmp(:);

  ! if ((periodic.ne.1).and.(rank.eq.0)) then
  !   mut(:,0) = mutin(:)
  ! endif
  if ((periodic.ne.1).and.(rank.eq.px-1)) then
    mut(:,this%k1) = 2.*mut(:,this%kmax)-mut(:,this%kmax-1)
  endif

end subroutine
!   !************************!
!   !    Laminar routines    !
!   !************************!

subroutine init_laminar(this)
    class(Laminar_TurbModel) :: this
end subroutine init_laminar

subroutine set_bc_laminar(this,mu,rho,walldist,centerBC,periodic,rank,px)
      class(Laminar_TurbModel) :: this
      real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: rho,mu
      real(8),dimension(1:this%imax),        intent(IN) :: walldist
      integer,                               intent(IN) :: centerBC,periodic, rank, px
end subroutine set_bc_laminar

subroutine init_sol_laminar(this)
    class(Laminar_TurbModel) :: this
end subroutine init_sol_laminar

subroutine init_mem_laminar(this)
    class(Laminar_TurbModel) :: this
end subroutine init_mem_laminar

subroutine set_mut_laminar(this,u,w,rho,mu,mui,walldist,Rp,dRp,dru,dz,mut)
      class(Laminar_TurbModel) :: this
      real(8), dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,mu,mui
      real(8), dimension(1:this%imax),        intent(IN) :: walldist
      real(8), dimension(0:this%i1),          intent(IN) :: Rp,dRp,dru
      real(8),                                intent(IN) :: dz
      real(8), dimension(0:this%i1,0:this%k1),intent(OUT):: mut
      mut = 0
end subroutine set_mut_laminar

subroutine advance_laminar(this,u,w,rho,mu,mui,muk,mut,beta,temp,&
                             Ru,Rp,dru,drp,dz,walldist,              &
                             alpha1,alpha2,alpha3,                   &
                             modification,rank,centerBC,periodic,    &
                             residual1, residual2, residual3)
      class(Laminar_TurbModel) :: this
      real(8), dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,mu,mui,muk,mut,beta,temp
      real(8), dimension(0:this%i1),          intent(IN) :: Ru,Rp,dru,drp
      real(8), dimension(1:this%i1),          intent(IN) :: walldist
      real(8),                                intent(IN) :: dz,alpha1,alpha2,alpha3
      integer,                                intent(IN) :: modification,rank,centerBC,periodic
      real(8),                                intent(OUT):: residual1,residual2,residual3
end subroutine advance_laminar




end module mod_turbmodels
