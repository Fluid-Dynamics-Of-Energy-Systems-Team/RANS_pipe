module mod_turbmodels
  implicit none


  !************************!
  !     Abstract class     !
  !************************!
  
  type, abstract, public :: TurbModel
  integer i1,k1,imax,kmax
  real(8), dimension(:,:), allocatable :: nuSA,Pk,om,k,bF1,bF2,eps,v2

  contains
    procedure(init), deferred :: init
    procedure(set_mut), deferred :: set_mut
    procedure(advance_turb), deferred :: advance_turb
    procedure(set_bc), deferred :: set_bc

  end type TurbModel

  interface
    subroutine init(this)
      import :: TurbModel
      class(TurbModel) :: this
    end subroutine init
    subroutine set_mut(this,u,w,rho,mu,mui,walldist,dRp,dru,dz,mut)
      import :: TurbModel
      class(TurbModel) :: this
      real(8), dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,mu,mui
      real(8), dimension(1:this%imax),        intent(IN) :: walldist
      real(8), dimension(0:this%i1),          intent(IN) :: dRp,dru
      real(8),                                intent(IN) :: dz
      real(8), dimension(0:this%i1,0:this%k1),intent(OUT):: mut
    end subroutine set_mut
    subroutine advance_turb(this,u,w,rho,mu,mui,muk,mut,beta,temp,&
                             Ru,Rp,dru,drp,dz,walldist,           &
                             alpha1,alpha2,modification,          &
                             rank,centerBC,periodic,              &
                             residual1, residual2)
      import :: TurbModel
      class(TurbModel) :: this
      real(8), dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,mu,mui,muk,mut,beta,temp
      real(8), dimension(0:this%i1),          intent(IN) :: Ru,Rp,dru,drp
      real(8), dimension(1:this%i1),          intent(IN) :: walldist
      real(8),                                intent(IN) :: dz,alpha1,alpha2
      integer,                                intent(IN) :: modification,rank,centerBC,periodic
      real(8),                                intent(OUT):: residual1,residual2
    end subroutine advance_turb
    subroutine set_bc(this,periodic, rank, px)
      import :: TurbModel
      class(TurbModel) :: this
      integer, intent(IN) :: periodic, rank, px
    end subroutine set_bc

  end interface


  ! !************************!
  ! !      Laminar class     !
  ! !************************!
  
  ! type, extends(TurbModel), public :: Laminar_TurbModel
  ! contains
  !   procedure :: init => init_laminar
  !   procedure :: init_mem => init_mem_laminar
  !   procedure :: set_mut  => set_mut_laminar
  !   procedure :: advance_turb => advance_laminar
  ! end type Laminar_TurbModel










  ! !************************!
  ! !         VF class       !
  ! !************************!

  ! type, extends(KE_TurbModel), public :: V2F_TurbModel
  ! contains
  !   procedure :: init => init_V2F
  !   procedure :: init_mem => init_mem_V2F
  !   procedure :: set_mut => set_mut_V2F
  !   procedure :: advance => advance_V2F
  ! end type V2F_TurbModel





contains


!   !************************!
!   !    Laminar routines    !
!   !************************!

! subroutine init_laminar(this)
!     class(Laminar_TurbModel) :: this
! end subroutine init_laminar

! subroutine init_mem_laminar(this)
!     class(Laminar_TurbModel) :: this
! end subroutine init_mem_laminar

! subroutine set_mut_laminar(this, u, w, rho, mu, mut)
!     implicit none
!     class(Laminar_TurbModel) :: this
!     real(8), dimension(:,:), intent(IN) :: u, w, rho, mu
!     real(8), dimension(:,:), intent(OUT) :: mut
! end subroutine set_mut_laminar

! subroutine advance_laminar(this, u, w, rho, mu)
!     class(Laminar_TurbModel) :: this
!     real(8), dimension(:,:), intent(IN) :: u, w, rho, mu
! end subroutine advance_laminar




end module mod_turbmodels
