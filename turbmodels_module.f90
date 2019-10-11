module mod_turbmodels
  implicit none


  !************************!
  !     Abstract class     !
  !************************!
  
  type, abstract, public :: TurbModel
  integer i1, k1
  contains
    procedure(init), deferred :: init
    procedure(set_mut), deferred :: set_mut
    procedure(advance_turb), deferred :: advance_turb
  end type TurbModel

  interface
    subroutine init(this)
      import :: TurbModel
      class(TurbModel) :: this
    end subroutine init
    subroutine set_mut(this, u, w, rho, mu, mut)
      import :: TurbModel
      class(TurbModel) :: this
      real(8), dimension(:,:), intent(IN) :: u, w, rho, mu
      real(8), dimension(:,:), intent(OUT) :: mut
    end subroutine set_mut
    subroutine advance_turb( this, u, w, rho, mu )
      import :: TurbModel
      class(TurbModel) :: this
      real(8), dimension(:,:), intent(IN) :: u, w, rho, mu
    end subroutine advance_turb
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
  end type Laminar_TurbModel

  ! !************************!
  ! !         SA class       !
  ! !************************!
  
  ! type, extends(TurbModel), public :: SA_TurbModel
  ! contains
  !   procedure :: init => init_SA
  !   procedure :: init_mem => init_mem_SA
  !   procedure :: set_mut => set_mut_SA
  !   procedure :: advance => advance_SA
  ! end type SA_TurbModel


  ! !************************!
  ! !         KE class       !
  ! !************************!
  
  ! type,abstract,extends(TurbModel), public :: KE_TurbModel
  ! contains
  !   procedure(init), deferred :: init_KE
  !   procedure(set_mut), deferred :: set_mut_KE
  !   procedure(advance), deferred :: advance_KE
  ! end type KE_TurbModel

  ! interface
  !   subroutine init_KE(this)
  !     import :: KE_TurbModel
  !     class(KE_TurbModel) :: this
  !   end subroutine init_KE
  !   subroutine set_mut_KE(this, u, w, rho, mu, mut)
  !     import :: KE_TurbModel
  !     class(KE_TurbModel) :: this
  !     real(8), dimension(:,:), intent(IN) :: u, w, rho, mu
  !     real(8), dimension(:,:), intent(OUT) :: mut
  !   end subroutine set_mut_KE
  !   subroutine advance_KE( this, u, w, rho, mu )
  !     import :: KE_TurbModel
  !     class(KE_TurbModel) :: this
  !     real(8), dimension(:,:), intent(IN) :: u, w, rho, ekm
  !   end subroutine advance_KE
  ! end interface

  ! !************************!
  ! !         MK class       !
  ! !************************!

  ! type, extends(KE_TurbModel), public :: MK_TurbModel
  ! contains
  !   procedure :: init => init_MK
  !   procedure :: init_mem => init_mem_MK
  !   procedure :: set_mut => set_mut_MK
  !   procedure :: advance => advance_MK
  ! end type MK_TurbModel

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

  ! !************************!
  ! !         SST class      !
  ! !************************!

  ! type, extends(TurbModel), public :: SST_TurbModel
  ! contains
  !   procedure :: init => init_SST
  !   procedure :: init_mem => init_mem_SST
  !   procedure :: set_mut => set_mut_SST
  !   procedure :: advance => advance_SST
  ! end type SST_TurbModel



contains

subroutine init_laminar(this)
    class(Laminar_TurbModel) :: this
end subroutine init_laminar

subroutine init_mem_laminar(this)
    class(Laminar_TurbModel) :: this
end subroutine init_mem_laminar

subroutine set_mut_laminar(this, u, w, rho, mu, mut)
    implicit none
    class(Laminar_TurbModel) :: this
    real(8), dimension(:,:), intent(IN) :: u, w, rho, mu
    real(8), dimension(:,:), intent(OUT) :: mut
end subroutine set_mut_laminar

subroutine advance_laminar(this, u, w, rho, mu)
    class(Laminar_TurbModel) :: this
    real(8), dimension(:,:), intent(IN) :: u, w, rho, mu
end subroutine advance_laminar


end module mod_turbmodels
