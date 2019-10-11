module mod_turbmodels
  implicit none


  !************************!
  !     Abstract class     !
  !************************!
  
  type, abstract, public :: TurbModel
  integer i1, k1
  real(8), dimension(:,:), private, allocatable :: Pk
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

  !************************!
  !         SA class       !
  !************************!
  
  type, extends(TurbModel), public :: SA_TurbModel
  real(8), dimension(:,:), private, allocatable :: nuSAnew
  contains
    procedure :: init => init_SA
    procedure :: init_mem => init_mem_SA
    procedure :: set_mut => set_mut_SA
    procedure :: advance_turb => advance_SA
  end type SA_TurbModel

  !************************!
  !         SST class      !
  !************************!

  type, extends(TurbModel), public :: SST_TurbModel
  contains
    procedure :: init => init_SST
    procedure :: init_mem => init_mem_SST
    procedure :: set_mut => set_mut_SST
    procedure :: advance_turb => advance_SST
  end type SST_TurbModel



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





contains


  !************************!
  !    Laminar routines    !
  !************************!

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

  !************************!
  !      SA routines      !
  !************************!

subroutine init_SA(this)
    implicit none
    class(SA_TurbModel) :: this

end subroutine init_SA
subroutine init_mem_SA(this)
    implicit none
    class(SA_TurbModel) :: this
    allocate(this%Pk(0:this%i1, 0:this%k1),this%nuSAnew(0:this%i1, 0:this%k1))

end subroutine init_mem_SA

subroutine set_mut_SA(this, u, w, rho, mu, mut)
    implicit none
    class(SA_TurbModel) :: this
    real(8), dimension(:,:), intent(IN) :: u, w, rho, mu
    real(8), dimension(:,:), intent(OUT) :: mut

end subroutine set_mut_SA
subroutine advance_SA(this, u, w, rho, mu)
    implicit none
    class(SA_TurbModel) :: this
    real(8), dimension(:,:), intent(IN) :: u, w, rho, mu
end subroutine advance_SA


  !************************!
  !      SST routines      !
  !************************!

subroutine init_SST(this)
    implicit none
    class(SST_TurbModel) :: this

end subroutine init_SST
subroutine init_mem_SST(this)
    implicit none
    class(SST_TurbModel) :: this

end subroutine init_mem_SST

subroutine set_mut_SST(this, u, w, rho, mu, mut)
    implicit none
    class(SST_TurbModel) :: this
    real(8), dimension(:,:), intent(IN) :: u, w, rho, mu
    real(8), dimension(:,:), intent(OUT) :: mut

end subroutine set_mut_SST
subroutine advance_SST(this, u, w, rho, mu)
    implicit none
    class(SST_TurbModel) :: this
    real(8), dimension(:,:), intent(IN) :: u, w, rho, mu
end subroutine advance_SST
end module mod_turbmodels
