module mod_tdm
  implicit none

!****************************************************************************************

  !************************!
  !     Abstract class     !
  !************************!
  
  type, abstract, public :: TurbDiffModel
  integer i1,k1,imax,kmax
  character(len=3)                      :: name
  real(8), dimension(:,:), allocatable :: Pkt,kt,epst, yp
  real(8), dimension(:),   allocatable :: alphatin, Pktin
  contains
    procedure(init_tdm),       deferred :: init
    procedure(set_alphat_tdm), deferred :: set_alphat
  end type TurbDiffModel

  interface
    subroutine init_tdm(this)
      import :: TurbDiffModel
      class(TurbDiffModel) :: this
    end subroutine init_tdm

    subroutine set_alphat_tdm(this,u,w,rho,temp,mu,mui,lam_cp,mut,alphat)
      import :: TurbDiffModel
      class(TurbDiffModel) :: this
      real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,temp,mu,mui,lam_cp, mut
      real(8),dimension(0:this%i1,0:this%k1),intent(OUT):: alphat
    end subroutine set_alphat_tdm

  end interface

  

class(TurbDiffModel), allocatable :: turbdiff_model
!****************************************************************************************

  !************************!
  ! Constant Prandtl class !
  !************************!

  type, extends(TurbDiffModel), public :: CPrt_TurbDiffModel
  real(8) :: Prt
  contains
    procedure :: init => init_constprt
    procedure :: set_alphat  => set_alphat_constprt
  end type CPrt_TurbDiffModel

contains

  type(CPrt_TurbDiffModel) function init_CPrt_TurbDiffModel(i1,k1,imax,kmax,name, Prt)
    integer, intent(in) :: i1,k1,imax,kmax
    character(len=2), intent(IN) :: name
    real(8), intent(in) :: Prt
    init_CPrt_TurbDiffModel%name=name
    init_CPrt_TurbDiffModel%i1 = i1
    init_CPrt_TurbDiffModel%k1 = k1
    init_CPrt_TurbDiffModel%imax = imax
    init_CPrt_TurbDiffModel%kmax = kmax
    init_CPrt_TurbDiffModel%Prt = Prt
  end function init_CPrt_TurbDiffModel

  subroutine init_constprt(this)
    class(CPrt_TurbDiffModel) :: this
  end subroutine init_constprt

  subroutine set_alphat_constprt(this,u,w,rho,temp,mu,mui,lam_cp,mut,alphat)
    class(CPrt_TurbDiffModel) :: this
    real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,temp,mu,mui,lam_cp, mut
    real(8),dimension(0:this%i1,0:this%k1),intent(OUT):: alphat
    alphat = mut/this%Prt
  end subroutine set_alphat_constprt

end module mod_tdm
