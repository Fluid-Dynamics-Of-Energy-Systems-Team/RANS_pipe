module mod_tdm
  implicit none

!****************************************************************************************

  !************************!
  !     Abstract class     !
  !************************!
  
  type, abstract, public :: TurbDiffModel
  integer i1,k1,imax,kmax
  character(len=3)                      :: name
  real(8), dimension(:,:), allocatable :: Pkt,kt,epst
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

    subroutine set_alpha_tdm(this,u,w,rho,temp,mu,mui,ekh,alphat)
      import :: TurbDiffModel
      class(TurbDiffModel) :: this
      real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,temp,mu,mui,ekh
      real(8),dimension(0:this%i1,0:this%k1),intent(OUT):: alphat
    end subroutine set_alpha_tdm

!    subroutine set_alphat_tdm(this,mut,lam_cp,mu,alphat)
!      import :: TurbDiffModel
!      class(TurbDiffModel) :: this
!      real(8), dimension(0:this%i1,0:this%k1),intent(IN) :: mut, lam_cp, mu
!      real(8), dimension(0:this%i1,0:this%k1),intent(OUT):: alphat
!    end subroutine set_alphat_tdm


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
  subroutine init_constprt(this)
    class(CPrt_TurbDiffModel) :: this
  end subroutine init_constprt

  subroutine set_alphat_constprt(this,u,w,rho,temp,mu,mui,ekh,alphat)
    use mod_tm, only : turb_model
    class(CPrt_TurbDiffModel) :: this
    real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,temp,mu,mui,ekh
    real(8),dimension(0:this%i1,0:this%k1),intent(OUT):: alphat
    alphat = turb_model%mut/this%Prt
  end subroutine set_alphat_constprt

!  subroutine set_alphat_constprt(this,mut,lam_cp,mu, alphat)
!    class(CPrt_TurbDiffModel) :: this
!    real(8), dimension(0:this%i1,0:this%k1),intent(IN) :: mut, lam_cp,mu
!    real(8), dimension(0:this%i1,0:this%k1),intent(OUT):: alphat
!    alphat = mut/this%Prt
!  end subroutine set_alphat_constprt

end module mod_tdm
