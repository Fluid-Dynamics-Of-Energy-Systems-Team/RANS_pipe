module Prt_tdm
  use mod_tdm
  implicit none

!****************************************************************************************

  !************************!
  !    Var Prt class     !
  !************************!
  
  type,abstract,extends(TurbDiffModel), public :: VarPrt_KaysCrawford_TurbDiffModel
  real(8), dimension(:,:), allocatable :: Prt,Pr,mut_mu
  
  contains
    procedure :: init => init_VarPrt
    procedure :: set_alphat  => set_alphat_VarPrt
    procedure :: init_mem_VarPrt

  end type VarPrt_KaysCrawford_TurbDiffModel

contains
  subroutine init_VarPrt(this)
    implicit none
    class(VarPrt_KaysCrawford_TurbDiffModel) :: this  
    call this%init_mem_VarPrt()
  end subroutine init_VarPrt

subroutine init_mem_VarPrt(this)
    implicit none
    class(SA_TurbModel) :: this
    allocate(this%Prt(0:this%i1,0:this%k1),this%Pr(0:this%i1,0:this%k1),this%mut_mu(0:this%i1,0:this%k1))
end subroutine init_mem_VarPrt

subroutine set_alphat_VarPrt(this,mu,lam_cp,alphat)
  implicit none
  class(VarPrt_KaysCrawford_TurbDiffModel) :: this
  real(8), dimension(0:this%i1,0:this%k1), intent(IN) :: mu, lam_cp !lam_cp==ekh in the code
  real(8), dimension(0:this%i1,0:this%k1), intent(OUT):: alphat

  integer  i,k
  real(8) :: c1,c2,c3,c4

  ! constant for kays
  c1 = 0.5882
  c2 = 0.228
  c3 =-0.0441
  c4 =-5.165

  do k=1,this%kmax

    do i=1,this%imax
      Pr(i,k)= mu(i,k)/lam_cp(i,k)
      mut_mu(i,k)= mut(i,k)/mu(i,k)
      
      ! Approximation of the turbulent Prandlt number (W. Keys Turb. Pr, Where are we? 1992)
      Prt(i,k)= 1/(c1+c2*mut_mu(i,k)+c3*(mut_mu(i,k)**2.0)*(1-exp(c4/mut_mu(i,k))))
      !!!!!!!!!!!!!!!!!
      alphat(i,k)= mut(i,k)/Prt(i,k)
    enddo
  enddo
end subroutine set_alphat_VarPrt





end module
