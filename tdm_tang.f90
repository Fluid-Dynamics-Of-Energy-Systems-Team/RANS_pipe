module Prt_tdm
  use mod_tdm
  implicit none

!****************************************************************************************

  !************************!
  !    Var Prt class     !
  !************************!
  
  type,abstract,extends(TurbDiffModel), public :: VarPrt_Tang_TurbDiffModel
  real(8), dimension(:,:), allocatable :: Prt,Pr,mut_mu
  
  contains
    procedure :: init => init_VarPrt
    procedure :: set_alphat  => set_alphat_VarPrt
    procedure :: init_mem_VarPrt

  end type VarPrt_Tang_TurbDiffModel

contains
  subroutine init_VarPrt(this)
    implicit none
    class(VarPrt_Tang_TurbDiffModel) :: this  
    call this%init_mem_VarPrt()
  end subroutine init_VarPrt

subroutine init_mem_VarPrt(this)
    implicit none
    class(SA_TurbModel) :: this
    allocate(this%Prt(0:this%i1,0:this%k1),this%Pr(0:this%i1,0:this%k1),this%mut_mu(0:this%i1,0:this%k1))
end subroutine init_mem_VarPrt

subroutine set_alphat_VarPrt(this,mu,lam_cp,alphat)
  implicit none
  class(VarPrt_Tang_TurbDiffModel) :: this
  real(8), dimension(0:this%i1,0:this%k1), intent(IN) :: mu, lam_cp !lam_cp==ekh in the code
  real(8), dimension(0:this%i1,0:this%k1), intent(OUT):: alphat

  integer  i,k
  real(8) :: Atang

  ! const 
  Atang = 15.0

  do k=1,this%kmax

    do i=1,this%imax
      Pr(i,k)= mu(i,k)/lam_cp(i,k)
      mut_mu(i,k)= mut(i,k)/mu(i,k)
      
      ! Approximation of the turbulent Prandlt number (Tang et al IJHMT 2016) 
      if (mut_mu(i,k).lt.0.2) then
        Prt(i,k) = 1.0 
      elseif (mut_mu(i,k).le.10.) then
        Prt(i,k) = 0.85 + Pr(i,k)/Atang
      else
        Prt(i,k) =0.85
      end
      !!!!!!!!!!!!!!!!!
      alphat(i,k)= mut(i,k)/Prt(i,k)
    enddo
  enddo
end subroutine set_alphat_VarPrt





end module
