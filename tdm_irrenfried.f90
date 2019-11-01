module Prt_tdm
  use mod_tdm
  implicit none

!****************************************************************************************

  !************************!
  !    Var Prt class     !
  !************************!
  
  type,abstract,extends(TurbDiffModel), public :: VarPrt_Irrenfried_TurbDiffModel
  real(8), dimension(:,:), allocatable :: Prt,Pr,mut_mu
  
  contains
    procedure :: init => init_VarPrt
    procedure :: set_alphat  => set_alphat_VarPrt
    procedure :: init_mem_VarPrt

  end type VarPrt_Irrenfried_TurbDiffModel

contains
  subroutine init_VarPrt(this)
    implicit none
    class(VarPrt_Irrenfried_TurbDiffModel) :: this  
    call this%init_mem_VarPrt()
  end subroutine init_VarPrt

subroutine init_mem_VarPrt(this)
    implicit none
    class(SA_TurbModel) :: this
    allocate(this%Prt(0:this%i1,0:this%k1),this%Pr(0:this%i1,0:this%k1),this%mut_mu(0:this%i1,0:this%k1))
end subroutine init_mem_VarPrt

subroutine set_alphat_VarPrt(this,mu,lam_cp,alphat)
  implicit none
  class(VarPrt_Irrenfried_TurbDiffModel) :: this
  real(8), dimension(0:this%i1,0:this%k1), intent(IN) :: mu, lam_cp !lam_cp==ekh in the code
  real(8), dimension(0:this%i1,0:this%k1), intent(OUT):: alphat

  integer  i,k
  real(8) :: gam,Agam,Prtinf,PeT

  ! for Irrenfried and steiner
  Prtinf = 1.0  !turbulent prandlt in the free stream (assumption)

  do k=1,this%kmax

    do i=1,this%imax
      Pr(i,k)= mu(i,k)/lam_cp(i,k)
      mut_mu(i,k)= mut(i,k)/mu(i,k)
      
      ! Approximation of the turbulent Prandlt number (C. Irrenfried, H. Steiner IJHFF 2017)
      PeT      = mut_mu(i,k)*Pr(i,k)
      gam      = 1.0/(Prtinf+0.1*(Pr(i,k)**0.83))
      Agam     = ((2/Prtinf)-2*gam)**0.5
      Prt(i,k) = (gam+3.0*PeT*Agam-((3.0*PeT)**2.0)*(1-exp(-Agam/(3.0*PeT))))**(-1.0)

      !!!!!!!!!!!!!!!!!
      alphat(i,k)= mut(i,k)/Prt(i,k)
    enddo
  enddo
end subroutine set_alphat_VarPrt





end module
