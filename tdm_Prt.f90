module Prt_tdm
  use mod_tdm
  implicit none

!****************************************************************************************

  !************************!
  !    Var Prt class     !
  !************************!
  
  type,abstract,extends(TurbDiffModel), public :: VarPrt_TurbDiffModel
  real(8), dimension(:,:), allocatable :: Prt,Pr,mut_mu
  
  contains
    procedure :: init => init_VarPrt
    procedure :: set_alphat  => set_alphat_VarPrt
    procedure :: init_mem_VarPrt

  end type VarPrt_TurbDiffModel

contains
  subroutine init_VarPrt(this)
    implicit none
    class(VarPrt_TurbDiffModel) :: this  
    call this%init_mem_VarPrt()
  end subroutine init_VarPrt

subroutine init_mem_VarPrt(this)
    implicit none
    class(VarPrt_TurbDiffModel) :: this
    allocate(this%Prt(0:this%i1,0:this%k1),this%Pr(0:this%i1,0:this%k1),this%mut_mu(0:this%i1,0:this%k1))
end subroutine init_mem_VarPrt

subroutine set_alphat_VarPrt(this,Prtmodel,mu,lam_cp,alphat)
  implicit none
  class(VarPrt_TurbDiffModel) :: this
  real(8), dimension(0:this%i1,0:this%k1), intent(IN) :: mu, lam_cp !lam_cp==ekh in the code
  real(8), dimension(0:this%i1,0:this%k1), intent(OUT):: alphat

  integer  i,k,Prtmodel
  real(8) :: c1,c2,c3,c4,gam,Agam,Prtinf,PeT,Atang

  ! constant for kays
  c1 = 0.5882
  c2 = 0.228
  c3 =-0.0441
  c4 =-5.165
  ! for Irrenfried and steiner
  Prtinf = 1.0  !turbulent prandlt in the free stream (assumption)
  ! const Tang
  Atang = 15.0

  do k=1,this%kmax

    do i=1,this%imax
      Pr(i,k)= mu(i,k)/lam_cp(i,k)
      mut_mu(i,k)= mut(i,k)/mu(i,k)
      
      ! Approximation of the turbulent Prandlt number 
      !W.M. Kays, M.E. Crawford, Convective Heat and Mass Transfer, McGraw-Hill Inc, New York, 1993.)
      if(Prtmodel.eq.1) then
        Prt(i,k)= 1/(c1+c2*mut_mu(i,k)+c3*(mut_mu(i,k)**2.0)*(1-exp(c4/mut_mu(i,k))))

      ! Approximation of the turbulent Prandlt number (C. Irrenfried, H. Steiner IJHFF 2017)
      elseif (Prtmodel.eq.2) then
        PeT      = mut_mu(i,k)*Pr(i,k)
        gam      = 1.0/(Prtinf+0.1*(Pr(i,k)**0.83))
        Agam     = ((2/Prtinf)-2*gam)**0.5
        Prt(i,k) = (gam+3.0*PeT*Agam-((3.0*PeT)**2.0)*(1-exp(-Agam/(3.0*PeT))))**(-1.0)

      ! Approximation of the turbulent Prandlt number (W. Kays Turb. Pr, Where are we? 1992)
      elseif (Prtmodel.eq.3) then
        PeT      = mut_mu(i,k)*Pr(i,k)
        if (mut_mu(i,k).lt.0.2) then
          Prt(i,k) = 1.07 
        else
          Prt(i,k) = (2.0/PeT)+0.85
        end

      ! Approximation of the turbulent Prandlt number (Tang et al IJHMT 2016) 
      elseif (Prtmodel.eq.4) then
        if (mut_mu(i,k).lt.0.2) then
          Prt(i,k) = 1.0 
        elseif (mut_mu(i,k).le.10.) then
          Prt(i,k) = 0.85 + Pr(i,k)/Atang
        else
          Prt(i,k) =0.85
        end
      end
      !!!!!!!!!!!!!!!!!
      alphat(i,k)= mut(i,k)/Prt(i,k)
    enddo
  enddo
end subroutine set_alphat_VarPrt





end module
