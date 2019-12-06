module nk_tm
  use ktet_tm
  implicit none

!****************************************************************************************

  !************************!
  !         NK class      !
  !************************!
  !!!!!!  Nagano, Y., and Kim, C. "A two-equation model for heat transport in 
  !         wall turbulent shear flows", Journal of Heat transfer (1988).
  type, extends(KtEt_TurbDiffModel), public :: NK_TurbDiffModel
  contains
    procedure :: set_constants => set_constants_NK
    procedure :: set_alphat => set_alphat_NK
    procedure :: rhs_epst_KtEt => rhs_epst_KtEt_NK
    procedure :: init_w_inflow => init_w_inflow_NK  !MISSING
  end type NK_TurbModel


contains
!****************************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************************


  !************************!
  !       NK routines      !
  !************************!

type(NK_TurbDiffModel) function init_NK_TurbDiffModel(i1,k1,imax,kmax,name)
  integer, intent(in) :: i1,k1,imax,kmax
  character(len=2), intent(IN) :: name
  init_NK_TurbDiffModel%name=name
  init_NK_TurbDiffModel%i1 = i1
  init_NK_TurbDiffModel%k1 = k1
  init_NK_TurbDiffModel%imax = imax
  init_NK_TurbDiffModel%kmax = kmax
end function init_NK_TurbDiffModel

subroutine set_constants_NK(this)
  class(NK_TurbDiffModel) :: this
  this%sigmakt = 1.0
  this%sigmaet = 1.0
  this%clambda = 0.11
  this%cp1 = 1.8
  this%cp2 = 0.72
  this%cd1 = 2.2
  this%cd2 = 0.8
end subroutine set_constants_NK

subroutine set_alphat_NK(this,u,w,rho,temp,mu,mui,ekh,alphat)
  use mod_tm, only : turb_model
  use mod_mesh, only : Rp,dRp,dRu,dz,walldist
  implicit none
  class(NK_TurbDiffModel) :: this
  real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,temp,mu,mui,ekh
  real(8),dimension(0:this%i1,0:this%k1),intent(OUT):: alphat
  integer  im,ip,km,kp,i,k 
  real(8) cfi, sti, utau
  real(8),dimension(0:this%k1) ::   tauw
  real(8),dimension(0:this%i1,0:this%k1) :: Ret, Reeps, yp
  real(8), dimension(0:this%i1,0:this%k1) :: kine, eps, Tt

  eps  = turb_model%eps
  kine = turb_model%k
  Tt   = turb_model%Tt

  do k=1,this%kmax
    km=k-1
    kp=k+1
    tauw(k) = mui(this%imax,k)*0.5*(w(this%imax,km)+w(this%imax,k))/walldist(this%imax)
    utau    = (tauw(k)/rho(imax,k))*0.5
    do i=1,this%imax
      im=i-1
      ip=i+1
      ! yplus hsould be an input so it can be changed from yplus to ystar
      this%yp(i,k) = sqrt(rho(i,k))/mu(i,k)*(walldist(i))*tauw(k)**0.5       

      cfi =  2*tauw(k)/rho(imax,k)/utau                     ! Skin friction
      sti =  Qwall/(rho(imax,k)*cp(imax,k)*utau*temp(imax,k))  ! Stanton number


      this%Ttemp(i,k)   = this%kt(i,k)/(this%epst(i,k)+1.0e-20)
      this%Tmix(i,k)    = (Tt(i,k) * this%Ttemp(i,k) )**0.5


      this%flambda(i,k) =(1 - exp(-(2*sti/cfi)*yp(i,k)*(Pr**0.5)/30.5))**2.0      
      alphat(i,k) = rho(i,k)*this%clambda*this%flambda(i,k)*k(i,k)*((Tt(i,k)*this%Ttemp(i,k))**0.5)

    enddo
  enddo
end subroutine set_alphat_NK

subroutine rhs_epst_KtEt_NK(this,putout,dimpl,temp,rho,mu,ekh,alphat)
  use mod_tm, only : turb_model
  implicit none
  class(NK_TurbDiffModel) :: this
  real(8), dimension(0:this%i1,0:this%k1), intent(IN) :: temp,rho,mu,ekh,alphat
  real(8), dimension(0:this%i1,0:this%k1), intent(OUT):: putout,dimpl
  integer ib,ie,kb,ke,i,k 
  real(8) d2Tdxdr
  real(8), dimension(0:this%i1,0:this%k1) :: Pk, kine, Tt
     ! ce2 is from the k-epsilon model
  ib = 1
  ie = this%i1-1
  Pk   = turb_model%Pk
  kine = turb_model%k
  Tt   = turb_mdel%Tt
  kb = 1
  ke = this%k1-1

  do k=kb,ke
    do i=ib,ie
      d2Tdxdr = (1-this%flambda(i,k))*rho(i,k)*alphat(i,k)*ekh(i,k) * &
             ((((temp(i,k)-temp(i,km))/dz)-((temp(im,k)-temp(im,km))/dz) )/dRu(i))**2

      putout(i,k) = putout(i,k) + (this%cp1*this%epst(i,k)/(this%kt(i,k)+1.0e-20)*this%Pkt(i,k) + &
                   this%cp2*this%epst(i,k)/ kine(i,k)* Pk(i,k) + d2Tdxdr)     /rho(i,k)
      dimpl(i,k)  = dimpl(i,k)  + this%cd1*this%epst(i,k)/(this%kt(i,k)+1.0e-20) + this%cd2/Tt(i,k)           

    enddo
  enddo
end subroutine rhs_epst_KtEt_NK


end module nk_tdm