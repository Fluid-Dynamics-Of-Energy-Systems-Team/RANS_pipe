module abe_tm
  use ke_tm
  implicit none

!****************************************************************************************

  !************************!
  !         Abe class       !
  !************************!

  type, extends(KE_TurbModel), public :: Abe_TurbModel
  contains
    procedure :: set_constants => set_constants_Abe
    procedure :: set_mut => set_mut_Abe
  end type Abe_TurbModel


contains
!****************************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************************

  !************************!
  !       Abe routines      !
  !************************!

type(Abe_TurbModel) function init_Abe_TurbModel(name)
  character(len=3), intent(IN) :: name
  init_Abe_TurbModel%name=name
end function init_Abe_TurbModel

subroutine set_constants_Abe(this)
  class(Abe_TurbModel) :: this
  this%cmu = 0.09
  this%sigmak = 1.4
  this%sigmae = 1.4
  this%ce1 = 1.5
  this%ce2 = 1.9
end subroutine set_constants_Abe

subroutine set_mut_Abe(this,u,w,rho,mu,mui,mut)
  use mod_param, only : i1,k1,imax,kmax
  use mod_mesh, only : walldist
  implicit none
  class(Abe_TurbModel) :: this
  real(8),dimension(0:i1,0:k1),intent(IN) :: u,w,rho,mu,mui
  real(8),dimension(0:i1,0:k1),intent(OUT):: mut
  integer  im,ip,km,kp,i,k
  real(8),dimension(0:k1) ::   tauw, utau
  real(8),dimension(0:i1,0:k1) :: Ret, yp, Reps
  real(8) :: nu_wall

  do k=1,kmax
    km=k-1
    kp=k+1
    tauw(k) = mui(imax,k)*0.5*(w(imax,km)+w(imax,k))/walldist(imax)
    utau(k) = (tauw(k)/(0.5*(rho(imax,k)+rho(i1,k))))**0.5
    do i=1,imax
      im=i-1
      ip=i+1
      nu_wall = mui(imax,k)/(0.5*(rho(imax,k)+rho(i1,k)))
      this%yp(i,k) = rho(i,k)*utau(k)*walldist(i)/mu(i,k)
      Ret(i,k)     = rho(i,k)*(this%k(i,k)**2.)/(mu(i,k)*this%eps(i,k))        ! not sure if r2 or r
      Reps(i,k)    = (walldist(i)*((mu(i,k)/rho(i,k))*this%eps(i,k))**.25)/(mu(i,k)/rho(i,k))
      this%fmu (i,k)= ((1.-exp(-Reps(i,k)/14.))**2) * (1.+(5./(Ret(i,k)**0.75))*exp(-(Ret(i,k)/200)**2))
      this%feps(i,k)= ((1.-exp(-Reps(i,k)/3.1))**2) * (1.-0.3                  *exp(-(Ret(i,k)/6.5)**2))
      mut(i,k) = min(1.,rho(i,k)*this%cmu*this%fmu(i,k)*this%k(i,k)**2./(this%eps(i,k)))
    enddo
  enddo
end subroutine set_mut_Abe

end module abe_tm