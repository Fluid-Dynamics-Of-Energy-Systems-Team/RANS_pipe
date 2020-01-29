module mk_tm
  use ke_tm
  implicit none

!****************************************************************************************

  !************************!
  !         MK class       !
  !************************!

  type, extends(KE_TurbModel), public :: MK_TurbModel
  contains
    procedure :: set_constants => set_constants_MK
    procedure :: set_mut => set_mut_MK
  end type MK_TurbModel

contains
!****************************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************************

  !************************!
  !       MK routines      !
  !************************!

type(MK_TurbModel) function init_MK_TurbModel(name)
  character(len=2), intent(IN) :: name
  init_MK_TurbModel%name=name
end function init_MK_TurbModel

subroutine set_constants_MK(this)
  class(MK_TurbModel) :: this
  this%sigmak = 1.4
  this%sigmae = 1.3
  this%cmu    = 0.09
  this%ce1    = 1.4
  this%ce2    = 1.8
end subroutine set_constants_MK

subroutine set_mut_MK(this,u,w,rho,mu,mui,mut)
  use mod_param, only : kmax,imax,k1,i1,k,i,Re, modifDiffTerm
  use mod_mesh,  only : walldist
  implicit none
  class(MK_TurbModel) :: this
  real(8),dimension(0:i1,0:k1),intent(IN) :: u,w,rho,mu,mui
  real(8),dimension(0:i1,0:k1),intent(OUT):: mut
  integer  im,ip,km,kp
  real(8),dimension(0:k1) ::   tauw, utau
  real(8),dimension(0:i1,0:k1) :: Ret, yp
  real(8) :: rho_wall, mu_wall

  do k=1,kmax
    km=k-1
    kp=k+1
    tauw(k) = mui(imax,k)*0.5*(w(imax,km)+w(imax,k))/walldist(imax)
    
    do i=1,imax
      im=i-1
      ip=i+1
      if (modifDiffTerm .eq. 1) then
        this%yp(i,k) = sqrt(rho(i,k))/mu(i,k)*walldist(i)*tauw(k)**0.5          ! ystar
     else
        this%yp(i,k) = sqrt(rho(imax,k))/mu(imax,k)*walldist(i)*tauw(k)**0.5    ! yplus
      endif
      Ret(i,k)     = rho(i,k)*(this%k(i,k)**2.)/(mu(i,k)*this%eps(i,k))        ! not sure if r2 or r
      this%fmu(i,k)= (1.-exp(-this%yp(i,k)/70.))*(1+3.45/Ret(i,k)**0.5)
      this%feps(i,k) = (1.-2./9.*exp(-(Ret(i,k)/6.)**2.))*(1.-exp(-this%yp(i,k)/5.))**2.0
      mut(i,k) = min(1.,rho(i,k)*this%cmu*this%fmu(i,k)*this%k(i,k)**2./(this%eps(i,k)))
    enddo
  enddo
end subroutine set_mut_MK

end module mk_tm