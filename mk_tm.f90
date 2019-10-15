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
    procedure :: advance_turb => advance_MK
    procedure :: set_bc => set_bc_MK
    procedure :: production_MK
  end type MK_TurbModel


contains
!****************************************************************************************

  !************************!
  !       MK routines      !
  !************************!

type(MK_TurbModel) function init_MK_TurbModel(i1,k1,imax,kmax)
  integer, intent(in) :: i1,k1,imax,kmax
  init_MK_TurbModel%i1 = i1
  init_MK_TurbModel%k1 = k1
  init_MK_TurbModel%imax = imax
  init_MK_TurbModel%kmax = kmax
end function init_MK_TurbModel

subroutine set_constants_MK(this)
  class(MK_TurbModel) :: this
  this%sigmak = 1.4
  this%sigmae = 1.3
  this%cmu = 0.09
  this%ce1 = 1.4
  this%ce2 = 1.8
end subroutine set_constants_MK

subroutine set_mut_MK(this,u,w,rho,mu,mui,walldist,Rp,dRp,dru,dz,mut)
  implicit none
  class(MK_TurbModel) :: this
  real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,mu,mui
  real(8),dimension(1:this%imax),        intent(IN) :: walldist
  real(8),dimension(0:this%i1),          intent(IN) :: Rp,dRp,dru
  real(8),                               intent(IN) :: dz
  real(8),dimension(0:this%i1,0:this%k1),intent(OUT):: mut
  integer  im,ip,km,kp,i,k
  real(8),dimension(0:this%k1) ::   tauw
  real(8),dimension(0:this%i1,0:this%k1) :: Ret, yp

  do k=1,this%kmax
    km=k-1
    kp=k+1
    tauw(k) = mui(this%imax,k)*0.5*(w(this%imax,km)+w(this%imax,k))/walldist(this%imax)
    do i=1,this%imax
      im=i-1
      ip=i+1
      yp(i,k)     = sqrt(rho(i,k))/mu(i,k)*(walldist(i))*tauw(k)**0.5           ! ystar
      Ret(i,k)    = rho(i,k)*(this%k(i,k)**2.)/(mu(i,k)*this%eps(i,k))        ! not sure if r2 or r
      this%fmu(i,k)     = (1.-exp(-yp(i,k)/70.))*(1+3.45/Ret(i,k)**0.5)
      this%f1(i,k)      = 1.
      this%f2(i,k)      = (1.-2./9.*exp(-(Ret(i,k)/6.)**2.))*(1.-exp(-yp(i,k)/5.))**2.0
      mut(i,k) = min(1.,rho(i,k)*this%cmu*this%fmu(i,k)*this%k(i,k)**2./(this%eps(i,k)))
    enddo
  enddo
end subroutine set_mut_MK

subroutine advance_MK(this,u,w,rho,mu,mui,muk,mut,beta,temp, &
                      Ru,Rp,dru,drp,dz,walldist,             &
                      alpha1,alpha2,alpha3,                  &
                      modification,rank,centerBC,periodic,   &
                      residual1, residual2, residual3)
  class(MK_TurbModel) :: this
  real(8), dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,mu,mui,muk,mut,beta,temp
  real(8), dimension(0:this%i1),          intent(IN) :: Ru,Rp,dru,drp
  real(8), dimension(1:this%i1),          intent(IN) :: walldist
  real(8),                                intent(IN) :: dz,alpha1,alpha2, alpha3
  integer,                                intent(IN) :: modification,rank,centerBC,periodic
  real(8),                                intent(OUT):: residual1,residual2, residual3
  real(8), dimension(0:this%i1,0:this%k1) :: rho_mod

  !1, our modification, 2, Aupoix modification
  if ((modification == 1) .or. (modification == 2)) then
    rho_mod = rho
  else
    rho_mod = 1.0
  endif

  call this%production_MK(u,w,temp,rho,mut,beta,Rp,Ru,dRu,dRp,dz)
  call this%solve_eps_KE(residual2,u,w,rho,mu,mui,muk,mut,rho_mod, &
                       Ru,Rp,dru,drp,dz, &
                       alpha2,modification,rank,centerBC,periodic)
  call this%solve_k_KE(residual1,u,w,rho,mu,mui,muk,mut,rho_mod, &
                       Ru,Rp,dru,drp,dz, &
                       alpha1,modification,rank,centerBC,periodic)
end

subroutine set_bc_MK(this,mu,rho,walldist,centerBC,periodic,rank,px)
  class(MK_TurbModel) :: this
  real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: rho,mu
  real(8),dimension(1:this%imax),        intent(IN) :: walldist
  integer,                               intent(IN) :: centerBC,periodic, rank, px
  real(8),dimension(0:this%k1) :: BCvalue
  real(8),dimension(0:this%i1) :: tmp

  this%k(this%i1,:) = -this%k(this%imax,:)
  BCvalue(:) = 2.0*mu(this%imax,:)/rho(this%imax,:)*this%k(this%imax,:)/walldist(this%imax)**2
  this%eps(this%i1,:) = 2.0*BCvalue(:) - this%eps(this%imax,:)

  ! channel
  if (centerBC.eq.-1) then
    this%k(0,:)  = -this%k(1,:)
    BCvalue(:)   = 2.0*mu(1,:)/rho(1,:)*this%k(1,:)/walldist(1)**2
    this%eps(0,:)= 2.0*BCvalue(:) - this%eps(1,:)
  endif

  ! pipe/BL
  if (centerBC.eq.1) then
    this%k  (0,:) = this%k  (1,:)
    this%eps(0,:) = this%eps(1,:)
  endif

  call shiftf(this%k,  tmp,rank); this%k  (:,0)      =tmp(:);
  call shiftf(this%eps,tmp,rank); this%eps(:,0)      =tmp(:);
  call shiftb(this%k,  tmp,rank); this%k  (:,this%k1)=tmp(:);
  call shiftb(this%eps,tmp,rank); this%eps(:,this%k1)=tmp(:);
  
  ! developing
  if (periodic.eq.1) return
  if (rank.eq.0) then
    this%k  (:,0) = kin(:)
    this%eps(:,0) = ein(:)
  endif
  if (rank.eq.px-1) then
    this%k  (:,this%k1)= 2.0*this%k  (:,this%kmax)-this%k  (:,this%kmax-1)
    this%eps(:,this%k1)= 2.0*this%eps(:,this%kmax)-this%eps(:,this%kmax-1)
  endif

end subroutine set_bc_MK

subroutine production_MK(this,u,w,temp,rho,mut,beta,Rp,Ru,dRu,dRp,dz)
  implicit none
  class(MK_TurbModel) :: this
  real(8), dimension(0:this%i1,0:this%k1), intent(IN) :: u,w,temp,rho,mut,beta
  real(8), dimension(0:this%i1),           intent(IN) :: Rp,Ru,dRu,dRp
  real(8),                                 intent(IN) :: dz
  real(8), dimension(0:this%i1,0:this%k1) :: div
  integer im,ip,jm,jp,km,kp,ib,ie,kb,ke,i,k
  real(8) :: Fr_1, ctheta

  Fr_1      = 0.0 !!!NOTE: this was originally in the param!!!!
  ctheta    = 0.3 !!!NOTE: this was originally in the param!!!!

  ib = 1
  ie = this%i1-1

  kb = 1
  ke = this%k1-1

  do k=kb,ke
    kp=k+1
    km=k-1
    do i=ib,ie
      ip=i+1
      im=i-1
            
      ! Production of turbulent kinetic energy
      this%Pk(i,k) = mut(i,k)*(  &
        2.*(((w(i,k)-w(i,km))/dz)**2. +  &
        ((u(i,k)-u(im,k))/dRu(i))**2. +  &
        ((u(i,k)+u(im,k))/(2.*Rp(i)))**2.) +  &
        (((w(ip,km)+w(ip,k)+w(i,km)+w(i,k))/4.-(w(im,km)+w(im,k)+w(i,km)+w(i,k))/4.)/dRu(i)  &
        +((u(i,kp)+u(im,kp)+u(i,k)+u(im,k))/4.-(u(im,km)+u(i,km)+u(im,k)+u(i,k))/4.)/(dz)  &
        )**2.)

      div(i,k) =(Ru(i)*u(i,k)-Ru(im)*u(im,k))/(Rp(i)*dru(i))  &
        +(      w(i,k) -      w(i,km))/dz

      this%Pk(i,k) = this%Pk(i,k) - 2./3.*(rho(i,k)*this%k(i,k)+mut(i,k)*(div(i,k)))*(div(i,k))

      ! turbulent time scale
      this%Tt(i,k)=this%k(i,k)/this%eps(i,k)

      ! Bouyancy production
      this%Gk(i,k)=-ctheta*beta(i,k)*Fr_1*this%Tt(i,k)  &
        *  (mut(i,k)*(((w(ip,km)+w(ip,k)+w(i,km)+w(i,k))/4.-(w(im,km)+w(im,k)+w(i,km)+w(i,k))/4.)/dRu(i)  &
                     +((u(i,kp)+u(im,kp)+u(i,k)+u(im,k))/4.-(u(im,km)+u(i,km)+u(im,k)+u(i,k))/4.)/(dz))*  &
        (temp(ip,k)-temp(im,k))/(dRp(i)+dRp(im))  )  &
        +(2.*mut(i,k)*((w(i,k)-w(i,km))/dz-2./3.*(rho(i,k)*this%k(i,k)))*(temp(i,kp)-temp(i,km))/(2.*dz)  &
        )

      this%Gk(i,k) = this%Gk(i,k) + ctheta*beta(i,k)*Fr_1*this%Tt(i,k)*2./3.*mut(i,k)*div(i,k)*(temp(i,kp)-temp(i,km))/(2.*dz)

    enddo
  enddo
end subroutine production_MK

end module mk_tm