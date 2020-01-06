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
    procedure :: advance_turb => advance_Abe
    procedure :: set_bc => set_bc_Abe
    procedure :: init_w_inflow => init_w_inflow_Abe
    procedure :: rhs_eps_KE => rhs_eps_Abe
    procedure :: rhs_k_KE => rhs_k_Abe
    procedure :: production_Abe
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

subroutine init_w_inflow_Abe(this,Re,systemsolve)
  use mod_param, only : i1,k1,imax,kmax
  class(Abe_TurbModel) :: this
  real(8), intent(IN) :: Re
  integer, intent(IN) :: systemsolve
  real(8), dimension(0:i1) :: dummy
  character(len=5)  :: Re_str
  integer           :: Re_int,i,k
  Re_int = int(Re)
  write(Re_str,'(I5.5)') Re_int
  if (systemsolve .eq. 1) open(29,file = 'pipe/Inflow_'//TRIM(this%name)//'_'//Re_str//'.dat',form='unformatted')
  if (systemsolve .eq. 2) open(29,file = 'channel/Inflow_'//TRIM(this%name)//'_'//Re_str//'.dat',form='unformatted')
  if (systemsolve .eq. 3) open(29,file = 'symchan/Inflow_'//TRIM(this%name)//'_'//Re_str//'.dat',form='unformatted')
  read(29) dummy(:),this%kin(:),this%epsin(:),dummy(:),dummy(:),dummy(:),this%mutin(:),dummy(:)
  close(29)
  do k=0,k1
    this%eps(:,k) = this%epsin(:)
    this%k(:,k) = this%kin(:)
  enddo
end subroutine init_w_inflow_Abe

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
      !this%yp(i,k) = sqrt(rho(i,k))/mu(i,k)*(walldist(i))*tauw(k)**0.5         ! ystar
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

subroutine set_bc_Abe(this,mu,rho,periodic,rank,px)
  use mod_param, only : i1,k1,imax,kmax,k
  use mod_mesh, only : top_bcvalue,bot_bcvalue,top_bcnovalue,bot_bcnovalue,walldist
  implicit none
  class(Abe_TurbModel) :: this
  real(8),dimension(0:i1,0:k1),intent(IN) :: rho,mu
  integer,                     intent(IN) :: periodic, rank, px
  real(8),dimension(0:i1) :: tmp
  real(8)                 :: topBCvalue, botBCvalue
  
  do k = 0,k1 
    this%k(0,k)  = bot_bcnovalue(k)*this%k(1,k)         !symmetry or 0 value
    this%k(i1,k) = top_bcnovalue(k)*this%k(imax,k)      !symmetry or 0 value
    botBCvalue   = 2.0*mu(1,k)/rho(1,k)*((this%k(1,k)**0.5)/walldist(1))**2                
    this%eps(0,k)= (1.-bot_bcvalue(k))*(2.0*botBCvalue-this%eps(1,k))      +bot_bcvalue(k)*this%eps(1,k)        !symmetry or bc value
    topBCvalue   = 2.0*mu(imax,k)/rho(imax,k)*((this%k(imax,k)**0.5)/walldist(imax))**2                          !bcvalue
    this%eps(i1,k) = (1.-top_bcvalue(k))*(2.0*topBCvalue-this%eps(imax,k)) +top_bcvalue(k)*this%eps(imax,k)!symmetry or bc value
  enddo

  call shiftf(this%k,  tmp,rank); this%k  (:,0) =tmp(:);
  call shiftf(this%eps,tmp,rank); this%eps(:,0) =tmp(:);
  call shiftb(this%k,  tmp,rank); this%k  (:,k1)=tmp(:);
  call shiftb(this%eps,tmp,rank); this%eps(:,k1)=tmp(:);
  
  ! developing
  if (periodic.eq.1) return
  if (rank.eq.0) then
    this%k  (:,0) = this%kin(:)
    this%eps(:,0) = this%epsin(:)
  endif
  if (rank.eq.px-1) then
    this%k  (:,k1)= 2.0*this%k  (:,kmax)-this%k  (:,kmax-1)
    this%eps(:,k1)= 2.0*this%eps(:,kmax)-this%eps(:,kmax-1)
  endif

end subroutine set_bc_Abe

subroutine rhs_eps_Abe(this,putout,dimpl,rho)
  use mod_param, only : i1,k1,imax,kmax,k,i
  implicit none
  class(Abe_TurbModel) :: this
  real(8), dimension(0:i1,0:k1), intent(IN) :: rho
  real(8), dimension(0:i1,0:k1), intent(OUT):: putout,dimpl
  do k=1,kmax
    do i=1,imax
      putout(i,k) = putout(i,k) +(this%ce1*(1/this%Tt(i,k))*this%Pk(i,k))/rho(i,k) !&
      dimpl(i,k)  = dimpl(i,k)  + this%ce2*this%feps(i,k)*(1./this%Tt(i,k))   ! note, ce2*f2*rho*epsilon/T/(rho*epsilon), set implicit and divided by density
    enddo
  enddo
end subroutine rhs_eps_Abe

subroutine rhs_k_Abe(this,putout,dimpl,rho)
  use mod_param, only : i1,k1,imax,kmax,i,k
  implicit none
  class(Abe_TurbModel) :: this
  real(8), dimension(0:i1,0:k1), intent(IN) :: rho
  real(8), dimension(0:i1,0:k1), intent(OUT):: putout,dimpl
  do k=1,kmax
    do i=1,imax
      putout(i,k) = putout(i,k)+ this%Pk(i,k)/rho(i,k)!+this%Gk(i,k))/rho(i,k)
      dimpl(i,k)  = dimpl(i,k) + this%eps(i,k)/this%k(i,k) ! note, rho*epsilon/(rho*k), set implicit and divided by density
    enddo
  enddo
end subroutine rhs_k_Abe

subroutine advance_Abe(this,u,w,rho,mu,mui,muk,mut,beta,temp, &
                      alpha1,alpha2,alpha3,                  &
                      modification,rank,periodic,   &
                      residual1, residual2, residual3)
  use mod_param, only : i1,k1,imax,kmax
  class(Abe_TurbModel) :: this
  real(8), dimension(0:i1,0:k1),intent(IN) :: u,w,rho,mu,mui,muk,mut,beta,temp
  real(8),                                intent(IN) :: alpha1,alpha2, alpha3
  integer,                                intent(IN) :: modification,rank,periodic
  real(8),                                intent(OUT):: residual1,residual2, residual3
  real(8), dimension(0:i1,0:k1) :: rho_mod

  !1, our modification, 2, Aupoix modification
  if ((modification == 1) .or. (modification == 2)) then
    rho_mod = rho
  else
    rho_mod = 1.0
  endif

  call this%production_Abe(u,w,temp,rho,mut,beta)
  call this%solve_eps_KE(residual2,u,w,rho,mu,mui,muk,mut,rho_mod, &
                         alpha2,modification,rank,periodic)
  call this%solve_k_KE(residual1,u,w,rho,mu,mui,muk,mut,rho_mod, &
                       alpha1,modification,rank,periodic)
end

subroutine production_Abe(this,u,w,temp,rho,mut,beta)
  use mod_param, only : i1,k1,imax,kmax,i,k
  use mod_mesh, only : Rp,Ru,dRu,dRp,dzw,dzp
  implicit none
  class(Abe_TurbModel) :: this
  real(8), dimension(0:i1,0:k1), intent(IN) :: u,w,temp,rho,mut,beta
  real(8), dimension(0:i1,0:k1) :: div
  integer im,ip,jm,jp,km,kp
  real(8) :: Fr_1, ctheta

  Fr_1      = 0.0 !!!NOTE: this was originally in the param!!!!
  ctheta    = 0.3 !!!NOTE: this was originally in the param!!!!

  do k=1,kmax
    kp=k+1
    km=k-1
    do i=1,imax
      ip=i+1
      im=i-1
      ! Production of turbulent kinetic energy
      this%Pk(i,k) = mut(i,k)*(  &
        2.*(((w(i,k)-w(i,km))/dzw(k))**2.      +  &
            ((u(i,k)-u(im,k))/dRu(i))**2.     +  &
            ((u(i,k)+u(im,k))/(2.*Rp(i)))**2.) +  &
        (((w(ip,km)+w(ip,k)+w(i,km)+w(i,k))/4.-(w(im,km)+w(im,k)+w(i,km)+w(i,k))/4.)/dRu(i)  &
        +((u(i,kp)+u(im,kp)+u(i,k)+u(im,k))/4.-(u(im,km)+u(i,km)+u(im,k)+u(i,k))/4.)/dzw(k)  &
        )**2.)

      div(i,k) =(Ru(i)*u(i,k)-Ru(im)*u(im,k))/(Rp(i)*dru(i))  &
               +(      w(i,k) -      w(i,km))/dzw(k)

      this%Pk(i,k) = this%Pk(i,k) - 2./3.*(rho(i,k)*this%k(i,k)+mut(i,k)*(div(i,k)))*(div(i,k))

      ! turbulent time scale
      this%Tt(i,k)=this%k(i,k)/this%eps(i,k)

      this%Gk(i,k)=-ctheta*beta(i,k)*Fr_1*this%Tt(i,k)  &
        *  (mut(i,k)*(((w(ip,km)+w(ip,k)+w(i,km)+w(i,k))/4.-(w(im,km)+w(im,k)+w(i,km)+w(i,k))/4.)/dRu(i)  &
                     +((u(i,kp)+u(im,kp)+u(i,k)+u(im,k))/4.-(u(im,km)+u(i,km)+u(im,k)+u(i,k))/4.)/dzw(k))*  &
        (temp(ip,k)-temp(im,k))/(dRp(i)+dRp(im))  )  &
        +(2.*mut(i,k)*((w(i,k)-w(i,km))/dzw(k)-2./3.*(rho(i,k)*this%k(i,k)))*(temp(i,kp)-temp(i,km))/(dzp(k)+dzp(km)))
      
      this%Gk(i,k) = this%Gk(i,k) + ctheta*beta(i,k)*Fr_1*this%Tt(i,k) &
                     *2./3.*mut(i,k)*div(i,k)*(temp(i,kp)-temp(i,km))/(dzp(k)+dzp(km))

    enddo
  enddo
end subroutine production_Abe

end module abe_tm