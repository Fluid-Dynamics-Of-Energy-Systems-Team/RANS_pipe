module vf_tm
  use ke_tm
  implicit none

!****************************************************************************************

  !************************!
  !         VF class       !
  !************************!

  type, extends(KE_TurbModel), public :: VF_TurbModel
  contains
    procedure :: set_constants => set_constants_VF
    procedure :: set_mut => set_mut_VF
    procedure :: advance_KE => advance_VF
    procedure :: set_bc => set_bc_VF
    procedure :: init_w_inflow => init_w_inflow_VF
    procedure :: production_KE => production_VF
    procedure :: rhs_v2_VF
    procedure :: solve_v2_VF
    procedure :: fillhem
  end type VF_TurbModel

contains

!****************************************************************************************

  !************************!
  !      VF routines       !
  !************************!

type(VF_TurbModel) function init_VF_TurbModel(name)
  character(len=2), intent(IN) :: name
  init_VF_TurbModel%name=name
end function init_VF_TurbModel

subroutine set_constants_VF(this)
  implicit none
  class(VF_TurbModel) :: this
  this%sigmak = 1.0
  this%sigmae = 1.3
  this%cmu    = 0.22
  this%ce1    = 1.4
  this%ce2    = 1.9
end subroutine

subroutine init_w_inflow_VF(this,nuSAin,pkin,kin,epsin,omin,mutin,v2in)
  use mod_param, only : i1,k1,k
  class(VF_TurbModel) :: this
  real(8), dimension(0:i1), intent(IN) :: nuSAin,pkin,kin,epsin,omin,mutin,v2in
    this%epsin = epsin
    this%kin  = kin
    this%v2in = v2in
    this%pkin = pkin
    do k=0,k1
      this%eps(:,k) = this%epsin(:)
      this%k(:,k) = this%kin(:)
      this%v2(:,k) = this%v2in(:)
      this%Pk(:,k) = this%pkin(:)
    enddo
end subroutine init_w_inflow_VF

subroutine set_mut_VF(this,u,w,rho,mu,mui,mut)
  use mod_param, only :k1,i1,imax,kmax,i,k
  use mod_mesh, only : mesh,dzw,dzp,dru, walldist,rp
  implicit none
  class(VF_TurbModel) :: this
  real(8), dimension(0:i1,0:k1), intent(IN) :: u, w, rho, mu, mui
  real(8), dimension(0:i1,0:k1), intent(OUT):: mut
  real(8), dimension(0:k1) :: tauw
  real(8) :: StR
  real(8) :: dz
  integer :: im,ip,km,kp

  do k=1,kmax
    km=k-1
    kp=k+1
    tauw(k) = mui(imax,k)*0.5*(w(imax,km)+w(imax,k))/walldist(imax)
    do i=1,imax
      im=i-1
      ip=i+1
      this%yp(i,k) = sqrt(rho(i,k))/mu(i,k)*(walldist(i))*tauw(k)**0.5           ! ystar
      StR= (2.*(((w(i,k)-w(i,km))/dz)**2.          + &
                ((u(i,k)-u(im,k))/dru(i))**2.      + &
                ((u(i,k)+u(im,k))/(2.*Rp(i)))**2.) +  &
        (((w(ip,km)+w(ip,k)+w(i,km)+w(i,k))/4.-(w(im,km)+w(im,k)+w(i,km)+w(i,k))/4.)/dru(i) &
        +((u(i,kp)+u(im,kp)+u(i,k)+u(im,k))/4.-(u(im,km)+u(i,km)+u(im,k)+u(i,k))/4.)/(dz)  )**2.)
      
      ! StR= (2.*(((w(i,k)-w(i,km))/dzw(k))**2.      + &
      !           ((u(i,k)-u(im,k))/dru(i))**2.      + &
      !           ((u(i,k)+u(im,k))/(2.*Rp(i)))**2.) + &
      !   (((w(ip,km)+w(ip,k)+w(i,km)+w(i,k))/4.-(w(im,km)+w(im,k)+w(i,km)+w(i,k))/4.)/dru(i) &
      !   +((u(i,kp)+u(im,kp)+u(i,k)+u(im,k))/4.-(u(im,km)+u(i,km)+u(im,k)+u(i,k))/4.)/dzw(k)))**2.)
      

      this%Tt(i,k) = max(this%k(i,k)/this%eps(i,k),6.0*(mu(i,k)/(rho(i,k)*this%eps(i,k)))**0.5)

      ! Srsq(i,k) = Str*rNew(i,k)*0.5
      ! if (modVF.eq.1) then
      !   ! Modifications: Lien&Kalitzin 2001 "Computations of transonic flow with the v2f turbulence model"
      !   Tt(i,k)   = max(Tt(i,k), 1.0e-8)
      !   Tt(i,k)   = min(Tt(i,k),0.6*kNew(i,k)/(3.**0.5*v2New(i,k)*cmu*(2.*Srsq(i,k))**0.5))
      ! endif
      this%fmu(i,k) = this%v2(i,k)*this%Tt(i,k)/(this%k(i,k)**2./this%eps(i,k))
      this%f1(i,k)  = 1.0 + 0.045*(this%k(i,k)/this%v2(i,k))**0.5
      this%f2(i,k)  = 1.0
      mut(i,k) = min(1.,rho(i,k)*this%cmu*this%v2(i,k)*this%Tt(i,k))
    enddo
  enddo
end subroutine set_mut_VF

subroutine set_bc_VF(this,mu,rho,periodic,rank,px)
  use mod_mesh, only : mesh,top_bcnovalue,bot_bcnovalue,top_bcvalue,bot_bcvalue,walldist
  use mod_param, only : k1,i1,kmax,imax,i,k
  implicit none
  class(VF_TurbModel) :: this
  real(8),dimension(0:i1,0:k1),intent(IN) :: rho,mu
  integer,                     intent(IN) :: periodic, rank, px
  real(8),dimension(0:k1) :: BCvalue
  real(8),dimension(0:i1) :: tmp
  real(8) :: topBCvalue, botBCvalue

  do k = 0,k1 
    this%k(0,k)  = bot_bcnovalue(k)*this%k(1,k)        !symmetry or 0 value
    this%k(i1,k) = top_bcnovalue(k)*this%k(imax,k)     !symmetry or 0 value
    this%v2(0,k) = bot_bcnovalue(k)*this%v2(1,k)       !symmetry or 0 value
    this%v2(i1,k)= top_bcnovalue(k)*this%v2(imax,k)    !symmetry or 0 value
    botBCvalue = 2.0*mu(1,k)/rho(1,k)*this%k(1,k)/walldist(1)**2                                            !bcvalue
    this%eps(0,k)       = (1.-bot_bcvalue(k))*(2.0*botBCvalue-this%eps(1,k)) + bot_bcvalue(k)*this%eps(1,k) !symmetry or bc value
    topBCvalue = 2.0*mu(imax,k)/rho(imax,k)*this%k(imax,k)/walldist(imax)**2                                !bcvalue
    this%eps(i1,k) = (1.-top_bcvalue(k))*(2.0*topBCvalue-this%eps(imax,k)) + top_bcvalue(k)*this%eps(imax,k)!symmetry or bc value
  enddo

  call shiftf(this%k,  tmp,rank); this%k  (:,0)      =tmp(:);
  call shiftf(this%eps,tmp,rank); this%eps(:,0)      =tmp(:);
  call shiftf(this%v2 ,tmp,rank); this%v2 (:,0)      =tmp(:);
  call shiftb(this%k,  tmp,rank); this%k  (:,k1)=tmp(:);
  call shiftb(this%eps,tmp,rank); this%eps(:,k1)=tmp(:);
  call shiftb(this%v2 ,tmp,rank); this%v2 (:,k1)=tmp(:);

  ! developing
  if (periodic.eq.1) return
  if (rank.eq.0) then
    this%k  (:,0) = this%kin(:)
    this%eps(:,0) = this%epsin(:)
    this%v2  (:,0)= this%v2in(:)
  endif
  if (rank.eq.px-1) then
    this%k  (:,k1)= 2.0*this%k  (:,kmax)-this%k  (:,kmax-1)
    this%eps(:,k1)= 2.0*this%eps(:,kmax)-this%eps(:,kmax-1)
    this%v2 (:,k1)= 2.0*this%v2 (:,kmax)-this%v2 (:,kmax-1)
  endif

end subroutine set_bc_VF

subroutine advance_VF(this,u,w,rho,mu,mui,muk,mut,beta,temp, &
                      alpha1,alpha2,alpha3,                  &
                      modification,rank,periodic,   &
                      residual1, residual2, residual3)
  use mod_param, only : k1,i1,kmax,imax,i,k
  implicit none
  class(VF_TurbModel) :: this
  real(8), dimension(0:i1,0:k1),intent(IN) :: u,w,rho,mu,mui,muk,mut,beta,temp
  real(8),                                intent(IN) :: alpha1,alpha2,alpha3
  integer,                                intent(IN) :: modification,rank,periodic
  real(8),                                intent(OUT):: residual1,residual2,residual3
  real(8), dimension(0:i1,0:k1) :: rho_mod

  !1, our modification, 2, Aupoix modification
  if ((modification == 1) .or. (modification == 2)) then
    rho_mod = rho
  else
    rho_mod = 1.0
  endif

  call this%fillhem(mu,rho)
  ! call SOLVEhelm(this%fv2,rank,this%Lh,centerBC)
  call this%production_KE(u,w,temp,rho,mu,mut,beta)
  call this%solve_eps_KE(residual2,u,w,rho,mu,mui,muk,mut,rho_mod, &
                       alpha2,modification,rank,periodic)
  call this%solve_k_KE(residual1,u,w,rho,mu,mui,muk,mut,rho_mod, &
                       alpha1,modification,rank,periodic)
  call this%solve_v2_VF(residual3,u,w,rho,mu,mui,muk,mut,rho_mod, &
                       alpha3,modification,rank,periodic)
end subroutine advance_VF

subroutine solve_v2_VF(this,resV2,u,w,rho,mu,mui,muk,mut,rho_mod, &
                       alphav2,modification,rank,periodic)
  use mod_param, only : k1,i1,kmax,imax,i,k
  use mod_math
  use mod_mesh, only : mesh,top_bcnovalue,bot_bcnovalue,dru,ru,rp,drp
  implicit none
  class(VF_TurbModel) :: this
  real(8),dimension(0:i1,0:k1), intent(IN) :: u, w, rho,mu,mui,muk,mut,rho_mod
  real(8),                                intent(IN) :: alphav2
  integer,                                intent(IN) :: modification,rank,periodic
  real(8),                                intent(OUT):: resV2
  real(8), dimension(0:i1,0:k1) :: dnew,dimpl
  real(8), dimension(imax)           :: a,b,c,rhs
  real(8) dz

  resV2=0.0; dnew=0.0; dimpl = 0.0;

  call advecc(dnew,dimpl,this%v2,u,w,rank,periodic,.true.)
  call this%rhs_v2_VF(dnew,dimpl)    
  call diffc(dnew,this%v2,mu,mui,muk,mut,this%sigmak,rho,modification)

  do k=1,kmax
    do i=1,imax
      
      if ((modification == 0) .or. (modification == 1)) then
        a(i) = (mui(i-1,k)+0.5*(mut(i,k)+mut(i-1,k))/this%sigmak)/((0.5*(rho_mod(i-1,k)+rho_mod(i,k)))**0.5)
        a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/rho(i,k)/(rho_mod(i,k)**0.5)

        c(i) = (mui(i  ,k)+0.5*(mut(i,k)+mut(i+1,k))/this%sigmak)/((0.5*(rho_mod(i+1,k)+rho_mod(i,k)))**0.5)
        c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/rho(i,k)/(rho_mod(i,k)**0.5)
      else if (modification == 2) then
        a(i) = (mui(i-1,k)+0.5*(mut(i,k)+mut(i-1,k))/this%sigmak)/(0.5*(rho_mod(i-1,k)+rho_mod(i,k)))
        a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/rho(i,k)

        c(i) = (mui(i  ,k)+0.5*(mut(i,k)+mut(i+1,k))/this%sigmak)/(0.5*(rho_mod(i+1,k)+rho_mod(i,k)))
        c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/rho(i,k)
      endif

      b(i) = (rho_mod(i,k)*(-a(i)-c(i)) + dimpl(i,k))
      a(i) = a(i)*rho_mod(i-1,k)
      c(i) = c(i)*rho_mod(i+1,k)
      rhs(i) =  dnew(i,k) + ((1-alphav2)/alphav2)*b(i)*this%v2(i,k)
    enddo

    i=1
    b(i)=b(i)+bot_bcnovalue(k)*a(i)
    rhs(i) = dnew(i,k) + ((1-alphav2)/alphav2)*b(i)*this%v2(i,k)

    i=imax
    b(i)=b(i)+top_bcnovalue(k)*c(i)
    rhs(i) = dnew(i,k) + ((1-alphav2)/alphav2)*b(i)*this%v2(i,k)

    call matrixIdir(imax,a,b/alphav2,c,rhs)

    do i=1,imax
      resV2 = resV2 + ((this%v2(i,k) - rhs(i))/(this%v2(i,k)+1.0e-20))**2.0
      this%v2(i,k) = min(2.0/3.0*this%k(i,k), max(rhs(i), 1.0e-8))
    enddo
  enddo
end subroutine solve_v2_VF

subroutine production_VF(this,u,w,temp,rho,mu,mut,beta)
  use mod_param, only :k1,i1,imax,kmax,i,k
  use mod_mesh, only : rp,ru,dru,drp,dzw,dzp
  implicit none
  class(VF_TurbModel) :: this
  real(8), dimension(0:i1,0:k1), intent(IN) :: u,w,temp,rho,mu,mut,beta
  real(8), dimension(0:i1,0:k1) :: div
  integer im,ip,km,kp
  real(8) :: Fr_1,ctheta,StR
  real(8) :: dz

  do k=1,kmax
    kp=k+1
    km=k-1
    do i=1,imax
      ip=i+1
      im=i-1

      !Production of turbulent kinetic energy
      this%Pk(i,k) = mut(i,k)*(  &
        2.*(((w(i,k)-w(i,km))/dz)**2.          +  &
            ((u(i,k)-u(im,k))/dRu(i))**2.      +  &
            ((u(i,k)+u(im,k))/(2.*Rp(i)))**2.) +  &
        (((w(ip,km)+w(ip,k)+w(i,km)+w(i,k))/4.-(w(im,km)+w(im,k)+w(i,km)+w(i,k))/4.)/dRu(i)  &
        +((u(i,kp)+u(im,kp)+u(i,k)+u(im,k))/4.-(u(im,km)+u(i,km)+u(im,k)+u(i,k))/4.)/(dz)  &
        )**2.)

      ! this%Pk(i,k) = mut(i,k)*(  &
      !   2.*(((w(i,k)-w(i,km))/dzw(k))**2.      +  &
      !       ((u(i,k)-u(im,k))/dRu(i))**2.      +  &
      !       ((u(i,k)+u(im,k))/(2.*Rp(i)))**2.) +  &
      !   (((w(ip,km)+w(ip,k)+w(i,km)+w(i,k))/4.-(w(im,km)+w(im,k)+w(i,km)+w(i,k))/4.)/dRu(i)  &
      !   +((u(i,kp)+u(im,kp)+u(i,k)+u(im,k))/4.-(u(im,km)+u(i,km)+u(im,k)+u(i,k))/4.)/dzw(k)  &
      !   )**2.)

      div(i,k) =(Ru(i)*u(i,k)-Ru(im)*u(im,k))/(Rp(i)*dru(i))  &
               +(      w(i,k) -      w(i,km))/dz

      ! div(i,k) =(Ru(i)*u(i,k)-Ru(im)*u(im,k))/(Rp(i)*dru(i))  &
      !          +(      w(i,k) -      w(i,km))/dzw(k)
      
      this%Pk(i,k) = this%Pk(i,k) - 2./3.*(rho(i,k)*this%k(i,k)+mut(i,k)*(div(i,k)))*(div(i,k))

      ! turbulent time scale
      this%Tt(i,k)   = max(this%k(i,k)/this%eps(i,k), 6.0*(mu(i,k)/(rho(i,k)*this%eps(i,k)))**0.5)
            
      ! if (modVF.eq.1) then
      !   StR = (2.*(((W(i,k)-W(i,km))/dz)**2. + &
      !     ((U(i,k)-U(im,k))/(Ru(i)-Ru(im)))**2. + &
      !     ((U(i,k)+U(im,k))/(2.*Rp(i)))**2.) + &
      !     (((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4. &
      !      -(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/dRu(i) &
      !     +((U(i,kp) +U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz)  )**2.)
    
      !   !               Srsq(i,k) = Pk(i,k)*rho(i,k)/(2.*ekmt(i,k))
      !   Srsq(i,k) = Str*rho(i,k)*0.5
      !   ! Modifications: Lien&Kalitzin 2001 "Computations of transonic flow with the v2f turbulence model"
      !   Tt(i,k)   = max(Tt(i,k), 1.0e-8)
      !   Tt(i,k)   = min(Tt(i,k),0.6*putink(i,k)/(3.**0.5*putinv2(i,k)*cmu*(2.*Srsq(i,k))**0.5))
      ! endif

      ! Bouyancy production
      this%Gk(i,k)=-ctheta*beta(i,k)*Fr_1*this%Tt(i,k)  &
        *  (mut(i,k)*(((w(ip,km)+w(ip,k)+w(i,km)+w(i,k))/4.-(w(im,km)+w(im,k)+w(i,km)+w(i,k))/4.)/dRu(i)  &
                     +((u(i,kp)+u(im,kp)+u(i,k)+u(im,k))/4.-(u(im,km)+u(i,km)+u(im,k)+u(i,k))/4.)/(dz))*  &
        (temp(ip,k)-temp(im,k))/(dRp(i)+dRp(im))  )  &
        +(2.*mut(i,k)*((w(i,k)-w(i,km))/dz-2./3.*(rho(i,k)*this%k(i,k)))*(temp(i,kp)-temp(i,km))/(2.*dz)  &
        )

      ! this%Gk(i,k)=-ctheta*beta(i,k)*Fr_1*this%Tt(i,k)  &
      !   *  (mut(i,k)*(((w(ip,km)+w(ip,k)+w(i,km)+w(i,k))/4.-(w(im,km)+w(im,k)+w(i,km)+w(i,k))/4.)/dRu(i)  &
      !                +((u(i,kp)+u(im,kp)+u(i,k)+u(im,k))/4.-(u(im,km)+u(i,km)+u(im,k)+u(i,k))/4.)/dzw(k))*  &
      !   (temp(ip,k)-temp(im,k))/(dRp(i)+dRp(im))  )  &
      !   +(2.*mut(i,k)*((w(i,k)-w(i,km))/dzw(k)-2./3.*(rho(i,k)*this%k(i,k)))*(temp(i,kp)-temp(i,km))/(2.*dzp(k))  &
      !   )


      this%Gk(i,k) = this%Gk(i,k) + ctheta*beta(i,k)*Fr_1*this%Tt(i,k)*2./3.*mut(i,k)*div(i,k)*(temp(i,kp)-temp(i,km))/(2.*dz)
      ! this%Gk(i,k) = this%Gk(i,k) + ctheta*beta(i,k)*Fr_1*this%Tt(i,k)*2./3.*mut(i,k)*div(i,k)*(temp(i,kp)-temp(i,km))/(2.*dzp(k))

    enddo
  enddo
end subroutine production_VF

subroutine rhs_v2_VF(this,putout,dimpl)
  use mod_param, only : k1,i1,imax,kmax,i,k
  implicit none
  class(VF_TurbModel) :: this
  real(8), dimension(0:i1,0:k1), intent(OUT):: putout,dimpl
      
  do k=1,kmax
    do i=1,imax
      !v'2 equation
      putout(i,k) = putout(i,k) + this%k(i,k)*this%fv2(i,k)       ! note, source is rho*k*f/rho
      dimpl(i,k)  = dimpl(i,k)  + 6.*this%eps(i,k)/this%k(i,k)    ! note, 6*rho*v'2*epsilon/k/(rho*v'2), set implicit and divided by density
    enddo
  enddo

end subroutine rhs_v2_VF

subroutine fillhem(this,mu,rho)
  use mod_param, only : k1,i1,kmax,imax,i,k
  implicit none
  class(VF_TurbModel) :: this
  real(8), dimension(0:i1,0:k1), intent(IN) :: mu, rho
  real(8), dimension(0:i1,0:k1) :: Srsq

  do  k=1,kmax
    do i=1,imax         
      ! time scale
      this%Tt(i,k)   = max(this%k(i,k)/this%eps(i,k), 6.0*(mu(i,k)/(rho(i,k)*this%eps(i,k)))**0.5)           
      ! lenght scale
      this%Lh(i,k)=0.23*max(this%k(i,k)**1.5/this%eps(i,k),70.*((mu(i,k)/rho(i,k))**3./this%eps(i,k))**0.25)
            
      ! if (modVF.eq.1) then
      !   ! Modifications: Lien&Kalitzin 2001 "Computations of transonic flow with the v2f turbulence model"
      !   Tt(i,k)   = max(Tt(i,k), 1.0e-8)
      !   Tt(i,k)   = min(Tt(i,k),0.6*knew(i,k)/(3.**0.5*v2new(i,k)*cmu*(2.*Srsq(i,k))**0.5))

      !   Lh(i,k)=min(knew(i,k)**1.5/enew(i,k),knew(i,k)**1.5/(3.**0.5*v2new(i,k)*cmu*(2.*Srsq(i,k))**0.5))
      !   Lh(i,k)=0.23*max(Lh(i,k),70.*((ekm(i,k)/rnew(i,k))**3./enew(i,k))**0.25)
      ! endif

      ! f-equation also has Gk: Kenjeres et al 2005 "Contribution to elliptic relaxation modelling of turbulent natural and mixed convection"
      this%fv2(i,k)= - (1.4-1.)*(2./3.-this%v2(i,k)/this%k(i,k))/this%Tt(i,k) &
                - 0.3*(this%Pk(i,k)+this%Gk(i,k))/(rho(i,k)*this%k(i,k))-5.*this%v2(i,k)/(this%k(i,k)*this%Tt(i,k))
      this%fv2(i,k) = this%fv2(i,k)/this%Lh(i,k)**2.0
    enddo
  enddo
end subroutine fillhem


end module
