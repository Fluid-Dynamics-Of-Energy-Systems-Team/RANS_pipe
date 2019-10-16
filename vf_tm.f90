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
    procedure :: advance_turb => advance_VF
    procedure :: set_bc => set_bc_VF
    procedure :: production_VF
    procedure :: rhs_v2_VF
    procedure :: solve_v2_VF
    procedure :: fillhem
  end type VF_TurbModel

contains

!****************************************************************************************

  !************************!
  !      VF routines       !
  !************************!

subroutine set_constants_VF(this)
  implicit none
  class(VF_TurbModel) :: this
  this%sigmak = 1.0
  this%sigmae = 1.3
  this%cmu    = 0.22
  this%ce1    = 1.4
  this%ce2    = 1.9
end subroutine

subroutine set_mut_VF(this,u,w,rho,mu,mui,walldist,Rp,dRp,dru,dz,mut)
  implicit none
  class(VF_TurbModel) :: this
  real(8), dimension(0:this%i1,0:this%k1), intent(IN) :: u, w, rho, mu, mui
  real(8), dimension(1:this%imax),         intent(IN) :: walldist
  real(8), dimension(0:this%i1),           intent(IN) :: Rp,dRp, dru
  real(8),                                 intent(IN) :: dz
  real(8), dimension(0:this%i1,0:this%k1), intent(OUT):: mut
  real(8), dimension(0:this%k1) :: tauw
  real(8) :: StR
  integer :: im,ip,km,kp,i,k

  do k=1,this%kmax
    km=k-1
    kp=k+1
    tauw(k) = mui(this%imax,k)*0.5*(w(this%imax,km)+w(this%imax,k))/walldist(this%imax)
    do i=1,this%imax
      im=i-1
      ip=i+1
      StR= (2.*(((w(i,k)-w(i,km))/dz)**2. + &
        ((u(i,k)-u(im,k))/dru(i))**2. + &
        ((u(i,k)+u(im,k))/(2.*Rp(i)))**2.) +  &
        (((w(ip,km)+w(ip,k)+w(i,km)+w(i,k))/4.-(w(im,km)+w(im,k)+w(i,km)+w(i,k))/4.)/dru(i) &
        +((u(i,kp)+u(im,kp)+u(i,k)+u(im,k))/4.-(u(im,km)+u(i,km)+u(im,k)+u(i,k))/4.)/(dz)  )**2.)
      
      this%Tt(i,k)   = max(this%k(i,k)/this%eps(i,k),6.0*(mu(i,k)/(rho(i,k)*this%eps(i,k)))**0.5)

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

subroutine advance_VF(this,u,w,rho,mu,mui,muk,mut,beta,temp, &
                      Ru,Rp,dru,drp,dz,walldist,             &
                      alpha1,alpha2,alpha3,                  &
                      modification,rank,centerBC,periodic,   &
                      residual1, residual2, residual3)
  implicit none
  class(VF_TurbModel) :: this
  real(8), dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,mu,mui,muk,mut,beta,temp
  real(8), dimension(0:this%i1),          intent(IN) :: Ru,Rp,dru,drp
  real(8), dimension(1:this%i1),          intent(IN) :: walldist
  real(8),                                intent(IN) :: dz,alpha1,alpha2,alpha3
  integer,                                intent(IN) :: modification,rank,centerBC,periodic
  real(8),                                intent(OUT):: residual1,residual2,residual3
  real(8), dimension(0:this%i1,0:this%k1) :: rho_mod

  !1, our modification, 2, Aupoix modification
  if ((modification == 1) .or. (modification == 2)) then
    rho_mod = rho
  else
    rho_mod = 1.0
  endif

  call this%fillhem(mu,rho)
  call SOLVEhelm(this%fv2,Ru,Rp,dRu,dRp,dz,rank,this%Lh,centerBC)
  call this%production_VF(u,w,temp,rho,mu,mut,beta,Rp,Ru,dRu,dRp,dz)
  call this%solve_eps_KE(residual2,u,w,rho,mu,mui,muk,mut,rho_mod, &
                       Ru,Rp,dru,drp,dz, &
                       alpha2,modification,rank,centerBC,periodic)
  call this%solve_k_KE(residual1,u,w,rho,mu,mui,muk,mut,rho_mod, &
                       Ru,Rp,dru,drp,dz, &
                       alpha1,modification,rank,centerBC,periodic)
  call this%solve_v2_VF(residual3,u,w,rho,mu,mui,muk,mut,rho_mod, &
                       Ru,Rp,dru,drp,dz, &
                       alpha3,modification,rank,centerBC,periodic)
end

type(VF_TurbModel) function init_VF_TurbModel(i1,k1,imax,kmax,name)
  integer,          intent(IN) :: i1,k1,imax,kmax
  character(len=2), intent(IN) :: name
  init_VF_TurbModel%name=name
  init_VF_TurbModel%i1 = i1
  init_VF_TurbModel%k1 = k1
  init_VF_TurbModel%imax = imax
  init_VF_TurbModel%kmax = kmax
end function init_VF_TurbModel

subroutine set_bc_VF(this,mu,rho,walldist,centerBC,periodic,rank,px)
  class(VF_TurbModel) :: this
  real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: rho,mu
  real(8),dimension(1:this%imax),        intent(IN) :: walldist
  integer,                               intent(IN) :: centerBC,periodic, rank, px
  real(8),dimension(0:this%k1) :: BCvalue
  real(8),dimension(0:this%i1) :: tmp
  

  this%k (this%i1,:) = -this%k (this%imax,:)
  this%v2(this%i1,:) = -this%v2(this%imax,:)
  BCvalue(:)  = 2.0*mu(this%imax,:)/rho(this%imax,:)*this%k(this%imax,:)/walldist(this%imax)**2
  this%eps(this%i1,:)  = 2.0*BCvalue(:) - this%eps(this%imax,:)

  ! channel
  if (centerBC.eq.-1) then
    this%k(0,:)  = -this%k(1,:)
    this%v2(0,:) = -this%v2(1,:)
    BCvalue(:)   = 2.0*mu(1,:)/rho(1,:)*this%k(1,:)/walldist(1)**2
    this%eps(0,:)= 2.0*BCvalue(:) - this%eps(1,:)
  endif

  ! pipe/BL
  if (centerBC.eq.1) then
    this%k  (0,:)= this%k  (1,:)
    this%eps(0,:)= this%eps(1,:)
    this%v2 (0,:)= this%v2 (1,:)
  endif  

  call shiftf(this%k,  tmp,rank); this%k  (:,0)      =tmp(:);
  call shiftf(this%eps,tmp,rank); this%eps(:,0)      =tmp(:);
  call shiftf(this%v2 ,tmp,rank); this%v2 (:,0)      =tmp(:);
  call shiftb(this%k,  tmp,rank); this%k  (:,this%k1)=tmp(:);
  call shiftb(this%eps,tmp,rank); this%eps(:,this%k1)=tmp(:);
  call shiftb(this%v2 ,tmp,rank); this%v2 (:,this%k1)=tmp(:);

  ! developing
  if (periodic.eq.1) return
  ! if (rank.eq.0) then
  !   this%k  (:,0) = kin(:)
  !   this%eps(:,0) = ein(:)
  !   this%v2  (:0) = v2in(:)
  ! endif
  ! if (rank.eq.px-1) then
  !   this%k  (:,this%k1)= 2.0*this%k  (:,this%kmax)-this%k  (:,this%kmax-1)
  !   this%eps(:,this%k1)= 2.0*this%eps(:,this%kmax)-this%eps(:,this%kmax-1)
  !   this%v2 (:,this%k1)= 2.0*this%v2 (:,this%kmax)-this%v2 (:,this%kmax-1)
  ! endif

end subroutine set_bc_VF

subroutine solve_v2_VF(this,resV2,u,w,rho,mu,mui,muk,mut,rho_mod, &
                       Ru,Rp,dru,dRp,dz, &
                       alphav2,modification,rank,centerBC,periodic)
  implicit none
  class(VF_TurbModel) :: this
  real(8),dimension(0:this%i1,0:this%k1), intent(IN) :: u, w, rho,mu,mui,muk,mut,rho_mod
  real(8),dimension(0:this%i1),           intent(IN) :: dRu,ru,rp,dRp
  real(8),                                intent(IN) :: dz, alphav2
  integer,                                intent(IN) :: modification,rank,centerBC,periodic
  real(8),                                intent(OUT):: resV2
  real(8), dimension(0:this%i1,0:this%k1) :: dnew,dimpl
  real(8), dimension(this%imax)           :: a,b,c,rhs
  integer :: i,k
      
  dnew  = 0.0; dimpl = 0.0;

  call advecc(dnew,dimpl,this%v2,u,w,Ru,Rp,dru,dz,this%i1,this%k1,rank,periodic,.true.)
  call this%rhs_v2_VF(dnew,dimpl)    
  call diffc(dnew,this%v2,mu,mui,muk,mut,this%sigmak,rho,Ru,Rp,dru,dz,rank,modification)

  do k=1,this%kmax
    do i=1,this%imax
      
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

      b(i) = (rho_mod(i,k)*(-a(i)-c(i)) + dimpl(i,k))/alphav2

      a(i) = a(i)*rho_mod(i-1,k)
      c(i) = c(i)*rho_mod(i+1,k)

      rhs(i) =  dnew(i,k) + (1-alphav2)*b(i)*this%v2(i,k)
    enddo

    i=1
    b(i)=b(i)+centerBC*a(i)

    i=this%imax
    b(i) = b(i) - (c(i) /alphav2)
    rhs(i) = dnew(i,k) + (1-alphav2)*b(i)*this%v2(i,k)

    call matrixIdir(this%imax,a,b,c,rhs)

    do i=1,this%imax
      resV2 = resV2 + ((this%v2(i,k) - rhs(i))/(this%v2(i,k)+1.0e-20))**2.0
      this%v2(i,k) = min(2.0/3.0*this%k(i,k), max(rhs(i), 1.0e-8))
    enddo
  enddo
end subroutine solve_v2_VF

subroutine production_VF(this,u,w,temp,rho,mu,mut,beta,Rp,Ru,dRu,dRp,dz)
  implicit none
  class(VF_TurbModel) :: this
  real(8), dimension(0:this%i1,0:this%k1), intent(IN) :: u,w,temp,rho,mu,mut,beta
  real(8), dimension(0:this%i1),           intent(IN) :: Rp,Ru,dRu,dRp
  real(8),                                 intent(IN) :: dz
  real(8), dimension(0:this%i1,0:this%k1) :: div
  integer im,ip,jm,jp,km,kp,ib,ie,kb,ke,i,k
  real(8) :: Fr_1,ctheta,StR

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
      this%Tt(i,k)   = max(this%k(i,k)/this%eps(i,k), 6.0*(mu(i,k)/(rho(i,k)*this%eps(i,k)))**0.5)
            
      ! if (modVF.eq.1) then
      !   StR = (2.*(((W(i,k)-W(i,km))/dz)**2. + &
      !     ((U(i,k)-U(im,k))/(Ru(i)-Ru(im)))**2. + &
      !     ((U(i,k)+U(im,k))/(2.*Rp(i)))**2.) + &
      !     (((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4. &
      !     -(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/dRu(i) &
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

      this%Gk(i,k) = this%Gk(i,k) + ctheta*beta(i,k)*Fr_1*this%Tt(i,k)*2./3.*mut(i,k)*div(i,k)*(temp(i,kp)-temp(i,km))/(2.*dz)
    enddo
  enddo
end subroutine production_VF

subroutine rhs_v2_VF(this,putout,dimpl)
  implicit none
  class(VF_TurbModel) :: this
  real(8), dimension(0:this%i1,0:this%k1), intent(OUT):: putout,dimpl
  integer ib,ie,kb,ke,i,k

  ib = 1
  ie = this%i1-1

  kb = 1
  ke = this%k1-1
      
  do k=kb,ke
    do i=ib,ie
      !v'2 equation
      putout(i,k) = putout(i,k) + this%k(i,k)*this%fv2(i,k)       ! note, source is rho*k*f/rho
      dimpl(i,k)  = dimpl(i,k)  + 6.*this%eps(i,k)/this%k(i,k)    ! note, 6*rho*v'2*epsilon/k/(rho*v'2), set implicit and divided by density
    enddo
  enddo

end subroutine rhs_v2_VF

subroutine fillhem(this,mu,rho)
  implicit none
  class(VF_TurbModel) :: this
  real(8), dimension(0:this%i1,0:this%k1), intent(IN) :: mu, rho
  real(8), dimension(0:this%i1,0:this%k1) :: Srsq
  integer i,k

  do  k=1,this%kmax
    do i=1,this%imax
         
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
