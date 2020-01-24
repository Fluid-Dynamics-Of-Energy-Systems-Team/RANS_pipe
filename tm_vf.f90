module vf_tm
  use ke_tm
  implicit none

!****************************************************************************************

  !************************!
  !         VF class       !
  !************************!

  type, extends(KE_TurbModel), public :: VF_TurbModel
  real(8), allocatable :: c1,c2,cl,ceta
  contains
    procedure :: set_constants => set_constants_VF
    procedure :: set_mut => set_mut_VF
    procedure :: advance_turb => advance_VF
    procedure :: set_bc => set_bc_VF
    procedure :: init_w_inflow => init_w_inflow_VF
    ! procedure :: init_sol_KE => init_sol_VF
    procedure :: rhs_v2_VF
    procedure :: rhs_eps_KE => rhs_eps_VF

    procedure :: solve_v2_VF
    procedure :: fillhelm
    procedure :: calc_turbulent_timescale => calc_turbulent_timescale_VF
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
  this%cmu    = 0.22
  this%sigmak = 1.0
  this%sigmae = 1.3
  this%ce2    = 1.9
  this%cl     = 0.23
  this%ceta    = 70.
  this%ce1    = 1.4
  this%c1     = 1.4
  this%c2     = 0.3
  this%fv2 = 0.
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
  use mod_param, only :k1,i1,imax,kmax,i,k, modifDiffTerm,Re
  use mod_mesh, only : mesh,dzw,dzp,dru, walldist,rp
  implicit none
  class(VF_TurbModel) :: this
  real(8), dimension(0:i1,0:k1), intent(IN) :: u, w, rho, mu, mui
  real(8), dimension(0:i1,0:k1), intent(OUT):: mut
  real(8), dimension(0:k1) :: tauw, utau
  real(8) :: StR
  real(8) :: dz, mu_wall, rho_wall
  integer :: im,ip,km,kp

  do k=1,kmax
    km=k-1
    kp=k+1
    tauw(k) = mui(imax,k)*0.5*(w(imax,km)+w(imax,k))/walldist(imax)
    utau(k) = sqrt(tauw(k)/(0.5*(rho(imax,k)+rho(i1,k))))
    rho_wall = 0.5*(rho(i1,k)+rho(imax,k))
    mu_wall = 0.5*(mu(i1,k)+mu(imax,k))
    do i=1,imax
      im=i-1
      ip=i+1
      if (modifDiffTerm .eq. 1) then
        this%yp(i,k) = walldist(i)*sqrt(rho(i,k)/rho_wall)*(mu_wall/mu(i,k))*Re*utau(k)         ! ystar
      else
        this%yp(i,k) = walldist(i)*Re*utau(k)
      endif
      this%f1(i,k)  = 1.0 + 0.045*(this%k(i,k)/this%v2(i,k))**0.5
      mut(i,k) = min(1.,rho(i,k)*this%cmu*this%v2(i,k)*this%Tt(i,k))
    enddo
  enddo
end subroutine set_mut_VF

subroutine rhs_eps_VF(this,putout,dimpl,rho)
  use mod_param, only :k1,i1,imax,kmax,k,i
  implicit none
  class(VF_TurbModel) :: this
  real(8), dimension(0:i1,0:k1), intent(IN) :: rho
  real(8), dimension(0:i1,0:k1), intent(OUT):: putout,dimpl
  
  do k=1,kmax
    do i=1,imax
      putout(i,k) = putout(i,k) +(this%ce1*this%f1(i,k)*(this%Pk(i,k)+this%Gk(i,k)))/(this%Tt(i,k)*rho(i,k)) 
      dimpl(i,k)  = dimpl(i,k)  + this%ce2/this%Tt(i,k)   ! note, ce2*f2*rho*epsilon/T/(rho*epsilon), set implicit and divided by density
    enddo
  enddo
end subroutine rhs_eps_VF

subroutine rhs_v2_VF(this,putout,dimpl)
  use mod_param, only : k1,i1,imax,kmax,i,k
  implicit none
  class(VF_TurbModel) :: this
  real(8), dimension(0:i1,0:k1), intent(OUT):: putout,dimpl
      
  do k=1,kmax
    do i=1,imax
      putout(i,k) = putout(i,k) + this%k(i,k)*this%fv2(i,k)       ! note, source is rho*k*f/rho
      dimpl(i,k)  = dimpl(i,k)  + 6.*this%eps(i,k)/this%k(i,k)    ! note, 6*rho*v'2*epsilon/k/(rho*v'2), set implicit and divided by density
    enddo
  enddo

end subroutine rhs_v2_VF

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
  use mod_mesh, only : mesh,ru,rp,dru,drp,dz
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
  
  call this%production_KE(u,w,temp,rho,mu,mut,beta)
  call this%fillhelm(mu,rho)
  call SOLVEhelm(this%fv2,Ru,Rp,dRu,dRp,dz,rank,this%Lh,mesh%centerBC)
  ! call solvepois_cr(this%fv2,0,rank,mesh%centerBC)
  call this%solve_eps_KE(residual2,u,w,rho,mu,mui,muk,mut,beta,temp,rho_mod, &
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

subroutine calc_turbulent_timescale_VF(this,rho,mu)
  use mod_param, only : kmax,imax,i,k, i1,k1
  implicit none
  class(VF_TurbModel) :: this
  real(8), dimension(0:i1,0:k1), intent(IN) :: rho, mu
  do  k=1,kmax
    do i=1,imax  
      this%Tt(i,k) = max(this%k(i,k)/this%eps(i,k), 6.0*(mu(i,k)/(rho(i,k)*this%eps(i,k)))**0.5)           
    enddo
  enddo
end subroutine


subroutine fillhelm(this,mu,rho)
  use mod_param, only : k1,i1,kmax,imax,i,k
  implicit none
  class(VF_TurbModel) :: this
  real(8), dimension(0:i1,0:k1), intent(IN) :: mu, rho
  real(8), dimension(0:i1,0:k1) :: Srsq
  real(8), dimension(imax,kmax) :: rhs_fv2

  do  k=1,kmax
    do i=1,imax         

      this%Lh(i,k)=0.23*max(this%k(i,k)**1.5/this%eps(i,k),70.*((mu(i,k)/rho(i,k))**3./this%eps(i,k))**0.25)
      this%fv2(i,k)= - (1.4-1.)*(2./3.-this%v2(i,k)/this%k(i,k))/this%Tt(i,k) &
                - 0.3*(this%Pk(i,k)+this%Gk(i,k))/(rho(i,k)*this%k(i,k))-5.*this%v2(i,k)/(this%k(i,k)*this%Tt(i,k))
      this%fv2(i,k) = this%fv2(i,k)/this%Lh(i,k)**2.0 
    enddo
  enddo
end subroutine fillhelm


end module
