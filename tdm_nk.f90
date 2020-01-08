module nk_tdm
  use ktet_tdm, only : KtEt_TurbDiffModel
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
    procedure :: rhs_kt_KtEt => rhs_kt_KtEt_NK
    procedure :: set_bc => set_bc_NK
  end type NK_TurbDiffModel


contains
!****************************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************************


  !************************!
  !       NK routines      !
  !************************!

type(NK_TurbDiffModel) function init_NK_TurbDiffModel(name)
  character(len=2), intent(IN) :: name
  init_NK_TurbDiffModel%name=name
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
  this%A = 30.5
end subroutine set_constants_NK

subroutine set_BC_NK(this,ekh,rho,periodic,rank,px)
  use mod_param, only : k1,i1,kmax,imax,k
  use mod_mesh,  only : top_bcvalue,bot_bcvalue, top_bcnovalue,bot_bcnovalue,walldist
  use mod_param, only : isothermalBC
  implicit none
  class(NK_TurbDiffModel) :: this
  real(8),dimension(0:i1,0:k1),intent(IN) :: rho,ekh
  integer,                     intent(IN) :: periodic, rank, px
  real(8),dimension(0:i1) :: tmp
  real(8) :: topBCvalue, botBCvalue
  
  !isothermal
  if (isothermalBC.eq.1) then
    do k = 0,k1 
      this%kt(0,k)    = bot_bcnovalue(k)*this%kt(1,k)         !dkt/dy = 0 (1) | or kt=0 (-1) 
      this%kt(i1,k)   = top_bcnovalue(k)*this%kt(imax,k)      !dkt/dy = 0 (1) | or kt=0 (-1)
      this%epst(0,k)  = bot_bcnovalue(k)*this%epst(1,k)       !depst/dy =0 (symmetry:1) | or epst =0 (wall:-1)
      this%epst(i1,k) = top_bcnovalue(k)*this%epst(imax,k)    !depst/dy =0 (symmetry:1) | or epst =0 (wall:-1)
    enddo
  !isoflux
  else
    do k = 0,k1 
      this%kt(0,k)        =this%kt(1,k)   ! dkt/dy = 0 
      this%kt(i1,k)  =this%kt(imax,k)     ! dkt/dy = 0 
      this%epst(0,k)      =this%epst(1,k) ! depst/dy = 0 
      this%epst(i1,k)=this%epst(imax,k)   ! depst/dy = 0 
    enddo
  endif

  call shiftf(this%kt,  tmp,rank); this%kt  (:,0) =tmp(:);
  call shiftf(this%epst,tmp,rank); this%epst(:,0) =tmp(:);
  call shiftb(this%kt,  tmp,rank); this%kt  (:,k1)=tmp(:);
  call shiftb(this%epst,tmp,rank); this%epst(:,k1)=tmp(:);
  
  ! developing
  if (periodic.eq.1) return
  if (rank.eq.0) then
    this%kt  (:,0) = this%ktin(:)
    this%epst(:,0) = this%epstin(:)
    this%Pkt(:,0) = this%Pktin(:)
  endif
  if (rank.eq.px-1) then
    this%kt  (:,k1)= 2.0*this%kt  (:,kmax)-this%kt  (:,kmax-1)
    this%epst(:,k1)= 2.0*this%epst(:,kmax)-this%epst(:,kmax-1)
    this%Pkt(:,k1) = 2.0*this%Pkt (:,kmax)-this%Pkt (:,kmax-1)
  endif

end subroutine set_BC_NK

subroutine set_alphat_NK(this,u,w,rho,temp,mu,mui,lam_cp,mut,alphat)
  use mod_param, only : k1,i1,kmax,imax,i,k,Pr,Re
  use mod_tm,    only : turb_model
  use mod_mesh,  only : dzp,walldist,drp,dzp
  use mod_common,only : cp, cnew, ekhi,cpi
  implicit none
  class(NK_TurbDiffModel) :: this
  real(8),dimension(0:i1,0:k1),intent(IN) :: u,w,rho,temp,mu,mui,lam_cp,mut
  real(8),dimension(0:i1,0:k1),intent(OUT):: alphat
  integer  im,ip,km,kp
  real(8) cf, st, utau, rhowall,twall,qwall,nu,phi,Reh
  real(8),dimension(0:k1) ::   tauw, Qwall_vec
  real(8),dimension(0:i1,0:k1) :: Ret, Reeps, yp
  real(8),dimension(0:i1,0:k1) :: kine, eps, Tt
  
  eps  = turb_model%eps
  kine = turb_model%k
  Tt   = turb_model%Tt

  do k=1,kmax
    km=k-1
    kp=k+1
    tauw(k) = mui(imax,k)*0.5*(w(imax,km)+w(imax,k))/walldist(imax)
    utau    = (tauw(k)/rho(imax,k))*0.5
    do i=1,imax
      im=i-1
      ip=i+1
      ! yplus hsould be an input so it can be changed from yplus to ystar
      ! this%yp(i,k) = sqrt(rho(i,k))/mu(i,k)*(walldist(i))*tauw(k)**0.5       
      nu = mu(i,k)/rho(i,k)
      this%yp(i,k) = walldist(i)*utau/nu

      rhowall = 0.5*(rho(i1,k)+rho(imax,k))
      twall=0.5*(temp(i1,k)+temp(imax,k))
      
      !skin friction: C_f = tau_w/(0.5*rho*U_r**2)
      !cfi =  2*tauw(k)/rho(imax,k)/utau                     ! Skin friction
      cf = tauw(k)/(0.5*rhowall*utau**2)
      !Stanton number: qwall/(rho*cp*U_r*(T_w-T_r))
      qwall= (ekhi(imax,k)/rhowall)*cpi(imax,k)*(temp(i1,k)-temp(imax,k))/drp(imax)
      st =  qwall/(rhowall*cpi(imax,k)*utau*twall)  
      
      phi = (2/cf)*st
      this%Ttemp(i,k)   = this%kt(i,k)/(this%epst(i,k)+1.0e-20)
      this%Tmix(i,k)    = (Tt(i,k) * this%Ttemp(i,k) )**0.5
      Reh = kine(i,k)*this%Tmix(i,k)

      this%flambda(i,k) =(1 - exp(-(Pr**0.5/this%A)*(2/cf)*st*this%yp(i,k)))**2.0      
      alphat(i,k) = rho(i,k)*this%clambda*this%flambda(i,k)*Reh
    enddo
  enddo

end subroutine set_alphat_NK

subroutine rhs_epst_KtEt_NK(this,putout,dimpl,temp,rho,mu,lam_cp,alphat)
  use mod_param,only : k1,i1,kmax,imax,i,k
  use mod_tm,   only : turb_model
  use mod_mesh, only : drp, dzp
  implicit none
  class(NK_TurbDiffModel) :: this
  real(8), dimension(0:i1,0:k1), intent(IN) :: temp,rho,mu,lam_cp,alphat
  real(8), dimension(0:i1,0:k1), intent(OUT):: putout,dimpl
  integer ip,im,kp,km
  real(8) d2Tdxdr
  real(8), dimension(0:i1,0:k1) :: Pk, kine, Tt

     ! ce2 is from the k-epsilon model
  Pk   = turb_model%Pk
  kine = turb_model%k
  Tt   = turb_model%Tt
  
  do k=1,kmax
    kp=k+1
    km=k-1
    do i=1,imax
      ip=i+1
      im=i-1
      ! d2Tdxdr = (1-this%flambda(i,k))*rho(i,k)*alphat(i,k)*ekh(i,k) * &
      !        ((((temp(i,k)-temp(i,km))/dz)-((temp(im,k)-temp(im,km))/dz) )/dRu(i))**2
      d2Tdxdr = (1-this%flambda(i,k))*rho(i,k)*alphat(i,k)*lam_cp(i,k) * &
             ((((temp(i,k)-temp(i,km))/dzp(km))-((temp(im,k)-temp(im,km))/dzp(km)) )/dRp(im))**2

      putout(i,k) = putout(i,k) + (this%cp1*this%epst(i,k)/(this%kt(i,k)+1.0e-20)*this%Pkt(i,k) + &
                   this%cp2*this%epst(i,k)/ kine(i,k)* Pk(i,k) + d2Tdxdr)     /rho(i,k)
      dimpl(i,k)  = dimpl(i,k)  + this%cd1*this%epst(i,k)/(this%kt(i,k)+1.0e-20) + this%cd2/Tt(i,k)           

    enddo
  enddo
end subroutine rhs_epst_KtEt_NK

subroutine rhs_kt_KtEt_NK(this,putout,dimpl,rho)
   use mod_param, only : k1,i1,kmax,imax,k,i
   use mod_mesh, only : drp,dzp
   use mod_common, only : ekh
  implicit none
  class(NK_TurbDiffModel) :: this
  real(8), dimension(0:i1,0:k1), intent(IN) :: rho
  real(8), dimension(0:i1,0:k1), intent(OUT):: putout,dimpl
  real(8) :: dsqrttdx,dsqrttdy
  integer :: im,km,ip,kp
  !kt equation
  do k=1,kmax
    kp = k+1
    km = k-1
    do i=1,imax
      ip = i+1
      im = i-1
      dsqrttdx = (sqrt(this%kt(i,kp))-sqrt(this%kt(i,km)))/(dzp(k)+dzp(km))
      dsqrttdy = (sqrt(this%kt(ip,k))-sqrt(this%kt(im,k)))/(drp(i)+drp(im))
      
      putout(i,k) = putout(i,k)+2.0*this%Pkt(i,k)/rho(i,k)-2*ekh(i,k)*(dsqrttdy**2+dsqrttdx**2) 
      dimpl(i,k)  = dimpl(i,k) +2.0*this%epst(i,k)/(this%kt(i,k)+1.0e-20) ! note, rho*epsilon/(rho*k), set implicit and divided by density
    enddo
  enddo
end subroutine rhs_kt_KtEt_NK


end module nk_tdm