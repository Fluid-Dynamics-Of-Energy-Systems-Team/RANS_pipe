module dwx_tdm
  use ktet_tdm, only : KtEt_TurbDiffModel
  implicit none

!****************************************************************************************

  !************************!
  !         DWX class      !
  !************************!
  !!!!!!  Deng, B., Wu, W., and Xi, S. "A near-wall two-equation heat transfer model 
  !         for wall turbulent flows", International Journal of Heat and mass transfer (2000).
  type, extends(KtEt_TurbDiffModel), public :: DWX_TurbDiffModel
  contains
    procedure :: set_constants => set_constants_DWX
    procedure :: set_alphat => set_alphat_DWX
    procedure :: rhs_epst_KtEt => rhs_epst_KtEt_DWX
    procedure :: set_bc => set_bc_DWX
  end type DWX_TurbDiffModel


contains
!****************************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************************

  !************************!
  !       DWX routines     !
  !************************!

type(DWX_TurbDiffModel) function init_DWX_TurbDiffModel(i1,k1,imax,kmax,name)
  integer, intent(in) :: i1,k1,imax,kmax
  character(len=3), intent(IN) :: name
  init_DWX_TurbDiffModel%name=name
  init_DWX_TurbDiffModel%i1 = i1
  init_DWX_TurbDiffModel%k1 = k1
  init_DWX_TurbDiffModel%imax = imax
  init_DWX_TurbDiffModel%kmax = kmax
end function init_DWX_TurbDiffModel

subroutine set_constants_DWX(this)
  class(DWX_TurbDiffModel) :: this
  this%sigmakt = 1.0
  this%sigmaet = 1.0
  this%clambda = 0.1
  this%cp1 = 2.34
  this%cp2 = 1.0
  this%cd1 = 2.0
  this%cd2 = 0.9
end subroutine set_constants_DWX

subroutine set_bc_DWX(this,ekh,rho,periodic,rank,px)
  use mod_mesh, only : mesh
  use mod_param, only : isothermalBC
  use mod_common, only : cpi, ekhi
  implicit none
  class(DWX_TurbDiffModel) :: this
  real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: rho,ekh
  integer,                               intent(IN) :: periodic, rank, px
  real(8),dimension(0:this%k1) ::  top_bcvalue, bot_bcvalue,top_bcnovalue, bot_bcnovalue
  real(8), dimension(1:this%imax) :: walldist
  real(8),dimension(0:this%i1) :: tmp
  real(8) :: topBCvalue, botBCvalue
  integer :: k

  walldist = mesh%walldist
  top_bcvalue = mesh%top_bcvalue
  bot_bcvalue = mesh%bot_bcvalue
  top_bcnovalue = mesh%top_bcnovalue
  bot_bcnovalue = mesh%bot_bcnovalue
  !isothermal
  if (isothermalBC.eq.1) then
    do k = 0,this%k1 
      this%kt(0,k)        = bot_bcnovalue(k)*this%kt(1,k)         !dkt/dy = 0 (1) | or kt=0 (-1) 
      this%kt(this%i1,k)  = top_bcnovalue(k)*this%kt(this%imax,k) !dkt/dy = 0 (1) | or kt=0 (-1)
      this%Pkt(0,k)       = bot_bcnovalue(k)*this%Pkt(1,k)
      this%Pkt(this%i1,k) = top_bcnovalue(k)*this%Pkt(this%imax,k)
      !botBCvalue = 2.0*ekh(1,k)/rho(1,k)*this%kt(1,k)/walldist(1)**2                                                          !bcvalue
      botBCvalue = ekhi(1,k)/(0.5*(rho(0,k)+rho(1,k))) &
                    *((this%kt(1,k)**0.5)/walldist(1))**2                                                    !NOTE: CHANGE BY STEPHAN
      this%epst(0,k)       = (1.-bot_bcvalue(k))*(2.0*botBCvalue-this%epst(1,k))         +bot_bcvalue(k)*this%epst(1,k)        !symmetry or bc value
      

      !topBCvalue = 2.0*ekh(this%imax,k)/rho(this%imax,k)*this%kt(this%imax,k)/walldist(this%imax)**2                          !bcvalue
      topBCvalue = ekhi(this%imax,k)/(.5*(rho(this%i1,k)+rho(this%imax,k))) &   
                    *((this%kt(this%imax,k)**0.5)/walldist(this%imax))**2
      this%epst(this%i1,k) = (1.-top_bcvalue(k))*(2.0*topBCvalue-this%epst(this%imax,k)) +top_bcvalue(k)*this%epst(this%imax,k)!symmetry or bc value
    enddo
  !isoflux
  else
    do k= 0,this%k1 
      ! this%kt(0,k)        =this%kt(1,k)         ! dkt/dy = 0 
      ! this%kt(this%i1,k)  =this%kt(this%imax,k) ! dkt/dy = 0 
      ! this%epst(0,k)      =this%epst(1,k)         ! depst/dy = 0 
      this%kt(0,k)        = this%kt(1,k)         !dkt/dy = 0 (1) | or kt=0 (-1) 
      this%kt(this%i1,k)  = this%kt(this%imax,k) !dkt/dy = 0 (1) | or kt=0 (-1)
      ! this%kt(0,k)        = this%kt(1,k)         !dkt/dy = 0 (1) | or kt=0 (-1) 
      ! this%kt(this%i1,k)  = this%kt(this%imax,k) !dkt/dy = 0 (1) | or kt=0 (-1)
      botBCvalue = 2.0*ekh(1,k)/rho(1,k)*(this%kt(1,k)**0.5/walldist(1))**2                                                    !NOTE: CHANGE BY STEPHAN
      this%epst(0,k)       = (2.0*botBCvalue-this%epst(1,k))              !symmetry or bc value
      topBCvalue = 2.0*ekh(this%imax,k)/rho(this%imax,k)*(this%kt(this%imax,k)**0.5/walldist(this%imax))**2
      this%epst(this%i1,k) = (2.0*topBCvalue-this%epst(this%imax,k))!symmetry or bc value

      this%Pkt(0,k)       =this%Pkt(1,k)
      this%Pkt(this%i1,k) =this%Pkt(this%imax,k)
    enddo
  endif

  call shiftf(this%kt,  tmp,rank); this%kt  (:,0)      =tmp(:);
  call shiftf(this%epst,tmp,rank); this%epst(:,0)      =tmp(:);
  call shiftf(this%Pkt, tmp,rank); this%Pkt (:,0)      =tmp(:);
  call shiftb(this%kt,  tmp,rank); this%kt  (:,this%k1)=tmp(:);
  call shiftb(this%epst,tmp,rank); this%epst(:,this%k1)=tmp(:);
  call shiftb(this%Pkt, tmp,rank); this%Pkt (:,this%k1)=tmp(:);
  ! developing
  if (periodic.eq.1) return
  if (rank.eq.0) then
    this%kt  (:,0) = this%ktin(:)
    this%epst(:,0) = this%epstin(:)
    this%Pkt(:,0)  = this%Pktin(:)
  endif
  if (rank.eq.px-1) then
    this%kt  (:,this%k1)=2.0*this%kt  (:,this%kmax)-this%kt  (:,this%kmax-1)
    this%epst(:,this%k1)=2.0*this%epst(:,this%kmax)-this%epst(:,this%kmax-1)
    this%Pkt (:,this%k1)=2.0*this%Pkt (:,this%kmax)-this%Pkt (:,this%kmax-1)
  endif
 
end subroutine set_bc_DWX


subroutine set_alphat_DWX(this,u,w,rho,temp,mu,mui,lam_cp,mut,alphat)
  use mod_tm, only : turb_model
  use mod_mesh, only : mesh
  use mod_common, only : cpi, ekhi
  implicit none
  class(DWX_TurbDiffModel) :: this
  real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,temp,mu,mui,lam_cp, mut
  real(8),dimension(0:this%i1,0:this%k1),intent(OUT):: alphat
  integer  im,ip,km,kp,i,k
  real(8),dimension(0:this%i1,0:this%k1) :: Ret, Reeps, yp
  real(8), dimension(0:this%i1,0:this%k1) :: kine, eps, Tt
  real(8), dimension(1:this%imax) :: walldist
  real(8) :: nu

  walldist = mesh%walldist

  eps  = turb_model%eps
  kine = turb_model%k
  Tt = turb_model%Tt

  do k=1,this%kmax
    km=k-1
    kp=k+1
    do i=1,this%imax
      im=i-1
      ip=i+1
      nu = mu(i,k)/rho(i,k)
      Ret(i,k)     = (kine(i,k)**2.)/(nu*eps(i,k))                          !k^2/(eps*nu)
      Reeps(i,k)   = (walldist(i)*(nu*eps(i,k))**0.25)/nu                   !y*(nu*eps)^(1/4)/nu
      this%Ttemp(i,k)   = this%kt(i,k)/(this%epst(i,k)+1.0e-20)             !kt/epst      
      this%flambda(i,k) =((1 - exp(-Reeps(i,k)/16.))**2.0)*(1+(3/(Ret(i,k)**0.75)))     !f_lambda=(1-exp(Reps/16))^2 * (1+3/Rt^(3/4))
      alphat(i,k) = rho(i,k)*this%clambda*this%flambda(i,k)*(kine(i,k)**2/eps(i,k))*(2.0*this%Ttemp(i,k)/Tt(i,k))**0.5              
    enddo
  enddo

end subroutine set_alphat_DWX

subroutine rhs_epst_KtEt_DWX(this,putout,dimpl,temp,rho,mu,lam_cp,alphat)
  use mod_tm, only : turb_model
  use mod_mesh, only : mesh
  implicit none
  class(DWX_TurbDiffModel) :: this
  real(8), dimension(0:this%i1,0:this%k1), intent(IN) :: rho,mu,temp,lam_cp,alphat
  real(8), dimension(0:this%i1,0:this%k1), intent(OUT):: putout,dimpl
  real(8),dimension(0:this%i1,0:this%k1) :: Reeps,Ret
  integer ib,ie,kb,ke,i,k 
  real(8) ce2,fd1,feps,fd2
  real(8), dimension(0:this%i1,0:this%k1) :: kine, eps,Tt
  real(8), dimension(1:this%imax) :: walldist
  real(8) :: nu
  walldist = mesh%walldist

  eps  = turb_model%eps
  kine = turb_model%k
  ce2  = turb_model%ce2
  Tt   = turb_model%Tt
  ib = 1
  ie = this%i1-1

  kb = 1
  ke = this%k1-1

  do k=kb,ke
    do i=ib,ie
      nu = mu(i,k)/rho(i,k)
      Ret(i,k)     = (kine(i,k)**2.)/(nu*eps(i,k))              !k^2/(eps*nu)
      Reeps(i,k)   = (walldist(i)*(nu*eps(i,k))**0.25)/nu       !y*(nu*eps)^(1/4)/nu
      fd1  = 1 - (exp(-Reeps(i,k)/1.7))**2.0                    ! (1-exp(-R_eps/1.7))^2
      feps = 1 - 0.3*exp(-(Ret(i,k)/6.5)**2.0)                  ! (1-0.3*exp(-(Ret/6.5)^2))
      fd2  = (1/this%cd2)*(ce2*feps-1.0)*(1 - (exp(-Reeps(i,k)/5.8))**2.0)  
      !putout(i,k) = putout(i,k) + (this%cp1*this%Pkt(i,k)/this%Tmix(i,k))/rho(i,k)
      this%Tmix(i,k)    = Tt(i,k) * this%Ttemp(i,k)                                       !tau_u*tau_t = k*kt/(epst*eps)
      putout(i,k) = putout(i,k) + ((this%cp1*(1/this%Tmix(i,k))**0.5)*this%Pkt(i,k))/rho(i,k) !NOTE: CHANGED BY STEPHAN ! -cp1*sqrt(epst*eps/kt*k)*alphat*dTdxi*dTdxi
      dimpl(i,k)  = dimpl(i,k)  + this%cd1*fd1*(1/this%Ttemp(i,k)) + this%cd2*fd2*(1/Tt(i,k))
    enddo
  enddo
end subroutine rhs_epst_KtEt_DWX

end module dwx_tdm