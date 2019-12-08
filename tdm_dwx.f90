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

subroutine set_alphat_DWX(this,u,w,rho,temp,mu,mui,lam_cp,mut,alphat)
  use mod_tm, only : turb_model
  use mod_mesh, only : mesh
  implicit none
  class(DWX_TurbDiffModel) :: this
  real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,temp,mu,mui,lam_cp, mut
  real(8),dimension(0:this%i1,0:this%k1),intent(OUT):: alphat
  integer  im,ip,km,kp,i,k
  real(8),dimension(0:this%k1) ::   tauw
  real(8),dimension(0:this%i1,0:this%k1) :: Ret, Reeps, yp
  real(8), dimension(0:this%i1,0:this%k1) :: kine, eps, Tt
  real(8), dimension(1:this%imax) :: walldist

  walldist = mesh%walldist

  eps  = turb_model%eps
  kine = turb_model%k
  Tt = turb_model%Tt

  do k=1,this%kmax
    km=k-1
    kp=k+1
    tauw(k) = mui(this%imax,k)*0.5*(w(this%imax,km)+w(this%imax,k))/walldist(this%imax)
    do i=1,this%imax
      im=i-1
      ip=i+1
      ! yplus hsould be an input so it can be changed from yplus to ystar
      this%yp(i,k) = sqrt(rho(i,k))/mu(i,k)*(walldist(i))*tauw(k)**0.5       
  
      Ret(i,k)     = rho(i,k)*(kine(i,k)**2.)/(mu(i,k)*eps(i,k))        
      Reeps(i,k)   = walldist(i)*((mu(i,k)*eps(i,k)/rho(i,k))**0.25)*rho(i,k)/mu(i,k)   

      this%Ttemp(i,k)   = this%kt(i,k)/(this%epst(i,k)+1.0e-20)
      this%Tmix(i,k)    = (Tt(i,k) * this%Ttemp(i,k) )**0.5

      this%flambda(i,k) =((1 - exp(-Reeps(i,k)/16))**2.0)*(1+(3/(Ret(i,k)**0.75)))
                
      alphat(i,k) = rho(i,k)*this%clambda*this%flambda(i,k)*kine(i,k)*Tt(i,k)*(2.0*this%Ttemp(i,k)/Tt(i,k))**0.5
                 

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
      Ret(i,k)     = rho(i,k)*(kine(i,k)**2.)/(mu(i,k)*eps(i,k)) 
      Reeps(i,k)   = walldist(i)*((mu(i,k)*eps(i,k)/rho(i,k))**0.25)*rho(i,k)/mu(i,k)   

      fd1  = 1 - (exp(-Reeps(i,k)/1.7))**2.0
      feps = 1 - 0.3*exp(-((Ret(i,k)/6.5)**2.0))   
      fd2  = (1/0.9)*(ce2*feps-1.0)*(1 - (exp(-Reeps(i,k)/5.8))**2.0)  
      
      putout(i,k) = putout(i,k) + (this%cp1*this%Pkt(i,k)/this%Tmix(i,k))/rho(i,k)
      dimpl(i,k)  = dimpl(i,k)  + this%cd1*fd1/this%Ttemp(i,k) + this%cd2*fd2/Tt(i,k)
    enddo
  enddo
end subroutine rhs_epst_KtEt_DWX

end module dwx_tdm