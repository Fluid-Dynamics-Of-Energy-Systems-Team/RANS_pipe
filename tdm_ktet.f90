module ktet_tdm
  use mod_tdm, only : TurbDiffModel
  implicit none

!****************************************************************************************

  !************************!
  !         KtEt class       !
  !************************!
  
  type,abstract,extends(TurbDiffModel), public :: KtEt_TurbDiffModel
  real(8), dimension(:,:), allocatable :: flambda, resKt, resEt
  real(8), dimension(:),   allocatable :: epstin, ktin
  real(8) :: sigmakt,sigmaet,clambda,cp1,cp2,cd1,cd2,A
  contains
    procedure(set_alpha_KtEt), deferred :: set_alphat
    procedure(rhs_epst_KtEt), deferred :: rhs_epst_KtEt
    procedure(rhs_kt_KtEt), deferred :: rhs_kt_KtEt
    procedure(set_constants), deferred :: set_constants
    procedure :: advance_turbdiff => advance_KtEt
    procedure :: get_sol => get_sol_KtEt
    procedure :: init_w_inflow => init_w_inflow_KtEt
    procedure :: init => init_KtEt
    procedure :: diffusion_epst_KtEt
    procedure :: solve_epst_KtEt
    procedure :: production_KtEt
    procedure :: solve_kt_KtEt
    procedure :: init_mem_KtEt
    procedure :: init_sol => init_sol_KtEt
    procedure :: get_profile => get_profile_KtEt
  end type KtEt_TurbDiffModel

  interface

    subroutine set_constants(this)
      import :: KtEt_TurbDiffModel
      class(KtEt_TurbDiffModel) :: this
    end subroutine set_constants

    subroutine set_alpha_KtEt(this,u,w,rho,temp,mu,mui,lam_cp,mut,alphat)
      use mod_param, only : k1,i1
      import :: KtEt_TurbDiffModel
      implicit none
      class(KtEt_TurbDiffModel) :: this
      real(8),dimension(0:i1,0:k1),intent(IN) :: u,w,rho,temp,mu,mui,lam_cp,mut
      real(8),dimension(0:i1,0:k1),intent(OUT):: alphat
    end subroutine set_alpha_KtEt

    subroutine rhs_epst_KtEt(this,putout,dimpl,temp,rho,mu,lam_cp,alphat)
      use mod_param, only : k1,i1
      import :: KtEt_TurbDiffModel
      implicit none
      class(KtEt_TurbDiffModel) :: this
      real(8), dimension(0:i1,0:k1), intent(IN) :: rho,mu,temp,lam_cp,alphat
      real(8), dimension(0:i1,0:k1), intent(OUT):: putout,dimpl
    end subroutine rhs_epst_KtEt
      subroutine rhs_kt_KtEt(this,putout,dimpl,rho)
      use mod_param, only : k1,i1
      import :: KtEt_TurbDiffModel
      implicit none
      class(KtEt_TurbDiffModel) :: this
      real(8), dimension(0:i1,0:k1), intent(IN) :: rho
      real(8), dimension(0:i1,0:k1), intent(OUT):: putout,dimpl
    end subroutine rhs_kt_KtEt
  end interface

contains
!****************************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************************

  !************************!
  !       KE routines      !
  !************************!

subroutine init_KtEt(this)   
  class(KtEt_TurbDiffModel) :: this
  call this%init_mem_KtEt()
  call this%set_constants()
  call this%init_sol()
end subroutine init_KtEt

subroutine init_sol_KtEt(this)   
  use mod_param, only : i1,i
  implicit none
  class(KtEt_TurbDiffModel) :: this
  do i=0,i1
    this%kt(i,:)  =.01
    this%epst(i,:)=1.0
    this%Pkt(i,:) = 0.1 
    this%ktin(i) = 0.0
    this%epstin(i) = 0.0
    this%Pktin(i) = 0 
    this%ResKt(i,:) = 0.
    this%ResEt(i,:) = 0.
  enddo
end subroutine init_sol_KtEt

subroutine init_mem_KtEt(this)  
  use mod_param, only : k1,i1
  class(KtEt_TurbDiffModel) :: this
  allocate(this%epst(0:i1,0:k1),this%kt (0:i1,0:k1),     &
           this%Pkt(0:i1,0:k1), this%flambda(0:i1,0:k1), &
           this%Ttemp (0:i1,0:k1),this%Tmix(0:i1,0:k1),  &
           this%yp(0:i1,0:k1),this%resEt(0:i1,0:k1),     &
           this%resKt(0:i1,0:k1))
           
  allocate(this%Pktin (0:i1), &
           this%epstin(0:i1),this%ktin(0:i1),this%alphatin(0:i1))
end subroutine init_mem_KtEt

subroutine get_sol_KtEt(this,Prt,epst,kt, Pkt, resKt, resEt)
  use mod_param, only : k1,i1
  class(KtEt_TurbDiffModel) :: this
  real(8),dimension(0:i1,0:k1), intent(OUT):: Prt,epst,kt, Pkt,resKt,resEt
  Prt  =0
  epst =this%epst
  kt   =this%kt
  Pkt  = this%Pkt
  resKt = this%resKt
  resEt = this%resEt
end subroutine get_sol_KtEt

subroutine init_w_inflow_KtEt(this,Prtin,alphatin,ktin,epstin,pktin)
  use mod_param, only : i1,k1,k
  implicit none
  class(KtEt_TurbDiffModel) :: this
  real(8), dimension(0:i1), intent(IN) :: Prtin,alphatin,ktin,epstin,pktin
  this%epstin(:)=epstin
  this%ktin(:)=ktin
  this%Pktin(:)=Pktin
  do k=0,k1
    this%epst(:,k) = this%epstin(:)
    this%kt(:,k) = this%ktin(:)
    this%Pkt(:,k) = this%Pktin(:)
  enddo
  this%ResKt = 0.
  this%ResEt = 0.
end subroutine init_w_inflow_KtEt

subroutine get_profile_KtEt(this,p_prt,p_kt,p_epst,p_Pkt,k)
  use mod_param, only : i1
  class(KtEt_TurbDiffModel) :: this
  integer,                               intent(IN) :: k
  real(8),dimension(0:i1),          intent(OUT):: p_prt, p_kt,p_epst,p_Pkt
  p_prt(:) = 0.0
  p_kt(:)   =this%kt(:,k)
  p_epst(:) =this%epst(:,k)
  p_Pkt(:)  =this%Pkt(:,k)
end subroutine get_profile_KtEt

subroutine advance_KtEt(this,u,w,c,temp,rho,mu,ekh,ekhi,ekhk,alphat, &
                        alpha1,alpha2,rank,periodic,residual1, residual2)
  use mod_param,  only : k1,i1
  use mod_common, only : cp
  class(KtEt_TurbDiffModel) :: this
  real(8), dimension(0:i1,0:k1),intent(IN) :: u,w,c,temp,rho,mu,ekh,ekhi,ekhk,alphat
  real(8),                                intent(IN) :: alpha1,alpha2
  integer,                                intent(IN) :: rank,periodic
  real(8),                                intent(OUT):: residual1,residual2
  
  call this%production_KtEt(c,temp,rho,cp,alphat)
  call this%solve_epst_KtEt(residual2,u,w,temp,rho,mu,ekh,ekhi,ekhk,alphat,alpha2,rank,periodic)
  call this%solve_kt_KtEt(residual1,u,w,rho,ekh,ekhi,ekhk,alphat,alpha1,rank,periodic)
end

subroutine solve_kt_KtEt(this,resKt,u,w,rho,ekh,ekhi,ekhk,alphat, &
                       alphakt,rank,periodic)
  use mod_param,only : k1,i1,kmax,imax,k,i,isothermalBC
  use mod_math, only : matrixIdir
  use mod_mesh, only : ru,rp,dz,dru,drp,top_bcnovalue, bot_bcnovalue
  implicit none
  class(KtEt_TurbDiffModel) :: this
  real(8),dimension(0:i1,0:k1), intent(IN) :: u, w, rho,ekh,ekhi,ekhk,alphat
  real(8),                                intent(IN) :: alphakt
  integer,                                intent(IN) :: rank,periodic
  real(8),                                intent(OUT):: resKt
  real(8), dimension(0:i1,0:k1) :: dnew,dimpl
  real(8), dimension(imax)           :: a,b,c,rhs

  resKt  = 0.0;  dnew  = 0.0; dimpl = 0.0;

  call advecc(dnew,dimpl,this%kt,u,w,rank,periodic,.true.)
  call this%rhs_kt_KtEt(dnew,dimpl,rho) 
  ! call diffc(dnew,this%kt,ekh,ekhi,ekhk,alphat,this%sigmakt,rho,0)
  call this%diffusion_epst_KtEt(dnew,this%kt,ekhk,alphat,this%sigmakt,rho)

  do k=1,kmax
    do i=1,imax
      a(i) = (ekhi(i-1,k)+0.25*(alphat(i,k)+alphat(i-1,k))*(rho(i,k)+rho(i-1,k))/this%sigmakt)
      a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/rho(i,k)
      c(i) = (ekhi(i  ,k)+0.25*(alphat(i,k)+alphat(i+1,k))*(rho(i,k)+rho(i+1,k))/this%sigmakt)
      c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/rho(i,k)
      b(i) = (-a(i)-c(i)) + dimpl(i,k)

      rhs(i) = dnew(i,k) + ((1-alphakt)/alphakt)*b(i)*this%kt(i,k)
    enddo
    if (isothermalBC.eq.1) then

      i=1
      b(i) = b(i) + bot_bcnovalue(k)*a(i) !symmetry = -1 ; wall = 1 
      rhs(i) = dnew(i,k) + ((1-alphakt)/alphakt)*b(i)*this%kt(i,k)
      i=imax
      b(i) = b(i) + top_bcnovalue(k)*c(i)
      rhs(i) = dnew(i,k) + ((1-alphakt)/alphakt)*b(i)*this%kt(i,k)

    else
      i=1
      b(i) = b(i) + a(i) !symmetry = -1 ; wall = 1 
      rhs(i) = dnew(i,k) + ((1-alphakt)/alphakt)*b(i)*this%kt(i,k)
      i=imax
      b(i) = b(i) + c(i)
      rhs(i) = dnew(i,k) + ((1-alphakt)/alphakt)*b(i)*this%kt(i,k)
    endif
    call matrixIdir(imax,a,b/alphakt,c,rhs)

    do i=1,imax
      this%resKt(i,k) = ((this%kt(i,k) - rhs(i)))!/(this%kt(i,k)+1.0e-20))**2.0
      resKt = resKt + ((this%kt(i,k) - rhs(i))/(this%kt(i,k)+1.0e-20))**2.0
      this%kt(i,k) = max(rhs(i), 1.0e-8)
    enddo
  enddo
end

subroutine solve_epst_KtEt(this,resEt,u,w,temp,rho,mu,ekh,ekhi,ekhk,alphat, &
                       alphaet,rank,periodic)
  use mod_param,only : k1,i1,kmax,imax,k,i,isothermalBC
  use mod_math, only : matrixIdir
  use mod_mesh, only : dru,drp,ru,rp,dz,bot_bcvalue,top_bcvalue,top_bcnovalue,bot_bcnovalue
  implicit none
  class(KtEt_TurbDiffModel) :: this
  real(8),dimension(0:i1,0:k1), intent(IN) :: u, w,temp,rho,mu,ekh,ekhi,ekhk,alphat
  real(8),                                intent(IN) :: alphaet
  integer,                                intent(IN) :: rank,periodic
  real(8),                                intent(OUT):: resEt
  real(8), dimension(0:i1,0:k1) :: dnew,dimpl
  real(8), dimension(imax)           :: a,b,c,rhs
  

  resEt  = 0.0; dnew  = 0.0; dimpl = 0.0;

  call advecc(dnew,dimpl,this%epst,u,w,rank,periodic,.true.)
  call this%rhs_epst_KtEt(dnew,dimpl,temp,rho,mu,ekh,alphat) 
  call this%diffusion_epst_KtEt(dnew,this%epst,ekhk,alphat,this%sigmaet,rho)

  do k=1,kmax
    do i=1,imax
      a(i) = (ekhi(i-1,k)+0.25*(alphat(i,k)+alphat(i-1,k))*(rho(i,k)+rho(i-1,k))/this%sigmaet)
      a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/rho(i,k)
      c(i) = (ekhi(i  ,k)+0.25*(alphat(i,k)+alphat(i+1,k))*(rho(i,k)+rho(i+1,k))/this%sigmaet)
      c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/rho(i,k)
      b(i) = ((-a(i)-c(i))+ dimpl(i,k)  )  
      rhs(i) = dnew(i,k) + ((1-alphaet)/alphaet)*b(i)*this%epst(i,k)
    enddo


    if (isothermalBC.eq.1) then
      i=1
      b(i) = b(i)+bot_bcvalue(k)*a(i)
      rhs(i) = dnew(i,k) - (1.-bot_bcvalue(k))*a(i)*this%epst(i-1,k) + ((1-alphaet)/alphaet)*b(i)*this%epst(i,k)   !wall with value

      i=imax
      b(i) = b(i)+top_bcvalue(k)*c(i)
      rhs(i) = dnew(i,k) - (1.-top_bcvalue(k))*c(i)*this%epst(i+1,k) + ((1-alphaet)/alphaet)*b(i)*this%epst(i,k)   !wall with value
    else
      i=1
      b(i) = b(i) + a(i) !symmetry = -1 ; wall = 1 
      rhs(i) = dnew(i,k) + ((1-alphaet)/alphaet)*b(i)*this%epst(i,k)
      
      i=imax
      b(i) = b(i) + c(i)
      rhs(i) = dnew(i,k) + ((1-alphaet)/alphaet)*b(i)*this%epst(i,k)
    endif



    call matrixIdir(imax,a,b/alphaet,c,rhs)
  
    do i=1,imax
      this%resEt(i,k) = ((this%epst(i,k) - rhs(i))/(this%epst(i,k)+1.0e-20))**2.0
      resEt = resEt + ((this%epst(i,k) - rhs(i))/(this%epst(i,k)+1.0e-20))**2.0
      this%epst(i,k) = max(rhs(i), 1.0e-8)
    enddo
  enddo
end subroutine solve_epst_KtEt

subroutine production_KtEt(this,c,temp,rho,cp,alphat)
  use mod_param,only : k1,i1,kmax,imax,k,i
  use mod_mesh, only : dzp,drp
  implicit none
  class(KtEt_TurbDiffModel) :: this
  real(8), dimension(0:i1,0:k1), intent(IN) :: c,temp,rho,cp,alphat
  real(8), dimension(0:i1,0:k1) :: div
  integer im,ip,km,kp
  real(8) :: dTdx,dTdy


  do k=1,kmax
    kp=k+1
    km=k-1
    do i=1,imax
      ip=i+1
      im=i-1
      ! Production of temperature fluctuations  Pkt= <u'j T'> dTdxj= (lambda_t/cp)/cp dCdxj dTdxj = = (lambda_t/cp)dTdxj^2
     dTdx = (temp(i,kp)-temp(i,km))/(dzp(k)+dzp(km))
     dTdy = (temp(ip,k)-temp(im,k))/(drp(i)+drp(im))
     this%Pkt(i,k) = rho(i,k)*alphat(i,k)*(dTdx*dTdx + dTdy*dTdy) !- 2*this%epst(i,k) this part is put implict !NOTE: CHANGED BY STEPHAN
    enddo
  enddo
end subroutine production_KtEt

! subroutine rhs_kt_KtEt(this,putout,dimpl,rho)
!   use mod_param, only : k1,i1,kmax,imax,k,i
!   implicit none
!   class(KtEt_TurbDiffModel) :: this
!   real(8), dimension(0:i1,0:k1), intent(IN) :: rho
!   real(8), dimension(0:i1,0:k1), intent(OUT):: putout,dimpl
  
!   !kt equation
!   do k=1,kmax
!     do i=1,imax
!       putout(i,k) = putout(i,k)+2.0*this%Pkt(i,k)/rho(i,k) !-2.0*(1/rho(i,k))*this%epst(i,k)
!       dimpl(i,k)  = dimpl(i,k) +2.0*(1/rho(i,k))*this%epst(i,k)/(this%kt(i,k)+1.0e-20) ! note, rho*epsilon/(rho*k), set implicit and divided by density
!       !dimpl(i,k)  = dimpl(i,k) +2.0*this%epst(i,k)/(this%kt(i,k)+1.0e-20) ! note, rho*epsilon/(rho*k), set implicit and divided by density
!     enddo
!   enddo
! end subroutine rhs_kt_KtEt

subroutine diffusion_epst_KtEt(this,putout,putin,ekhk,alphat,sigma,rho)
  use mod_param, only : k1,i1,kmax,imax,k,i
  use mod_mesh, only : dzw,dzp
  implicit none
  class(KtEt_TurbDiffModel) :: this
  real(8), dimension(0:i1, 0:k1), intent(IN) :: putin, ekhk, alphat, rho
  real(8),                        intent(IN) :: sigma
  real(8), dimension(0:i1, 0:k1), intent(OUT):: putout
  integer :: kp,km

  do k=1,kmax
    kp=k+1
    km=k-1
    do i=1,imax
      putout(i,k) = putout(i,k) + 1.0/rho(i,k)*( &
        ( (ekhk(i,k ) + 0.25*(alphat(i,k)+alphat(i,kp))*(rho(i,k)+rho(i,kp))/sigma)*(putin(i,kp)-putin(i,k ))/dzp(k) &
        - (ekhk(i,km) + 0.25*(alphat(i,k)+alphat(i,km))*(rho(i,k)+rho(i,km))/sigma)*(putin(i,k )-putin(i,km))/dzp(km) &
        )/dzw(k)   )
    enddo
  enddo
  
end subroutine diffusion_epst_KtEt




end module