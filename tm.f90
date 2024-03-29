module mod_tm
  implicit none

!****************************************************************************************

  !************************!
  !     Abstract class     !
  !************************!
  
  type, abstract, public :: TurbModel
  character(len=3)                     :: name
  real(8), dimension(:,:), allocatable :: nuSA,Pk,om,k,bF1,bF2,eps,v2,yp, Tt
  real(8), dimension(:),   allocatable :: mutin, Pkin
  real(8), allocatable :: ce1,ce2
  contains
    procedure(init_tm), deferred :: init
    procedure(set_mut_tm), deferred :: set_mut
    procedure(advance_turb_tm), deferred :: advance_turb
    procedure(set_bc_tm), deferred :: set_bc
    procedure(init_sol_tm), deferred :: init_sol
    procedure(get_profile_tm), deferred :: get_profile
    procedure(init_w_inflow_tm), deferred :: init_w_inflow
    procedure(get_sol_tm), deferred :: get_sol
    procedure :: set_mut_bc
  end type TurbModel

  interface
    subroutine init_tm(this)
      import :: TurbModel
      class(TurbModel) :: this
    end subroutine init_tm
    subroutine init_sol_tm(this)
      import :: TurbModel
      class(TurbModel) :: this
    end subroutine init_sol_tm
    subroutine set_mut_tm(this,u,w,rho,mu,mui,mut)
      use mod_param, only : i1,k1
      import :: TurbModel
      class(TurbModel) :: this
      real(8), dimension(0:i1,0:k1),intent(IN) :: u,w,rho,mu,mui
      real(8), dimension(0:i1,0:k1),intent(OUT):: mut
    end subroutine set_mut_tm
    subroutine advance_turb_tm(this,u,w,rho,mu,mui,muk,mut,beta,temp,&
                             alpha1,alpha2,alpha3,                   &
                             modification,rank,periodic,    &
                             residual1, residual2, residual3)
      use mod_param, only : i1,k1
      import :: TurbModel
      class(TurbModel) :: this
      real(8), dimension(0:i1,0:k1),intent(IN) :: u,w,rho,mu,mui,muk,mut,beta,temp
      real(8),                                intent(IN) :: alpha1,alpha2,alpha3
      integer,                                intent(IN) :: modification,rank,periodic
      real(8),                                intent(OUT):: residual1,residual2,residual3
    end subroutine advance_turb_tm
    subroutine set_bc_tm(this,mu,rho,periodic,rank,px)
      use mod_param, only : i1,k1
      import :: TurbModel
      class(TurbModel) :: this
      real(8),dimension(0:i1,0:k1),intent(IN) :: rho,mu
      integer,                     intent(IN) :: periodic, rank, px
    end subroutine set_bc_tm
    subroutine get_profile_tm(this,p_nuSA,p_k,p_eps,p_om,p_v2,p_Pk,p_bF1,p_bF2,yp,k)
      use mod_param, only : i1,imax
      import :: TurbModel
      class(TurbModel) :: this
      integer,                               intent(IN) :: k
      real(8),dimension(0:i1),          intent(OUT):: p_nuSA,p_k,p_eps,p_om,p_v2,p_Pk,p_bF1,yp
      real(8),dimension(1:imax),        intent(OUT):: p_bF2
    end subroutine get_profile_tm
    subroutine get_sol_tm(this,nuSA,k,eps,om,v2,pk, gk,yp)
      use mod_param, only : i1,k1
      import :: TurbModel
      class(TurbModel) :: this
      real(8),dimension(0:i1,0:k1),intent(OUT):: nuSA,k,eps,om,v2,yp,pk,gk
    end subroutine get_sol_tm
    subroutine init_w_inflow_tm(this,nuSAin,pkin,kin,epsin,omin,mutin,v2in)
      use mod_param, only : i1
      import :: TurbModel
      class(TurbModel) :: this
      real(8), dimension(0:i1), intent(IN) :: nuSAin,pkin,kin,epsin,omin,mutin,v2in
    end subroutine

  end interface

class(TurbModel), allocatable :: turb_model
!****************************************************************************************

  !************************!
  !      Laminar class     !
  !************************!
  
  type, extends(TurbModel), public :: Laminar_TurbModel
  contains
    procedure :: init => init_laminar
    procedure :: init_mem => init_mem_laminar
    procedure :: set_mut  => set_mut_laminar
    procedure :: advance_turb => advance_laminar
    procedure :: init_sol => init_sol_laminar
    procedure :: set_bc => set_bc_laminar
    procedure :: get_profile => get_profile_laminar
    procedure :: init_w_inflow => init_w_inflow_laminar
    procedure :: get_sol => get_sol_laminar
  end type Laminar_TurbModel

contains
!****************************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************************

   !************************!
   !    Abstract routines    !
   !************************!

subroutine set_mut_bc(this,mut,periodic,px,rank)
  use mod_param, only : i1,k1,imax,kmax
  use mod_mesh, only : top_bcnovalue, bot_bcnovalue
  class(TurbModel) :: this
  integer,                       intent(IN) :: periodic,px,rank
  real(8), dimension(0:i1,0:k1), intent(OUT):: mut
  real(8), dimension(0:i1)                  :: tmp
 
  mut(i1,:) = top_bcnovalue(:)*mut(imax,:)
  mut(0,:)       = bot_bcnovalue(:)*mut(1,:)

  call shiftf(mut,tmp,rank); mut(:,0)      =tmp(:);
  call shiftb(mut,tmp,rank); mut(:,k1)=tmp(:);

  if ((periodic.ne.1).and.(rank.eq.0)) then
    mut(:,0) = this%mutin(:)
  endif
  if ((periodic.ne.1).and.(rank.eq.px-1)) then
    mut(:,k1) = 2.*mut(:,kmax)-mut(:,kmax-1)
  endif

end subroutine

!****************************************************************************************

   !************************!
   !    Laminar routines    !
   !************************!

subroutine init_laminar(this)
  class(Laminar_TurbModel) :: this
  this%name='lam'
  call this%init_mem()
end subroutine init_laminar

subroutine set_bc_laminar(this,mu,rho,periodic,rank,px)
  use mod_param, only : i1,k1,imax,kmax
  class(Laminar_TurbModel) :: this
  real(8),dimension(0:i1,0:k1),intent(IN) :: rho,mu
  integer,intent(IN) :: periodic, rank, px
end subroutine set_bc_laminar

subroutine init_sol_laminar(this)
  class(Laminar_TurbModel) :: this
end subroutine init_sol_laminar

subroutine init_w_inflow_laminar(this,nuSAin,pkin,kin,epsin,omin,mutin,v2in)
  use mod_param, only : i1
  class(Laminar_TurbModel) :: this
  real(8), dimension(0:i1), intent(IN) :: nuSAin,pkin,kin,epsin,omin,mutin,v2in
  this%mutin = mutin
end subroutine init_w_inflow_laminar

subroutine init_mem_laminar(this)
  use mod_param, only : i1,k1
  class(Laminar_TurbModel) :: this
  allocate(this%yp(0:i1,0:k1))
  allocate(this%mutin(0:i1))
end subroutine init_mem_laminar

subroutine get_profile_laminar(this,p_nuSA,p_k,p_eps,p_om,p_v2,p_Pk,p_bF1,p_bF2,yp,k)
  use mod_param, only : i1,imax
  class(Laminar_TurbModel) :: this
  integer,                               intent(IN) :: k
  real(8),dimension(0:i1),          intent(OUT):: p_nuSA,p_k,p_eps,p_om,p_v2,p_Pk,p_bF1,yp
  real(8),dimension(1:imax),        intent(OUT):: p_bF2
  p_nuSA(:)=0
  p_k(:)   =0
  p_eps(:) =0
  p_v2(:)  =0
  p_om(:)  =0
  p_Pk(:)  =0
  p_bF1(:) =0
  p_bF2(:) =0
  yp(:) = this%yp(:,k)
end subroutine get_profile_laminar

subroutine get_sol_laminar(this,nuSA,k,eps,om,v2,pk, gk,yp)
  use mod_param, only : i1,k1
  class(Laminar_TurbModel) :: this
  real(8),dimension(0:i1,0:k1), intent(OUT):: nuSA,k,eps,om,v2,yp, pk,gk
  nuSA=0
  k   =0    
  eps =0
  v2  =0
  om  =0
  yp  = this%yp
  pk=0.
  gk=0.
end subroutine get_sol_laminar


subroutine set_mut_laminar(this,u,w,rho,mu,mui,mut)
  use mod_param, only : i1,k1,imax,kmax,i,k
  use mod_mesh, only : walldist
  class(Laminar_TurbModel) :: this
  real(8), dimension(0:i1,0:k1),intent(IN) :: u,w,rho,mu,mui
  real(8), dimension(0:i1,0:k1),intent(OUT):: mut
  real(8), dimension(0:k1)                 :: tauw
  do k=1,kmax
    tauw(k) = mui(imax,k)*0.5*(w(imax,k)+w(imax,k))/walldist(imax)
    do i=1,imax
      this%yp(i,k) = sqrt(rho(i,k))/mu(i,k)*(walldist(i))*tauw(k)**0.5   ! ystar
    enddo
  enddo
  mut = 0
end subroutine set_mut_laminar

subroutine advance_laminar(this,u,w,rho,mu,mui,muk,mut,beta,temp,&
                             alpha1,alpha2,alpha3,               &
                             modification,rank,periodic,&
                             residual1, residual2, residual3)
  use mod_param, only : i1,k1,imax,kmax
  class(Laminar_TurbModel) :: this
  real(8), dimension(0:i1,0:k1),intent(IN) :: u,w,rho,mu,mui,muk,mut,beta,temp
  real(8),                      intent(IN) :: alpha1,alpha2,alpha3
  integer,                      intent(IN) :: modification,rank,periodic
  real(8),                      intent(OUT):: residual1,residual2,residual3
end subroutine advance_laminar




end module mod_tm
