module mod_tm
  implicit none

!****************************************************************************************

  !************************!
  !     Abstract class     !
  !************************!
  
  type, abstract, public :: TurbModel
  integer i1,k1,imax,kmax
  character(len=3)                     :: name
  real(8), dimension(:,:), allocatable :: nuSA,Pk,om,k,bF1,bF2,eps,v2,yp
  real(8), dimension(:),   allocatable :: mutin, Pkin
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
    subroutine set_mut_tm(this,u,w,rho,mu,mui,walldist,Rp,dRp,dru,dz,mut)
      import :: TurbModel
      class(TurbModel) :: this
      real(8), dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,mu,mui
      real(8), dimension(1:this%imax),        intent(IN) :: walldist
      real(8), dimension(0:this%i1),          intent(IN) :: Rp,dRp,dru
      real(8),                                intent(IN) :: dz
      real(8), dimension(0:this%i1,0:this%k1),intent(OUT):: mut
    end subroutine set_mut_tm
    subroutine advance_turb_tm(this,u,w,rho,mu,mui,muk,mut,beta,temp,&
                             Ru,Rp,dru,drp,dz,walldist,              &
                             alpha1,alpha2,alpha3,                   &
                             modification,rank,centerBC,periodic,    &
                             residual1, residual2, residual3)
      import :: TurbModel
      class(TurbModel) :: this
      real(8), dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,mu,mui,muk,mut,beta,temp
      real(8), dimension(0:this%i1),          intent(IN) :: Ru,Rp,dru,drp
      real(8), dimension(1:this%imax),        intent(IN) :: walldist
      real(8),                                intent(IN) :: dz,alpha1,alpha2,alpha3
      integer,                                intent(IN) :: modification,rank,centerBC,periodic
      real(8),                                intent(OUT):: residual1,residual2,residual3
    end subroutine advance_turb_tm
    subroutine set_bc_tm(this,mu,rho,periodic,rank,px)
      import :: TurbModel
      class(TurbModel) :: this
      real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: rho,mu
      integer,                               intent(IN) :: periodic, rank, px
    end subroutine set_bc_tm
    subroutine get_profile_tm(this,p_nuSA,p_k,p_eps,p_om,p_v2,p_Pk,p_bF1,p_bF2,yp,k)
      import :: TurbModel
      class(TurbModel) :: this
      integer,                               intent(IN) :: k
      real(8),dimension(0:this%i1),          intent(OUT):: p_nuSA,p_k,p_eps,p_om,p_v2,p_Pk,p_bF1,yp
      real(8),dimension(1:this%imax),        intent(OUT):: p_bF2
    end subroutine get_profile_tm
    subroutine get_sol_tm(this,nuSA,k,eps,om,v2,yp)
      import :: TurbModel
      class(TurbModel) :: this
      real(8),dimension(0:this%i1,0:this%k1),intent(OUT):: nuSA,k,eps,om,v2,yp
    end subroutine get_sol_tm
    subroutine init_w_inflow_tm(this,Re, systemsolve)
      import :: TurbModel
      class(TurbModel) :: this
      real(8), intent(IN) :: Re
      integer, intent(IN) :: systemsolve
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
  use mod_mesh, only : mesh
  class(TurbModel) :: this
  integer,                                 intent(IN) :: periodic,px,rank
  real(8), dimension(0:this%i1,0:this%k1), intent(OUT):: mut
  real(8), dimension(0:this%i1) :: tmp
 
  mut(this%i1,:) = mesh%top_bcnovalue(:)*mut(this%imax,:)
  mut(0,:)       = mesh%bot_bcnovalue(:)*mut(1,:)

  call shiftf(mut,tmp,rank); mut(:,0)      =tmp(:);
  call shiftb(mut,tmp,rank); mut(:,this%k1)=tmp(:);

  if ((periodic.ne.1).and.(rank.eq.0)) then
    mut(:,0) = this%mutin(:)
  endif
  if ((periodic.ne.1).and.(rank.eq.px-1)) then
    mut(:,this%k1) = 2.*mut(:,this%kmax)-mut(:,this%kmax-1)
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
      class(Laminar_TurbModel) :: this
      real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: rho,mu
      integer,                               intent(IN) :: periodic, rank, px
end subroutine set_bc_laminar

subroutine init_sol_laminar(this)
    class(Laminar_TurbModel) :: this
end subroutine init_sol_laminar

subroutine init_w_inflow_laminar(this, Re, systemsolve)
    class(Laminar_TurbModel) :: this
    real(8), intent(IN) :: Re
    integer, intent(IN) :: systemsolve
    real(8), dimension(0:this%i1) :: dummy
    character(len=5)  :: Re_str
    integer           :: Re_int
    Re_int = int(Re)
    write(Re_str,'(I5.5)') Re_int
    if (systemsolve .eq. 1) open(29,file = 'pipe/Inflow_'//this%name//'_'//Re_str//'.dat',form='unformatted')
    if (systemsolve .eq. 2) open(29,file = 'channel/Inflow_'//this%name//'_'//Re_str//'.dat',form='unformatted')
    if (systemsolve .eq. 3) open(29,file = 'symchan/Inflow_'//this%name//'_'//Re_str//'.dat',form='unformatted')
    
    read(29) dummy(:),dummy(:),dummy(:),dummy(:),dummy(:), &
             dummy(:),this%mutin(:),dummy(:)
    close(29)
end subroutine init_w_inflow_laminar

subroutine init_mem_laminar(this)
    class(Laminar_TurbModel) :: this
    allocate(this%yp(0:this%i1,0:this%k1))
    allocate(this%mutin(0:this%i1))
end subroutine init_mem_laminar

subroutine get_profile_laminar(this,p_nuSA,p_k,p_eps,p_om,p_v2,p_Pk,p_bF1,p_bF2,yp,k)
    class(Laminar_TurbModel) :: this
    integer,                               intent(IN) :: k
    real(8),dimension(0:this%i1),          intent(OUT):: p_nuSA,p_k,p_eps,p_om,p_v2,p_Pk,p_bF1,yp
    real(8),dimension(1:this%imax),        intent(OUT):: p_bF2
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

subroutine get_sol_laminar(this,nuSA,k,eps,om,v2,yp)
    class(Laminar_TurbModel) :: this
    real(8),dimension(0:this%i1,0:this%k1), intent(OUT):: nuSA,k,eps,om,v2,yp
    nuSA=0
    k   =0    
    eps =0
    v2  =0
    om  =0
    yp  = this%yp
end subroutine get_sol_laminar


subroutine set_mut_laminar(this,u,w,rho,mu,mui,walldist,Rp,dRp,dru,dz,mut)
      class(Laminar_TurbModel) :: this
      real(8), dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,mu,mui
      real(8), dimension(1:this%imax),        intent(IN) :: walldist
      real(8), dimension(0:this%i1),          intent(IN) :: Rp,dRp,dru
      real(8),                                intent(IN) :: dz
      real(8), dimension(0:this%i1,0:this%k1),intent(OUT):: mut
      real(8), dimension(0:this%k1) :: tauw
      integer :: i,k,km
      do k=1,this%kmax
        km=k-1
        tauw(k) = mui(this%imax,k)*0.5*(w(this%imax,k)+w(this%imax,k))/walldist(this%imax)
        do i=1,this%imax
          this%yp(i,k) = sqrt(rho(i,k))/mu(i,k)*(walldist(i))*tauw(k)**0.5   ! ystar
        enddo
      enddo
      mut = 0
end subroutine set_mut_laminar

subroutine advance_laminar(this,u,w,rho,mu,mui,muk,mut,beta,temp,&
                             Ru,Rp,dru,drp,dz,walldist,              &
                             alpha1,alpha2,alpha3,                   &
                             modification,rank,centerBC,periodic,    &
                             residual1, residual2, residual3)
      class(Laminar_TurbModel) :: this
      real(8), dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,mu,mui,muk,mut,beta,temp
      real(8), dimension(0:this%i1),          intent(IN) :: Ru,Rp,dru,drp
      real(8), dimension(1:this%i1),          intent(IN) :: walldist
      real(8),                                intent(IN) :: dz,alpha1,alpha2,alpha3
      integer,                                intent(IN) :: modification,rank,centerBC,periodic
      real(8),                                intent(OUT):: residual1,residual2,residual3
end subroutine advance_laminar




end module mod_tm
