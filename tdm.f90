module mod_tdm
  implicit none

!****************************************************************************************

  !************************!
  !     Abstract class     !
  !************************!
  
  type, abstract, public :: TurbDiffModel
  integer i1,k1,imax,kmax
  character(len=4)                      :: name
  real(8), dimension(:,:), allocatable :: Pkt,kt,epst, yp
  real(8), dimension(:),   allocatable :: alphatin, Pktin
  contains
    procedure(init_tdm),       deferred :: init
    procedure(set_alphat_tdm), deferred :: set_alphat
    procedure(get_sol_tdm),    deferred :: get_sol
    ! procedure(init_w_inflow_tdm), deferred :: init_w_inflow
    procedure(get_profile_tdm), deferred :: get_profile
    procedure :: advance_turbdiff => advance_turbdiff_tdm
    procedure :: set_alphat_bc
    procedure :: set_bc => set_bc_tdm

  end type TurbDiffModel

  interface
    subroutine init_tdm(this)
      import :: TurbDiffModel
      class(TurbDiffModel) :: this
    end subroutine init_tdm

    subroutine set_alphat_tdm(this,u,w,rho,temp,mu,mui,lam_cp,mut,alphat)
      import :: TurbDiffModel
      class(TurbDiffModel) :: this
      real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,temp,mu,mui,lam_cp, mut
      real(8),dimension(0:this%i1,0:this%k1),intent(OUT):: alphat
    end subroutine set_alphat_tdm

    subroutine get_profile_tdm(this,p_prt,p_kt,p_epst,p_Pkt,k)
      import :: TurbDiffModel
      class(TurbDiffModel) :: this
      integer,                               intent(IN) :: k
      real(8),dimension(0:this%i1),          intent(OUT):: p_prt,p_kt,p_epst,p_Pkt
    end subroutine get_profile_tdm

    subroutine get_sol_tdm(this,Prt,epst,kt)
      import :: TurbDiffModel
      class(TurbDiffModel) :: this
      real(8),dimension(0:this%i1,0:this%k1), intent(OUT):: Prt,epst,kt
    end subroutine get_sol_tdm

  end interface

  

class(TurbDiffModel), allocatable :: turbdiff_model
!****************************************************************************************

  !************************!
  ! Constant Prandtl class !
  !************************!

  type, extends(TurbDiffModel), public :: CPrt_TurbDiffModel
  real(8) :: Prt
  contains
    procedure :: init => init_constprt
    procedure :: set_alphat  => set_alphat_constprt
    procedure :: get_sol => get_sol_constprt
    procedure :: get_profile => get_profile_constprt
  end type CPrt_TurbDiffModel

contains

!****************************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************************
  !*************************!
  !    Abstract routines    !
  !*************************!

  subroutine advance_turbdiff_tdm(this,u,w,c,temp,rho,mu,ekh,ekhi,ekhk,alphat, &
                      alpha1,alpha2,                   &
                      modification,rank,periodic,    &
                      residual1, residual2)
  class(TurbDiffModel) :: this
  real(8), dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,c,temp,rho,mu,ekh,ekhi,ekhk,alphat
  real(8),                                intent(IN) :: alpha1,alpha2
  integer,                                intent(IN) :: modification,rank,periodic
  real(8),                                intent(OUT):: residual1,residual2
  end subroutine advance_turbdiff_tdm

  subroutine set_alphat_bc(this,alphat,periodic,px,rank)
    use mod_mesh, only : mesh
    class(TurbDiffModel) :: this
    integer,                                 intent(IN) :: periodic,px,rank
    real(8), dimension(0:this%i1,0:this%k1), intent(OUT):: alphat
    real(8), dimension(0:this%i1) :: tmp
   
    alphat(this%i1,:) = mesh%top_bcnovalue(:)*alphat(this%imax,:)
    alphat(0,:)       = mesh%bot_bcnovalue(:)*alphat(1,:)

    call shiftf(alphat,tmp,rank); alphat(:,0)      =tmp(:);
    call shiftb(alphat,tmp,rank); alphat(:,this%k1)=tmp(:);

    if ((periodic.ne.1).and.(rank.eq.0)) then
      alphat(:,0) = this%alphatin(:)
    endif
    if ((periodic.ne.1).and.(rank.eq.px-1)) then
      alphat(:,this%k1) = 2.*alphat(:,this%kmax)-alphat(:,this%kmax-1)
    endif
  end subroutine

  subroutine set_bc_tdm(this,ekh,rho,periodic,rank,px)
    implicit none
    class(TurbDiffModel) :: this
    real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: rho,ekh
    integer,                               intent(IN) :: periodic, rank, px
  end subroutine set_bc_tdm

  !*************************!
  !Constant Prandtl routines!
  !*************************!

  type(CPrt_TurbDiffModel) function init_CPrt_TurbDiffModel(i1,k1,imax,kmax,name, Prt)
    integer, intent(in) :: i1,k1,imax,kmax
    character(len=3), intent(IN) :: name
    real(8), intent(in) :: Prt
    init_CPrt_TurbDiffModel%name=name
    init_CPrt_TurbDiffModel%i1 = i1
    init_CPrt_TurbDiffModel%k1 = k1
    init_CPrt_TurbDiffModel%imax = imax
    init_CPrt_TurbDiffModel%kmax = kmax
    init_CPrt_TurbDiffModel%Prt = Prt
  end function init_CPrt_TurbDiffModel

  subroutine init_constprt(this)
    class(CPrt_TurbDiffModel) :: this
  end subroutine init_constprt

  subroutine set_alphat_constprt(this,u,w,rho,temp,mu,mui,lam_cp,mut,alphat)
    class(CPrt_TurbDiffModel) :: this
    real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,temp,mu,mui,lam_cp, mut
    real(8),dimension(0:this%i1,0:this%k1),intent(OUT):: alphat
    alphat = mut/this%Prt
  end subroutine set_alphat_constprt

  subroutine get_sol_constprt(this,Prt,epst,kt)
    class(CPrt_TurbDiffModel) :: this
    real(8),dimension(0:this%i1,0:this%k1), intent(OUT):: Prt,epst,kt
    Prt  =this%Prt
    epst =0.    
    kt   =0.
  end subroutine get_sol_constprt

  subroutine get_profile_constprt(this,p_prt,p_kt,p_epst,p_Pkt,k)
    class(CPrt_TurbDiffModel) :: this
    integer,                               intent(IN) :: k
    real(8),dimension(0:this%i1),          intent(OUT):: p_prt,p_kt,p_epst,p_Pkt
    p_prt = this%Prt
    p_kt = 0
    p_epst = 0
    p_pkt = 0
  end subroutine get_profile_constprt



end module mod_tdm
