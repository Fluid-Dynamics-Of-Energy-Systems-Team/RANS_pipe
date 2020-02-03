module mod_tdm
  implicit none

!****************************************************************************************

  !************************!
  !     Abstract class     !
  !************************!
  
  type, abstract, public :: TurbDiffModel
  character(len=4)                      :: name
  real(8), dimension(:,:), allocatable :: Pkt,kt,epst, yp,Tmix,Ttemp
  real(8), dimension(:),   allocatable :: alphatin, Pktin,Prtin
  contains
    procedure(init_tdm),          deferred :: init
    procedure(set_alphat_tdm),    deferred :: set_alphat
    procedure(get_sol_tdm),       deferred :: get_sol
    procedure(init_w_inflow_tdm), deferred :: init_w_inflow
    procedure(get_profile_tdm),   deferred :: get_profile
    procedure :: advance_turbdiff => advance_turbdiff_tdm
    procedure :: set_bc           => set_bc_tdm
    procedure :: set_alphat_bc
  end type TurbDiffModel

  interface
    subroutine init_tdm(this)
      import :: TurbDiffModel
      class(TurbDiffModel) :: this
    end subroutine init_tdm

    subroutine set_alphat_tdm(this,u,w,rho,temp,mu,mui,lam_cp,mut,alphat)
      use mod_param, only : k1,i1
      import :: TurbDiffModel
      class(TurbDiffModel) :: this
      real(8),dimension(0:i1,0:k1),intent(IN) :: u,w,rho,temp,mu,mui,lam_cp, mut
      real(8),dimension(0:i1,0:k1),intent(OUT):: alphat
    end subroutine set_alphat_tdm

    subroutine get_profile_tdm(this,p_prt,p_kt,p_epst,p_Pkt,k)
      use mod_param, only : k1,i1
      import :: TurbDiffModel
      class(TurbDiffModel) :: this
      integer,                               intent(IN) :: k
      real(8),dimension(0:i1),          intent(OUT):: p_prt,p_kt,p_epst,p_Pkt
    end subroutine get_profile_tdm

    subroutine get_sol_tdm(this,Prt,epst,kt, Pkt, resKt, resEt)
      use mod_param, only : k1,i1
      import :: TurbDiffModel
      class(TurbDiffModel) :: this
      real(8),dimension(0:i1,0:k1), intent(OUT):: Prt,epst,kt,Pkt,resKt,resEt
    end subroutine get_sol_tdm

    subroutine init_w_inflow_tdm(this,Prtin,alphatin,ktin,epstin,pktin)
      use mod_param, only : i1
      import :: TurbDiffModel
      class(TurbDiffModel) :: this
      real(8), dimension(0:i1), intent(IN) :: Prtin,alphatin,ktin,epstin,pktin
    end subroutine

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
    procedure :: init_w_inflow => init_w_inflow_constprt
  end type CPrt_TurbDiffModel

contains

!****************************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************************
  !*************************!
  !    Abstract routines    !
  !*************************!

  subroutine advance_turbdiff_tdm(this,u,w,c,temp,rho,mu,ekh,ekhi,ekhk,alphat, &
                                  alpha1,alpha2,rank,periodic,residual1,residual2)
  use mod_param, only : k1,i1
  class(TurbDiffModel) :: this
  real(8), dimension(0:i1,0:k1),intent(IN) :: u,w,c,temp,rho,mu,ekh,ekhi,ekhk,alphat
  real(8),                                intent(IN) :: alpha1,alpha2
  integer,                                intent(IN) :: rank,periodic
  real(8),                                intent(OUT):: residual1,residual2
  end subroutine advance_turbdiff_tdm

  subroutine set_alphat_bc(this,alphat,periodic,px,rank)
    use mod_param, only : k1,i1,kmax,imax
    use mod_mesh, only : top_bcnovalue, bot_bcnovalue
    implicit none
    class(TurbDiffModel) :: this
    integer,                       intent(IN) :: periodic,px,rank
    real(8), dimension(0:i1,0:k1), intent(OUT):: alphat
    real(8), dimension(0:i1) :: tmp
   
    alphat(i1,:) = top_bcnovalue(:)*alphat(imax,:)
    alphat(0,:)  = bot_bcnovalue(:)*alphat(1,:)

    call shiftf(alphat,tmp,rank); alphat(:,0) =tmp(:);
    call shiftb(alphat,tmp,rank); alphat(:,k1)=tmp(:);

    if ((periodic.ne.1).and.(rank.eq.0)) then
      alphat(:,0) = this%alphatin(:)
    endif
    if ((periodic.ne.1).and.(rank.eq.px-1)) then
      alphat(:,k1) = 2.*alphat(:,kmax)-alphat(:,kmax-1)
    endif
  end subroutine

  subroutine set_bc_tdm(this,ekh,rho,periodic,rank,px)
    use mod_param, only : k1,i1
    implicit none
    class(TurbDiffModel) :: this
    real(8),dimension(0:i1,0:k1),intent(IN) :: rho,ekh
    integer,                               intent(IN) :: periodic, rank, px
  end subroutine set_bc_tdm

  !*************************!
  !Constant Prandtl routines!
  !*************************!

  type(CPrt_TurbDiffModel) function init_CPrt_TurbDiffModel(name, Prt)
    character(len=3), intent(IN) :: name
    real(8), intent(in) :: Prt
    init_CPrt_TurbDiffModel%name=name
    init_CPrt_TurbDiffModel%Prt = Prt
  end function init_CPrt_TurbDiffModel

  subroutine init_constprt(this)
    use mod_param, only : i1
    class(CPrt_TurbDiffModel) :: this
    allocate(this%alphatin(0:i1))
  end subroutine init_constprt

  subroutine set_alphat_constprt(this,u,w,rho,temp,mu,mui,lam_cp,mut,alphat)
    use mod_param, only : k1,i1
    class(CPrt_TurbDiffModel) :: this
    real(8),dimension(0:i1,0:k1),intent(IN) :: u,w,rho,temp,mu,mui,lam_cp, mut
    real(8),dimension(0:i1,0:k1),intent(OUT):: alphat
    alphat = mut/(this%Prt*rho)
  end subroutine set_alphat_constprt

  subroutine get_sol_constprt(this,Prt,epst,kt, Pkt, resKt, resEt)
    use mod_param, only : k1,i1
    class(CPrt_TurbDiffModel) :: this
    real(8),dimension(0:i1,0:k1), intent(OUT):: Prt,epst,kt, Pkt, resKt, resEt
    Prt  =this%Prt
    epst =0.    
    kt   =0.
    Pkt = 0.
    resEt = 0.
    resKt = 0.
  end subroutine get_sol_constprt

  subroutine get_profile_constprt(this,p_prt,p_kt,p_epst,p_Pkt,k)
    use mod_param, only : i1
    class(CPrt_TurbDiffModel) :: this
    integer,                               intent(IN) :: k
    real(8),dimension(0:i1),          intent(OUT):: p_prt,p_kt,p_epst,p_Pkt
    p_prt = this%Prt
    p_kt = 0
    p_epst = 0
    p_pkt = 0
  end subroutine get_profile_constprt

  subroutine init_w_inflow_constprt(this,Prtin,alphatin,ktin,epstin,pktin)
    use mod_param, only : i1
    implicit none
    class(Cprt_TurbDiffModel) :: this
    real(8), dimension(0:i1), intent(IN) :: Prtin,alphatin,ktin,epstin,pktin
    this%Prt = Prtin(0)
  end subroutine

end module mod_tdm
