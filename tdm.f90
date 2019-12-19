module mod_tdm
  implicit none

!****************************************************************************************

  !************************!
  !     Abstract class     !
  !************************!
  
  type, abstract, public :: TurbDiffModel
  integer i1,k1,imax,kmax
  character(len=4)                      :: name
  real(8), dimension(:,:), allocatable :: Pkt,kt,epst, yp,Tmix,Ttemp
  real(8), dimension(:),   allocatable :: alphatin, Pktin,Prtin
  contains
    procedure(init_tdm),       deferred :: init
    procedure(set_alphat_tdm), deferred :: set_alphat
    procedure(get_sol_tdm),    deferred :: get_sol
    procedure(init_w_inflow_tdm), deferred :: init_w_inflow
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

    subroutine get_sol_tdm(this,Prt,epst,kt, Pkt, resKt, resEt)
      import :: TurbDiffModel
      class(TurbDiffModel) :: this
      real(8),dimension(0:this%i1,0:this%k1), intent(OUT):: Prt,epst,kt,Pkt,resKt,resEt
    end subroutine get_sol_tdm

    subroutine init_w_inflow_tdm(this,Re,systemsolve)
      import :: TurbDiffModel
      class(TurbDiffModel) :: this
      real(8), intent(IN) :: Re
      integer, intent(IN) :: systemsolve
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
    allocate(this%alphatin(0:this%i1))
  end subroutine init_constprt

  subroutine set_alphat_constprt(this,u,w,rho,temp,mu,mui,lam_cp,mut,alphat)
    class(CPrt_TurbDiffModel) :: this
    real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,temp,mu,mui,lam_cp, mut
    real(8),dimension(0:this%i1,0:this%k1),intent(OUT):: alphat
    alphat = mut/this%Prt
  end subroutine set_alphat_constprt

  subroutine get_sol_constprt(this,Prt,epst,kt, Pkt, resKt, resEt)
    class(CPrt_TurbDiffModel) :: this
    real(8),dimension(0:this%i1,0:this%k1), intent(OUT):: Prt,epst,kt, Pkt, resKt, resEt
    Prt  =this%Prt
    epst =0.    
    kt   =0.
    Pkt = 0.
    resEt = 0.
    resKt = 0.
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

  subroutine init_w_inflow_constprt(this,Re,systemsolve)
    use mod_tm, only : turb_model
    implicit none
    class(Cprt_TurbDiffModel) :: this
    real(8), intent(IN) :: Re
    integer, intent(IN) :: systemsolve
    real(8), dimension(0:this%i1) :: dummy, Prtin
    character(len=5)  :: Re_str
    character(len=100) :: fname
    integer           :: Re_int
    Re_int = int(Re)
    write(Re_str,'(I5.5)') Re_int
    fname = 'Inflow_'//trim(turb_model%name)//'_'//trim(this%name)//'_'//Re_str//'.dat'
    if (systemsolve .eq. 1) open(29,file = 'pipe/'//trim(fname),form='unformatted')
    if (systemsolve .eq. 2) open(29,file = 'channel/'//trim(fname),form='unformatted')
    if (systemsolve .eq. 3) open(29,file = 'symchan/'//trim(fname),form='unformatted')
    read(29) dummy(:),dummy(:),dummy(:),dummy(:),dummy(:), &
             dummy(:),dummy(:),dummy(:),this%alphatin(:), Prtin(:), &
             dummy(:),dummy(:),dummy(:)
    close(29)
    this%Prt = Prtin(0)
  end subroutine init_w_inflow_constprt



end module mod_tdm
