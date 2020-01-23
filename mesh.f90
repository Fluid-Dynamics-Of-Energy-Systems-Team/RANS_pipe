module mod_mesh

  implicit none
!****************************************************************************************
  
  !*************************!
  !     Abstract class      !
  !*************************!

  type, abstract, public :: AbstractMesh
  character(len=10)                  :: name
  real(8), allocatable               :: dpdz,start
  integer, allocatable               :: centerBC,numDomain           
  
  contains
    procedure :: init_mem               => init_mem
    procedure :: discretize_streamwise  => discretize_streamwise
    procedure :: discretize_streamwise2 => discretize_streamwise2    
    procedure :: discretize_wall_normal => discretize_wall_normal
    procedure :: calc_walldist          => calc_walldist
    procedure :: set_carthesian         => set_carthesian
    procedure(init), deferred     :: init
    procedure(set_bc), deferred   :: set_bc
  end type AbstractMesh

  interface

    subroutine init(this,LoD,K_start_heat,x_start_heat,rank,px)
      import :: AbstractMesh
      class(AbstractMesh) :: this
      real(8), intent(IN) :: LoD, x_start_heat
      integer, intent(IN) :: px, K_start_heat,rank
    end subroutine init

    subroutine set_bc(this,K_start_heat,x_start_heat, rank)
      import :: AbstractMesh
      class(AbstractMesh) :: this
      integer, intent(IN) :: K_start_heat, rank
      real(8), intent(IN) :: x_start_heat
    end subroutine set_bc
  end interface

class(AbstractMesh),   allocatable :: mesh
real(8), dimension(:), allocatable :: top_bcnovalue,bot_bcnovalue,top_bcvalue1,bot_bcvalue1, &
                                      top_bcvalue,bot_bcvalue,ubot_bcvalue, &
                                      drp,dru,rp,ru,dzp,dzw,zp,zw,y_fa,y_cv,walldistu,walldist
real(8) :: dz

!****************************************************************************************

  !*************************!
  !    Pipe Mesh Class      !
  !*************************!
  
  type, extends(AbstractMesh), public :: Pipe_Mesh
  contains
    procedure :: init => init_pipe
    procedure :: set_bc => set_bc_pipe
    ! procedure ::   calc_walldist => calc_walldist
  end type Pipe_Mesh

  !****************************************************************************************

  !*************************!
  !   Channel Mesh Class    !
  !*************************!
  
  type, extends(AbstractMesh), public :: Channel_Mesh
  contains
    procedure :: init => init_channel
    procedure :: set_bc => set_bc_channel
    procedure :: calc_walldist => calc_walldist_twowall
  end type Channel_Mesh

!****************************************************************************************

  !*************************!
  ! Sym. Channel Mesh Class !
  !*************************!
  
  type, extends(AbstractMesh), public :: SymChannel_Mesh
  contains
    procedure :: init => init_symchannel
    procedure :: set_bc => set_bc_symchannel
  end type SymChannel_Mesh

!****************************************************************************************

  !*************************!
  !Boundary Layer Mesh Class!
  !*************************!
  
  type, extends(AbstractMesh), public :: BLayer_Mesh
  contains
    procedure :: init => init_bl
    procedure :: set_bc => set_bc_blayer
  end type BLayer_Mesh

contains
!****************************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************************

  !************************!
  !   Abstract routines    !
  !************************!

  subroutine init_mem(this)
    use mod_param, only : i1,k1,imax
    implicit none
    class(AbstractMesh) :: this
    allocate(Ru  (0:i1),Rp  (0:i1), &
             dru (0:i1),drp (0:i1), &
             y_fa(0:i1),y_cv(0:i1))
    allocate(zw  (0:k1),zp  (0:k1), &
             dzw (0:k1),dzp (0:k1))
    allocate(wallDist(1:imax), wallDistu(0:i1))
    allocate(top_bcvalue  (0:k1),bot_bcvalue  (0:k1), &
             top_bcvalue1 (0:k1),bot_bcvalue1 (0:k1), &
             top_bcnovalue(0:k1),bot_bcnovalue(0:k1), &
             ubot_bcvalue (0:k1))
  end subroutine init_mem

  subroutine set_carthesian(this)
    use mod_param, only : i1,k1
    implicit none
    class(AbstractMesh) :: this
    integer :: i
    do i=0,i1
      ru(i)=1.0
      rp(i)=1.0
    enddo
  end subroutine

  subroutine discretize_wall_normal(this, fA, fB, gridSize)
    use mod_param, only : i1,k1,imax,kmax
    implicit None
    class(AbstractMesh) :: this
    real(8), intent(IN) :: fB, fA, gridsize
    real(8) :: fact
    integer :: i

    ru(0) = 0

    !apply stretching
    do i = 1,imax
      fact       = (i-0.)/(imax-0.)
      ru(i) = (1.-tanh(fB*(fA-fact))/tanh(fA*fB))
      ru(i) = ru(i)/(1.-tanh(fB*(fA-1.))/tanh(fA*fB))
    enddo

    !normalize with the gridsize
    do i=0,imax
      ru(i)=ru(i)/ru(imax)*gridSize
    enddo

    !calculate the cell centers and differences
    do i =1,imax
      Rp(i)  = (Ru(i)+Ru(i-1))/2.0
      dru(i) = (Ru(i)-Ru(i-1))
    enddo

    dru(i1) = dru(imax)
    Ru(i1)  = Ru(imax) + dru(i1)
    Rp(i1)  = Ru(imax) + dru(i1)/2.0

    dru(0) = dru(1)
    Rp(0)  = Ru(0) - dru(0)/2.0

    do i = 0,imax
      drp(i) = Rp(i+1) - Rp(i)
    enddo

    do i=0,i1
      y_cv(i)=rp(i)
      y_fa(i)=ru(i)
    enddo

  end subroutine discretize_wall_normal

  subroutine discretize_streamwise(this, LoD, px)
    use mod_param, only : i1,kmax
    implicit None
    class(AbstractMesh) :: this
    real(8), intent(IN) :: LoD
    integer, intent(IN) :: px 
    dz    = 1.0*LoD/(kmax*px)
  end subroutine discretize_streamwise

  subroutine discretize_streamwise2(this, LoD,rank, px)
    use mod_param, only : i1,k1,imax,kmax,K_start_heat
    use mod_math, only : splint, spline
    implicit None
    include 'mpif.h'
    class(AbstractMesh) :: this
    real(8), intent(IN) :: LoD
    integer, intent(IN) :: rank, px 
    real(8) :: L,a,c,H,value,tmp
    real(8), dimension(0:2000) :: y, x2tab, x, ys
    integer :: i,nelem, k, ierr
    integer :: tabkhi,tabklo = 0 
    character*5 cha
    
    a = 10.
    L = 0.1
    c = 0.8
    H = 0.1
    
    !create function to interpolate on
    nelem=2000
    do i=0,nelem
      x(i) = (i+0.)/(nelem+0)
      y(i) = H*((tanh(a*(x(i)/L-0.5))+1.)/2.)+c*x(i)
    enddo
    do i=0,nelem
      ys(i) = (y(i)-y(0))/(y(nelem)-y(0))
    enddo
    call spline(ys,x,  nelem+1,x2tab)

    !calculate the cell faces
    do i = 0,kmax!*px
      call splint(ys,x, x2tab,nelem+1,(i+rank*kmax+0.)/(kmax*px),zw(i),tabkhi,tabklo) 
      zw(i) = zw(i)*LoD
    enddo
    call shiftv_b(zw, k1, value, rank); zw(k1) = value;
    if (rank .eq. px-1) zw(k1) =  zw(kmax)+(zw(kmax)-zw(kmax-1))

    !calculate the cell centers
    do k =1,k1
      zp(k)  = (zw(k)+zw(k-1))/2.0
    enddo
    call shiftv_f(zp, k1, value, rank); zp(0) = value;
    if (rank .eq. 0) then
      zp(0) =  -zp(1)
      tmp = zw(K_start_heat)
    endif

    call MPI_Bcast( tmp, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr);
    this%start = tmp
    !calculate the differences
    do i = 0,kmax
      dzp(i) = (zp(i+1)+zp(i))/2.
    enddo
    if (rank .eq. 0)   dzp(0) = dzp(1)
    do i = 1,k1
      dzw(i) = zw(i)-zw(i-1)
    enddo
    dzw(0) = dzw(1)

    dz    = 1.0*LoD/(kmax*px)
    do k=0,k1
       dzw(k) = dz 
       dzp(k) = dz
       zw(k)  = (k+kmax*rank)*dz
       zp(k)  = (k+kmax*rank)*dz - (0.5)*dz
    enddo

  end subroutine discretize_streamwise2

  subroutine calc_walldist(this, gridSize)
    use mod_param, only : i1,imax
    implicit None
    class(AbstractMesh) :: this
    real(8), intent(IN) :: gridSize
    integer :: i
    do i = 1,imax
      wallDist(i) = gridSize - rp(i)
    enddo
    do i = 0,i1
      wallDistu(i) = gridSize - ru(i)
    enddo
  end subroutine

!****************************************************************************************

  !*************************!
  !    Pipe Mesh routines   !
  !*************************!

  type(Pipe_Mesh) function init_PipeMesh(name)
    character(len=4), intent(IN) :: name
    init_PipeMesh%name=name
  end function init_PipeMesh

  subroutine init_pipe(this, LoD, K_start_heat, x_start_heat, rank,px)
    use mod_param, only : bulkmod
    implicit none
    class(Pipe_Mesh) :: this
    real(8), intent(IN) :: LoD, x_start_heat
    integer, intent(IN) :: px, K_start_heat,rank
    real(8) :: pi, gridsize, fA, fB
    this%numDomain = 1
    this%centerBC = 1
    gridSize  = 0.5
    fA        = 0.12
    fB        = 2.4
    ! this%dpdz = 4.0
    if (bulkmod .eq. 1) then
      this%dpdz = 0.00005
    else
      this%dpdz = 4.0
    endif

    call this%init_mem()
    call this%discretize_streamwise(LoD, px)
    call this%discretize_wall_normal(fA,fB,gridSize)
    call this%calc_walldist(gridsize)
    call this%set_bc(K_start_heat, x_start_heat, rank)

  end subroutine init_pipe
  subroutine set_bc_pipe(this, K_start_heat, x_start_heat,rank)
    use mod_param, only :k1,k
    implicit none
    class(Pipe_Mesh) :: this
    real(8), intent(IN) :: x_start_heat
    integer, intent(IN) :: K_start_heat, rank
    !bc for the momentum and turbulent scalars
    bot_bcnovalue(:) = 1 ! symmetry
    top_bcnovalue(:) =-1 ! wall
    bot_bcvalue(:)   = 1 ! symmetry
    top_bcvalue(:)   = 0 ! wall
    ubot_bcvalue(:)  = 1 ! zero vertical velocity
    
    ! bc for the temperature
    bot_bcvalue1(:)  = 1 ! symmetry
    do k=0,k1
      if ((rank.eq.0) .and. (k.lt.K_start_heat)) then
        top_bcvalue1(k) =1 ! no heat flux (symmetry)
      else
        top_bcvalue1(k) =0 ! heat flux or isothermal
      endif
    enddo

  end subroutine


!****************************************************************************************

  !*************************!
  !  Channel Mesh routines  !
  !*************************!

  type(Channel_Mesh) function init_ChannelMesh(name)
    character(len=7), intent(IN) :: name
    init_ChannelMesh%name=name
  end function init_ChannelMesh

  subroutine init_channel(this, LoD, K_start_heat, x_start_heat,rank, px)
    use mod_param, only : bulkmod
    implicit none
    class(Channel_Mesh) :: this
    real(8), intent(IN) :: LoD, x_start_heat
    integer, intent(IN) :: px, K_start_heat,rank
    real(8) :: pi, gridsize, fA, fB
    this%numDomain = -1
    this%centerBC = -1
    gridSize  = 2.0
    fA        = 0.5
    fB        = 4.6

    if (bulkmod .eq. 1) then
      this%dpdz = 0.0005
    else
      this%dpdz = 1.0
    endif

    call this%init_mem()
    call this%discretize_streamwise(LoD, px)
    call this%discretize_wall_normal(fA,fB,gridSize)
    call this%set_bc(K_start_heat, x_start_heat,rank)
    call this%calc_walldist(gridsize)
    call this%set_carthesian()

  end subroutine init_channel


  subroutine set_bc_channel(this, K_start_heat, x_start_heat,rank)
    use mod_param, only :k1,k
    implicit none
    class(Channel_Mesh) :: this
    real(8), intent(IN) :: x_start_heat
    integer, intent(IN) :: K_start_heat, rank

    !bc for the momentum and turbulent scalars
    bot_bcnovalue(:) =-1 ! wall
    top_bcnovalue(:) =-1 ! wall
    bot_bcvalue(:)   = 0 ! wall
    top_bcvalue(:)   = 0 ! wall
    ubot_bcvalue(:)  = 1 ! zero vertical velocity

    ! bc for the temperature
    do k=0,k1
      if ((rank.eq.0) .and. (k.lt.K_start_heat)) then
        top_bcvalue1(k) = 1 ! no heat flux (symmetry)
        bot_bcvalue1(k) = 1 ! no heat flux (symmetry)
      else
        top_bcvalue1(k) = 0 ! heat flux or isothermal
        bot_bcvalue1(k) = 0 ! heat flux or isothermal
      endif
    enddo

  end subroutine

  subroutine calc_walldist_twowall(this, gridSize)
    use mod_param, only : imax,i1
    implicit None
    class(Channel_Mesh) :: this
    real(8), intent(IN) :: gridSize
    integer :: i

    do i = 1,imax
      if (rp(i).le.1) then
        wallDist(i) = rp(i)
      else
        wallDist(i) = gridSize-rp(i)
      endif
    enddo
    do i = 0, i1 
      if (ru(i).le.1) then
        wallDistu(i) = ru(i)
      else
        wallDistu(i) = gridSize-ru(i)
      endif
    enddo
  end subroutine

!****************************************************************************************

  !**************************!
  !Sym. Channel Mesh routines!
  !**************************!
  type(SymChannel_Mesh) function init_SymChannelMesh(name)
    character(len=7), intent(IN) :: name
    init_SymChannelMesh%name=name
  end function init_SymChannelMesh

  subroutine init_symchannel(this, LoD, K_start_heat, x_start_heat,rank, px)
    implicit none
    class(SymChannel_Mesh) :: this
    real(8), intent(IN) :: LoD, x_start_heat
    integer, intent(IN) :: px, K_start_heat,rank
    real(8) :: pi, gridsize, fA, fB
    this%numDomain = -1
    this%centerBC = 1
    gridSize  = 1.0
    fA        = 0.12
    fB        = 2.4
    this%dpdz = 1.0

    call this%init_mem()
    call this%discretize_streamwise(LoD, px)
    call this%discretize_wall_normal(fA,fB,gridSize)
    call this%set_bc(K_start_heat, x_start_heat,rank)
    call this%calc_walldist(gridsize)
    call this%set_carthesian()

  end subroutine init_symchannel

  subroutine set_bc_symchannel(this, K_start_heat, x_start_heat,rank)
    use mod_param, only :k1,k
    class(SymChannel_Mesh) :: this
    real(8), intent(IN) :: x_start_heat
    integer, intent(IN) :: K_start_heat, rank
    
    !bc for the momentum and turbulent scalars
    bot_bcnovalue(:) = 1 ! symmetry
    bot_bcvalue(:)   = 1 ! symmetry
    ubot_bcvalue(:)  = 1 ! zero vertical velocity
    top_bcnovalue(:) =-1 ! wall
    top_bcvalue(:)   = 0 ! wall
    
    ! bc for the temperature
    bot_bcvalue1(:)  = 1 ! symmetry
    do k=0,k1
      if ((rank.eq.0) .and. (k.lt.K_start_heat)) then
        top_bcvalue1(k) = 1 ! no heat flux (symmetry)
      else
        top_bcvalue1(k) = 0 ! heat flux or isothermal
      endif
    enddo
    
  end subroutine

!****************************************************************************************

  !****************************!
  !Boundary Layer Mesh routines!
  !****************************!

  type(BLayer_Mesh) function init_BLayerMesh(name)
    character(len=2), intent(IN) :: name
    init_BLayerMesh%name=name
  end function init_BLayerMesh

  subroutine init_bl(this, LoD, K_start_heat, x_start_heat,rank, px)
    implicit none
    class(BLayer_Mesh) :: this
    real(8), intent(IN) :: LoD, x_start_heat
    integer, intent(IN) :: px, K_start_heat,rank
    real(8) :: pi, gridsize, fA, fB
    this%numDomain = -1
    this%centerBC = 1
    gridSize  = 1
    fA        = 0.12
    fB        = 2.4
    this%dpdz = 1.0

    call this%init_mem()
    call this%discretize_streamwise(LoD, px)
    call this%discretize_wall_normal(fA,fB,gridSize)
    call this%set_bc(K_start_heat, x_start_heat,rank)
    call this%calc_walldist(gridsize)
    call this%set_carthesian()

  end subroutine init_bl

  subroutine set_bc_blayer(this, K_start_heat, x_start_heat,rank)
    use mod_param, only : k1,i,k
    class(BLayer_Mesh) :: this
    real(8), intent(IN) :: x_start_heat
    integer, intent(IN) :: K_start_heat, rank
    
    !bc for the momentum and turbulent scalars
    bot_bcvalue(:)   = 1 ! symmetry
    bot_bcnovalue(:) = 1 ! symmetry
    ubot_bcvalue(:)  = 0 ! 0: set the wall to du/dy =0        
    do k=0,k1
      if ((rank.eq.0) .and. (k.lt.K_start_heat)) then
        top_bcnovalue(k) = 1 !symmetry
        top_bcvalue(k)   = 1 !symmetry
      else
        top_bcnovalue(k) =-1 !wall
        top_bcvalue(k)   = 0 !wall
      endif
    enddo
    
    ! bc for the temperature
    bot_bcvalue1(:)      = 1 ! symmetry
    do k=0,k1
      if ((rank.eq.0) .and. (k.lt.K_start_heat)) then
        top_bcvalue1(k)  = 1 ! no heat flux (symmetry)
      else
        top_bcvalue1(k)  = 0 ! heat flux or isothermal
      endif
    enddo
  end subroutine
end module mod_mesh

