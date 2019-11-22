module module_mesh

  implicit none
!****************************************************************************************
  
  !*************************!
  !     Abstract class      !
  !*************************!

  type, abstract, public :: AbstractMesh
  integer                            :: i1, k1, imax, kmax
  real(8), dimension(:), allocatable :: Ru,Rp,dru,drp,y_fa,y_cv     !0:i1
  real(8), allocatable               :: dz,dpdz
  real(8), dimension(:), allocatable :: wallDist                    !1:imax
  integer, allocatable               :: centerBC,numDomain           
  real(8), dimension(:), allocatable :: top_bcvalue,  bot_bcvalue,   &
                                        top_bcnovalue,bot_bcnovalue, &
                                        top_bcvalue1, bot_bcvalue1,  &
                                        ubot_bcvalue
  contains
    procedure :: init_mem               => init_mem
    procedure :: discretize_streamwise  => discretize_streamwise
    procedure :: discretize_wall_normal => discretize_wall_normal
    procedure :: calc_walldist          => calc_walldist
    procedure :: set_carthesian         => set_carthesian
    procedure(init), deferred     :: init
    procedure(set_bc), deferred   :: set_bc
  end type AbstractMesh

  interface

    subroutine init(this, LoD, K_start_heat, x_start_heat, rank,px)
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

class(AbstractMesh),  allocatable :: mesh
!****************************************************************************************

  !*************************!
  !    Pipe Mesh Class      !
  !*************************!
  
  type, extends(AbstractMesh), public :: Pipe_Mesh
  contains
    procedure :: init => init_pipe
    procedure :: set_bc => set_bc_pipe
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
    implicit none
    class(AbstractMesh) :: this
    allocate(this%Ru  (0:this%i1),this%Rp  (0:this%i1), &
             this%dru (0:this%i1),this%drp (0:this%i1), &
             this%y_fa(0:this%i1),this%y_cv(0:this%i1))
    allocate(this%wallDist(1:this%imax))
    allocate(this%top_bcvalue  (0:this%k1), this%bot_bcvalue  (0:this%k1), &
             this%top_bcvalue1 (0:this%k1), this%bot_bcvalue1 (0:this%k1), &
             this%top_bcnovalue(0:this%k1), this%bot_bcnovalue(0:this%k1), &
             this%ubot_bcvalue (0:this%k1))
  end subroutine init_mem

  subroutine set_carthesian(this)
    implicit none
    class(AbstractMesh) :: this
    integer :: i
    do i=0,this%i1
      this%ru(i)=1.0
      this%rp(i)=1.0
    enddo
  end subroutine

  subroutine discretize_wall_normal(this, fB, fA, gridSize)
    implicit None
    class(AbstractMesh) :: this
    real(8), intent(IN) :: fB, fA, gridsize
    real(8) :: fact
    integer :: i

    this%ru(0) = 0

    !apply stretching
    do i = 1,this%imax
      fact       = (i-0.)/(this%imax-0.)
      this%ru(i) = (1.-tanh(fB*(fA-fact))/tanh(fA*fB))
      this%ru(i) = this%ru(i)/(1.-tanh(fB*(fA-1.))/tanh(fA*fB))
    enddo

    !normalize with the gridsize
    do i=0,this%imax
      this%ru(i)=this%ru(i)/this%ru(this%imax)*gridSize
    enddo

    !calculate the cell centers and differences
    do i =1,this%imax
      this%Rp(i)  = (this%Ru(i)+this%Ru(i-1))/2.0
      this%dru(i) = (this%Ru(i)-this%Ru(i-1))
    enddo

    this%dru(this%i1) = this%dru(this%imax)
    this%Ru(this%i1)  = this%Ru(this%imax) + this%dru(this%i1)
    this%Rp(this%i1)  = this%Ru(this%imax) + this%dru(this%i1)/2.0

    this%dru(0) = this%dru(1)
    this%Rp(0)  = this%Ru(0) - this%dru(0)/2.0

    do i = 0,this%imax
      this%drp(i) = this%Rp(i+1) - this%Rp(i)
    enddo

    do i=0,this%i1
      this%y_cv(i)=this%rp(i)
      this%y_fa(i)=this%ru(i)
    enddo

  end subroutine discretize_wall_normal

  subroutine discretize_streamwise(this, LoD, px)
    implicit None
    class(AbstractMesh) :: this
    real(8), intent(IN) :: LoD
    integer, intent(IN) :: px 
    this%dz    = 1.0*LoD/(this%kmax*px)
  end subroutine discretize_streamwise

  subroutine calc_walldist(this, gridSize)
    implicit None
    class(AbstractMesh) :: this
    real(8), intent(IN) :: gridSize
    integer :: i
    do i = 1,this%imax
      this%wallDist(i) = gridSize - this%rp(i)
    enddo
  end subroutine

!****************************************************************************************

  !*************************!
  !    Pipe Mesh routines   !
  !*************************!


  subroutine init_pipe(this, LoD, K_start_heat, x_start_heat, rank,px)
    implicit none
    class(Pipe_Mesh) :: this
    real(8), intent(IN) :: LoD, x_start_heat
    integer, intent(IN) :: px, K_start_heat,rank
    real(8) :: pi, gridsize, fA, fB
    gridSize  = 0.5
    fA        = 0.12
    fB        = 2.4
    this%dpdz = 4.0

    call this%init_mem()
    call this%discretize_streamwise(LoD, px)
    call this%discretize_wall_normal(fA,fB,gridSize)
    call this%calc_walldist(gridsize)
    call this%set_bc(K_start_heat, x_start_heat, rank)

  end subroutine init_pipe
  subroutine set_bc_pipe(this, K_start_heat, x_start_heat,rank)
    class(Pipe_Mesh) :: this
    real(8), intent(IN) :: x_start_heat
    integer, intent(IN) :: K_start_heat, rank
    integer :: i,k
    !bc for the momentum and turbulent scalars
    this%bot_bcnovalue(:) = 1 ! symmetry
    this%bot_bcvalue(:)   = 1 ! symmetry
    this%ubot_bcvalue(:)  = 1 ! zero vertical velocity
    this%top_bcnovalue(:) =-1 ! wall
    this%top_bcvalue(:)   = 0 ! wall
    
    ! bc for the temperature
    this%bot_bcvalue1(:)  = 1 ! symmetry
    do k=0,this%k1
      if ( (k+rank*this%kmax)*this%dz.lt. x_start_heat) then
        this%top_bcvalue1(k) =1 ! no heat flux (symmetry)
      else
        this%top_bcvalue1(k) =0 ! heat flux or isothermal
      endif
    enddo

  end subroutine


!****************************************************************************************

  !*************************!
  !  Channel Mesh routines  !
  !*************************!

  subroutine init_channel(this, LoD, K_start_heat, x_start_heat,rank, px)
    implicit none
    class(Channel_Mesh) :: this
    real(8), intent(IN) :: LoD, x_start_heat
    integer, intent(IN) :: px, K_start_heat,rank
    real(8) :: pi, gridsize, fA, fB
    gridSize  = 2.0
    fA        = 0.5
    fB        = 4.6
    this%dpdz = 1.0

    call this%init_mem()
    call this%discretize_streamwise(LoD, px)
    call this%discretize_wall_normal(fA,fB,gridSize)
    call this%set_bc(K_start_heat, x_start_heat,rank)
    call this%calc_walldist(gridsize)
    call this%set_carthesian()

  end subroutine init_channel


  subroutine set_bc_channel(this, K_start_heat, x_start_heat,rank)
    class(Channel_Mesh) :: this
    real(8), intent(IN) :: x_start_heat
    integer, intent(IN) :: K_start_heat, rank
    integer :: i,k

    !bc for the momentum and turbulent scalars
    this%bot_bcnovalue(:) =-1 ! wall
    this%bot_bcvalue(:)   = 0 ! wall
    this%ubot_bcvalue(:)  = 1 ! zero vertical velocity
    this%top_bcnovalue(:) =-1 ! wall
    this%top_bcvalue(:)   = 0 ! wall

    ! bc for the temperature
    do k=0,this%k1
      if ((k+rank*this%kmax)*this%dz.lt.x_start_heat) then
        this%top_bcvalue1(k) = 1 ! no heat flux (symmetry)
        this%bot_bcvalue1(k) = 1 ! no heat flux (symmetry)
      else
        this%top_bcvalue1(k) = 0 ! heat flux or isothermal
        this%bot_bcvalue1(k) = 0 ! heat flux or isothermal
      endif
    enddo

  end subroutine

  subroutine calc_walldist_twowall(this, gridSize)
    implicit None
    class(Channel_Mesh) :: this
    real(8), intent(IN) :: gridSize
    integer :: i
    do i = 1,this%imax
      if (this%rp(i).le.1) then
        this%wallDist(i) = this%rp(i)
      else
        this%wallDist(i) = gridSize-this%rp(i)
      endif
    enddo
  end subroutine

!****************************************************************************************

  !**************************!
  !Sym. Channel Mesh routines!
  !**************************!

  subroutine init_symchannel(this, LoD, K_start_heat, x_start_heat,rank, px)
    implicit none
    class(SymChannel_Mesh) :: this
    real(8), intent(IN) :: LoD, x_start_heat
    integer, intent(IN) :: px, K_start_heat,rank
    real(8) :: pi, gridsize, fA, fB
    
    gridSize  = 1.0
    fA        = 0.12
    fB        = 2.4
    this%dpdz      = 1.0

    call this%init_mem()
    call this%discretize_streamwise(LoD, px)
    call this%discretize_wall_normal(fA,fB,gridSize)
    call this%set_bc(K_start_heat, x_start_heat,rank)
    call this%calc_walldist(gridsize)
    call this%set_carthesian()

  end subroutine init_symchannel

  subroutine set_bc_symchannel(this, K_start_heat, x_start_heat,rank)
    class(SymChannel_Mesh) :: this
    real(8), intent(IN) :: x_start_heat
    integer, intent(IN) :: K_start_heat, rank
    integer :: i,k
    
    !bc for the momentum and turbulent scalars
    this%bot_bcnovalue(:) = 1 ! symmetry
    this%bot_bcvalue(:)   = 1 ! symmetry
    this%ubot_bcvalue(:)  = 1 ! zero vertical velocity
    this%top_bcnovalue(:) =-1 ! wall
    this%top_bcvalue(:)   = 0 ! wall
    
    ! bc for the temperature
    this%bot_bcvalue1(:)  = 1 ! symmetry
    do k=0,this%k1
      if ((k+rank*this%kmax)*this%dz.lt.x_start_heat) then
        this%top_bcvalue1(k) = 1 ! no heat flux (symmetry)
      else
        this%top_bcvalue1(k) = 0 ! heat flux or isothermal
      endif
    enddo
    
  end subroutine

!****************************************************************************************

  !****************************!
  !Boundary Layer Mesh routines!
  !****************************!

  subroutine init_bl(this, LoD, K_start_heat, x_start_heat,rank, px)
    implicit none
    class(BLayer_Mesh) :: this
    real(8), intent(IN) :: LoD, x_start_heat
    integer, intent(IN) :: px, K_start_heat,rank
    real(8) :: pi, gridsize, fA, fB
    
    gridSize  = 1
    fA        = 0.12
    fB        = 2.4
    this%dpdz      = 1.0

    call this%init_mem()
    call this%discretize_streamwise(LoD, px)
    call this%discretize_wall_normal(fA,fB,gridSize)
    call this%set_bc(K_start_heat, x_start_heat,rank)
    call this%calc_walldist(gridsize)
    call this%set_carthesian()

  end subroutine init_bl

  subroutine set_bc_blayer(this, K_start_heat, x_start_heat,rank)
    class(BLayer_Mesh) :: this
    real(8), intent(IN) :: x_start_heat
    integer, intent(IN) :: K_start_heat, rank
    integer :: i,k

    !bc for the momentum and turbulent scalars
    this%bot_bcvalue(:)   = 1 ! symmetry
    this%bot_bcnovalue(:) = 1 ! symmetry
    this%ubot_bcvalue(:)  = 0 ! 0: set the wall to du/dy =0        
    do k=0,this%k1
      if ((rank.eq.0) .and. (k.le.K_start_heat)) then
        this%top_bcnovalue(k) = 1 !symmetry
        this%top_bcvalue(k)   = 1 !symmetry
      else
        this%top_bcnovalue(k) =-1 !wall
        this%top_bcvalue(k)   = 0 !wall
      endif
    enddo
    
    ! bc for the temperature
    this%bot_bcvalue1(:)      = 1 ! symmetry
    do k=0,this%k1
      if ((k+rank*this%kmax)*this%dz.lt.x_start_heat) then
        this%top_bcvalue1(k)  = 1 ! no heat flux (symmetry)
      else
        this%top_bcvalue1(k)  = 0 ! heat flux or isothermal
      endif
    enddo
  end subroutine
end module module_mesh

