module mod_eosmodels
  implicit none
!****************************************************************************************
  
  !************************!
  !     Abstract class     !
  !************************!
  
  type, abstract, public :: EOSModel
  real(8) Re, Pr
  contains
    procedure(init), deferred :: init
    procedure(set_w_enth), deferred :: set_w_enth
    procedure(set_w_temp), deferred :: set_w_temp
  end type EOSModel

  interface
    subroutine init(this)
      import :: EOSModel
      class(EOSModel) :: this
    end subroutine init
    subroutine set_w_enth( this, enth, prop, output )
      import :: EOSModel
      class(EOSModel) :: this
      real(8), intent(IN) :: enth
      character(len=1), intent(IN) :: prop
      real(8), intent(OUT) :: output
    end subroutine set_w_enth
    subroutine set_w_temp( this, temp, prop, output )
      import :: EOSModel
      class(EOSModel) :: this
      real(8), intent(IN) :: temp
      character(len=1), intent(IN) :: prop
      real(8), intent(OUT) :: output
    end subroutine set_w_temp
  end interface

!****************************************************************************************

  !************************!
  !     Ideal Gas class    !
  !************************!
  
  type, extends(EOSModel), public :: IG_EOSModel
  contains
    procedure :: init => initialize_ig
    procedure :: set_w_enth => set_w_enth_ig
    procedure :: set_w_temp => set_w_temp_ig

  end type IG_EOSModel

!****************************************************************************************

  !************************!
  !       Table class      !
  !************************!
  
  type, extends(EOSModel), public :: Table_EOSModel
  integer ntab
  character(len=20) filename
  real(8), dimension(:), private, allocatable :: tempTab,rhoTab,betaTab, muTab,lamTab, &
                                                 cpTab,enthTab,lamocpTab, temp2Tab,    &
                                                 rho2Tab,beta2Tab, mu2Tab,lam2Tab,     &
                                                 cp2Tab,enth2Tab,lamocp2Tab
  contains
    procedure :: init => initialize_table
    procedure :: set_w_enth => set_w_enth_table
    procedure :: set_w_temp => set_w_temp_table
    procedure, private :: allocate_mem => allocate_mem
    procedure, private :: read_table => read_table
    procedure, private :: calc_interp_coeff => calc_interp_coeff
  end type Table_EOSModel


contains
!****************************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************************

  !************************!
  !   Ideal Gas routines   !
  !************************!

  subroutine initialize_ig(this)
    implicit none
    class(IG_EOSModel) :: this
  end subroutine initialize_ig

  subroutine set_w_enth_ig(this, enth, prop, output)
    implicit none
    class(IG_EOSModel) :: this
    real(8), intent(IN) :: enth
    character(len=1), intent(IN) :: prop
    real(8), intent(OUT) :: output
    select case (prop)
      case ("T")
        output = enth + 1.0
      case ("D")
        output = 1./(enth+1.)
      case ("L")
        output = 1./(this%Re*this%Pr)
      case ("C")
        output = 1.
      case ("V")
        output = 1./this%Re
      case ("B")
        output = 1./(enth+1)
      case default
        write(*,*) "Property doesn't exist!!!"
    end select
  end subroutine set_w_enth_ig
  
  subroutine set_w_temp_ig(this, temp, prop, output)
    implicit none
    class(IG_EOSModel) :: this
    real(8), intent(IN) :: temp
    character(len=1), intent(IN) :: prop
    real(8), intent(OUT) :: output
    select case (prop)
      case ("H")
        output = temp - 1.0
      case default
        write(*,*) "Property doesn't exist!!!"
    end select
  end subroutine set_w_temp_ig
  
!****************************************************************************************

  !************************!
  !     Table routines     !
  !************************!

  subroutine initialize_table(this)
    implicit none
    class(Table_EOSModel) :: this
    call this%allocate_mem()
    call this%read_table()
    call this%calc_interp_coeff()
  end subroutine initialize_table

  subroutine set_w_enth_table(this, enth, prop, output)
    use mod_math
    implicit none
    class(Table_EOSModel) :: this
    real(8), intent(IN) :: enth
    character(len=1), intent(IN) :: prop
    real(8), intent(OUT):: output
    integer :: tabkhi,tabklo = 0 
    select case (prop)
      case ("T")
        call splint(this%enthTab,this%tempTab,  this%temp2Tab,  this%nTab,enth,output,tabkhi,tabklo)
      case ("D")
        call splint(this%enthTab,this%rhoTab,   this%rho2Tab,   this%nTab,enth,output,tabkhi,tabklo)
      case ("L")
        call splint(this%enthTab,this%lamocpTab,this%lamocp2Tab,this%nTab,enth,output,tabkhi,tabklo)
        output = output/(this%Re*this%Pr)
      case ("C")
        call splint(this%enthTab,this%cpTab,    this%cp2Tab,    this%nTab,enth,output,tabkhi,tabklo)
      case ("V")
        call splint(this%enthTab,this%muTab,    this%mu2Tab,    this%nTab,enth,output,tabkhi,tabklo)
        output = output/this%Re
      case ("B")
        call splint(this%enthTab,this%betaTab,  this%beta2Tab,  this%nTab,enth,output,tabkhi,tabklo)
      case default
        write(*,*) "Property doesn't exist!!!"
    end select
  end subroutine set_w_enth_table

  subroutine set_w_temp_table(this, temp, prop, output)
    use mod_math
    implicit none
    class(Table_EOSModel) :: this
    real(8), intent(IN) :: temp
    character(len=1), intent(IN) :: prop
    real(8), intent(OUT) :: output
    integer :: tabkhi,tabklo = 0 
    select case (prop)
      case ("H")
        call splint(this%tempTab,this%enthTab,this%enth2Tab,this%nTab,temp,output,tabkhi,tabklo)
      case default
        write(*,*) "Property doesn't exist!!!"
    end select
  end subroutine set_w_temp_table
  
  subroutine allocate_mem(this) 
    implicit none
    class(Table_EOSModel) :: this
    allocate(this%tempTab(1:this%nTab),this%rhoTab(1:this%nTab),   this%betaTab(1:this%nTab),  &
             this%muTab(1:this%nTab),  this%lamTab(1:this%nTab),   this%cpTab(1:this%nTab),    &
             this%enthTab(1:this%nTab),this%lamocpTab(1:this%nTab),this%temp2Tab(1:this%nTab), &
             this%rho2Tab(1:this%nTab),this%beta2Tab(1:this%nTab), this%mu2Tab(1:this%nTab),   &
             this%lam2Tab(1:this%nTab),this%cp2Tab(1:this%nTab),   this%enth2Tab(1:this%nTab), &
             this%lamocp2Tab(1:this%nTab))
  end subroutine allocate_mem

  subroutine read_table(this)
    implicit none
    class(Table_EOSModel) :: this
    integer i, ierr
    open(27,file=this%filename)
    do i=1,this%nTab
      read (27,*) this%tempTab(i),this%rhoTab(i), this%muTab(i),this%lamTab(i), &
                  this%cpTab(i),  this%enthTab(i),this%betaTab(i)
      this%lamocpTab(i) = this%lamTab(i)/this%cpTab(i)
    enddo
    close(27)
  end subroutine read_table

  subroutine calc_interp_coeff(this)
    use mod_math
    implicit none
    class(Table_EOSModel) :: this
    call spline(this%tempTab,this%enthTab,  this%nTab,this%enth2Tab)
    call spline(this%enthTab,this%rhoTab,   this%nTab,this%rho2Tab)
    call spline(this%enthTab,this%muTab,    this%nTab,this%mu2Tab)
    call spline(this%enthTab,this%lamTab,   this%nTab,this%lam2Tab)
    call spline(this%enthTab,this%cpTab,    this%nTab,this%cp2Tab)
    call spline(this%enthTab,this%lamocpTab,this%nTab,this%lamocp2Tab)
    call spline(this%enthTab,this%tempTab,  this%nTab,this%temp2Tab)
    call spline(this%enthTab,this%betaTab,  this%nTab,this%beta2Tab)
  end subroutine calc_interp_coeff

end module mod_eosmodels
