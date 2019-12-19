module mod_eos
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
    procedure(set_enth_w_qwall), deferred :: set_enth_w_qwall
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
    subroutine set_enth_w_qwall( this, qwall, enth_in, drp, enth_out )
      import :: EOSModel
      class(EOSModel) :: this
      real(8), intent(IN) :: qwall,drp
      real(8), intent(OUT) :: enth_out, enth_in
    end subroutine set_enth_w_qwall
  end interface

class(EOSModel),  allocatable :: eos_model
!****************************************************************************************


  !************************!
  !     Ideal Gas class    !
  !************************!
  
  type, extends(EOSModel), public :: IG_EOSModel
  contains
    procedure :: init => initialize_ig
    procedure :: set_w_enth => set_w_enth_ig
    procedure :: set_w_temp => set_w_temp_ig
    procedure :: set_enth_w_qwall => set_enth_w_qwall_ig

  end type IG_EOSModel

!****************************************************************************************

  !************************!
  !       Table class      !
  !************************!
  
  type, extends(EOSModel), public :: Table_EOSModel
  integer ntab
  character(len=40) filename
  real(8), dimension(:), private, allocatable :: tempTab,rhoTab,betaTab, muTab,lamTab, &
                                                 cpTab,enthTab,lamocpTab, temp2Tab,    &
                                                 rho2Tab,beta2Tab, mu2Tab,lam2Tab,     &
                                                 cp2Tab,enth2Tab,lamocp2Tab
  contains
    procedure :: init => initialize_table
    procedure :: set_w_enth => set_w_enth_table
    procedure :: set_w_temp => set_w_temp_table
    procedure :: set_enth_w_qwall => set_enth_w_qwall_table
    procedure, private :: allocate_mem => allocate_mem
    procedure, private :: read_table => read_table
    procedure, private :: calc_interp_coeff => calc_interp_coeff
    procedure :: check_qwall
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
        output = 1.!./(enth+1.)
      case ("L")
        output = 1./(this%Re*this%Pr)
      case ("C")
        output = 1.
      case ("V")
        output = 1./this%Re
      case ("B")
        output = 1.!1./(enth+1)
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

  subroutine set_enth_w_qwall_ig(this, qwall, enth_in, drp, enth_out)
    implicit none
    class(IG_EOSModel) :: this
    real(8), intent(IN) :: qwall, drp
    real(8), intent(OUT) :: enth_out,enth_in
    real(8) ::  ekh_imax
    call this%set_w_enth(enth_in,"L",ekh_imax)
    enth_out = enth_in + drp*qwall/(ekh_imax*this%Re*this%Pr) 
  end subroutine set_enth_w_qwall_ig
  
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

  subroutine set_enth_w_qwall_table(this, qwall, enth_in, drp, enth_out)
    implicit none
    class(Table_EOSModel) :: this
    real(8), intent(IN)   :: qwall,drp
    real(8), intent(OUT)  :: enth_out,enth_in
    integer :: niter,success
    real(8) :: error, error1
    success = 1
    error = 1000.0
    niter = 0

    if (enth_in.gt.2)    enth_in =  2.0
    if (enth_in.lt.-0.1) enth_in = -0.1
    enth_out = enth_in

    do while (abs(error).gt.1.0e-10)
      call this%check_qwall(qwall,enth_in, drp, enth_out, error)
      call this%check_qwall(qwall,enth_in, drp, enth_out+1.0e-8, error1)
      enth_out = enth_out - error/((error1-error)/1.0e-8)
      if (niter.gt.200) then
        error = 0.0
        success = 0
      endif
      niter = niter + 1
    enddo

    if (success.eq.0) then
      write (*,*) 'newton didnt converge, enthimax= ',enth_in,', ', niter, ', ', enth_out
    stop
    endif
  end subroutine set_enth_w_qwall_table

  subroutine check_qwall(this, qwall, enth_in, drp, enth_found, error)
    implicit none
    class(Table_EOSModel) :: this
    real(8), intent(IN)   :: qwall, enth_in, enth_found, drp
    real(8), intent(OUT)  :: error
    real(8) :: lamocpinter
    call this%set_w_enth(0.5*(enth_in+enth_found), 'L', lamocpinter)
    lamocpinter = lamocpinter*(this%Re*this%Pr)
    error = enth_found - enth_in - drp*qwall/lamocpinter 
  end subroutine check_qwall

end module mod_eos
