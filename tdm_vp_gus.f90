module vp_tdm
  use mod_tdm
  !************************!
  !        VPrt class      !
  !************************!
  
  type, abstract, extends(TurbDiffModel) :: VPrt_TurbDiffModel
  real(8), dimension(:,:), allocatable :: Prt,Pr,mut_mu
  contains
    procedure(set_alphat_vp_tdm), deferred :: set_alphat
    procedure :: init => init_mem_vp_tdm
  end type VPrt_TurbDiffModel
  interface

    subroutine set_alphat_vp_tdm(this,u,w,rho,temp,mu,mui,ekh,alphat)
      import :: VPrt_TurbDiffModel
      class(VPrt_TurbDiffModel) :: this
      real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,temp,mu,mui,ekh
      real(8),dimension(0:this%i1,0:this%k1),intent(OUT):: alphat
    end subroutine set_alphat_vp_tdm


!    subroutine set_alphat_vp_tdm(this,mut,lam_cp,mu,alphat)
!      import :: VPrt_TurbDiffModel
!      class(VPrt_TurbDiffModel) :: this
!      real(8), dimension(0:this%i1,0:this%k1),intent(IN) :: mut, lam_cp, mu
!      real(8), dimension(0:this%i1,0:this%k1),intent(OUT):: alphat
!    end subroutine set_alphat_vp_tdm

  end interface

  !************************!
  !   Kays Crawford class  !
  !************************!

  type,extends(VPrt_TurbDiffModel), public :: KaysCrawford_TurbDiffModel
  contains
    procedure :: set_alphat  => set_alphat_kc
  end type KaysCrawford_TurbDiffModel

  !************************!
  !        Kays class      !
  !************************!

  type,extends(VPrt_TurbDiffModel), public :: Kays_TurbDiffModel
  contains
    procedure :: set_alphat  => set_alphat_kays
  end type Kays_TurbDiffModel

  !************************!
  !       Tang class       !
  !************************!
  
  type,extends(VPrt_TurbDiffModel), public :: Tang_TurbDiffModel
  contains
    procedure :: set_alphat  => set_alphat_tang
  end type Tang_TurbDiffModel

  !************************!
  !   Irrendfried class    !
  !************************!
  
  type,extends(VPrt_TurbDiffModel), public :: Irrenfried_TurbDiffModel
  contains
    procedure :: set_alphat  => set_alphat_irrenfried
  end type Irrenfried_TurbDiffModel

contains
!****************************************************************************************

  !************************!
  !     VPrt routines      !
  !************************!

  subroutine init_mem_vp_tdm(this)
      implicit none
      class(VPrt_TurbDiffModel) :: this
      allocate(this%Prt(0:this%i1,0:this%k1),this%Pr(0:this%i1,0:this%k1),this%mut_mu(0:this%i1,0:this%k1))
  end subroutine init_mem_vp_tdm

  !************************!
  ! Kays Crawford routines !
  !************************!

  subroutine set_alphat_kc(this,u,w,rho,temp,mu,mui,ekh,alphat)
    use mod_tm, only : turb_model
    implicit none
    class(KaysCrawford_TurbDiffModel) :: this
    real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,temp,mu,mui,ekh
    real(8),dimension(0:this%i1,0:this%k1),intent(OUT):: alphat
    real(8), dimension(0:this%i1,0:this%k1) :: mut
    integer :: i,k
    real(8) :: c1,c2,c3,c4

    mut = turb_model%mut

    ! constant for kays
    c1 = 0.5882
    c2 = 0.228
    c3 =-0.0441
    c4 =-5.165

    do k=1,this%kmax
      do i=1,this%imax
        this%Pr(i,k)= mu(i,k)/ekh(i,k)
        this%mut_mu(i,k)= mut(i,k)/mu(i,k)
        ! Approximation of the turbulent Prandlt number (W. Keys Turb. Pr, Where are we? 1992)
        this%Prt(i,k)= 1/(c1+c2*this%mut_mu(i,k)+c3*(this%mut_mu(i,k)**2.0)*(1-exp(c4/this%mut_mu(i,k))))
        alphat(i,k)= mut(i,k)/this%Prt(i,k)
      enddo
    enddo
  end subroutine set_alphat_kc

  !************************!
  !     Kays routines      !
  !************************!

  subroutine set_alphat_kays(this,u,w,rho,temp,mu,mui,ekh,alphat)
    use mod_tm, only : turb_model
    implicit none
    class(Kays_TurbDiffModel) :: this
    real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,temp,mu,mui,ekh
    real(8),dimension(0:this%i1,0:this%k1),intent(OUT):: alphat
    real(8), dimension(0:this%i1,0:this%k1) :: mut
    integer :: i,k
    real(8) :: PeT

    mut = turb_model%mut

    do k=1,this%kmax
      do i=1,this%imax
        this%Pr(i,k)= mu(i,k)/ekh(i,k)
        this%mut_mu(i,k)= mut(i,k)/mu(i,k)
        ! Approximation of the turbulent Prandlt number (W. Keys Turb. Pr, Where are we? 1992)
        PeT      = this%mut_mu(i,k)*this%Pr(i,k)
        if (this%mut_mu(i,k).lt.0.2) then
          this%Prt(i,k) = 1.07 
        else
          this%Prt(i,k) = (2.0/PeT)+0.85
        endif
        alphat(i,k)= mut(i,k)/this%Prt(i,k)
      enddo
    enddo
  end subroutine set_alphat_kays

  !************************!
  !    Tang routines       !
  !************************!

  subroutine set_alphat_tang(this,u,w,rho,temp,mu,mui,ekh,alphat)
    use mod_tm, only : turb_model
    implicit none
    class(Tang_TurbDiffModel) :: this
    real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,temp,mu,mui,ekh
    real(8),dimension(0:this%i1,0:this%k1),intent(OUT):: alphat
    real(8), dimension(0:this%i1,0:this%k1) :: mut
    integer :: i,k
    real(8) :: Atang

    mut = turb_model%mut

    Atang = 15.0

    do k=1,this%kmax
      do i=1,this%imax
        this%Pr(i,k)= mu(i,k)/ekh(i,k)
        this%mut_mu(i,k)= mut(i,k)/mu(i,k)
        ! Approximation of the turbulent Prandlt number (Tang et al IJHMT 2016) 
        if (this%mut_mu(i,k).lt.0.2) then
          this%Prt(i,k) = 1.0 
        elseif (this%mut_mu(i,k).le.10.) then
          this%Prt(i,k) = 0.85 + this%Pr(i,k)/Atang
        else
          this%Prt(i,k) =0.85
        endif
        alphat(i,k)= mut(i,k)/this%Prt(i,k)
      enddo
    enddo
  end subroutine set_alphat_tang

  !************************!
  !   Irrenfried routines  !
  !************************!

  subroutine set_alphat_irrenfried(this,u,w,rho,temp,mu,mui,ekh,alphat)
    use mod_tm, only : turb_model
    implicit none
    class(Irrenfried_TurbDiffModel) :: this
    real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,temp,mu,mui,ekh
    real(8),dimension(0:this%i1,0:this%k1),intent(OUT):: alphat
    real(8), dimension(0:this%i1,0:this%k1) :: mut
    integer :: i,k
    real(8) :: gam,Agam,Prtinf,PeT

    mut = turb_model%mut

    ! for Irrenfried and steiner
    Prtinf = 1.0  !turbulent prandlt in the free stream (assumption)
    do k=1,this%kmax
      do i=1,this%imax
        this%Pr(i,k)= mu(i,k)/ekh(i,k)
        this%mut_mu(i,k)= mut(i,k)/mu(i,k)
        ! Approximation of the turbulent Prandlt number (C. Irrenfried, H. Steiner IJHFF 2017)
        PeT      = this%mut_mu(i,k)*this%Pr(i,k)
        gam      = 1.0/(Prtinf+0.1*(this%Pr(i,k)**0.83))
        Agam     = ((2/Prtinf)-2*gam)**0.5
        this%Prt(i,k) = (gam+3.0*PeT*Agam-((3.0*PeT)**2.0)*(1-exp(-Agam/(3.0*PeT))))**(-1.0)
        alphat(i,k)= mut(i,k)/this%Prt(i,k)
      enddo
    enddo
  end subroutine set_alphat_irrenfried

end module vp_tdm
