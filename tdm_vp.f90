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
    procedure :: get_sol => get_sol_vp_tdm
    procedure :: init_w_inflow => init_w_inflow_vp_tdm
    procedure :: get_profile => get_profile_vp_tdm
  end type VPrt_TurbDiffModel
  interface
    subroutine set_alphat_vp_tdm(this,u,w,rho,temp,mu,mui,lam_cp,mut,alphat)
      use mod_param, only : k1,i1
      import :: VPrt_TurbDiffModel
      class(VPrt_TurbDiffModel) :: this
      real(8),dimension(0:i1,0:k1),intent(IN) :: u,w,rho,temp,mu,mui,lam_cp, mut
      real(8),dimension(0:i1,0:k1),intent(OUT):: alphat
    end subroutine set_alphat_vp_tdm
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

  !************************!
  !       Bae class        !
  !************************!
  
  type,extends(VPrt_TurbDiffModel), public :: Bae_TurbDiffModel
  real(8), dimension(:,:), allocatable :: yplus
  real(8), allocatable :: Aplus, B
  contains
    procedure :: set_alphat  => set_alphat_bae
    procedure :: init  => init_mem_bae_tdm

  end type Bae_TurbDiffModel


contains
!****************************************************************************************

  !************************!
  !     VPrt routines      !
  !************************!

  subroutine init_mem_vp_tdm(this)
    use mod_param, only : k1,i1
    implicit none
    class(VPrt_TurbDiffModel) :: this
    allocate(this%Prt(0:i1,0:k1),this%Pr(0:i1,0:k1),this%mut_mu(0:i1,0:k1))
    allocate(this%Prtin(0:i1), this%alphatin(0:i1))
  end subroutine init_mem_vp_tdm

  subroutine get_sol_vp_tdm(this,Prt,epst,kt, Pkt, resKt, resEt)
    use mod_param, only : k1,i1
    class(VPrt_TurbDiffModel) :: this
    real(8),dimension(0:i1,0:k1), intent(OUT):: Prt,epst,kt, Pkt, resKt, resEt
    Prt  =this%Prt
    epst =0.    
    kt   =0.
    Pkt = 0.
    resEt = 0.
    resKt = 0
  end subroutine get_sol_vp_tdm

  subroutine get_profile_vp_tdm(this,p_prt,p_kt,p_epst,p_Pkt,k)
    use mod_param, only : i1
    class(VPrt_TurbDiffModel) :: this
    integer,                               intent(IN) :: k
    real(8),dimension(0:i1),          intent(OUT):: p_prt,p_kt,p_epst,p_Pkt
    p_prt = this%Prt(:,k)
    p_kt = 0.
    p_epst = 0.
    p_pkt = 0.
  end subroutine get_profile_vp_tdm

subroutine init_w_inflow_vp_tdm(this,Prtin,alphatin,ktin,epstin,pktin)
  use mod_param, only : i1,k1,k
  implicit none
  class(VPrt_TurbDiffModel) :: this
  real(8), dimension(0:i1), intent(IN) :: Prtin,alphatin,ktin,epstin,pktin
  this%Prtin(:)=Prtin
  this%alphatin(:)=alphatin
  do k=0,k1
    this%Prt(:,k) = this%Prtin(:)
  enddo
end subroutine init_w_inflow_vp_tdm

  !************************!
  ! Kays Crawford routines !
  !************************!

  subroutine set_alphat_kc(this,u,w,rho,temp,mu,mui,lam_cp,mut,alphat)
    use mod_param, only : k1,i1,kmax,imax,k,i
    implicit none
    class(KaysCrawford_TurbDiffModel) :: this
    real(8),dimension(0:i1,0:k1),intent(IN) :: u,w,rho,temp,mu,mui,lam_cp, mut
    real(8),dimension(0:i1,0:k1),intent(OUT):: alphat
    real(8) :: c1,c2,c3,c4

    ! constant for kays
    c1 = 0.5882
    c2 = 0.228
    c3 =-0.0441
    c4 =-5.165

    do k=1,kmax
      do i=1,imax
        this%Pr(i,k)= mu(i,k)/lam_cp(i,k)
        this%mut_mu(i,k)= mut(i,k)/mu(i,k)
        ! Approximation of the turbulent Prandlt number (W. Keys Turb. Pr, Where are we? 1992)
        this%Prt(i,k)= 1/(c1+c2*this%mut_mu(i,k)+c3*(this%mut_mu(i,k)**2.0)*(1-exp(c4/this%mut_mu(i,k))))
        alphat(i,k)= mut(i,k)/(this%Prt(i,k)*rho(i,k))
      enddo
    enddo
  end subroutine set_alphat_kc

  !************************!
  !     Kays routines      !
  !************************!

  subroutine set_alphat_kays(this,u,w,rho,temp,mu,mui,lam_cp,mut,alphat)
    use mod_param, only : k1,i1,kmax,imax,k,i
    implicit none
    class(Kays_TurbDiffModel) :: this
    real(8),dimension(0:i1,0:k1),intent(IN) :: u,w,rho,temp,mu,mui,lam_cp, mut
    real(8),dimension(0:i1,0:k1),intent(OUT):: alphat
    real(8) :: PeT

    do k=1,kmax
      do i=1,imax
        this%Pr(i,k)= mu(i,k)/lam_cp(i,k)
        this%mut_mu(i,k)= mut(i,k)/mu(i,k)
        ! Approximation of the turbulent Prandlt number (W. Keys Turb. Pr, Where are we? 1992)
        PeT = this%mut_mu(i,k)*this%Pr(i,k)
        if (this%mut_mu(i,k).lt.0.2) then
          this%Prt(i,k) = 1.07 
        else
          this%Prt(i,k) = (2.0/PeT)+0.85
        endif
        alphat(i,k)= mut(i,k)/(this%Prt(i,k)*rho(i,k))
      enddo
    enddo
  end subroutine set_alphat_kays

  !************************!
  !    Tang routines       !
  !************************!

  subroutine set_alphat_tang(this,u,w,rho,temp,mu,mui,lam_cp,mut,alphat)
    use mod_param, only : k1,i1,kmax,imax,k,i
    implicit none
    class(Tang_TurbDiffModel) :: this
    real(8),dimension(0:i1,0:k1),intent(IN) :: u,w,rho,temp,mu,mui,lam_cp, mut
    real(8),dimension(0:i1,0:k1),intent(OUT):: alphat
    real(8) :: Atang

    Atang = 15.0

    do k=1,kmax
      do i=1,imax
        this%Pr(i,k)= mu(i,k)/lam_cp(i,k)
        this%mut_mu(i,k)= mut(i,k)/mu(i,k)
        ! Approximation of the turbulent Prandlt number (Tang et al IJHMT 2016) 
        if (this%mut_mu(i,k).lt.0.2) then
          this%Prt(i,k) = 1.0 
        elseif (this%mut_mu(i,k).le.10.) then
          this%Prt(i,k) = 0.85 + this%Pr(i,k)/Atang
        else
          this%Prt(i,k) =0.85
        endif
        alphat(i,k)= mut(i,k)/(this%Prt(i,k)*rho(i,k))
      enddo
    enddo
  end subroutine set_alphat_tang

  !************************!
  !   Irrenfried routines  !
  !************************!

  subroutine set_alphat_irrenfried(this,u,w,rho,temp,mu,mui,lam_cp,mut,alphat)
    use mod_param, only : k1,i1,kmax,imax,k,i
    implicit none
    class(Irrenfried_TurbDiffModel) :: this
    real(8),dimension(0:i1,0:k1),intent(IN) :: u,w,rho,temp,mu,mui,lam_cp, mut
    real(8),dimension(0:i1,0:k1),intent(OUT):: alphat
    real(8) :: gam,Agam,Prtinf,PeT

    ! for Irrenfried and steiner
    Prtinf = 1.0  !turbulent prandlt in the free stream (assumption)
    do k=1,kmax
      do i=1,imax
        this%Pr(i,k)= mu(i,k)/lam_cp(i,k)
        this%mut_mu(i,k)= mut(i,k)/mu(i,k)
        ! Approximation of the turbulent Prandlt number (C. Irrenfried, H. Steiner IJHFF 2017)
        PeT      = this%mut_mu(i,k)*this%Pr(i,k)
        gam      = 1.0/(Prtinf+0.1*(this%Pr(i,k)**0.83))
        Agam     = ((2/Prtinf)-2*gam)**0.5
        this%Prt(i,k) = (gam+3.0*PeT*Agam-((3.0*PeT)**2.0)*(1-exp(-Agam/(3.0*PeT))))**(-1.0)
        alphat(i,k)= mut(i,k)/(this%Prt(i,k)*rho(i,k))
      enddo
    enddo
  end subroutine set_alphat_irrenfried

  !************************!
  !       Bae routines     !
  !************************!
  type(Bae_TurbDiffModel) function init_Bae_TurbDiffModel(name, Aplus, B)
    character(len=3),         intent(IN) :: name
    real(8) ,                 intent(IN) :: Aplus, B
    
    init_Bae_TurbDiffModel%name=name
    init_Bae_TurbDiffModel%Aplus = Aplus
    init_Bae_TurbDiffModel%B = B
  end function init_Bae_TurbDiffModel


  subroutine init_mem_bae_tdm(this)
    use mod_param, only : k1,i1,kmax,imax,k,i
    implicit none
    class(Bae_TurbDiffModel) :: this
    allocate(this%Prt(0:i1,0:k1),this%Pr(0:i1,0:k1),this%mut_mu(0:i1,0:k1),&
             this%yplus(0:i1,0:k1))
    allocate(this%Prtin(0:i1), this%alphatin(0:i1))
  end subroutine init_mem_bae_tdm


  subroutine set_alphat_bae(this,u,w,rho,temp,mu,mui,lam_cp,mut,alphat)
    use mod_param, only : k1,i1,kmax,imax,k,i
    use mod_mesh,  only : walldistu, walldist
    use mod_common,only : cp
    implicit none
    class(Bae_TurbDiffModel) :: this
    real(8),dimension(0:i1,0:k1),intent(IN) :: u,w,rho,temp,mu,mui,lam_cp, mut
    real(8),dimension(0:i1,0:k1),intent(OUT):: alphat
    integer :: km,kp,ip,im
    real(8),dimension(0:k1) ::   tauw
    real(8) :: dwdy, drhody, dcpdy, dTdy, wcenter, sigma_t,f1, f2, Prt0
    sigma_t = 0.9
        
    do k=1,kmax
      km = k-1
      kp = k+1
      tauw(k) = mui(imax,k)*0.5*(w(imax,km)+w(imax,k))/walldist(imax)

      do i=1,imax
        ip = i+1
        im = i-1

        this%yplus(i,k) = sqrt(rho(i,k))/mu(i,k)*(walldist(i))*tauw(k)**0.5       


        wcenter = (w(i,k)+w(i,km))/2. !velocity at cell center
        dwdy = ( &
                (w(ip,k)+w(ip,km))/2. &
               -(w(im,k)+w(im,km))/2. &
               )/(walldist(ip)-walldist(im))
        
        drhody = ( (rho(ip,k) + rho(i ,k))/2.0 &
                  -(rho(i, k) + rho(im,k))/2.0 &
                 )/(walldistu(ip)-walldistu(i))

        dTdy = (   (temp(ip,k) + temp(i ,k))/2.0 &
                  -(temp(i, k) + temp(im,k))/2.0 &
                 )/(walldistu(ip)-walldistu(i))

        dcpdy = (  (cp(ip,k) + cp(i ,k))/2.0 &
                  -(cp(i, k) + cp(im,k))/2.0 &
                 )/(walldistu(ip)-walldistu(i))

        f1 = 1-exp(this%yplus(i,k)/this%Aplus)
        f2 = 0.5*(1+tanh((this%B-this%yplus(i,k))/10.))
        !Prt0 = ( 1 + (w/rho)*abs(drho/dy / dw/dy) )/ (1 + T/rho * abs(drho/dy / dT/dy)) + T/cp * abs(dCp/dy / dT/dy ) 
        Prt0 = (1. + (wcenter/rho(i,k)) * abs(drhody/dwdy) ) &
              /(1+ (temp(i,k)/rho(i,k))*abs(drhody/(dTdy+1e-20)) + (temp(i,k)/cp(i,k))*abs(dcpdy/(dTdy+1e-20)))

        this%Prt(i,k) = sigma_t-f1*f2*(sigma_t-Prt0)
        alphat(i,k)= mut(i,k)/(this%Prt(i,k)*rho(i,k))

      enddo
    enddo
  end subroutine set_alphat_bae

end module vp_tdm
