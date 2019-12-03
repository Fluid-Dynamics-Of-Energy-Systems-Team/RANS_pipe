module sst_tm
  use mod_tm
  implicit none
  
!****************************************************************************************

  !************************!
  !         SST class      !
  !************************!

  type, extends(TurbModel), public :: SST_TurbModel
  real(8), dimension(:,:), allocatable :: Gk,cdKOM,Tt
  real(8), dimension(:),   allocatable :: omin,kin,bF1in
  contains
    procedure :: init => init_SST
    procedure :: init_sol => init_sol_SST
    procedure :: init_mem_SST
    procedure :: set_mut => set_mut_SST
    procedure :: advance_turb => advance_SST
    procedure :: set_bc => set_bc_SST
    procedure :: get_profile => get_profile_SST
    procedure :: get_sol => get_sol_SST
    procedure :: init_w_inflow => init_w_inflow_SST
    procedure :: solve_k_SST
    procedure :: solve_om_sst
    procedure :: diffusion_k_SST
    procedure :: diffusion_om_SST
    procedure :: production_SST
    procedure :: rhs_k_SST
    procedure :: rhs_om_SST
  end type SST_TurbModel

contains
!****************************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************************

   !************************!
   !      SST routines      !
   !************************!

subroutine init_SST(this)
    implicit none
    class(SST_TurbModel) :: this
    this%name='SST'
    call this%init_mem_SST()
    call this%init_sol()
end subroutine init_SST

subroutine init_sol_SST(this)
  class(SST_TurbModel) :: this
  integer i

  do i=0,this%i1
    this%k(i,:)  = 0.0!0.1
    this%om(i,:)  = 1.0
    this%bF1(i,:)  = 1.0
    this%kin = 0.0
    this%omin = 1.0
    this%bF1in = 1.0
    this%mutin = this%kin/this%omin
  enddo
end subroutine init_sol_SST

subroutine init_mem_SST(this)
    implicit none
    class(SST_TurbModel) :: this
    allocate(this%om (0:this%i1,0:this%k1),this%k(0:this%i1,0:this%k1),    &
             this%bF1(0:this%i1,0:this%k1),this%bF2(this%imax,this%kmax),  &
             this%Gk (0:this%i1,0:this%k1),this%Pk (0:this%i1,0:this%k1),  &
             this%Tt (0:this%i1,0:this%k1),this%cdKOM(this%imax,this%kmax),&
             this%yp (0:this%i1,0:this%k1))
    allocate(this%mutin(0:this%i1),this%Pkin (0:this%i1), &
             this%bF1in(0:this%i1),                       &
             this%omin (0:this%i1),this%kin  (0:this%i1))
end subroutine init_mem_SST

subroutine init_w_inflow_SST(this,Re, systemsolve)
    implicit none
    class(SST_TurbModel) :: this
    real(8), intent(IN) :: Re
    integer, intent(IN) :: systemsolve
    real(8), dimension(0:this%i1) :: dummy
    character(len=5)  :: Re_str
    integer           :: Re_int,k
    Re_int = int(Re)
    write(Re_str,'(I5.5)') Re_int
    if (systemsolve .eq. 1) open(29,file = 'pipe/Inflow_'//this%name//'_'//Re_str//'.dat',form='unformatted')
    if (systemsolve .eq. 2) open(29,file = 'channel/Inflow_'//this%name//'_'//Re_str//'.dat',form='unformatted')
    if (systemsolve .eq. 3) open(29,file = 'symchan/Inflow_'//this%name//'_'//Re_str//'.dat',form='unformatted')

    read(29) dummy(:),this%kin(:),dummy(:),dummy(:),this%omin(:), &
         dummy(:),this%mutin(:),dummy(:)
    close(29)
    do k=0,this%k1
      this%om(:,k) = this%omin(:)
      this%k(:,k) = this%kin(:)
    enddo
end subroutine init_w_inflow_SST



subroutine set_mut_SST(this,u,w,rho,mu,mui,walldist,Rp,dRp,dru,dz,mut)
  use module_mesh, only : mesh
  implicit none
  class(SST_TurbModel) :: this
  real(8), dimension(0:this%i1,0:this%k1), intent(IN) :: u, w, rho, mu, mui
  real(8), dimension(1:this%imax),         intent(IN) :: walldist
  real(8), dimension(0:this%i1),           intent(IN) :: Rp,dRp, dru
  real(8),                                 intent(IN) :: dz
  real(8), dimension(0:this%i1,0:this%k1), intent(OUT):: mut
  real(8), dimension(this%k1) :: tauw
  integer  im,ip,km,kp,i,k
  real(8)  sigma_om2,betaStar,gradkom,gamma1,gamma2,gamma3,gammaSST,zetaSST,StR, wallD
  real(8), dimension(0:this%k1) :: dzw, dzp

  dzw = mesh%dzw
  dzp = mesh%dzp

  !constants
  sigma_om2 = 0.856
  betaStar  = 0.09
      
  do k=1,this%kmax
    km=k-1
    kp=k+1
    tauw(k) = mui(this%imax,k)*0.5*(w(this%imax,km)+w(this%imax,k))/walldist(this%imax)
  
    do i=1,this%imax
      im=i-1
      ip=i+1
      this%yp(i,k)     = sqrt(rho(i,k))/mu(i,k)*(walldist(i))*tauw(k)**0.5   ! ystar
      wallD     = walldist(i)
      ! Vorticity rate
      StR = ( & 
            ( -( (w(ip,km)+w(ip,k)+w(i,km)+w(i ,k))/4.-(w(im,km)+w(im,k)+w(i,km)+w(i,k))/4. )/dru(i) &
          + ( (u(i,kp) +u(im,kp)+u(i,k)+u(im,k))/4.-(u(im,km)+u(i,km)+u(im,k)+u(i,k))/4.    )/(dz) &
            )**2.)

      ! StR = ( & 
      !       ( -( (w(ip,km)+w(ip,k)+w(i,km)+w(i ,k))/4.-(w(im,km)+w(im,k)+w(i,km)+w(i,k))/4. )/dru(i) &
      !     + ( (u(i,kp) +u(im,kp)+u(i,k)+u(im,k))/4.-(u(im,km)+u(i,km)+u(im,k)+u(i,k))/4.    )/dzw(k) &
      !       )**2.)


      StR = StR**0.5
               
      gradkom = ((this%k(ip,k) - this%k(im,k))/(dRp(i)+dRp(im))) *((this%om(ip,k) - this%om(im,k))/(dRp(i)+dRp(im))) &
               +((this%k(i,kp) - this%k(i,km))/(2.0*dz))         *((this%om(i,kp) - this%om(i,km))/(2.0*dz))

      ! gradkom = ((this%k(ip,k) - this%k(im,k))/(dRp(i)+dRp(im))) *((this%om(ip,k) - this%om(im,k))/(dRp(i)+dRp(im))) &
      !          +((this%k(i,kp) - this%k(i,km))/(2.0*dzp(k)))     *((this%om(i,kp) - this%om(i,km))/(2.0*dzp(k)))


      this%cdKOM(i,k) = 2.0*sigma_om2*rho(i,k)/this%om(i,k)*gradkom;
              
      gamma1   = 500.0*mu(i,k)/(rho(i,k)*this%om(i,k)*wallD**2.0);
      gamma2   = 4.0*sigma_om2*rho(i,k)*this%k(i,k)/(wallD*wallD*max(this%cdKOM(i,k), 1.0e-20));
      gamma3   = (this%k(i,k)**0.5)/(betaStar*this%om(i,k)*wallD)

      gammaSST = min(max(gamma1, gamma3), gamma2)
      this%bF1(i,k) = tanh(gammaSST**4.0)

      gammaSST = max(2.0*gamma3, gamma1)
      this%bF2(i,k) = tanh(gammaSST**2.0)

      zetaSST  = max(0.31*this%om(i,k), this%bF2(i,k)*StR)
      mut(i,k) = rho(i,k)*this%k(i,k)/zetaSST !!! NOTE this is the correct one !!!!
      ! mut(i,k) = rho(i,k)*this%k(i,k)/this%om(i,k)

    enddo
  enddo
end subroutine set_mut_SST

subroutine advance_SST(this,u,w,rho,mu,mui,muk,mut,beta,temp,&
                       Ru,Rp,dru,drp,dz,walldist,            &
                       alpha1,alpha2,alpha3,                 &
                       modification,rank,centerBC,periodic,  &
                       residual1, residual2, residual3)
  implicit none
  class(SST_TurbModel) :: this
  real(8),dimension(0:this%i1,0:this%k1), intent(IN) :: u, w, rho,mu,mui,muk,mut,beta,temp
  real(8),dimension(0:this%i1),           intent(IN) :: dru,ru,rp,drp
  real(8),dimension(1:this%imax),         intent(IN) :: walldist
  real(8),                                intent(IN) :: dz, alpha1, alpha2, alpha3
  integer,                                intent(IN) :: modification,rank,centerBC,periodic
  real(8),                                intent(OUT):: residual1, residual2, residual3
  real(8),dimension(0:this%i1,0:this%k1)             :: rho_mod

  !modification: 1, our modification | 2, Aupoix modification
  if ((modification == 1) .or. (modification == 2)) then
    rho_mod = rho
  else
    rho_mod = 1.0
  endif

  call this%production_SST(u,w,rho,temp,mut,beta,Rp,Ru,dRu,dRp,dz)
  call this%solve_k_SST(residual1,u,w,rho,mu,mui,muk,mut,rho_mod, &
                        Ru,Rp,dru,drp,dz, &
                        alpha1,modification,rank,centerBC,periodic)
  call this%solve_om_SST(residual2,u,w,rho,mu,mui,muk,mut,beta,temp,rho_mod, &
                        Ru,Rp,dru,drp,dz, &
                        alpha2,modification,rank,centerBC,periodic)
end subroutine advance_SST

subroutine set_bc_SST(this,mu,rho,walldist,centerBC,periodic,rank,px)
  use mod_mesh, only : top_bcvalue, bot_bcvalue,top_bcnovalue, bot_bcnovalue
  implicit none
  class(SST_TurbModel) :: this
  real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: rho,mu
  real(8),dimension(1:this%imax),        intent(IN) :: walldist
  integer,                               intent(IN) :: centerBC,periodic, rank, px
  real(8),dimension(0:this%k1) :: BCvalue
  real(8),dimension(0:this%i1) :: tmp
  real(8) :: botBCvalue, topBCvalue
  integer :: k

  do k = 0,this%k1
    this%k(0,k)         = bot_bcnovalue(k)*this%k(1,k)         !symmetry or 0 value
    this%k(this%i1,k)   = top_bcnovalue(k)*this%k(this%imax,k) !symmetry or 0 value
    this%bF1(0,k)= bot_bcnovalue(k)*this%bF1(1,k)
    this%bF1(this%i1,k)= top_bcnovalue(k)*this%bF1(this%imax,k)
    botBCvalue = 60.0/0.075*mu(1,k)/rho(1,k)/walldist(1)**2                                                           !bcvalue bot
    this%om(0,k)       = (1.-bot_bcvalue(k))*(2.0*botBCvalue-this%om(1,k))         + bot_bcvalue(k)*this%om(1,k)         !symmetry or bc value
    topBCvalue = 60.0/0.075*mu(this%imax,k)/rho(this%imax,k)/walldist(this%imax)**2                                   !bcvalue top
    this%om(this%i1,k) = (1.-top_bcvalue(k))*(2.0*topBCvalue-this%om(this%imax,k)) + top_bcvalue(k)*this%om(this%imax,k) !symmetry or bc value
  enddo

  ! this%k  (this%i1,:)= -this%k(this%imax,:)
  ! BCvalue(:)         = 60.0/0.075*mu(this%imax,:)/rho(this%imax,:)/walldist(this%imax)**2
  ! this%om (this%i1,:)= 2.0*BCvalue(:) - this%om(this%imax,:)
  ! this%bF1(this%i1,:)=  this%bF1(this%imax,:)
  
  ! !channel
  ! if (centerBC.eq.-1) then
  !   this%k(0,:) = -this%k(1,:)
  !   BCvalue(:)  = 60.0/0.075*mu(1,:)/rho(1,:)/walldist(1)**2
  !   this%om(0,:)= 2.0*BCvalue(:) - this%om(1,:)
  ! endif

  ! !pipe/BL
  ! if (centerBC.eq.1) then
  !   this%k (0,:) = this%k  (1,:)
  !   this%om(0,:) = this%om (1,:)
  !   this%bF1(0,:)= this%bF1(1,:)
  ! endif

  call shiftf(this%k,  tmp,rank); this%k  (:,0)      =tmp(:);
  call shiftf(this%om ,tmp,rank); this%om (:,0)      =tmp(:);
  call shiftf(this%bF1,tmp,rank); this%bF1(:,0)      =tmp(:)
  call shiftb(this%k  ,tmp,rank); this%k  (:,this%k1)=tmp(:);
  call shiftb(this%om ,tmp,rank); this%om (:,this%k1)=tmp(:);
  call shiftb(this%bF1,tmp,rank); this%bF1(:,this%k1)=tmp(:)

  !developing
  if (periodic.eq.1) return
  if (rank.eq.0) then
    this%k  (:,0) = this%kin(:)
    this%om (:,0) = this%omin(:)
    !this%bF1(:,0) = this%bF1(:,1)  ! ATTENTION (THIS WAS THE ORIGINAL)
    this%bF1(:,0) = this%bF1in(:) 
  endif

  if (rank.eq.px-1) then
    this%k  (:,this%k1) = 2.0*this%k  (:,this%kmax)-this%k  (:,this%kmax-1)
    this%om (:,this%k1) = 2.0*this%om (:,this%kmax)-this%om (:,this%kmax-1)
    this%bF1(:,this%k1) = 2.0*this%bF1(:,this%kmax)-this%bF1(:,this%kmax-1)
  endif
  
end subroutine set_bc_SST

subroutine get_profile_SST(this,p_nuSA,p_k,p_eps,p_om,p_v2,p_Pk,p_bF1,p_bF2,yp,k)
    class(SST_TurbModel) :: this
    integer,                               intent(IN) :: k
    real(8),dimension(0:this%i1),          intent(OUT):: p_nuSA,p_k,p_eps,p_om,p_v2,p_Pk,p_bF1,yp
    real(8),dimension(1:this%imax),        intent(OUT):: p_bF2

    p_nuSA(:)=0
    p_k(:)   =this%k(:,k)
    p_eps(:) =0
    p_om(:)  =this%om(:,k)
    p_Pk(:)  =this%Pk(:,k)
    p_bF1(:) = this%bF1(:,k)
    p_bF2(:) = this%bF2(:,k)
    p_v2(:)  = 0
    yp(:)    = this%yp(:,k)

end subroutine get_profile_SST

subroutine get_sol_SST(this,nuSA,k,eps,om,v2,yp)
    class(SST_TurbModel) :: this
    real(8),dimension(0:this%i1,0:this%k1), intent(OUT):: nuSA,k,eps,om,v2,yp
    nuSA=0
    k   =this%k    
    eps =0
    v2  =0
    om  =this%om
    yp  = this%yp
end subroutine get_sol_SST



subroutine production_SST(this,u,w,temp,rho,mut,beta,Rp,Ru,dRu,dRp,dz)
  use module_mesh, only : mesh
  implicit none
  class(SST_TurbModel) :: this
  real(8), dimension(0:this%i1,0:this%k1), intent(IN) :: u,w,temp,rho,mut,beta
  real(8), dimension(0:this%i1),           intent(IN) :: Rp,Ru,dRu,dRp
  real(8),                                 intent(IN) :: dz
  real(8), dimension(0:this%i1,0:this%k1) :: div
  integer                       :: im,ip,km,kp,ib,ie,kb,ke,i,k 
  real(8)                       :: sigma_om1,sigma_om2, &
                                   beta_1,beta_2,betaStar, &
                                   alfa_1,alfa_2, ctheta,Fr_1
  real(8), dimension(0:this%k1) :: dzw, dzp

  dzw = mesh%dzw
  dzp = mesh%dzp

  Fr_1      = 0.0 !!!NOTE: this was originally in the param!!!!
  ctheta    = 0.3 !!!NOTE: this was originally in the param!!!!
  sigma_om1 = 0.5
  sigma_om2 = 0.856
  beta_1    = 0.075
  beta_2    = 0.0828
  betaStar  = 0.09
  alfa_1    = beta_1/betaStar - sigma_om1*(0.41**2.0)/(betaStar**0.5)
  alfa_2    = beta_2/betaStar - sigma_om2*(0.41**2.0)/(betaStar**0.5)

  ib = 1
  ie = this%i1-1
  kb = 1
  ke = this%k1-1

  do k=kb,ke
    kp=k+1
    km=k-1
    do i=ib,ie
      ip=i+1
      im=i-1
      ! Production of turbulent kinetic energy
      this%Pk(i,k) = mut(i,k)*( &
        2.*(((w(i,k)-w(i,km))/dz)**2. +((u(i,k)-u(im,k))/dRu(i))**2.+((u(i,k)+u(im,k))/(2.*Rp(i)))**2.)+ &
        ( &
          ((w(ip,km)+w(ip,k)+w(i,km)+w(i,k))/4.-(w(im,km)+w(im,k)+w(i,km)+w(i,k))/4.)/dRu(i) + &
          ((u(i,kp)+u(im,kp)+u(i,k)+u(im,k))/4.-(u(im,km)+u(i,km)+u(im,k)+u(i,k))/4.)/dz &
        )**2.)

      ! this%Pk(i,k) = mut(i,k)*( &
      !   2.*(((w(i,k)-w(i,km))/dzw(k))**2. +((u(i,k)-u(im,k))/dRu(i))**2.+((u(i,k)+u(im,k))/(2.*Rp(i)))**2.)+ &
      !   ( &
      !     ((w(ip,km)+w(ip,k)+w(i,km)+w(i,k))/4.-(w(im,km)+w(im,k)+w(i,km)+w(i,k))/4.)/dRu(i) + &
      !     ((u(i,kp)+u(im,kp)+u(i,k)+u(im,k))/4.-(u(im,km)+u(i,km)+u(im,k)+u(i,k))/4.)/dzw(k) &
      !   )**2.)


      div(i,k) =(Ru(i)*u(i,k)-Ru(im)*u(im,k))/(Rp(i)*dRu(i)) +(w(i,k)-w(i,km))/dz

      ! div(i,k) =(Ru(i)*u(i,k)-Ru(im)*u(im,k))/(Rp(i)*dRu(i)) +(w(i,k)-w(i,km))/dzw(k)



      this%Pk(i,k) = this%Pk(i,k) - 2./3.*(rho(i,k)*this%k(i,k)+mut(i,k)*(div(i,k)))*(div(i,k))

      ! turbulent time scale
      this%Tt(i,k)   = 1.0/this%om(i,k)   ! 0.31 cmu/omega
               
      this%Gk(i,k)=-ctheta*beta(i,k)*Fr_1*this%Tt(i,k) &
             *(mut(i,k)*( &
                        ((w(ip,km)+w(ip,k)+w(i,km)+w(i,k))/4.-(w(im,km)+w(im,k)+w(i,km)+w(i,k))/4.)/dRu(i) + &
                        ((u(i,kp)+u(im,kp)+u(i,k)+u(im,k))/4.-(u(im,km)+u(i,km)+u(im,k)+u(i,k))/4.)/(dz) &
                        )* &
              (temp(ip,k)-temp(im,k))/(dRp(i)+dRp(im))  &
              ) &
             +(2.*mut(i,k)*((w(i,k)-w(i,km))/dz-2./3.*(rho(i,k)*this%k(i,k)))*(temp(i,kp)-temp(i,km))/(2.*dz))




      ! this%Gk(i,k)=-ctheta*beta(i,k)*Fr_1*this%Tt(i,k) &
      !        *(mut(i,k)*( &
      !                   ((w(ip,km)+w(ip,k)+w(i,km)+w(i,k))/4.-(w(im,km)+w(im,k)+w(i,km)+w(i,k))/4.)/dRu(i) + &
      !                   ((u(i,kp)+u(im,kp)+u(i,k)+u(im,k))/4.-(u(im,km)+u(i,km)+u(im,k)+u(i,k))/4.)/(dzw(k)) &
      !                   )* &
      !         (temp(ip,k)-temp(im,k))/(dRp(i)+dRp(im))  &
      !         ) &
      !        +(2.*mut(i,k)*((w(i,k)-w(i,km))/dzw(k)-2./3.*(rho(i,k)*this%k(i,k)))*(temp(i,kp)-temp(i,km))/(2.*dzp(k)))




      ! this%Gk(i,k) = this%Gk(i,k) + ctheta*beta(i,k)*Fr_1*this%Tt(i,k)*2./3.*mut(i,k)*div(i,k)*(temp(i,kp)-temp(i,km))/(2.*dz)

      this%Gk(i,k) = this%Gk(i,k) + ctheta*beta(i,k)*Fr_1*this%Tt(i,k)*2./3.*mut(i,k)*div(i,k)*(temp(i,kp)-temp(i,km))/(2.*dzp(k))


    enddo
  enddo
end subroutine production_SST

subroutine solve_k_SST(this,resK,u,w,rho,mu,mui,muk,mut,rho_mod, &
                       Ru,Rp,dru,drp,dz, &
                       alphak,modification,rank,centerBC,periodic)
  use mod_math
  use mod_mesh, only : top_bcnovalue, bot_bcnovalue
  implicit none
  class(SST_TurbModel) :: this
  real(8),dimension(0:this%i1,0:this%k1), intent(IN) :: u, w, rho,mu,mui,muk,mut,rho_mod
  real(8),dimension(0:this%i1),           intent(IN) :: dru,ru,rp,drp
  real(8),                                intent(IN) :: dz, alphak
  integer,                                intent(IN) :: modification,rank,centerBC,periodic
  real(8),                                intent(OUT):: resK
  real(8), dimension(0:this%i1,0:this%k1) :: dnew,dimpl,sigmakSST
  real(8), dimension(this%imax) :: a,b,c,rhs
  integer i,k
  
  resK  = 0.0; dnew=0.0; dimpl = 0.0;

  call advecc(dnew,dimpl,this%k,u,w,Ru,Rp,dru,dz,this%i1,this%k1,rank,periodic,.true.)
  call this%rhs_k_SST(dnew,dimpl,rho)

  ! calculating constant with blending function factor
  sigmakSST = 0.85*this%bF1 + 1.0*(1.0 - this%bF1)
  sigmakSST = 1.0/sigmakSST
  
  call this%diffusion_k_sst(dnew,this%k,mu,mui,muk,mut,sigmakSST,rho,dz,modification)

  do k=1,this%kmax
    do i=1,this%imax
      if ((modification == 0) .or. (modification == 1)) then
        a(i) = (mui(i-1,k)+(mut(i,k)+mut(i-1,k))/(sigmakSST(i,k)+sigmakSST(i-1,k)))/(0.5*(rho_mod(i-1,k)+rho_mod(i,k)))**0.5
        a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/rho(i,k)/rho_mod(i,k)**0.5

        c(i) = (mui(i  ,k)+(mut(i,k)+mut(i+1,k))/(sigmakSST(i,k)+sigmakSST(i+1,k)))/(0.5*(rho_mod(i+1,k)+rho_mod(i,k)))**0.5
        c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/rho(i,k)/rho_mod(i,k)**0.5
      else if (modification == 2) then
        a(i) = (mui(i-1,k)+(mut(i,k)+mut(i-1,k))/(sigmakSST(i,k)+sigmakSST(i-1,k)))/(0.5*(rho_mod(i-1,k)+rho_mod(i,k)))
        a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/rho(i,k)

        c(i) = (mui(i  ,k)+(mut(i,k)+mut(i+1,k))/(sigmakSST(i,k)+sigmakSST(i+1,k)))/(0.5*(rho_mod(i+1,k)+rho_mod(i,k)))
        c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/rho(i,k)
      endif

      b(i) = (rho_mod(i,k)*(-a(i)-c(i)) + dimpl(i,k))
      a(i) = a(i)*rho_mod(i-1,k)
      c(i) = c(i)*rho_mod(i+1,k)
      rhs(i) = dnew(i,k)  + ((1-alphak)/alphak)*b(i)*this%k(i,k)
    enddo
    
    i=1
    b(i) = b(i)+bot_bcnovalue(k)*a(i)
    rhs(i) = dnew(i,k)  + ((1-alphak)/alphak)*b(i)*this%k(i,k) 

    i=this%imax
    b(i) = b(i)+top_bcnovalue(k)*c(i)
    rhs(i) = dnew(i,k)  + ((1-alphak)/alphak)*b(i)*this%k(i,k)

    call matrixIdir(this%imax,a,b/alphak,c,rhs)

    do i=1,this%imax
      resK = resK + ((this%k(i,k) - rhs(i))/(this%k(i,k)+1.0e-20))**2.0
      this%k(i,k) = max(rhs(i), 1.0e-8)
    enddo
  enddo
end subroutine solve_k_SST

subroutine rhs_k_SST(this,putout,dimpl,rho)
  implicit none
  class(SST_TurbModel) :: this
  real(8), dimension(0:this%i1,0:this%k1), intent(IN) :: rho
  real(8), dimension(0:this%i1,0:this%k1), intent(OUT):: putout, dimpl

  integer ib,ie,kb,ke,i,k 
  real(8)  betaStar

  betaStar  = 0.09
  ib = 1
  ie = this%i1-1

  kb = 1
  ke = this%k1-1

  do k=kb,ke
    do i=ib,ie
      ! k- equation
      putout(i,k) = putout(i,k) + ( this%Pk(i,k) + this%Gk(i,k) )/rho(i,k)
      dimpl(i,k)  = dimpl(i,k)  + betaStar*this%om(i,k)            ! note, betaStar*rho*k*omega/(rho*k), set implicit and divided by density
    enddo
  enddo
end subroutine rhs_k_SST

subroutine diffusion_k_sst(this,putout,putin,ek,eki,ekk,ekmt,sigma,rho,dz,modification)
  implicit none
  class(SST_TurbModel) :: this
  real(8), dimension(0:this%i1,0:this%k1), intent(IN) :: putin,ek,eki,ekk,ekmt,sigma,rho
  real(8)                                , intent(IN) :: dz
  integer                                , intent(IN) :: modification
  real(8), dimension(0:this%i1,0:this%k1), intent(OUT):: putout

  integer   km,kp,i,k

  if (modification == 1) then       ! Inverse SLS
    do k=1,this%kmax
      kp=k+1
      km=k-1
      do i=1,this%imax
        putout(i,k) = putout(i,k) + 1.0/rho(i,k)/sqrt(rho(i,k))*( &
          ( (ekk(i,k ) + (ekmt(i,k)+ekmt(i,kp))/(sigma(i,k)+sigma(i,kp)))/ &
          sqrt(0.5*(rho(i,k)+rho(i,kp)))*(rho(i,kp)*putin(i,kp)-rho(i,k )*putin(i,k )) &
          -(ekk(i,km) + (ekmt(i,k)+ekmt(i,km))/(sigma(i,k)+sigma(i,km)))/ &
          sqrt(0.5*(rho(i,k)+rho(i,km)))*(rho(i,k )*putin(i,k )-rho(i,km)*putin(i,km)) &
          )/(dz*dz)   )
        ! putout(i,k) = putout(i,k) + 1.0/rho(i,k)/sqrt(rho(i,k))*( &
        !   ( (ekk(i,k ) + (ekmt(i,k)+ekmt(i,kp))/(sigma(i,k)+sigma(i,kp)))/ &
        !   sqrt(0.5*(rho(i,k)+rho(i,kp)))*(rho(i,kp)*putin(i,kp)-rho(i,k )*putin(i,k ))/dzp(k) &
        !   -(ekk(i,km) + (ekmt(i,k)+ekmt(i,km))/(sigma(i,k)+sigma(i,km)))/ &
        !   sqrt(0.5*(rho(i,k)+rho(i,km)))*(rho(i,k )*putin(i,k )-rho(i,km)*putin(i,km))/dzp(km) &
        !   )/dzw(k))
      enddo
    enddo
  elseif (modification == 2) then   ! Aupoix
    do k=1,this%kmax
      kp=k+1
      km=k-1
      do i=1,this%imax
        putout(i,k) = putout(i,k) + 1.0/rho(i,k)*( &
          ( (ekk(i,k ) + (ekmt(i,k)+ekmt(i,kp))/(sigma(i,k)+sigma(i,kp)))/ &
          (0.5*(rho(i,k)+rho(i,kp)))*(rho(i,kp)*putin(i,kp)-rho(i,k )*putin(i,k )) &
          -(ekk(i,km) + (ekmt(i,k)+ekmt(i,km))/(sigma(i,k)+sigma(i,km)))/ &
          (0.5*(rho(i,k)+rho(i,km)))*(rho(i,k )*putin(i,k )-rho(i,km)*putin(i,km)) &
          )/(dz*dz)   )
        ! putout(i,k) = putout(i,k) + 1.0/rho(i,k)*( &
        !   ( (ekk(i,k ) + (ekmt(i,k)+ekmt(i,kp))/(sigma(i,k)+sigma(i,kp)))/ &
        !   (0.5*(rho(i,k)+rho(i,kp)))*(rho(i,kp)*putin(i,kp)-rho(i,k )*putin(i,k ))/dzp(k) &
        !   -(ekk(i,km) + (ekmt(i,k)+ekmt(i,km))/(sigma(i,k)+sigma(i,km)))/ &
        !   (0.5*(rho(i,k)+rho(i,km)))*(rho(i,k )*putin(i,k )-rho(i,km)*putin(i,km))/dzp(km) &
        !   )/dzw(k))

      enddo
    enddo
  else                               ! Standard
    do k=1,this%kmax
      kp=k+1
      km=k-1
      do i=1,this%imax
        putout(i,k) = putout(i,k) + 1.0/rho(i,k)*( &
          ( (ekk(i,k ) + (ekmt(i,k)+ekmt(i,kp))/(sigma(i,k)+sigma(i,kp)))*(putin(i,kp)-putin(i,k )) &
           -(ekk(i,km) + (ekmt(i,k)+ekmt(i,km))/(sigma(i,k)+sigma(i,km)))*(putin(i,k )-putin(i,km))  &
          )/(dz*dz)   )
        ! putout(i,k) = putout(i,k) + 1.0/rho(i,k)*( &
        !   ( (ekk(i,k ) + (ekmt(i,k)+ekmt(i,kp))/(sigma(i,k)+sigma(i,kp)))*(putin(i,kp)-putin(i,k ))/dzp(k ) &
        !    -(ekk(i,km) + (ekmt(i,k)+ekmt(i,km))/(sigma(i,k)+sigma(i,km)))*(putin(i,k )-putin(i,km))/dzp(km) &
        !   )/dzw(k))
      enddo
    enddo
  endif
end subroutine diffusion_k_sst


subroutine solve_om_sst(this,resOm,u,w,rho,mu,mui,muk,mut,beta,temp,rho_mod, &
                       Ru,Rp,dru,drp,dz, &
                       alphae,modification,rank,centerBC,periodic)
  use mod_math
  use mod_mesh, only : top_bcvalue, bot_bcvalue
  implicit none
  class(SST_TurbModel) :: this
  real(8),dimension(0:this%i1,0:this%k1), intent(IN) :: u, w, rho,mu,mui,muk,mut,beta,temp,rho_mod
  real(8),dimension(0:this%i1),           intent(IN) :: dru,ru,rp,drp
  real(8),                                intent(IN) :: dz, alphae
  integer,                                intent(IN) :: modification,rank,centerBC,periodic
  real(8),                                intent(OUT):: resOm
  real(8), dimension(0:this%i1,0:this%k1) :: dnew,dimpl,sigmakSST
  real(8), dimension(this%imax) :: a,b,c,rhs
  integer i,k
  resOm = 0.0

  dnew=0.0; dimpl = 0.0;
  call advecc(dnew,dimpl,this%om,u,w,Ru,Rp,dru,dz,this%i1,this%k1,rank,periodic,.true.)
  call this%rhs_om_sst(dnew,dimpl,this%k,u,w,temp,rho,beta,Rp,Ru,dRu,dRp,dz)


  ! calculating constant with blending function factor
  sigmakSST = 0.5*this%bF1 + 0.856*(1.0 - this%bF1)
  sigmakSST = 1.0/sigmakSST
  call this%diffusion_om_SST(dnew,this%om,mu,mui,muk,mut,sigmakSST,rho,Ru,Rp,dru,dz,modification)
  do k=1,this%kmax
    do i=1,this%imax
  
      a(i) = (mui(i-1,k)+(mut(i,k)+mut(i-1,k))/(sigmakSST(i,k)+sigmakSST(i-1,k)))/(0.5*(rho_mod(i-1,k)+rho_mod(i,k)))**0.5
      a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/rho(i,k)/rho_mod(i,k)**0.5
  
      c(i) = (mui(i  ,k)+(mut(i,k)+mut(i+1,k))/(sigmakSST(i,k)+sigmakSST(i+1,k)))/(0.5*(rho_mod(i+1,k)+rho_mod(i,k)))**0.5
      c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/rho(i,k)/rho_mod(i,k)**0.5
  
      b(i) = ((-a(i)-c(i))*rho_mod(i,k)**0.5 + dimpl(i,k))
  
      a(i) = a(i)*rho_mod(i-1,k)**0.5
      c(i) = c(i)*rho_mod(i+1,k)**0.5
  
      rhs(i) = dnew(i,k)  + ((1-alphae)/alphae)*b(i)*this%om(i,k)
    enddo
  
    i=1
    b(i)=b(i)+bot_bcvalue(k)*a(i)
    rhs(i) = dnew(i,k) - (1-bot_bcvalue(k))*a(i)*this%om(i-1,k) + ((1-alphae)/alphae)*b(i)*this%om(i,k)  !wall
    
    ! if (centerBC.eq.-1) then
    !   rhs(i) = dnew(i,k) - a(i)*this%om(i-1,k) + ((1-alphae)/alphae)*b(i)*this%om(i,k)  !wall
    ! else
    !   b(i)=b(i)+a(i) !symmetry
    ! endif

    i = this%imax
    b(i)=b(i)+top_bcvalue(k)*c(i)
    rhs(i) = dnew(i,k) - (1-top_bcvalue(k))*c(i)*this%om(i+1,k) + ((1-alphae)/alphae)*b(i)*this%om(i,k)  !wall
    ! rhs(i) = dnew(i,k) - c(i)*this%om(i+1,k) + ((1-alphae)/alphae)*(b(i)*this%om(i,k)) !wall
  
    call matrixIdir(this%imax,a,b/alphae,c,rhs)
  
    do i=1,this%imax
      resOm = resOm + ((this%om(i,k) - rhs(i))/(this%om(i,k)+1.0e-20))**2.0
      this%om(i,k) = max(rhs(i), 1.0e-8)
    enddo
  enddo
end subroutine solve_om_sst

subroutine rhs_om_sst(this,putout,dimpl,putink,u,w,temp,rho,beta,Rp,Ru,dRu,dRp,dz)
  use module_mesh, only : mesh
  implicit none
  class(SST_TurbModel) :: this
  real(8), dimension(0:this%i1,0:this%k1), intent(IN) :: u,w,temp,rho,beta, putink
  real(8), dimension(0:this%i1),           intent(IN) :: Rp,Ru,dRu,dRp
  real(8),                                 intent(IN) :: dz
  real(8), dimension(0:this%i1,0:this%k1), intent(OUT):: putout, dimpl

  integer km,kp,im,ip,ib,ie,kb,ke,i,k
  real(8), dimension(0:this%i1,0:this%k1) :: div
  real(8) sigma_om1,sigma_om2,beta_1,beta_2,betaStar,alfa_1,alfa_2,alfaSST,betaSST,StR,GtR,ctheta,Fr_1
  real(8), dimension(0:this%k1) :: dzw, dzp

  dzw = mesh%dzw
  dzp = mesh%dzp

  Fr_1      = 0.0 !!!NOTE: this was originally in the param!!!!
  ctheta    = 0.3 !!!NOTE: this was originally in the param!!!!
  sigma_om1 = 0.5
  sigma_om2 = 0.856
  beta_1    = 0.075
  beta_2    = 0.0828
  betaStar  = 0.09
  alfa_1    = beta_1/betaStar - sigma_om1*(0.41**2.0)/(betaStar**0.5)
  alfa_2    = beta_2/betaStar - sigma_om2*(0.41**2.0)/(betaStar**0.5)

  ib = 1
  ie = this%i1-1

  kb = 1
  ke = this%k1-1

  do k=kb,ke
    kp=k+1
    km=k-1
    do i=ib,ie
      ip=i+1
      im=i-1
      ! omega- equation
      alfaSST   = alfa_1*this%bF1(i,k) + alfa_2*(1.0-this%bF1(i,k))
      betaSST   = beta_1*this%bF1(i,k) + beta_2*(1.0-this%bF1(i,k))

      StR = (2.*(((w(i,k)-w(i,km))/dz)**2. + &
        ((u(i,k)-u(im,k))/dRu(i))**2. + &
        ((u(i,k)+u(im,k))/(2.*Rp(i)))**2.) + &
        (((w(ip,km)+w(ip,k)+w(i,km)+w(i,k))/4. -(w(im,km)+w(im,k)+w(i,km)+w(i,k))/4.)/dRu(i) &
        +((u(i,kp) +u(im,kp)+u(i,k)+u(im,k))/4.-(u(im,km)+u(i,km)+u(im,k)+u(i,k))/4.)/(dz)  )**2.)

      ! StR = (2.*(((w(i,k)-w(i,km))/dzw(k))**2. + &
      !   ((u(i,k)-u(im,k))/dRu(i))**2. + &
      !   ((u(i,k)+u(im,k))/(2.*Rp(i)))**2.) + &
      !   (((w(ip,km)+w(ip,k)+w(i,km)+w(i,k))/4. -(w(im,km)+w(im,k)+w(i,km)+w(i,k))/4.)/dRu(i) &
      !   +((u(i,kp) +u(im,kp)+u(i,k)+u(im,k))/4.-(u(im,km)+u(i,km)+u(im,k)+u(i,k))/4.)/dzw(k)  )**2.)

      div(i,k) =(Ru(i)*u(i,k)-Ru(im)*u(im,k))/(Rp(i)*dru(i)) &
        +(      w(i,k) -      w(i,km))/dz

      ! div(i,k) =(Ru(i)*u(i,k)-Ru(im)*u(im,k))/(Rp(i)*dru(i)) &
      !   +(      w(i,k) -      w(i,km))/dzw(k)

         ! Bouyancy prodution divided by mut
      GtR=-ctheta*beta(i,k)*Fr_1*this%Tt(i,k) &
        *  ((((w(ip,km)+w(ip,k)+w(i,km)+w(i,k))/4.-(w(im,km)+w(im,k)+w(i,km)+w(i,k))/4.)/dRu(i) &
        +    ((u(i,kp)+u(im,kp)+u(i,k)+u(im,k))/4.-(u(im,km)+u(i,km)+u(im,k)+u(i,k))/4.)/(dz) )* &
        (temp(ip,k)-temp(im,k))/(dRp(i)+dRp(im))  ) &
        +(2*((w(i,k)-w(i,km))/dz-2./3.*(rho(i,k)*putink(i,k)))*(temp(i,kp)-temp(i,km))/(2.*dz) &
        )

      ! GtR=-ctheta*beta(i,k)*Fr_1*this%Tt(i,k) &
      !   *  ((((w(ip,km)+w(ip,k)+w(i,km)+w(i,k))/4.-(w(im,km)+w(im,k)+w(i,km)+w(i,k))/4.)/dRu(i) &
      !   +    ((u(i,kp)+u(im,kp)+u(i,k)+u(im,k))/4.-(u(im,km)+u(i,km)+u(im,k)+u(i,k))/4.)/dzw(k) )* &
      !   (temp(ip,k)-temp(im,k))/(dRp(i)+dRp(im))  ) &
      !   +(2*((w(i,k)-w(i,km))/dzw(k)-2./3.*(rho(i,k)*putink(i,k)))*(temp(i,kp)-temp(i,km))/(2.*dzp(k)) &
      !   )


      GtR = GtR + ctheta*beta(i,k)*Fr_1*this%Tt(i,k)*2./3.*div(i,k)*(temp(i,kp)-temp(i,km))/(2.*dz)

      ! GtR = GtR + ctheta*beta(i,k)*Fr_1*this%Tt(i,k)*2./3.*div(i,k)*(temp(i,kp)-temp(i,km))/(2.*dzp(k))


      putout(i,k) = putout(i,k) + (alfaSST*StR*rho(i,k) + alfaSST*GtR*rho(i,k) + (1.0-this%bF1(i,k))*this%cdKOM(i,k) ) /rho(i,k)
      dimpl(i,k)  = dimpl(i,k)  + betaSST*this%om(i,k) ! note, beta*rho*omega^2/(rho*omega), set implicit and divided by density

    enddo
  enddo
end subroutine rhs_om_sst

subroutine diffusion_om_SST(this, putout,putin,ek,eki,ekk,ekmt,sigma,rho,Ru,Rp,dru,dz,modification)
  use module_mesh, only : mesh
  implicit none
  class(SST_TurbModel) :: this
  real(8), dimension(0:this%i1,0:this%k1), intent(IN) :: putin,ek,eki,ekk,ekmt,sigma,rho
  real(8), dimension(0:this%i1),           intent(IN) :: Rp,Ru,dRu
  real(8),                                 intent(IN) :: dz
  integer                                , intent(IN) :: modification
  real(8), dimension(0:this%i1,0:this%k1), intent(OUT):: putout
  integer  km,kp,i,k
  real(8), dimension(0:this%k1) :: dzw, dzp

  dzw = mesh%dzw
  dzp = mesh%dzp

  if ((modification == 1) .or. (modification == 2)) then       ! Inverse SLS & Aupoix
    do k=1,this%kmax
      kp=k+1
      km=k-1
      do i=1,this%imax
        putout(i,k) = putout(i,k) + 1.0/rho(i,k)*( &
          ( (ekk(i,k ) + (ekmt(i,k)+ekmt(i,kp))/(sigma(i,k)+sigma(i,kp)))/sqrt(0.5*(rho(i,k)+rho(i,kp)))* &
          (putin(i,kp)*sqrt(rho(i,kp)) - putin(i,k )*sqrt(rho(i,k ))) &
          - (ekk(i,km) + (ekmt(i,k)+ekmt(i,km))/(sigma(i,k)+sigma(i,km)))/sqrt(0.5*(rho(i,k)+rho(i,km)))* &
          (putin(i,k )*sqrt(rho(i,k )) - putin(i,km)*sqrt(rho(i,km))) &
          )/(dz*dz)   )

        ! putout(i,k) = putout(i,k) + 1.0/rho(i,k)*( &
        !   ( (ekk(i,k ) + (ekmt(i,k)+ekmt(i,kp))/(sigma(i,k)+sigma(i,kp)))/sqrt(0.5*(rho(i,k)+rho(i,kp)))* &
        !   (putin(i,kp)*sqrt(rho(i,kp)) - putin(i,k )*sqrt(rho(i,k )))/dzp(k) &
        !   - (ekk(i,km) + (ekmt(i,k)+ekmt(i,km))/(sigma(i,k)+sigma(i,km)))/sqrt(0.5*(rho(i,k)+rho(i,km)))* &
        !   (putin(i,k )*sqrt(rho(i,k )) - putin(i,km)*sqrt(rho(i,km)))/dzp(km) &
        !   )/dzw(k)   )
      enddo
    enddo
  else                                                        ! Standard
    do k=1,this%kmax
      kp=k+1
      km=k-1
      do i=1,this%imax
        putout(i,k) = putout(i,k) + 1.0/rho(i,k)*( &
          ( (ekk(i,k ) + (ekmt(i,k)+ekmt(i,kp))/(sigma(i,k)+sigma(i,kp)))*(putin(i,kp)-putin(i,k )) &
           -(ekk(i,km) + (ekmt(i,k)+ekmt(i,km))/(sigma(i,k)+sigma(i,km)))*(putin(i,k )-putin(i,km)) &
          )/(dz*dz))
        ! putout(i,k) = putout(i,k) + 1.0/rho(i,k)*( &
        !   ( (ekk(i,k ) + (ekmt(i,k)+ekmt(i,kp))/(sigma(i,k)+sigma(i,kp)))*(putin(i,kp)-putin(i,k ))/dzp(k) &
        !    -(ekk(i,km) + (ekmt(i,k)+ekmt(i,km))/(sigma(i,k)+sigma(i,km)))*(putin(i,k )-putin(i,km))/dzp(km) &  
        !   )/dzw(k))

      enddo
    enddo
  endif
end subroutine diffusion_om_SST



end module sst_tm
