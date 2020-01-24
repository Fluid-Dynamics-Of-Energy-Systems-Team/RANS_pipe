module sst_tm
  use ke_tm
  implicit none
  
!****************************************************************************************

  !************************!
  !         SST class      !
  !************************!

  type, extends(KE_TurbModel), public :: SST_TurbModel
  real(8), dimension(:,:), allocatable :: cdKOM
  real(8), dimension(:),   allocatable :: omin,bF1in
  contains
    procedure :: init => init_SST
    procedure :: init_sol => init_sol_SST
    procedure :: set_mut => set_mut_SST
    procedure :: set_bc => set_bc_SST
    procedure :: get_profile => get_profile_SST
    procedure :: get_sol => get_sol_SST
    procedure :: init_w_inflow => init_w_inflow_SST
    procedure :: solve_k_KE => solve_k_SST
    procedure :: solve_eps_KE => solve_om_SST
    procedure :: calc_turbulent_timescale => calc_turbulent_timescale_SST
    ! procedure :: production_KE => production_SST
    procedure :: init_mem_SST
    procedure :: diffusion_k_SST
    procedure :: diffusion_om_SST
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


type(SST_TurbModel) function init_SST_TurbModel(name)
  character(len=3), intent(IN) :: name
  init_SST_TurbModel%name=name
end function init_SST_TurbModel

subroutine init_SST(this)
    implicit none
    class(SST_TurbModel) :: this
    this%name='SST'
    call this%init_mem_SST()
    call this%init_sol()
end subroutine init_SST

subroutine init_sol_SST(this)
  use mod_param, only : i1,i
  implicit none
  class(SST_TurbModel) :: this

  do i=0,i1
    this%k(i,:)  = 0.1!!0.1
    this%om(i,:)  = 1.0
    this%bF1(i,:)  = 1.0
    this%kin = 0.1
    this%omin = 1.0
    this%bF1in = 1.0
    this%mutin = this%kin/this%omin
  enddo
end subroutine init_sol_SST

subroutine init_mem_SST(this)
  use mod_param, only : k1,i1,kmax,imax  
  implicit none
  class(SST_TurbModel) :: this
  allocate(this%om (0:i1,0:k1),this%k(0:i1,0:k1),   this%bF1(0:i1,0:k1),this%bF2(imax,kmax),    &
           this%Gk (0:i1,0:k1),this%Pk (0:i1,0:k1), this%Tt (0:i1,0:k1),this%cdKOM(imax,kmax),  &
           this%yp (0:i1,0:k1), this%div(0:i1,0:k1))
  allocate(this%mutin(0:i1),this%Pkin (0:i1), this%bF1in(0:i1),this%omin (0:i1),this%kin(0:i1))
end subroutine init_mem_SST

subroutine init_w_inflow_SST(this,nuSAin,pkin,kin,epsin,omin,mutin,v2in)
  use mod_param, only : i1,k1,k
  class(SST_TurbModel) :: this
  real(8), dimension(0:i1), intent(IN) :: nuSAin,pkin,kin,epsin,omin,mutin,v2in
    this%omin = omin
    this%kin = kin
    do k=0,k1
      this%om(:,k) = this%omin(:)
      this%k(:,k) = this%kin(:)
    enddo
end subroutine init_w_inflow_SST

subroutine set_mut_SST(this,u,w,rho,mu,mui,mut)
  use mod_param, only : k1,i1,kmax,imax,k,i
  use mod_mesh, only : mesh,dzw,dzp,dru,drp,rp,ru,walldist,ru
  implicit none
  class(SST_TurbModel) :: this
  real(8), dimension(0:i1,0:k1), intent(IN) :: u, w, rho, mu, mui
  real(8), dimension(0:i1,0:k1), intent(OUT):: mut
  real(8), dimension(k1) :: tauw
  integer  im,ip,km,kp
  real(8)  sigma_om2,betaStar,gradkom,gamma1,gamma2,gamma3,gammaSST,zetaSST,StR, wallD

  !constants
  sigma_om2 = 0.856
  betaStar  = 0.09
      
  do k=1,kmax
    km=k-1
    kp=k+1
    tauw(k) = mui(imax,k)*0.5*(w(imax,km)+w(imax,k))/walldist(imax)
  
    do i=1,imax
      im=i-1
      ip=i+1
      this%yp(i,k)     = sqrt(rho(i,k))/mu(i,k)*(walldist(i))*tauw(k)**0.5   ! ystar
      wallD     = walldist(i)
      ! Vorticity rate
      StR = ( & 
            ( -( (w(ip,km)+w(ip,k)+w(i,km)+w(i ,k))/4.-(w(im,km)+w(im,k)+w(i,km)+w(i,k))/4. )/dru(i) &
          + (    (u(i,kp) +u(im,kp)+u(i,k)+u(im,k))/4.-(u(im,km)+u(i,km)+u(im,k)+u(i,k))/4. )/dzw(k) &
            )**2.)

      StR = StR**0.5
               
      gradkom = ((this%k(ip,k) - this%k(im,k))/(dRp(i)+dRp(im))) *((this%om(ip,k) - this%om(im,k))/(dRp(i)+dRp(im))) &
               +((this%k(i,kp) - this%k(i,km))/(dzp(k)+dzp(km))) *((this%om(i,kp) - this%om(i,km))/(dzp(k)+dzp(km)))

      this%cdKOM(i,k) = 2.0*sigma_om2*rho(i,k)/this%om(i,k)*gradkom;
              
      gamma1   = 500.0*mu(i,k)/(rho(i,k)*this%om(i,k)*wallD**2.0);
      gamma2   = 4.0*sigma_om2*rho(i,k)*this%k(i,k)/(wallD*wallD*max(this%cdKOM(i,k), 1.0e-20));
      gamma3   = (this%k(i,k)**0.5)/(betaStar*this%om(i,k)*wallD)

      gammaSST = min(max(gamma1, gamma3), gamma2)
      this%bF1(i,k) = tanh(gammaSST**4.0)

      gammaSST = max(2.0*gamma3, gamma1)
      this%bF2(i,k) = tanh(gammaSST**2.0)

      ! zetaSST  = max(0.31*this%om(i,k), this%bF2(i,k)*StR)
      zetaSST  = min(1.0/this%om(i,k), 0.31/(this%bF2(i,k)*StR))

      mut(i,k) = rho(i,k)*this%k(i,k)*zetaSST !!! NOTE this is the correct one !!!!
      ! mut(i,k) = rho(i,k)*this%k(i,k)/this%om(i,k)
      mut(i,k) = min(max(mut(i,k),0.0),100.0);

    enddo
  enddo
end subroutine set_mut_SST

subroutine set_bc_SST(this,mu,rho,periodic,rank,px)
  use mod_param,only : k1,i1,kmax,imax,k
  use mod_mesh, only : walldist,top_bcnovalue,bot_bcnovalue,top_bcvalue,bot_bcvalue
  implicit none
  class(SST_TurbModel) :: this
  real(8),dimension(0:i1,0:k1),intent(IN) :: rho,mu
  integer,                               intent(IN) :: periodic, rank, px
  real(8),dimension(0:k1) :: BCvalue
  real(8),dimension(0:i1) :: tmp
  real(8) :: botBCvalue, topBCvalue
  
  do k = 0,k1
    this%k(0,k)   = bot_bcnovalue(k)*this%k(1,k)      !symmetry or 0 value
    this%k(i1,k)  = top_bcnovalue(k)*this%k(imax,k)   !symmetry or 0 value
    this%bF1(0,k) = bot_bcnovalue(k)*this%bF1(1,k)
    this%bF1(i1,k)= top_bcnovalue(k)*this%bF1(imax,k)
    botBCvalue    = (60.0/0.075)*mu(1,k)/rho(1,k)/walldist(1)**2                                            !bcvalue bot
    this%om(0,k)  = (1.-bot_bcvalue(k))*(2.0*botBCvalue-this%om(1,k))    + bot_bcvalue(k)*this%om(1,k)    !symmetry or bc value
    topBCvalue    = (60.0/0.075)*mu(imax,k)/rho(imax,k)/walldist(imax)**2                                   !bcvalue top
    this%om(i1,k) = (1.-top_bcvalue(k))*(2.0*topBCvalue-this%om(imax,k)) + top_bcvalue(k)*this%om(imax,k) !symmetry or bc value
  enddo

  call shiftf(this%k,  tmp,rank); this%k  (:,0)      =tmp(:);
  call shiftf(this%om ,tmp,rank); this%om (:,0)      =tmp(:);
  call shiftf(this%bF1,tmp,rank); this%bF1(:,0)      =tmp(:)
  call shiftb(this%k  ,tmp,rank); this%k  (:,k1)=tmp(:);
  call shiftb(this%om ,tmp,rank); this%om (:,k1)=tmp(:);
  call shiftb(this%bF1,tmp,rank); this%bF1(:,k1)=tmp(:)

  !developing
  if (periodic.eq.1) return
  if (rank.eq.0) then
    this%k  (:,0) = this%kin(:)
    this%om (:,0) = this%omin(:)
    !this%bF1(:,0) = this%bF1(:,1)  ! ATTENTION (THIS WAS THE ORIGINAL)
    this%bF1(:,0) = this%bF1in(:) 
  endif

  if (rank.eq.px-1) then
    this%k  (:,k1) = 2.0*this%k  (:,kmax)-this%k  (:,kmax-1)
    this%om (:,k1) = 2.0*this%om (:,kmax)-this%om (:,kmax-1)
    this%bF1(:,k1) = 2.0*this%bF1(:,kmax)-this%bF1(:,kmax-1)
  endif
  
end subroutine set_bc_SST

subroutine get_profile_SST(this,p_nuSA,p_k,p_eps,p_om,p_v2,p_Pk,p_bF1,p_bF2,yp,k)
  use mod_param, only : i1,imax  
  class(SST_TurbModel) :: this
  integer,                               intent(IN) :: k
  real(8),dimension(0:i1),          intent(OUT):: p_nuSA,p_k,p_eps,p_om,p_v2,p_Pk,p_bF1,yp
  real(8),dimension(1:imax),        intent(OUT):: p_bF2
  p_nuSA(:)=0
  p_k(:)   =this%k(:,k)
  p_eps(:) =0
  p_om(:)  =this%om(:,k)
  p_Pk(:)  =this%Pk(:,k)
  p_bF1(:) =this%bF1(:,k)
  p_bF2(:) =this%bF2(:,k)
  p_v2(:)  =0
  yp(:)    =this%yp(:,k)
end subroutine get_profile_SST

subroutine get_sol_SST(this,nuSA,k,eps,om,v2,pk, gk,yp)
  use mod_param, only : k1,i1
  class(SST_TurbModel) :: this
  real(8),dimension(0:i1,0:k1), intent(OUT):: nuSA,k,eps,om,v2,yp,pk,gk
  nuSA=0
  k   =this%k    
  eps =0
  v2  =0
  om  =this%om
  pk  = this%pk
  gk = this%gk
  yp  =this%yp
end subroutine get_sol_SST

subroutine solve_k_SST(this,resK,u,w,rho,mu,mui,muk,mut,rho_mod, &
                       alphak,modification,rank,periodic)
  use mod_param, only : k1,i1,kmax,imax,i,k
  use mod_mesh,  only : drp,dru,ru,rp,drp,top_bcnovalue,bot_bcnovalue
  use mod_math
  implicit none
  class(SST_TurbModel) :: this
  real(8),dimension(0:i1,0:k1), intent(IN) :: u, w, rho,mu,mui,muk,mut,rho_mod
  real(8),                      intent(IN) :: alphak
  integer,                      intent(IN) :: modification,rank,periodic
  real(8),                      intent(OUT):: resK
  real(8), dimension(0:i1,0:k1) :: dnew,dimpl,sigmakSST
  real(8), dimension(imax)      :: a,b,c,rhs
  real(8) :: dz
  
  resK  = 0.0; dnew=0.0; dimpl = 0.0;

  call advecc(dnew,dimpl,this%k,u,w,rank,periodic,.true.)
  call this%rhs_k_SST(dnew,dimpl,rho)

  ! calculating constant with blending function factor
  sigmakSST = 0.85*this%bF1 + 1.0*(1.0 - this%bF1)
  sigmakSST = 1.0/sigmakSST
  
  call this%diffusion_k_sst(dnew,this%k,mu,mui,muk,mut,sigmakSST,rho,modification)

  do k=1,kmax
    do i=1,imax
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

    i=imax
    b(i) = b(i)+top_bcnovalue(k)*c(i)
    rhs(i) = dnew(i,k)  + ((1-alphak)/alphak)*b(i)*this%k(i,k)

    call matrixIdir(imax,a,b/alphak,c,rhs)

    do i=1,imax
      resK = resK + ((this%k(i,k) - rhs(i))/(this%k(i,k)+1.0e-20))**2.0
      this%k(i,k) = max(rhs(i), 1.0e-8)
    enddo
  enddo
end subroutine solve_k_SST

subroutine solve_om_sst(this,rese,u,w,rho,mu,mui,muk,mut,beta,temp,rho_mod, &
                        alphae,modification,rank,periodic)
  use mod_param, only : k1,i1,kmax,imax,i,k
  use mod_math
  use mod_mesh, only : top_bcvalue,bot_bcvalue,ru,rp,dru,drp
  implicit none
  class(SST_TurbModel) :: this
  real(8), dimension(0:i1,0:k1), intent(IN) :: u, w, rho,mu,mui,muk,mut,beta,temp,rho_mod
  real(8),                       intent(IN) :: alphae
  integer,                       intent(IN) :: modification,rank,periodic
  real(8),                       intent(OUT):: rese
  real(8), dimension(0:i1,0:k1) :: dnew,dimpl,sigmakSST
  real(8), dimension(imax)      :: a,b,c,rhs

  resE = 0.0
  dnew=0.0; dimpl = 0.0;
  call advecc(dnew,dimpl,this%om,u,w,rank,periodic,.true.)
  call this%rhs_om_sst(dnew,dimpl,this%k,u,w,temp,rho,beta,mut)


  ! calculating constant with blending function factor
  sigmakSST = 0.5*this%bF1 + 0.856*(1.0 - this%bF1)
  sigmakSST = 1.0/sigmakSST
  call this%diffusion_om_SST(dnew,this%om,mu,mui,muk,mut,sigmakSST,rho,modification)
  do k=1,kmax
    do i=1,imax
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

    i = imax
    b(i)=b(i)+top_bcvalue(k)*c(i)
    rhs(i) = dnew(i,k) - (1-top_bcvalue(k))*c(i)*this%om(i+1,k) + ((1-alphae)/alphae)*b(i)*this%om(i,k)  !wall
  
    call matrixIdir(imax,a,b/alphae,c,rhs)
  
    do i=1,imax
      resE = resE + ((this%om(i,k) - rhs(i))/(this%om(i,k)+1.0e-20))**2.0
      this%om(i,k) = max(rhs(i), 1.0e-8)
    enddo
  enddo
end subroutine solve_om_sst

subroutine rhs_k_SST(this,putout,dimpl,rho)
  use mod_param, only : k1,i1,kmax,imax,i,k
  implicit none
  class(SST_TurbModel) :: this
  real(8), dimension(0:i1,0:k1), intent(IN) :: rho
  real(8), dimension(0:i1,0:k1), intent(OUT):: putout, dimpl
  real(8)  betaStar

  betaStar  = 0.09

  do k=1,kmax
    do i=1,imax
      putout(i,k) = putout(i,k) + ( this%Pk(i,k) + this%Gk(i,k) )/rho(i,k)
      dimpl(i,k)  = dimpl(i,k)  + betaStar*this%om(i,k)            ! note, betaStar*rho*k*omega/(rho*k), set implicit and divided by density
    enddo
  enddo
end subroutine rhs_k_SST

subroutine diffusion_k_sst(this,putout,putin,ek,eki,ekk,ekmt,sigma,rho,modification)
  use mod_param, only : k1,i1,kmax,imax,i,k
  use mod_mesh,  only : dzw,dzp
  implicit none
  class(SST_TurbModel) :: this
  real(8), dimension(0:i1,0:k1), intent(IN) :: putin,ek,eki,ekk,ekmt,sigma,rho
  integer                      , intent(IN) :: modification
  real(8), dimension(0:i1,0:k1), intent(OUT):: putout
  integer   km,kp

  if (modification == 1) then       ! Inverse SLS
    do k=1,kmax
      kp=k+1
      km=k-1
      do i=1,imax
        putout(i,k) = putout(i,k) + 1.0/rho(i,k)/sqrt(rho(i,k))*( &
          ( (ekk(i,k ) + (ekmt(i,k)+ekmt(i,kp))/(sigma(i,k)+sigma(i,kp)))/ &
          sqrt(0.5*(rho(i,k)+rho(i,kp)))*(rho(i,kp)*putin(i,kp)-rho(i,k )*putin(i,k ))/dzp(k) &
          -(ekk(i,km) + (ekmt(i,k)+ekmt(i,km))/(sigma(i,k)+sigma(i,km)))/ &
          sqrt(0.5*(rho(i,k)+rho(i,km)))*(rho(i,k )*putin(i,k )-rho(i,km)*putin(i,km))/dzp(km) &
          )/dzw(k))
      enddo
    enddo
  elseif (modification == 2) then   ! Aupoix
    do k=1,kmax
      kp=k+1
      km=k-1
      do i=1,imax
        putout(i,k) = putout(i,k) + 1.0/rho(i,k)*( &
          ( (ekk(i,k ) + (ekmt(i,k)+ekmt(i,kp))/(sigma(i,k)+sigma(i,kp)))/ &
          (0.5*(rho(i,k)+rho(i,kp)))*(rho(i,kp)*putin(i,kp)-rho(i,k )*putin(i,k ))/dzp(k) &
          -(ekk(i,km) + (ekmt(i,k)+ekmt(i,km))/(sigma(i,k)+sigma(i,km)))/ &
          (0.5*(rho(i,k)+rho(i,km)))*(rho(i,k )*putin(i,k )-rho(i,km)*putin(i,km))/dzp(km) &
          )/dzw(k))
      enddo
    enddo
  else                               ! Standard
    do k=1,kmax
      kp=k+1
      km=k-1
      do i=1,imax
        putout(i,k) = putout(i,k) + 1.0/rho(i,k)*( &
          ( (ekk(i,k ) + (ekmt(i,k)+ekmt(i,kp))/(sigma(i,k)+sigma(i,kp)))*(putin(i,kp)-putin(i,k ))/dzp(k ) &
           -(ekk(i,km) + (ekmt(i,k)+ekmt(i,km))/(sigma(i,k)+sigma(i,km)))*(putin(i,k )-putin(i,km))/dzp(km) &
          )/dzw(k))
      enddo
    enddo
  endif
end subroutine diffusion_k_sst

subroutine rhs_om_sst(this,putout,dimpl,putink,u,w,temp,rho,beta,mut)
  use mod_param, only : k1,i1,kmax,imax,i,k
  use mod_mesh,  only : dzw,dzp,rp,ru,dru,drp
  implicit none
  class(SST_TurbModel) :: this
  real(8), dimension(0:i1,0:k1), intent(IN) :: u,w,temp,rho,beta, putink, mut
  real(8), dimension(0:i1,0:k1), intent(OUT):: putout, dimpl
  integer km,kp,im,ip
  real(8), dimension(0:i1,0:k1) :: div
  real(8) sigma_om1,sigma_om2,beta_1,beta_2,betaStar,alfa_1,alfa_2,alfaSST,betaSST,StR,GtR,ctheta,Fr_1

  
  sigma_om1 = 0.5
  sigma_om2 = 0.856
  beta_1    = 0.075
  beta_2    = 0.0828
  betaStar  = 0.09

  alfa_1    = beta_1/betaStar - sigma_om1*(0.41**2.0)/(betaStar**0.5)
  alfa_2    = beta_2/betaStar - sigma_om2*(0.41**2.0)/(betaStar**0.5)

  do k=1,kmax
    kp=k+1
    km=k-1
    do i=1,imax
      ip=i+1
      im=i-1
      ! omega- equation
      alfaSST   = alfa_1*this%bF1(i,k) + alfa_2*(1.0-this%bF1(i,k))
      betaSST   = beta_1*this%bF1(i,k) + beta_2*(1.0-this%bF1(i,k))

      putout(i,k) = putout(i,k) + ((alfaSST*rho(i,k)/mut(i,k))*(this%Pk(i,k) +this%Gk(i,k))  &
                                    + (1.0-this%bF1(i,k))*this%cdKOM(i,k) )/rho(i,k)
      dimpl(i,k)  = dimpl(i,k)  + betaSST*this%om(i,k) ! note, beta*rho*omega^2/(rho*omega), set implicit and divided by density
    enddo
  enddo
end subroutine rhs_om_sst

subroutine diffusion_om_SST(this, putout,putin,ek,eki,ekk,ekmt,sigma,rho,modification)
  use mod_param, only : k1,i1,kmax,imax,i,k
  use mod_mesh,  only : dzw,dzp,Ru,Rp,dru
  implicit none
  class(SST_TurbModel) :: this
  real(8), dimension(0:i1,0:k1), intent(IN) :: putin,ek,eki,ekk,ekmt,sigma,rho
  integer                      , intent(IN) :: modification
  real(8), dimension(0:i1,0:k1), intent(OUT):: putout
  integer :: km,kp

  if ((modification == 1) .or. (modification == 2)) then       ! Inverse SLS & Aupoix
    do k=1,kmax
      kp=k+1
      km=k-1
      do i=1,imax
        putout(i,k) = putout(i,k) + 1.0/rho(i,k)*( &
          ( (ekk(i,k ) + (ekmt(i,k)+ekmt(i,kp))/(sigma(i,k)+sigma(i,kp)))/sqrt(0.5*(rho(i,k)+rho(i,kp)))* &
          (putin(i,kp)*sqrt(rho(i,kp)) - putin(i,k )*sqrt(rho(i,k )))/dzp(k) &
          - (ekk(i,km) + (ekmt(i,k)+ekmt(i,km))/(sigma(i,k)+sigma(i,km)))/sqrt(0.5*(rho(i,k)+rho(i,km)))* &
          (putin(i,k )*sqrt(rho(i,k )) - putin(i,km)*sqrt(rho(i,km)))/dzp(km) &
          )/dzw(k)   )
      enddo
    enddo
  else                                                        ! Standard
    do k=1,kmax
      kp=k+1
      km=k-1
      do i=1,imax
        putout(i,k) = putout(i,k) + 1.0/rho(i,k)*( &
          ( (ekk(i,k ) + (ekmt(i,k)+ekmt(i,kp))/(sigma(i,k)+sigma(i,kp)))*(putin(i,kp)-putin(i,k ))/dzp(k) &
           -(ekk(i,km) + (ekmt(i,k)+ekmt(i,km))/(sigma(i,k)+sigma(i,km)))*(putin(i,k )-putin(i,km))/dzp(km) &  
          )/dzw(k))
      enddo
    enddo
  endif
end subroutine diffusion_om_SST

subroutine calc_turbulent_timescale_SST(this,rho,mu)
  use mod_param, only : k1,i1
  implicit none
  class(SST_TurbModel) :: this
  real(8), dimension(0:i1,0:k1), intent(IN) :: rho, mu
  this%Tt = 1/this%om
end subroutine

end module sst_tm
