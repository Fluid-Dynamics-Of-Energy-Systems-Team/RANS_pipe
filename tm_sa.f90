module sa_tm
  use mod_tm
  implicit none

!****************************************************************************************

  !************************!
  !         SA class       !
  !************************!

  type, extends(TurbModel), public :: SA_TurbModel
  real(8), dimension(:), allocatable :: nuSAin
  contains
    procedure :: init => init_SA
    procedure :: init_sol => init_sol_SA
    procedure :: set_mut => set_mut_SA
    procedure :: advance_turb => advance_SA
    procedure :: get_profile => get_profile_SA
    procedure :: get_sol => get_sol_SA 
    procedure :: set_bc => set_bc_SA
    procedure :: init_w_inflow => init_w_inflow_SA
    procedure :: init_mem_SA
    procedure :: solve_SA
    procedure :: production_SA
    procedure :: diffusion_SA
    procedure :: rhs_SA


  end type SA_TurbModel

contains

!****************************************************************************************

  !************************!
  !      SA routines      !
  !************************!

subroutine init_SA(this)
    implicit none
    class(SA_TurbModel) :: this  
    this%name='SA'
    call this%init_mem_SA()
    call this%init_sol()
end subroutine init_SA

subroutine init_sol_SA(this)
    use mod_param, only : i1
    implicit none
    class(SA_TurbModel) :: this
    integer i
    do i=0,i1
      this%nuSA(i,:) = 0.0
      this%nuSAin(i) = 0.0
      this%pkin(i)   = 0.0
      this%mutin(i)  = 0.0
    enddo
end subroutine init_sol_SA

subroutine init_w_inflow_SA(this,nuSAin,pkin,kin,epsin,omin,mutin,v2in)
  use mod_param, only : i1,k1,k
  class(SA_TurbModel) :: this
  real(8), dimension(0:i1), intent(IN) :: nuSAin,pkin,kin,epsin,omin,mutin,v2in
  this%nuSAin = nuSAin
  this%pkin   = pkin
  do k=0,k1
    this%nuSA(:,k) = this%nuSAin(:)
    this%pk(:,k)   = this%pkin(:)
  enddo
end subroutine init_w_inflow_SA

subroutine init_mem_SA(this)
  use mod_param, only : i1,k1
  implicit none
  class(SA_TurbModel) :: this
  allocate(this%nuSA(0:i1,0:k1),this%Pk(0:i1,0:k1),this%yp(0:i1,0:k1))
  allocate(this%mutin(0:i1),this%nuSAin(0:i1),this%Pkin(0:i1))
end subroutine init_mem_SA

subroutine set_mut_SA(this,u,w,rho,mu,mui,mut)
  use mod_param, only : i1,k1,imax,kmax
  use mod_mesh, only : walldist
  implicit none
  class(SA_TurbModel) :: this
  real(8), dimension(0:i1,0:k1), intent(IN) :: u, w, rho, mu, mui
  real(8), dimension(0:i1,0:k1), intent(OUT):: mut
  real(8), dimension(k1) :: tauw
  integer  im,ip,km,kp,i,k
  real*8   cv1_3,chi,fv1SA

  cv1_3 = (7.1)**3.0
  do k=1,kmax
    km=k-1
    kp=k+1
    tauw(k) = mui(imax,k)*0.5*(w(imax,km)+w(imax,k))/walldist(imax)
    do i=1,imax
      im=i-1
      ip=i+1
      this%yp(i,k) = sqrt(rho(i,k))/mu(i,k)*(walldist(i))*tauw(k)**0.5 !ystar
      chi     = this%nuSA(i,k)/(mu(i,k)/rho(i,k))
      fv1SA   = (chi**3.0)/(chi**3.0 + cv1_3);
      mut(i,k)= min(100.0,max(0.0,rho(i,k)*this%nuSA(i,k)*fv1SA))
    enddo
  enddo
end subroutine set_mut_SA

subroutine set_bc_SA(this,mu,rho,periodic,rank,px)
  use mod_param, only : k1,i1,kmax,imax,k
  use mod_mesh,  only : top_bcnovalue,bot_bcnovalue
  use mod_common, only : wnew
  implicit none
  class(SA_TurbModel) :: this
  real(8),dimension(0:i1,0:k1),intent(IN) :: rho,mu
  integer,                     intent(IN) :: periodic, rank, px
  real(8),dimension(0:i1)                 :: tmp
  
  do k = 0,kmax 
    this%nuSA(0,k)  = bot_bcnovalue(k)*this%nuSA(1,k)    !symmetry or 0 value
    this%nuSA(i1,k) = top_bcnovalue(k)*this%nuSA(imax,k) !symmetry or 0 value
  enddo  

  call shiftf(this%nuSA,tmp,rank); this%nuSA(:,0)  =tmp(:);
  call shiftb(this%nuSA,tmp,rank); this%nuSA(:,k1) =tmp(:);

  ! developing
  if (periodic.eq.1) return
  if (rank.eq.0) then
    this%nuSA(:,0) = this%nuSAin(:)
    this%Pk(:,0)   = this%pkin(:)
  endif
  if (rank.eq.px-1) then
    this%nuSA(:,k1) = 2.0*this%nuSA(:,kmax) - this%nuSA(:,kmax-1)
    this%pk(:,k1)   = 2.0*this%pk  (:,kmax) - this%pk  (:,kmax-1)
  endif

end subroutine set_bc_SA

subroutine get_profile_SA(this,p_nuSA,p_k,p_eps,p_om,p_v2,p_Pk,p_bF1,p_bF2,yp,k)
  use mod_param, only : i1,imax
  implicit  none
  class(SA_TurbModel) :: this
  integer,                          intent(IN) :: k
  real(8),dimension(0:i1),          intent(OUT):: p_nuSA,p_k,p_eps,p_om,p_v2,p_Pk,p_bF1,yp
  real(8),dimension(1:imax),        intent(OUT):: p_bF2
  p_nuSA(:)=this%nuSA(:,k)
  p_k(:)   =0
  p_eps(:) =0
  p_v2(:)  =0
  p_om(:)  =0
  p_Pk(:)  =this%Pk(:,k)
  p_bF1(:) =0
  p_bF2(:) =0
  yp(:)  = this%yp(:,k)
end subroutine get_profile_SA

subroutine get_sol_SA(this,nuSA,k,eps,om,v2,yp)
  use mod_param, only : i1,k1
  class(SA_TurbModel) :: this
  real(8),dimension(0:i1,0:k1), intent(OUT):: nuSA,k,eps,om,v2,yp
  nuSA=this%nuSA
  k   =0    
  eps =0
  v2  =0
  om  =0
  yp  = this%yp
end subroutine get_sol_SA

subroutine advance_SA(this,u,w,rho,mu,mui,muk,mut,beta,temp, &
                      alpha1,alpha2,alpha3,                  &
                      modification,rank,periodic,   &
                      residual1, residual2, residual3)
  use mod_param, only : i1,k1
  implicit none
  class(SA_TurbModel) :: this
  real(8),dimension(0:i1,0:k1), intent(IN) :: u,w,rho,mu,mui,muk,mut,beta,temp
  real(8),                      intent(IN) :: alpha1,alpha2,alpha3
  integer,                      intent(IN) :: modification,rank,periodic
  real(8),                      intent(OUT):: residual1, residual2,residual3
  real(8),dimension(0:i1,0:k1)             :: rho_mod

  !1, our modification, 2, Aupoix modification
  if ((modification == 1) .or. (modification == 2)) then
    rho_mod = rho
  else
    rho_mod = 1.0
  endif
  call this%production_SA(this%nuSA,u,w,rho,mu)
  call this%solve_SA(residual1,u,w,rho,mu,mui,muk,rho_mod,alpha1,modification,periodic,rank)
end subroutine advance_SA

subroutine solve_SA(this,resSA,u,w,rho,mu,mui,muk,rho_mod, &
                    alphak,modification,periodic,rank)
  use mod_param, only : i1,k1,imax,kmax,k,i
  use mod_mesh,  only : Ru,Rp,dru,drp,walldist,top_bcnovalue,bot_bcnovalue
  use mod_math
  use mod_common, only : res_nuSA
  implicit none
  class(SA_TurbModel) :: this
  real(8),dimension(0:i1,0:k1), intent(IN) :: u, w, rho,mu,mui,muk,rho_mod
  real(8),                      intent(IN) :: alphak
  integer,                      intent(IN) :: rank, modification, periodic
  real(8),                      intent(OUT):: resSA
  real(8) :: cb3
  real(8), dimension(0:i1,0:k1) :: dnew,tempArray,dimpl,eknu,eknui,eknuk
  real(8), dimension(imax)      :: a,b,c,rhs
  
  cb3 = 2.0/3.0
  resSA = 0.0
  dnew=0.0; dimpl = 0.0;
  call advecc(dnew,dimpl,this%nuSA,u,w,rank,periodic,.true.)
  call this%rhs_SA(dnew,dimpl,this%nuSA,rho,modification)

  do k=0,kmax+1
    do i=0,imax+1
      eknu (i,k) = mu(i,k)/rho(i,k)
      if (i.LE.imax) eknui(i,k) = mui(i,k)/(0.5*(rho(i,k)+rho(i+1,k)))
      if (k.LE.kmax) eknuk(i,k) = muk(i,k)/(0.5*(rho(i,k)+rho(i,k+1)))
    enddo
  enddo

  tempArray = 0.0
  call this%diffusion_SA(tempArray,this%nuSA,eknu,eknui,eknuk,1.0,rho_mod,modification)
  dnew = dnew + tempArray/cb3

  !> diffusion term in the r-direction, set implicit!
  !! conv.    d/dy (nu+nuSA) dnuSA/dy
  !! modified (1/rho) d/dy [sqrt(rho) (nu+nuSA) d(sqrt(rho)nuSA)/dy]
  do k=1,kmax
    do i=1,imax

      a(i) = (eknui(i-1,k)+0.5*(this%nuSA(i,k)+this%nuSA(i-1,k)))/cb3*((0.5*(rho_mod(i-1,k)+rho_mod(i,k)))**0.5)
      a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/rho_mod(i,k)
      c(i) = (eknui(i  ,k)+0.5*(this%nuSA(i,k)+this%nuSA(i+1,k)))/cb3*((0.5*(rho_mod(i+1,k)+rho_mod(i,k)))**0.5)
      c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/rho_mod(i,k)

      b(i) = ((-a(i)-c(i))*(rho_mod(i,k)**0.5)+dimpl(i,k))!/alphak

      a(i) = a(i)*(rho_mod(i-1,k)**0.5)
      c(i) = c(i)*(rho_mod(i+1,k)**0.5)

      rhs(i) = dnew(i,k) + ((1-alphak)/alphak)*b(i)*this%nuSA(i,k)
    enddo

    i=1
    b(i) = b(i)+bot_bcnovalue(k)*a(i)
    rhs(i) = dnew(i,k)+((1-alphak)/alphak)*b(i)*this%nuSA(i,k)
    
    i=imax
    b(i) = b(i)+top_bcnovalue(k)*c(i)
    rhs(i) = dnew(i,k)+((1-alphak)/alphak)*b(i)*this%nuSA(i,k)
  
    call matrixIdir(imax,a,b/alphak,c,rhs)
    do i=1,imax
      resSA = resSA + ((this%nuSA(i,k) - rhs(i))/(this%nuSA(i,k)+1.0e-20))**2.0
      res_nuSA(i,k) = ((this%nuSA(i,k) - rhs(i))/(this%nuSA(i,k)+1.0e-20))**2.0
      this%nuSA(i,k) = max(rhs(i), 1.0e-8)
    enddo
  enddo
end subroutine solve_SA

subroutine production_SA(this,nuSA,u,w,rho,mu)
  use mod_param, only : i1,k1,imax,kmax,i,k
  use mod_mesh,  only : dzw,dzp,dru,walldist
  implicit none
  class(SA_TurbModel) :: this
  real(8),dimension(0:i1,0:k1), intent(IN) :: u,w,rho,mu,nuSA
  real(8),dimension(0:i1,0:k1) :: Gk, Tt
  integer :: im,ip,km,kp
  real*8  :: cv1_3,cb1,cb2,cb3,cw1,cw2,cw3_6,inv_cb3,kappa_2,chi,fv1SA,fv2SA,shatSA,StR

  cv1_3     = (7.1)**3.0
  cb1       = 0.1355
  cb2       = 0.622
  cb3       = 2.0/3.0
  inv_cb3   = 1.0/cb3
  kappa_2   = (0.41)**2.0   ! von karman constant
  cw1       = (cb1/kappa_2)+(1.0+cb2)/cb3
  cw2       = 0.3
  cw3_6     = (2.0)**6.0

  do k=1,kmax
    kp=k+1
    km=k-1
    do i=1,imax
      ip=i+1
      im=i-1
      ! Bouyancy prodution and time scale, not defined for this model
      Gk(i,k)=0
      Tt(i,k)=1
      ! magnitude of rate of rotation: omega=sqrt(2*Wij*Wij), Wrz = 0.5*(dU/dz-dW/dr);  note, utheta=0 d/dtheta=0
      StR = ( ( -( (w(ip,km)+w(ip,k)+w(i,km)+w(i ,k))/4.-(w(im,km)+w(im,k)+w(i,km)+w(i,k))/4.)/dRu(i) &
           +(      (u(i,kp) +u(im,kp)+u(i,k)+u(im,k))/4.-(u(im,km)+u(i,km)+u(im,k)+u(i,k))/4.)/dzw(k) )**2.)
      StR = StR**0.5
      ! calculating shat from SA model
      chi         = this%nuSA(i,k)/(mu(i,k)/rho(i,k))
      fv1SA       = (chi**3.0)/((chi**3.0) + cv1_3)
      fv2SA       = 1.0 - (chi/(1.0 + (chi*fv1SA)))
      shatSA      = StR + fv2SA*this%nuSA(i,k)/(kappa_2*(walldist(i)**2.0))
      ! production term in SA model
      this%Pk(i,k) = cb1*this%nuSA(i,k)*shatSA
    enddo
  enddo
end subroutine production_SA

subroutine diffusion_SA(this,putout,ekmt,ek,eki,ekk,sigma,rho_mod,modification)
  use mod_param, only : k1,i1,kmax,imax,k,i
  use mod_mesh,  only : dzw,dzp,Ru,Rp,dru
  implicit none
  class(SA_TurbModel) :: this
  real(8), dimension(0:i1, 0:k1), intent(IN) :: rho_mod,ek,eki,ekk,ekmt
  real(8),                        intent(IN) :: sigma
  integer,                        intent(IN) :: modification
  real(8), dimension(0:i1, 0:k1), intent(OUT):: putout
  integer :: km,kp,im,ip
  
  !Important note: this function takes instead of ek, eki, ekk, ekmt: eknu, eknui, eknuk, nuSANew, respectively.
  ! For, Standard, Inverse SLS and Aupoix
  ! rho=1 for standard
  do k=1,kmax
    kp=k+1
    km=k-1
    do i=1,imax
      putout(i,k) = putout(i,k) + 1.0/rho_mod(i,k)*( &
        ( (ekk(i,k ) + 0.5*(ekmt(i,k)+ekmt(i,kp))/sigma)* &
        sqrt(0.5*(rho_mod(i,k)+rho_mod(i,kp)))*((rho_mod(i,kp)**0.5)*ekmt(i,kp)-(rho_mod(i,k )**0.5)*ekmt(i,k ))/dzp(k) &
        -(ekk(i,km) + 0.5*(ekmt(i,k)+ekmt(i,km))/sigma)* &
        sqrt(0.5*(rho_mod(i,k)+rho_mod(i,km)))*((rho_mod(i,k )**0.5)*ekmt(i,k )-(rho_mod(i,km)**0.5)*ekmt(i,km))/dzp(km) &
        )/dzw(k))   
    enddo
  enddo
  
  ! For Aupoix we need to substract the density gradient diffusion with molecular viscosity
  !-1/rho d/dx[nu*nusa/2 drhodx]
  if (modification == 2) then                                        
    do k=1,kmax     ! in the z-direction
      kp=k+1
      km=k-1
      do i=1,imax
        putout(i,k) = putout(i,k) - 1.0/rho_mod(i,k)*((                         &
           ekk(i,k )*0.5*(ekmt(i,k)+ekmt(i,kp))/2*(rho_mod(i,kp)-rho_mod(i,k ))/dzp(k) &
          -ekk(i,km)*0.5*(ekmt(i,k)+ekmt(i,km))/2*(rho_mod(i,k )-rho_mod(i,km))/dzp(km) &
          )/dzw(k))
      enddo
    enddo
    do i=1,imax     ! in the r-direction
      ip=i+1
      im=i-1
      do k=1,kmax
        putout(i,k) = putout(i,k) - 1.0/rho_mod(i,k)* &
          (  Ru(i )/((Rp(ip)-Rp(i ))*Rp(i)*dru(i))*(eki(i ,k)*0.5*(ekmt(i,k)+ekmt(ip,k))/2*(rho_mod(ip,k)-rho_mod(i ,k))) &
          -  Ru(im)/((Rp(i )-Rp(im))*Rp(i)*dru(i))*(eki(im,k)*0.5*(ekmt(i,k)+ekmt(im,k))/2*(rho_mod(i ,k)-rho_mod(im,k))) &
          )
      enddo
    enddo
  endif
end subroutine diffusion_SA

subroutine rhs_SA(this, putout,dimpl,nuSA,rho,modification)
  use mod_param, only : k1,i1,kmax,imax,k,i
  use mod_mesh,  only : dzw,dzp,drp,walldist
  implicit none
  class(SA_TurbModel) :: this
  real(8),dimension(0:i1,0:k1),intent(IN) :: rho,nuSA
  integer,                     intent(IN) :: modification
  real(8),dimension(0:i1,0:k1),intent(OUT):: putout,dimpl
  integer :: im,ip,km,kp
  real(8) :: cv1_3,cb1,cb2,cb3,cw1,cw2,cw3_6,inv_cb3,kappa_2,r_SA,g_SA,fw_SA,shatSA
  
  cv1_3     = (7.1)**3.0
  cb1       = 0.1355
  cb2       = 0.622
  cb3       = 2.0/3.0
  inv_cb3   = 1.0/cb3
  kappa_2   = (0.41)**2.0   ! von karman constant
  cw1       = (cb1/kappa_2) + (1.0+cb2)/cb3
  cw2       = 0.3
  cw3_6     = (2.0)**6.0

  if ((modification == 1) .or. (modification == 2)) then
    do k=1,kmax
      kp=k+1
      km=k-1
      do i=1,imax
        ip=i+1
        im=i-1
                
        shatSA = this%Pk(i,k)/(cb1*nuSA(i,k))
        ! destruction term in SA model
        r_SA         = min(nuSA(i,k)/(kappa_2*(walldist(i)**2.0)*shatSA), 10.0)
        g_SA         = r_SA + cw2*((r_SA**6.0) - r_SA)
        fw_SA        = g_SA*(((1.0 + cw3_6)/(g_SA**6.0 + cw3_6))**(1.0/6.0))
        ! gustavo: i think this is not correct
        !destrSA(i,k) = cw1/rho(i,k)*fw_SA*nuSAtmp(i,k)/(wallDist(i)**2)
        dimpl(i,k) = dimpl(i,k) + cw1*fw_SA*nuSA(i,k)/(walldist(i)**2.0)
        ! source term
        ! invSLS and Aupoix SA model=  advection + Pk + (1/rho)*cb2/cb3*(d(nuSA*sqrt(rho))/dr)^2 +(d(nuSA*sqrt(rho))/dz)^2
        putout(i,k) = putout(i,k) + this%Pk(i,k) + cb2*inv_cb3/rho(i,k) * ( &
             (((nuSA(ip,k)*(rho(ip,k)**0.5))- (nuSA(im,k)*(rho(im,k)**0.5)))/(dRp(i)+dRp(im)))**2.0 &
           + (((nuSA(i,kp)*(rho(i,kp)**0.5))- (nuSA(i,km)*(rho(i,km)**0.5)))/(dzp(k)+dzp(km)))**2.0  )

      enddo
    enddo
  else
    do k=1,kmax
      kp=k+1
      km=k-1
      do i=1,imax
        ip=i+1
        im=i-1

        shatSA = this%Pk(i,k)/(cb1*nuSA(i,k))
        ! destruction term in SA model
        r_SA         = min(nuSA(i,k)/(kappa_2*(walldist(i)**2.0)*shatSA), 10.0)
        g_SA         = r_SA + cw2*((r_SA**6.0) - r_SA)
        fw_SA        = g_SA*(((1.0 + cw3_6)/(g_SA**6.0 + cw3_6))**(1.0/6.0))
        ! gustavo: i think this is not correct
        !destrSA(i,k) = cw1/rho(i,k)*fw_SA*nuSAtmp(i,k)/(wallDist(i)**2)
        dimpl(i,k) = dimpl(i,k) + cw1*fw_SA*nuSA(i,k)/(walldist(i)**2.0)
        ! source term
        ! invSLS and Aupoix SA model=  advection + Pk + (1/rho)*cb2/cb3*(d(nuSA*sqrt(rho))/dr)^2 +(d(nuSA*sqrt(rho))/dz)^2
        !putout(i,k) = putout(i,k) + this%Pk(i,k) + cb2*inv_cb3* ( & !!ORIGINAL NOTE STEPHAN!!!
        putout(i,k) = putout(i,k) + this%Pk(i,k) + cb2*inv_cb3 * ( &
            ((nuSA(ip,k) - nuSA(im,k))/(dRp(i)+dRp(im)))**2.0 + ((nuSA(i,kp) - nuSA(i,km))/(dzp(k)+dzp(km)))**2.0  )
      enddo
    enddo

  endif
end subroutine rhs_SA


end module sa_tm
