module ktet_tdm
  use mod_tm
  implicit none

!****************************************************************************************

  !************************!
  !         KtEt class       !
  !************************!
  
  type,abstract,extends(TurbDiffModel), public :: KtEt_TurbDiffModel
  real(8), dimension(:,:), allocatable :: Ttemp,Tmix,flambda
  real(8), dimension(:),   allocatable :: epstin, ktin
  real(8) :: sigmakt,sigmaet,clambda,cp1,cp2,cd1,cd2
  contains
    procedure(set_alphat_KE), deferred :: set_alphat
    procedure(advance_KtEt), deferred :: advance_turbdiff
    
    procedure(set_constants), deferred :: set_constants
    procedure(init_w_inflow_KtEt), deferred :: init_w_inflow  ! MISSING
    procedure :: init => init_KtEt
    procedure :: rhs_kt_KtEt
    procedure :: rhs_epst_KtEt
    procedure :: diffusion_epst_KtEt
    procedure :: solve_epst_KtEt
    procedure :: set_bc_KtEt
    procedure :: production_KtEt
    procedure :: solve_kt_KtEt
    procedure :: init_mem_KtEt
    procedure :: init_sol => init_sol_KtEt
    procedure :: get_profile => get_profile_KtEt
    procedure :: get_sol => get_sol_KtEt

  end type KtEt_TurbDiffModel

  interface

    subroutine set_constants(this)
      import :: KtEt_TurbDiffModel
      class(KtEt_TurbDiffModel) :: this
    end subroutine set_constants

    subroutine set_alpha_KtEt(this,u,w,rho,temp,mu,mui,ekh,alphat)
      implicit none
      class(KtEt_TurbDiffModel) :: this
      real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,rho,temp,mu,mui,ekh
      real(8),dimension(0:this%i1,0:this%k1),intent(OUT):: alphat
    end subroutine set_alpha_KtEt

    subroutine rhs_epst_KtEt(this,putout,dimpl,temp,rho,mu,ekh,alphat)
      implicit none
      class(DWX_TurbDiffModel) :: this
      real(8), dimension(0:this%i1,0:this%k1), intent(IN) :: rho,mu,temp,ekh,alphat
      real(8), dimension(0:this%i1,0:this%k1), intent(OUT):: putout,dimpl
    end subroutine rhs_epst_KtEt  

    subroutine init_w_inflow_KtEt(this)
!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end subroutine init_w_inflow_KtEt
  end interface

contains
!****************************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!****************************************************************************************

  !************************!
  !       KE routines      !
  !************************!

subroutine init_KtEt(this)   
  class(KtEt_TurbDiffModel) :: this
  call this%init_mem_KtEt()
  call this%set_constants()
  call this%init_sol()
end subroutine init_KtEt

subroutine init_sol_KtEt(this)   
  class(KtEt_TurbDiffModel) :: this
  integer :: i
  do i=0,this%i1
    this%kt(i,:)  =0.0001
    this%epst(i,:)=0.0010
    this%Pkt(i,:) = 0 
    this%ktin(i) = 0.0
    this%epstin(i) = 0.0
    this%Pktin(i) = 0 
  enddo
end subroutine init_sol_KtEt

subroutine init_mem_KtEt(this)  
  class(KtEt_TurbDiffModel) :: this
  allocate(this%epst(0:this%i1,0:this%k1),this%kt (0:this%i1,0:this%k1),     &
           this%Pkt(0:this%i1,0:this%k1), this%flambda(0:this%i1,0:this%k1), &
           this%Ttemp (0:this%i1,0:this%k1),this%Tmix(0:this%i1,0:this%k1)
           
  allocate(this%alphatin(0:this%i1),this%Pktin (0:this%i1), &
           this%epstin(0:this%i1),this%ktin(0:this%i1))
end subroutine init_mem_KtEt

subroutine solve_kt_KtEt(this,resKt,u,w,rho,ekh,ekhi,ekhk,alphat,rho_mod, &
                       alphakt,modification,rank,centerBC,periodic)
  use mod_math
  use mod_mesh, only : bot_bcnovalue,top_bcnovalue,Rp,Ru,dRp,dRu,dz
  implicit none
  class(KtEt_TurbDiffModel) :: this
  real(8),dimension(0:this%i1,0:this%k1), intent(IN) :: u, w, rho,ekh,ekhi,ekhk,alphat
  real(8),                                intent(IN) :: alphakt
  integer,                                intent(IN) :: modification,rank,centerBC,periodic
  real(8),                                intent(OUT):: resKt
  real(8), dimension(0:this%i1,0:this%k1) :: dnew,dimpl
  real(8), dimension(this%imax)           :: a,b,c,rhs
  integer                                 :: i,k

  resKt  = 0.0;  dnew  = 0.0; dimpl = 0.0;

  call advecc(dnew,dimpl,this%kt,u,w,Ru,Rp,dRu,dz,this%i1,this%k1,rank,periodic,.true.)
  call this%rhs_kt_KtEt(dnew,dimpl,rho) 
  call diffc(dnew,this%kt,ekh,ekhi,ekhk,alphat,this%sigmakt,rho,Ru,Rp,dRu,dz,rank,modification)

  do k=1,this%kmax
    do i=1,this%imax
      if ((modification == 0) .or. (modification == 1)) then
      a(i) = (ekhi(i-1,k)+0.5*(alphat(i,k)+alphat(i-1,k))/this%sigmakt)
      a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/rho(i,k)

      c(i) = (ekhi(i  ,k)+0.5*(alphat(i,k)+alphat(i+1,k))/this%sigmakt)
      c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/rho(i,k)

      b(i) = (-a(i)-c(i)) + dimpl(i,k)

      rhs(i) = dnew(i,k) + ((1-alphakt)/alphakt)*b(i)*this%kt(i,k)
    enddo


     i=1
     b(i) = b(i) + bot_bcnovalue(k)*a(i) !symmetry = -1 ; wall = 1 
     rhs(i) = dnew(i,k) + ((1-alphakt)/alphakt)*b(i)*this%kt(i,k)
       
     i=this%imax
     b(i) = b(i) + top_bcnovalue(k)*c(i)
     rhs(i) = dnew(i,k) + ((1-alphakt)/alphakt)*b(i)*this%kt(i,k)

     call matrixIdir(this%imax,a,b/alphakt,c,rhs)

    do i=1,this%imax
      resKt = resKt + ((this%kt(i,k) - rhs(i))/(this%kt(i,k)+1.0e-20))**2.0
      this%kt(i,k) = max(rhs(i), 1.0e-8)
    enddo
  enddo
end
subroutine get_sol_KtEt(this,kt,epst,yp)
    class(KtEt_TurbDiffModel) :: this
    real(8),dimension(0:this%i1,0:this%k1), intent(OUT):: kt,epst,yp
    nuSA=0
    kt   =this%kt    
    epst =this%epst
    yp  = this%yp
end subroutine get_sol_KtEt

subroutine get_profile_KtEt(this,p_kt,p_epst,p_Pkt,yp,k)
  class(KtEt_TurbDiffModel) :: this
  integer,                               intent(IN) :: k
  real(8),dimension(0:this%i1),          intent(OUT):: p_kt,p_epst,p_Pkt,yp
  p_kt(:)   =this%kt(:,k)
  p_epst(:) =this%epst(:,k)
  p_Pkt(:)  =this%Pkt(:,k)
  yp(:)    =this%yp(:,k)
end subroutine get_profile_KtEt

subroutine advance_KtEt(this,u,w,c,temp,rho,mu,ekh,ekhi,ekhk,alphat, &
                      alpha1,alpha2,alpha3,                   &
                      modification,rank,centerBC,periodic,    &
                      residual1, residual2, residual3)
  use mod_mesh, only : Ru,Rp,dRp,dRu,dz,walldist
  class(DWX_TurbDiffModel) :: this
  real(8), dimension(0:this%i1,0:this%k1),intent(IN) :: u,w,c,temp,rho,mu,ekh,ekhi,ekhk,alphat
  real(8),                                intent(IN) :: alpha1,alpha2, alpha3
  integer,                                intent(IN) :: modification,rank,centerBC,periodic
  real(8),                                intent(OUT):: residual1,residual2, residual3
  real(8), dimension(0:this%i1,0:this%k1) :: rho_mod

  !1, our modification, 2, Aupoix modification
  if ((modification == 1) .or. (modification == 2)) then
    rho_mod = rho
  else
    rho_mod = 1.0
  endif

  call this%production_KtEt(this,c,temp,rho,cp,alphat,dRu,dRp,dz)
  call this%solve_epst_KtEt(residual2,u,w,temp,rho,mu,ekh,ekhi,ekhk,alphat,rho_mod, &
                       alpha2,modification,rank,centerBC,periodic)
  call this%solve_kt_KtEt(residual1,u,w,rho,ekh,ekhi,ekhk,alphat,rho_mod, &
                       alpha1,modification,rank,centerBC,periodic)
end

subroutine production_KtEt(this,c,temp,rho,cp,alphat)
  use module_mesh, only : mesh
  use mod_mesh, only : dRp,dRu,dz
  implicit none
  class(KtEt_TurbDiffModel) :: this
  real(8), dimension(0:this%i1,0:this%k1), intent(IN) :: c,temp,rho,cp,alphat
  real(8), dimension(0:this%i1,0:this%k1) :: div
  integer im,ip,jm,jp,km,kp,ib,ie,kb,ke,i,k
  real(8), dimension(0:this%k1) :: dzw, dzp

  dzw = mesh%dzw
  dzp = mesh%dzp

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
      ! Production of temperature fluctuations  Pkt= <u'j T'> dTdxj= (lambda_t/cp)/cp dCdxj dTdxj = = (lambda_t/cp)dTdxj^2
      this%Pkt(i,k) = alphat(i,k)/cp(i,k)                           *  &
                (((c(i,k)-c(i,km))/dz)*((temp(i,k)-temp(i,km))/dz)  +  &
               + ((c(i,k)-c(im,k))/dRu(i)*((temp(i,k)-temp(im,k))/dRu(i))))
!     this%Pkt(i,k) = alphat(i,k)*(((temp(i,k)-temp(i,km))/dz)**2.0+ ((temp(i,k)-temp(im,k))/dRu(i))**2.0)

    enddo
  enddo
end subroutine production_KtEt

subroutine solve_epst_KtEt(this,resEt,u,w,temp,rho,mu,ekh,ekhi,ekhk,alphat,rho_mod, &
                       alphaet,modification,rank,centerBC,periodic)
  use mod_math
  use mod_mesh, only : bot_bcvalue, top_bcvalue,Ru,Rp,dRp,dRu,dz
  implicit none
  class(TurbDiffModel) :: this
  real(8),dimension(0:this%i1,0:this%k1), intent(IN) :: u, w,temp,rho,mu,ekh,ekhi,ekhk,alphat
  real(8),                                intent(IN) :: alphaet
  integer,                                intent(IN) :: modification,rank,centerBC,periodic
  real(8),                                intent(OUT):: resEt
  real(8), dimension(0:this%i1,0:this%k1) :: dnew,dimpl
  real(8), dimension(this%imax)           :: a,b,c,rhs
  integer                                 :: i,k
  
  resEt  = 0.0; dnew  = 0.0; dimpl = 0.0;

  call advecc(dnew,dimpl,this%epst,u,w,Ru,Rp,dru,dz,this%i1,this%k1,rank,periodic,.true.)
  call this%rhs_epst_KtEt(dnew,dimpl,temp,rho,mu,ekh,alphat) 
  call this%diffusion_epst_KtEt(dnew,this%epst,ekhk,alphat,this%sigmaet,rho,dz,modification)

  do k=1,this%kmax
    do i=1,this%imax
                  
      a(i) = (ekhi(i-1,k)+0.5*(alphat(i,k)+alphat(i-1,k))/this%sigmaet)
      a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/rho(i,k)
  
      c(i) = (ekhi(i  ,k)+0.5*(alphat(i,k)+alphat(i+1,k))/this%sigmaet)
      c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/rho(i,k)
  
      b(i) = ((-a(i)-c(i))+ dimpl(i,k)  )  
      
      rhs(i) = dnew(i,k) + ((1-alphaet)/alphaet)*b(i)*this%epst(i,k)

    enddo

    i=1
    b(i) = b(i)+bot_bcvalue(k)*a(i)
    rhs(i) = dnew(i,k) - (1.-bot_bcvalue(k))*a(i)*this%epst(i-1,k) + ((1-alphaet)/alphaet)*b(i)*this%epst(i,k)   !wall with value

    i=this%imax
    b(i) = b(i)+top_bcvalue(k)*c(i)
    rhs(i) = dnew(i,k) - (1.-top_bcvalue(k))*c(i)*this%epst(i+1,k) + ((1-alphaet)/alphaet)*b(i)*this%epst(i,k)   !wall with value
      
    call matrixIdir(this%imax,a,b/alphaet,c,rhs)
  
    do i=1,this%imax
      resEt = resEt + ((this%epst(i,k) - rhs(i))/(this%epst(i,k)+1.0e-20))**2.0
      this%epst(i,k) = max(rhs(i), 1.0e-8)
  
    enddo
  enddo
end subroutine solve_epst_KtEt

subroutine rhs_kt_KtEt(this,putout,dimpl,rho)
  implicit none
  class(TurbDiffModel) :: this
  real(8), dimension(0:this%i1,0:this%k1), intent(IN) :: rho
  real(8), dimension(0:this%i1,0:this%k1), intent(OUT):: putout,dimpl
  integer ib,ie,kb,ke,i,k
  
  ib = 1
  ie = this%i1-1

  kb = 1
  ke = this%k1-1

  do k=kb,ke
    do i=ib,ie
      !kt equation
      putout(i,k) = putout(i,k)+2.0*this%Pkt(i,k)/rho(i,k)
      dimpl(i,k)  = dimpl(i,k) +2.0*this%epst(i,k)/(this%kt(i,k)+1.0e-20) ! note, rho*epsilon/(rho*k), set implicit and divided by density
    enddo
  enddo
end subroutine rhs_kt_KtEt

subroutine set_bc_KtEt(this,ekh,rho,centerBC,periodic,rank,px)
  use mod_mesh, only : top_bcvalue, bot_bcvalue,top_bcnovalue, bot_bcnovalue,walldist
  implicit none
  class(MK_TurbModel) :: this
  real(8),dimension(0:this%i1,0:this%k1),intent(IN) :: rho,ekh
  integer,                               intent(IN) :: centerBC,periodic, rank, px
  real(8),dimension(0:this%k1) :: BCvalue
  real(8),dimension(0:this%i1) :: tmp
  real(8) :: topBCvalue, botBCvalue
  integer :: k

  do k = 0,this%k1 
    this%kt(0,k)         = bot_bcnovalue(k)*this%kt(1,k)         !symmetry or 0 value
    this%kt(this%i1,k)   = top_bcnovalue(k)*this%kt(this%imax,k) !symmetry or 0 value
    botBCvalue = 2.0*ekh(1,k)/rho(1,k)*this%kt(1,k)/walldist(1)**2                                                          !bcvalue
    this%epst(0,k)       = (1.-bot_bcvalue(k))*(2.0*botBCvalue-this%epst(1,k))         +bot_bcvalue(k)*this%epst(1,k)        !symmetry or bc value
    topBCvalue = 2.0*ekh(this%imax,k)/rho(this%imax,k)*this%kt(this%imax,k)/walldist(this%imax)**2                          !bcvalue
    this%eps(this%i1,k) = (1.-top_bcvalue(k))*(2.0*topBCvalue-this%epst(this%imax,k)) +top_bcvalue(k)*this%epst(this%imax,k)!symmetry or bc value
  enddo


  ! this%kt(this%i1,:) = -this%kt(this%imax,:)
  ! BCvalue(:) = 2.0*ekh(this%imax,:)/rho(this%imax,:)*this%kt(this%imax,:)/walldist(this%imax)**2
  ! this%epst(this%i1,:) = 2.0*BCvalue(:) - this%epst(this%imax,:)

  ! ! channel
  ! if (centerBC.eq.-1) then
  !   this%kt(0,:)  = -this%kt(1,:)
  !   BCvalue(:)   = 2.0*ekh(1,:)/rho(1,:)*this%kt(1,:)/walldist(1)**2
  !   this%epst(0,:)= 2.0*BCvalue(:) - this%epst(1,:)
  ! endif

  ! ! pipe/BL
  ! if (centerBC.eq.1) then
  !   this%kt  (0,:) = this%kt  (1,:)
  !   this%epst(0,:) = this%epst(1,:)
  ! endif

  call shiftf(this%kt,  tmp,rank); this%kt  (:,0)      =tmp(:);
  call shiftf(this%epst,tmp,rank); this%epst(:,0)      =tmp(:);
  call shiftb(this%kt,  tmp,rank); this%kt  (:,this%k1)=tmp(:);
  call shiftb(this%epst,tmp,rank); this%epst(:,this%k1)=tmp(:);
  
  ! developing
  if (periodic.eq.1) return
  if (rank.eq.0) then
    this%kt  (:,0) = this%ktin(:)
    this%epst(:,0) = this%epstin(:)
  endif
  if (rank.eq.px-1) then
    this%kt  (:,this%k1)= 2.0*this%kt  (:,this%kmax)-this%kt  (:,this%kmax-1)
    this%epst(:,this%k1)= 2.0*this%epst(:,this%kmax)-this%epst(:,this%kmax-1)
  endif

end subroutine set_bc_KtEt

subroutine diffusion_epst_KtEt(this,putout,putin,ekhk,alphat,sigma,rho,modification)
  use module_mesh, only : mesh
  use mod_mesh, only : dz
  implicit none
  class(TurbDiffModel) :: this
  real(8), dimension(0:this%i1, 0:this%k1), intent(IN) :: putin, ekhk, alphat, rho
  real(8),                                  intent(IN) :: sigma
  real(8), dimension(0:this%i1, 0:this%k1), intent(OUT):: putout
  integer i,k,kp,km
  real(8), dimension(0:this%k1) :: dzw, dzp

  dzw = mesh%dzw
  dzp = mesh%dzp                              ! Standard
  do k=1,this%k1-1
    kp=k+1
    km=k-1
    do i=1,this%i1-1
      putout(i,k) = putout(i,k) + 1.0/rho(i,k)*( &
        ( (ekhk(i,k ) + 0.5*(alphat(i,k)+alphat(i,kp))/sigma)*(putin(i,kp)-putin(i,k )) &
        - (ekhk(i,km) + 0.5*(alphat(i,k)+alphat(i,km))/sigma)*(putin(i,k )-putin(i,km)))/(dz*dz)   )
      ! putout(i,k) = putout(i,k) + 1.0/rho(i,k)*( &
      !   ( (muk(i,k ) + 0.5*(mut(i,k)+mut(i,kp))/sigma)*(putin(i,kp)-putin(i,k ))/dzp(k) &
      !   - (muk(i,km) + 0.5*(mut(i,k)+mut(i,km))/sigma)*(putin(i,k )-putin(i,km))/dzp(km)&
      !   )/dzw(k)   )
    enddo
  enddo
  
end subroutine diffusion_epst_KtEt


end module