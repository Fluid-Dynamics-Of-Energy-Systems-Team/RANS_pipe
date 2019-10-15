module mod_common
  implicit none

  public

  integer istep,nsteps,inlet,centerBC,numDomain
  real(8) pi,dpdz


  !***************TURBULENCE
  real(8) sigmat,sigmak,sigmae,cmu,ce1,ce2,sigmah2
  ! 0:i1,0:k1
  real(8), dimension(:,:), allocatable :: fmu,f1,f2,ReT,yp,ReTauS,ypt,dterm,eterm,bF1,Pk,Gk,Tt
  ! imax,kmax
  real(8), dimension(:,:), allocatable :: fv2,Lh,LhT,bF2,cdKOM

  !***************GRID
  real(8) dz
  !0:i1
  real(8), dimension(:), allocatable :: Ru,Rp,y_fa,y_cv,dru,drp
  !0:k1
  real(8),  dimension(:), allocatable :: z1,z2
  !1:imax
  real(8),  dimension(:), allocatable :: wallDist

  !***************STATE PROPERTIES
  !0:i1,0:k1 
  real(8), dimension(:,:), allocatable :: ekm,ekmt,ekme,ekh,cp,temp,peclet,beta,ekhi,ekhk,ekmi,ekmk,cpi,cpk
  real(8) enth_wall

  real(8) dt,dtmax
  
  !***************EQUATION VARIABLES
  !0:i1,0:k1 
  real(8), dimension(:,:), allocatable :: Unew,Vnew,Wnew,rnew,v2new,h2new,Cnew,knew,enew,qcrit,nuSAnew,omNew
  !0:i1,0:k1 
  real(8), dimension(:,:), allocatable :: Uold,Wold,rold
  !0:i1,0:k1 
  real(8), dimension(:,:), allocatable :: dUdt,dVdt,dWdt
  !imax,kmax
  real(8), dimension(:,:), allocatable :: p
               
  !***************NUMERICAL CLUTTER                         
  !Nx
  integer, dimension(:), allocatable :: Xii
  !Nx,Mt
  integer, dimension(:,:), allocatable :: Xkk
  !Nt
  integer, dimension(:), allocatable :: Tkk
  !Nt,Mx
  integer, dimension(:,:), allocatable :: Tii
  !Mt,Nx
  real(8), dimension(:,:), allocatable :: W1, W2
  !Mx,Nt
  real(8), dimension(:,:), allocatable :: W1t, W2t

    
contains

  subroutine initMem()
    use mod_param
    implicit none

    !TURBULENCE
    allocate(fmu(0:i1,0:k1),f1(0:i1,0:k1),f2(0:i1,0:k1),ReT(0:i1,0:k1),yp(0:i1,0:k1),ReTauS(0:i1,0:k1),    &
             ypt(0:i1,0:k1),dterm(0:i1,0:k1),eterm(0:i1,0:k1),bF1(0:i1,0:k1),Pk(0:i1,0:k1),Gk(0:i1,0:k1),  &
             Tt(0:i1,0:k1))
    allocate(fv2(imax,kmax),Lh(imax,kmax),LhT(imax,kmax),bF2(imax,kmax),cdKOM(imax,kmax))
    !GRID
    allocate( Ru(0:i1),Rp(0:i1),y_fa(0:i1),y_cv(0:i1),dru(0:i1),drp(0:i1))
    allocate(z1(0:k1),z2(0:k1))
    allocate(wallDist(1:imax))
    !STATE VARIABLES
    allocate(ekm(0:i1,0:k1),ekmt(0:i1,0:k1),ekme(0:i1,0:k1),ekh(0:i1,0:k1),cp(0:i1,0:k1),temp(0:i1,0:k1),  &
             peclet(0:i1,0:k1),beta(0:i1,0:k1),ekhi(0:i1,0:k1), &
             ekhk(0:i1,0:k1),ekmi(0:i1,0:k1),ekmk(0:i1,0:k1),cpi(0:i1,0:k1),cpk(0:i1,0:k1))
    !EQUATION VARIABLES
    allocate(Unew(0:i1,0:k1),Vnew(0:i1,0:k1),Wnew(0:i1,0:k1),v2new(0:i1,0:k1),h2new(0:i1,0:k1),            &
             Cnew(0:i1,0:k1),knew(0:i1,0:k1),enew(0:i1,0:k1),nuSAnew(0:i1,0:k1),          &
             omNew(0:i1,0:k1),rnew(0:i1,0:k1))
    allocate(rold(0:i1,0:k1),Uold(0:i1,0:k1),Wold(0:i1,0:k1))
    allocate(dUdt(0:i1,0:k1),dVdt(0:i1,0:k1),dWdt(0:i1,0:k1))
    allocate(p(imax,kmax))

    !NUMERICAL STUFF
    allocate(qcrit(0:i1,0:k1))
    
    !NUMERICAL CLUTTER
    allocate(Xii(Nx))
    allocate(Xkk(Nx,Mt))
    allocate(Tkk(Nt))
    allocate(Tii(Nt,Mx))
    allocate(W1t(Mx,Nt),W2t(Mx,Nt))
    allocate(W1(Mt,Nx),W2(Mt,Nx))

  end

end module mod_common
      
