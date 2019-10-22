module mod_common
  implicit none

  public

  integer istep,nsteps,inlet,centerBC,numDomain
  real(8) pi


  !***************GRID
  real(8) dz,dpdz
  real(8), dimension(:), allocatable :: Ru,Rp,y_fa,y_cv,dru,drp !0:i1
  real(8),  dimension(:), allocatable :: z1,z2 !0:k1
  real(8),  dimension(:), allocatable :: wallDist !1:imax
  !***************STATE PROPERTIES
  real(8), dimension(:,:), allocatable :: ekm,ekmt,ekme,ekh,cp,temp,peclet,beta,ekhi,ekhk,ekmi,ekmk,cpi,cpk  !0:i1,0:k1 
  ! real(8) enth_wall
  real(8) dt,dtmax
  !***************EQUATION VARIABLES
  real(8), dimension(:,:), allocatable :: Unew,Vnew,Wnew,rnew,Cnew,Uold,Wold,rold,dUdt,dVdt,dWdt   !0:i1,0:k1 
  real(8), dimension(:,:), allocatable :: p,qcrit   !0:i1,0:k1 
               
  !***************NUMERICAL CLUTTER                         
  integer, dimension(:),   allocatable :: Xii       !Nx
  integer, dimension(:,:), allocatable :: Xkk       !Nx,Mt
  integer, dimension(:),   allocatable :: Tkk       !Nt
  integer, dimension(:,:), allocatable :: Tii       !Nt,Mx
  real(8), dimension(:,:), allocatable :: W1, W2    !Mt,Nx
  real(8), dimension(:,:), allocatable :: W1t, W2t  !Mx,Nt

    
contains

  subroutine initMem()
    use mod_param
    implicit none

    !GRID
    allocate( Ru(0:i1),Rp(0:i1),y_fa(0:i1),y_cv(0:i1),dru(0:i1),drp(0:i1))
    allocate(z1(0:k1),z2(0:k1))
    allocate(wallDist(1:imax))
    !STATE VARIABLES
    allocate(ekm(0:i1,0:k1),ekmt(0:i1,0:k1),ekme(0:i1,0:k1),ekh(0:i1,0:k1),cp(0:i1,0:k1),temp(0:i1,0:k1),  &
             peclet(0:i1,0:k1),beta(0:i1,0:k1),ekhi(0:i1,0:k1), &
             ekhk(0:i1,0:k1),ekmi(0:i1,0:k1),ekmk(0:i1,0:k1),cpi(0:i1,0:k1),cpk(0:i1,0:k1))
    !EQUATION VARIABLES
    allocate(Unew(0:i1,0:k1),Vnew(0:i1,0:k1),Wnew(0:i1,0:k1),Cnew(0:i1,0:k1), &
             rnew(0:i1,0:k1),rold(0:i1,0:k1),Uold(0:i1,0:k1),Wold(0:i1,0:k1), &
             dUdt(0:i1,0:k1),dVdt(0:i1,0:k1),dWdt(0:i1,0:k1))
    allocate(p(imax,kmax))

    !NUMERICAL STUFF
    allocate(qcrit(0:i1,0:k1))
    
    !NUMERICAL CLUTTER
    allocate(Xii(Nx),Xkk(Nx,Mt),Tkk(Nt),Tii(Nt,Mx),W1t(Mx,Nt),W2t(Mx,Nt),W1(Mt,Nx),W2(Mt,Nx))

  end

end module mod_common
      
