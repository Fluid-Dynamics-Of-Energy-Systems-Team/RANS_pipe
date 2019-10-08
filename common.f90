module mod_common
  implicit none

  public

  integer istep,nsteps,inlet,centerBC,numDomain
  real(8) pi,dpdz

  ! KOM SST
  real(8) sigmat,sigmak,sigmae,cmu,ce1,ce2,sigmah2
  ! 0:i1,0:k1
  real(8), dimension(:,:), allocatable :: fmu,f1,f2,ReT,yp,ReTauS,ypt,dterm,eterm,atmp,bF1
  ! imax,kmax
  real(8), dimension(:,:), allocatable :: fv2,Lh,at,LhT,bF2,cdKOM

  !0:i1
  real(8), dimension(:), allocatable :: Ru,Rp,y_fa,y_cv,dru,drp
  real(8) dz
  !0:k1
  real(8),  dimension(:), allocatable :: z1,z2
  !1:imax
  real(8),  dimension(:), allocatable :: wallDist

  !0:i1,0:k1 
  real(8), dimension(:,:), allocatable :: ekm,ekmt,ekme,ekh,cp,temp,tco,Pk,Gk,peclet,beta,ekhi,ekhk,ekmi,ekmk,cpi,cpk,Tt
  real(8) enth_wall

  real(8) dt,dtmax
  
  !0:i1,0:k1 
  real(8), dimension(:,:), allocatable :: Unew,Vnew,Wnew,v2new,h2new,Cnew,knew,enew,qcrit,nuSAnew,omNew

  !0:i1,0:k1 
  real(8), dimension(:,:), allocatable :: Uold,Vold,eold,v2old,Wold,Cold,kold,h2old,nuSAold,omOld

  !0:i1,0:k1 
  real(8), dimension(:,:), allocatable :: dUdt,dVdt,dWdt

  !0:i1,0:k1 
  real(8), dimension(:,:), allocatable :: rold,rnew

  !imax,kmax
  real(8), dimension(:,:), allocatable :: bdts,p,phi
  !iwork
  real(8), dimension(:), allocatable  :: work
  !isave
  real(8), dimension(:), allocatable  :: save
  !kmax
  real(8), dimension(:), allocatable  :: bdrs
  !imax
  real(8), dimension(:), allocatable  :: bdzs
  

  !1:nTab
  real(8), dimension(:), allocatable :: tempTab,rhoTab,betaTab, muTab,lamTab, cpTab,enthTab,lamocpTab, &
                                          temp2Tab,rho2Tab,beta2Tab, mu2Tab,lam2Tab, cp2Tab,enth2Tab,lamocp2Tab


                                          
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

    allocate(fmu(0:i1,0:k1),f1(0:i1,0:k1),f2(0:i1,0:k1),ReT(0:i1,0:k1),yp(0:i1,0:k1),ReTauS(0:i1,0:k1),    &
             ypt(0:i1,0:k1),dterm(0:i1,0:k1),eterm(0:i1,0:k1),atmp(0:i1,0:k1),bF1(0:i1,0:k1))
    allocate(fv2(imax,kmax),Lh(imax,kmax),at(imax,kmax),LhT(imax,kmax),bF2(imax,kmax),cdKOM(imax,kmax))

    allocate( Ru(0:i1),Rp(0:i1),y_fa(0:i1),y_cv(0:i1),dru(0:i1),drp(0:i1))
    allocate(z1(0:k1),z2(0:k1))
    allocate(wallDist(1:imax))

    allocate(ekm(0:i1,0:k1),ekmt(0:i1,0:k1),ekme(0:i1,0:k1),ekh(0:i1,0:k1),cp(0:i1,0:k1),temp(0:i1,0:k1),  &
             tco(0:i1,0:k1),Pk(0:i1,0:k1),Gk(0:i1,0:k1),peclet(0:i1,0:k1),beta(0:i1,0:k1),ekhi(0:i1,0:k1), &
             ekhk(0:i1,0:k1),ekmi(0:i1,0:k1),ekmk(0:i1,0:k1),cpi(0:i1,0:k1),cpk(0:i1,0:k1),Tt(0:i1,0:k1))

    allocate(Unew(0:i1,0:k1),Vnew(0:i1,0:k1),Wnew(0:i1,0:k1),v2new(0:i1,0:k1),h2new(0:i1,0:k1),            &
             Cnew(0:i1,0:k1),knew(0:i1,0:k1),enew(0:i1,0:k1),qcrit(0:i1,0:k1),nuSAnew(0:i1,0:k1),          &
             omNew(0:i1,0:k1))

    allocate(Uold(0:i1,0:k1),Vold(0:i1,0:k1),eold(0:i1,0:k1),v2old(0:i1,0:k1),Wold(0:i1,0:k1),             &
             Cold(0:i1,0:k1),kold(0:i1,0:k1),h2old(0:i1,0:k1),nuSAold(0:i1,0:k1),omOld(0:i1,0:k1))

    allocate(dUdt(0:i1,0:k1),dVdt(0:i1,0:k1),dWdt(0:i1,0:k1))

    allocate(rold(0:i1,0:k1),rnew(0:i1,0:k1))

    allocate(bdts(imax,kmax),p(imax,kmax),phi(imax,kmax))

    allocate(work(iwork))
    allocate(save(isave))
    allocate(bdrs(kmax))
    allocate(bdzs(imax))

    allocate(tempTab(1:nTab),rhoTab(1:nTab),betaTab(1:nTab),muTab(1:nTab),lamTab(1:nTab),cpTab(1:nTab), &
             enthTab(1:nTab),lamocpTab(1:nTab),temp2Tab(1:nTab),rho2Tab(1:nTab),beta2Tab(1:nTab),       &
             mu2Tab(1:nTab),lam2Tab(1:nTab),cp2Tab(1:nTab),enth2Tab(1:nTab),lamocp2Tab(1:nTab))


    allocate(Xii(Nx))
    allocate(Xkk(Nx,Mt))
    allocate(Tkk(Nt))
    allocate(Tii(Nt,Mx))
    allocate(W1t(Mx,Nt),W2t(Mx,Nt))
    allocate(W1(Mt,Nx),W2(Mt,Nx))

  end

end module mod_common
      
