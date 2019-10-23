module mod_common
  implicit none

  public

  integer istep,nsteps,inlet


  !***************STATE PROPERTIES
  !0:i1,0:k1 
  real(8), dimension(:,:), allocatable :: ekm,ekmt,ekme,ekh,cp,temp,peclet,beta,ekhi,ekhk,ekmi,ekmk,cpi,cpk
  real(8) dt,dtmax
  
  !***************EQUATION VARIABLES
  !0:i1,0:k1 
  real(8), dimension(:,:), allocatable :: Unew,Vnew,Wnew,rnew,cnew,qcrit,&!v2new,h2new,Cnew,knew,enew,,nuSAnew,omNew, &
                                          Uold,Wold,rold,dUdt,dVdt,dWdt
  real(8), dimension(:), allocatable   :: Win,ekmtin
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

    !STATE VARIABLES
    allocate(ekm(0:i1,0:k1),ekmt(0:i1,0:k1),ekme(0:i1,0:k1),ekh(0:i1,0:k1),cp(0:i1,0:k1),temp(0:i1,0:k1),  &
             peclet(0:i1,0:k1),beta(0:i1,0:k1),ekhi(0:i1,0:k1), &
             ekhk(0:i1,0:k1),ekmi(0:i1,0:k1),ekmk(0:i1,0:k1),cpi(0:i1,0:k1),cpk(0:i1,0:k1))
    !EQUATION VARIABLES
    allocate(Unew(0:i1,0:k1),Vnew(0:i1,0:k1),Wnew(0:i1,0:k1),&!,v2new(0:i1,0:k1),h2new(0:i1,0:k1),            &
             Cnew(0:i1,0:k1),&!knew(0:i1,0:k1),enew(0:i1,0:k1),nuSAnew(0:i1,0:k1), omNew(0:i1,0:k1),         &
             rnew(0:i1,0:k1),rold(0:i1,0:k1),Uold(0:i1,0:k1),Wold(0:i1,0:k1),             &
             dUdt(0:i1,0:k1),dVdt(0:i1,0:k1),dWdt(0:i1,0:k1))
    allocate(Win(0:i1),ekmtin(0:i1))
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
      
