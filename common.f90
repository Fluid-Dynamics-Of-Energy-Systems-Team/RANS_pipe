module mod_common
  implicit none

  public


  integer istep,nsteps,inlet
  real(8), dimension(:,:), allocatable :: ekmt,ekme,temp,peclet,beta, &
                                          ekh,ekhi,ekhk, & 
                                          ekm,ekmi,ekmk, &
                                          cp,cpi,cpk,alphat
  real(8) dt
  real(8), dimension(:,:), allocatable :: Unew,Wnew,rnew,cnew,qcrit, &
                                          Uold,Wold,rold,dUdt,dVdt,dWdt,p,res_nuSA
  real(8), dimension(:),   allocatable :: Win,ekmtin, dis


  real(8), dimension(:), allocatable :: wstress, sfriction, mom_thickness, dis_thickness, bl_thickness
  
  integer, dimension(:),   allocatable :: Xii     !Nx
  integer, dimension(:,:), allocatable :: Xkk     !Nx,Mt
  integer, dimension(:),   allocatable :: Tkk     !Nt
  integer, dimension(:,:), allocatable :: Tii     !Nt,Mx
  real(8), dimension(:,:), allocatable :: W1, W2  !Mt,Nx
  real(8), dimension(:,:), allocatable :: W1t, W2t!Mx,Nt
  
contains

  subroutine initMem()
    use mod_param
    implicit none
    !STATE VARIABLES
    allocate(ekmt(0:i1,0:k1),ekme(0:i1,0:k1),temp(0:i1,0:k1),peclet(0:i1,0:k1),beta(0:i1,0:k1), &
              ekh(0:i1,0:k1),ekhi(0:i1,0:k1),ekhk(0:i1,0:k1), &
              ekm(0:i1,0:k1),ekmi(0:i1,0:k1),ekmk(0:i1,0:k1), &
               cp(0:i1,0:k1), cpi(0:i1,0:k1), cpk(0:i1,0:k1), alphat(0:i1,0:k1))
    !EQUATION VARIABLES
    allocate(rnew(0:i1,0:k1),Unew(0:i1,0:k1),Wnew(0:i1,0:k1),Cnew(0:i1,0:k1), &
             rold(0:i1,0:k1),Uold(0:i1,0:k1),Wold(0:i1,0:k1),                 &
             dUdt(0:i1,0:k1),dVdt(0:i1,0:k1),dWdt(0:i1,0:k1),                 &
             res_nuSA(0:i1,0:k1))
    allocate(Win(0:i1),ekmtin(0:i1))
    allocate(dis(0:k1))
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

    allocate(wstress(0:k1),sfriction(0:k1),mom_thickness(0:k1),dis_thickness(0:k1),bl_thickness(0:k1))
  end

end module mod_common
      
