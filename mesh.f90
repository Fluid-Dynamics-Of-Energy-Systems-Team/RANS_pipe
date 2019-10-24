module mod_mesh
  !***************GRID
  real(8) :: dz,dpdz
  real(8), dimension(:), allocatable :: Ru,Rp,y_fa,y_cv,dru,drp !0:i1
  real(8), dimension(:), allocatable :: z1,z2 !0:k1
  real(8), dimension(:), allocatable :: wallDist !1:imax
  integer :: centerBC,numDomain

contains
!!********************************************************************
!!     mkgrid
!!********************************************************************
subroutine mkgrid(rank)
  use mod_param
  implicit none
  integer, intent(IN)       :: rank
  real(8)                   :: pi,fA,fB,fact,gridSize,Yplus,Rei
  real(8), dimension(imax)  :: delta
  
  !init mem
  allocate(Ru(0:i1),Rp(0:i1),y_fa(0:i1),y_cv(0:i1),dru(0:i1),drp(0:i1))
  allocate(z1(0:k1),z2(0:k1))
  allocate(wallDist(1:imax))

  
  pi    = 4.0*atan(1.0)
  Rei   = 1.0/Re
  dz    = 1.0*LoD/(kmax*px)
  ru(0) = 0

  !pipe
  if (systemSolve.eq.1) then
    if (rank.eq.0) print*,"************* SOLVING A PIPE FLOW *************!"
    numDomain = 1
    centerBC  = 1
    gridSize  = 0.5
    fA = 0.12
    fB = 2.4
    dpdz      = 4.0
  !channel
  elseif (systemSolve.eq.2) then
    if (rank.eq.0) print*,"************* SOLVING A CHANNEL FLOW *************!"
    numDomain = -1
    centerBC  = -1
    gridSize  = 2.0
    fA        = 0.5
    fB        = 4.6
    dpdz      = 1.0
  !bl
  elseif (systemSolve.eq.3) then
    if (rank.eq.0) print*,"************* SOLVING A BOUNDARY LAYER FLOW *************!"
    numDomain = -1
    centerBC  = 1
    gridSize  = 1.0
    fA = 0.12
    fB = 2.4
    dpdz      = 1.0
  else
    if (rank.eq.0) print '("systemSolve is ",i7," when it should be either 1 (pipe), 2(channel) or 3(BL)")', systemSolve
    stop
  endif

  !apply stretching
  do i = 1,imax
    fact = (i-0.)/(imax-0.)
    ru(i) = (1.-tanh(fB*(fA-fact))/tanh(fA*fB))
    ru(i) = ru(i)/(1.-tanh(fB*(fA-1.))/tanh(fA*fB))
    delta(i)=0.5*(ru(i)-ru(i-1))
  enddo

  !normalize with the gridsize
  do i=0,imax
    ru(i)=ru(i)/ru(imax)*gridSize
  enddo

  !calculate the cell centers and differences
  do i = 1 , imax
    Rp(i)  = (Ru(i)+Ru(i-1))/2.0
    dru(i) = (Ru(i)-Ru(i-1))
  enddo
  dru(i1) = dru(imax)
  Ru(i1) = Ru(imax) + dru(i1)
  Rp(i1) = Ru(imax) + dru(i1)/2.0
  dru(0)  = dru(1)
  Rp(0)  = Ru(0) - dru(0)/2.0
  do i = 0,imax
    drp(i) = Rp(i+1) - Rp(i)
  enddo

  !calculate wall distance
  !channel
  if (centerBC.eq.-1) then 
    do i = 1,imax
      if (rp(i).le.1) then
        wallDist(i) = rp(i)
      else
        wallDist(i) = gridSize-rp(i)
      endif
    enddo
  !bl/pipe
  else 
    do i = 1,imax
      wallDist(i) = gridSize - rp(i)
    enddo
  endif
  do i=0,i1
    y_cv(i)=rp(i)
    y_fa(i)=ru(i)
  enddo
 
  !conversion to carthesian coordinates
  if (numDomain.eq.-1) then
    do i=0,i1
      ru(i)=1.0
      rp(i)=1.0
    enddo
  endif

  !write the file
  if (rank.eq.0) then
    open(11,file = 'grid.txt')
    write(11,*) Re, imax
    do i=1,imax
      yplus = wallDist(i)*Re
      write(11,'(i5,4F12.6)') i,yplus,y_fa(i),y_cv(i),delta(max(1,i))
    enddo
    close(11)
  endif
  return
end
end module