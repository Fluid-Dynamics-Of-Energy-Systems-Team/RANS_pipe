
!>******************************************************************************************
!!      To calculate the rhs of the k equation  
!>******************************************************************************************
subroutine rhs_K(putout,dimpl,putink,putine,rho)
  implicit none
      include 'param.f90'
      include 'common.f90'

  integer ib,ie,kb,ke !< integers
  real*8, dimension(0:i1,0:k1) :: putout,rho,putink,putine,dimpl

  ib = 1
  ie = i1-1

  kb = 1
  ke = k1-1

  do k=kb,ke
    do i=ib,ie
      !k equation
      putout(i,k) = putout(i,k) + ( Pk(i,k) + Gk(i,k) )/rho(i,k)
      dimpl(i,k)  = dimpl(i,k) + putine(i,k)/putink(i,k)       ! note, rho*epsilon/(rho*k), set implicit and divided by density
    enddo
  enddo

end

!>******************************************************************************************
!!      To calculate the rhs of the epsilon equation 
!>******************************************************************************************
subroutine rhs_Epsilon(putout,dimpl,rho)
  implicit none
      include 'param.f90'
      include 'common.f90'

  integer ib,ie,kb,ke !< integers
  real*8, dimension(0:i1,0:k1) :: putout,rho,dimpl!,Tt

  ib = 1
  ie = i1-1

  kb = 1
  ke = k1-1

  do k=kb,ke
    do i=ib,ie
      !epsilon equation
      putout(i,k) = putout(i,k) +(ce1*f1(i,k)*Pk(i,k)/Tt(i,k) +  ce1*f1(i,k)*Gk(i,k)/Tt(i,k) )/rho(i,k)
      dimpl(i,k)  = dimpl(i,k)  + ce2*f2(i,k)/Tt(i,k)              ! note, ce2*f2*rho*epsilon/T/(rho*epsilon), set implicit and divided by density
    enddo
  enddo

end

!>************************************************************************************
!!
!!     Performes time integration with second order
!!     Adams-Bashforth scheme, i.e
!!
!!
!!     n+1     n
!!     dU     U     - U                               n
!!     ---- = ------------ = 1.5*( -ADV + DIFF + Force)     -
!!     dt        dt                                   n-1
!!     0.5*( -ADV + DIFF + Force)
!!
!!     This scheme is weakly instabel for pure advection,
!!     and therefore a very small amount of physical diffusion
!!     is necessary.
!!     The timestep is limited (see routine chkdt)
!!
!!************************************************************************************
subroutine advanceEpsilon(resE,Utmp,Wtmp,Rtmp,rho3,ftmp,rank)
  implicit none
      include 'param.f90'
      include 'common.f90'
  real*8 dnew(0:i1,0:k1),dimpl(0:i1,0:k1)
  real*8 Utmp(0:i1,0:k1),Wtmp(0:i1,0:k1),Rtmp(0:i1,0:k1),ftmp(imax,kmax)
  real*8 rho3(0:i1,0:k1)

  real*8     a  (imax)
  real*8     b  (imax)
  real*8     c  (imax)
  real*8     rhs(imax)

  real*8 resE

  resE  = 0.0
  dnew  = 0.0; dimpl = 0.0;

  call advecc(dnew,dimpl,eNew,utmp,wtmp,Ru,Rp,dru,dz,i1,k1,rank,periodic,.true.)
  call rhs_Epsilon(dnew,dimpl,Rtmp)    !new
  call diffEPS(dnew,eNew,ekm,ekmi,ekmk,ekmt,sigmae,Rtmp,Ru,Rp,dru,dz,rank,modifDiffTerm)

  if (centerBC.eq.-1) then
    do k=1,kmax
      do i=1,imax
                    
        a(i) = (ekmi(i-1,k)+0.5*(ekmt(i,k)+ekmt(i-1,k))/sigmae)/sqrt(0.5*(rho3(i-1,k)+rho3(i,k)))
        a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/Rtmp(i,k)/rho3(i,k)
    
        c(i) = (ekmi(i  ,k)+0.5*(ekmt(i,k)+ekmt(i+1,k))/sigmae)/sqrt(0.5*(rho3(i+1,k)+rho3(i,k)))
        c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/Rtmp(i,k)/rho3(i,k)
    
        b(i) = ((-a(i)-c(i))*(rho3(i,k)**1.5) + dimpl(i,k)  )/alphae
    
        a(i) = a(i)*(rho3(i-1,k)**1.5)
        c(i) = c(i)*(rho3(i+1,k)**1.5)
    
        rhs(i) = dnew(i,k) + (1-alphae)*b(i)*eNew(i,k)
      enddo
    
      i=1
      rhs(i) = dnew(i,k) - a(i)*eNew(i-1,k) + (1-alphae)*b(i)*eNew(i,k)
    
      i=imax
      rhs(i) = dnew(i,k) - c(i)*eNew(i+1,k) + (1-alphae)*b(i)*eNew(i,k)
    
      call matrixIdir(imax,a,b,c,rhs)
    
      do i=1,imax
        resE = resE + ((eNew(i,k) - rhs(i))/(eNew(i,k)+1.0e-20))**2.0
        eNew(i,k) = max(rhs(i), 1.0e-8)
    
      enddo
    enddo
  else
    do k=1,kmax
      do i=1,imax
                    
        a(i) = (ekmi(i-1,k)+0.5*(ekmt(i,k)+ekmt(i-1,k))/sigmae)/sqrt(0.5*(rho3(i-1,k)+rho3(i,k)))
        a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/Rtmp(i,k)/rho3(i,k)
    
        c(i) = (ekmi(i  ,k)+0.5*(ekmt(i,k)+ekmt(i+1,k))/sigmae)/sqrt(0.5*(rho3(i+1,k)+rho3(i,k)))
        c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/Rtmp(i,k)/rho3(i,k)
    
        b(i) = ((-a(i)-c(i))*(rho3(i,k)**1.5) + dimpl(i,k)  )/alphae
    
        a(i) = a(i)*(rho3(i-1,k)**1.5)
        c(i) = c(i)*(rho3(i+1,k)**1.5)
    
        rhs(i) = dnew(i,k) + (1-alphae)*b(i)*eNew(i,k)
      enddo
    
      i=1
      b(i)=b(i)+a(i)

      i=imax
      rhs(i) = dnew(i,k) - c(i)*eNew(i+1,k) + (1-alphae)*b(i)*eNew(i,k)
    
      call matrixIdir(imax,a,b,c,rhs)
    
      do i=1,imax
        resE = resE + ((eNew(i,k) - rhs(i))/(eNew(i,k)+1.0e-20))**2.0
        eNew(i,k) = max(rhs(i), 1.0e-8)
    
      enddo
    enddo
  endif
end
!>************************************************************************************
!!
!!     Performes time integration with second order
!!     Adams-Bashforth scheme, i.e
!!
!!
!!     n+1     n
!!     dU     U     - U                               n
!!     ---- = ------------ = 1.5*( -ADV + DIFF + Force)     -
!!     dt        dt                                   n-1
!!     0.5*( -ADV + DIFF + Force)
!!
!!     This scheme is weakly instabel for pure advection,
!!     and therefore a very small amount of physical diffusion
!!     is necessary.
!!     The timestep is limited (see routine chkdt)
!!
!!************************************************************************************
subroutine advanceK(resK,Utmp,Wtmp,Rtmp,rho3,ftmp,mrank)
  implicit none
      include 'param.f90'
      include 'common.f90'
  real*8 dnew(0:i1,0:k1),dimpl(0:i1,0:k1)
  real*8 Utmp(0:i1,0:k1),Wtmp(0:i1,0:k1),Rtmp(0:i1,0:k1),ftmp(imax,kmax)
  real*8 rho3(0:i1,0:k1)
 
  real*8     a  (imax)
  real*8     b  (imax)
  real*8     c  (imax)
  real*8     rhs(imax)

  real*8 resK
  integer mrank

  resK  = 0.0
  dnew  = 0.0; dimpl = 0.0;


  call advecc(dnew,dimpl,kNew,utmp,wtmp,Ru,Rp,dru,dz,i1,k1,mrank,periodic,.true.)
  call rhs_K(dnew,dimpl,kNew,eNew,Rtmp)    !new
  call diffc(dnew,kNew,ekm,ekmi,ekmk,ekmt,sigmak,Rtmp,Ru,Rp,dru,dz,mrank,modifDiffTerm)

  if ((modifDiffTerm == 0) .or. (modifDiffTerm == 1)) then
    do k=1,kmax
      do i=1,imax

        a(i) = (ekmi(i-1,k)+0.5*(ekmt(i,k)+ekmt(i-1,k))/sigmak)/((0.5*(rho3(i-1,k)+rho3(i,k)))**0.5)
        a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/Rtmp(i,k)/(rho3(i,k)**0.5)

        c(i) = (ekmi(i  ,k)+0.5*(ekmt(i,k)+ekmt(i+1,k))/sigmak)/((0.5*(rho3(i+1,k)+rho3(i,k)))**0.5)
        c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/Rtmp(i,k)/(rho3(i,k)**0.5)

        b(i) = (rho3(i,k)*(-a(i)-c(i)) + dimpl(i,k))/alphak

        a(i) = a(i)*rho3(i-1,k)
        c(i) = c(i)*rho3(i+1,k)

        rhs(i) = dnew(i,k) + (1-alphak)*b(i)*kNew(i,k)
      enddo

      i=1
      b(i) = b(i) + centerBC*a(i)
             
      i=imax
      b(i) = b(i) - (c(i) /alphak)
      !b(i) = (rho3(i,k)*(-(a(i)/rho3(i-1,k))-(c(i)/rho3(i+1,k))) - c(i) + dimpl(i,k) )/alphak
      !b(i) = (rho3(i,k)*(-a(i)-c(i)) - rho3(i+1,k)*c(i) + dimpl(i,k) )/alphak
      rhs(i) = dnew(i,k) + (1-alphak)*b(i)*kNew(i,k)

      call matrixIdir(imax,a,b,c,rhs)

      do i=1,imax
        resK = resK + ((kNew(i,k) - rhs(i))/(kNew(i,k)+1.0e-20))**2.0
        kNew(i,k) = max(rhs(i), 1.0e-8)
      enddo
    enddo
  else if (modifDiffTerm == 2) then
    do k=1,kmax
      do i=1,imax
        a(i) = (ekmi(i-1,k)+0.5*(ekmt(i,k)+ekmt(i-1,k))/sigmak)/(0.5*(rho3(i-1,k)+rho3(i,k)))
        a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/Rtmp(i,k)

        c(i) = (ekmi(i  ,k)+0.5*(ekmt(i,k)+ekmt(i+1,k))/sigmak)/(0.5*(rho3(i+1,k)+rho3(i,k)))
        c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/Rtmp(i,k)

        b(i) = (rho3(i,k)*(-a(i)-c(i)) + dimpl(i,k))/alphak

        a(i) = a(i)*rho3(i-1,k)
        c(i) = c(i)*rho3(i+1,k)

        rhs(i) = dnew(i,k) + (1-alphak)*b(i)*kNew(i,k)
      enddo

      i=1
      b(i) = b(i) + centerBC*a(i)
             
      i=imax
      b(i) = b(i) - (c(i) /alphak)
      !b(i) = (rho3(i,k)*(-(a(i)/rho3(i-1,k))-(c(i)/rho3(i+1,k))) - c(i) + dimpl(i,k) )/alphak
      !b(i) = (rho3(i,k)*(-a(i)-c(i)) - rho3(i+1,k)*c(i) + dimpl(i,k) )/alphak
      rhs(i) = dnew(i,k) + (1-alphak)*b(i)*kNew(i,k)

      call matrixIdir(imax,a,b,c,rhs)

      do i=1,imax
        resK = resK + ((kNew(i,k) - rhs(i))/(kNew(i,k)+1.0e-20))**2.0
        kNew(i,k) = max(rhs(i), 1.0e-8)
      enddo
    enddo
  endif

end


!>********************************************************************
!! diffusion term for epsilon in the z-direction, set as a source term...
!!********************************************************************
subroutine diffEPS(putout,putin,ek,eki,ekk,ekmt,sigma,rho,Ru,Rp,dru,dz,rank1,diffVersion)
  implicit none
      include 'param.f90'
  integer   km,kp, rank1,diffVersion
  real*8     putout(0:i1,0:k1),putin(0:i1,0:k1), &
    rho(0:i1,0:k1),ek(0:i1,0:k1),eki(0:i1,0:k1),ekk(0:i1,0:k1), &
    ekmt(0:i1,0:k1),dru(0:i1),dz,Ru(0:i1),Rp(0:i1),sigma, difcp, difcm

  if ((diffVersion == 1) .or. (diffVersion == 2)) then       ! Inverse SLS  and Aupoix
    do k=1,k1-1
      kp=k+1
      km=k-1
      do i=1,i1-1
        difcp = (ekk(i,k ) + 0.5*(ekmt(i,k)+ekmt(i,kp))/sigma)/sqrt(0.5*(rho(i,k)+rho(i,kp)))
        difcm = (ekk(i,km) + 0.5*(ekmt(i,k)+ekmt(i,km))/sigma)/sqrt(0.5*(rho(i,k)+rho(i,km)))
        putout(i,k) = putout(i,k) + 1.0/rho(i,k)/rho(i,k)*( &
          (     difcp * ((rho(i,kp)**1.5)*putin(i,kp)-(rho(i,k )**1.5)*putin(i,k )) &
          -difcm * ((rho(i,k )**1.5)*putin(i,k )-(rho(i,km)**1.5)*putin(i,km))  )/(dz*dz)   )
      enddo
    enddo
  else                               ! Standard
    do k=1,k1-1
      kp=k+1
      km=k-1
      do i=1,i1-1
        putout(i,k) = putout(i,k) + 1.0/rho(i,k)*( &
          ( (ekk(i,k ) + 0.5*(ekmt(i,k)+ekmt(i,kp))/sigma)*(putin(i,kp)-putin(i,k )) &
          -(ekk(i,km) + 0.5*(ekmt(i,k)+ekmt(i,km))/sigma)*(putin(i,k )-putin(i,km))   )/(dz*dz)   )
      enddo
    enddo
  endif


end



