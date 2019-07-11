!*************************************************************************************
!    to calculate the turbulent production
!*************************************************************************************
      subroutine calculate_Pk(tmpPk, tmpDiv, U,W,T,rho,i,im,ip,k,km,kp)
      implicit none
      include 'param.txt'
      include 'common.txt'

      integer im,ip,km,kp,ib,ie,kb,ke !< integers
      real*8, dimension(0:i1,0:k1) :: U,W,T,rho,div

      ! Production of turbulent kinetic energy
      tmpPk = ekmt(i,k)*(
     &         2.*(((W(i,k)-W(i,km))/dz)**2. +
     &             ((U(i,k)-U(im,k))/(Ru(i)-Ru(im)))**2. +
     &             ((U(i,k)+U(im,k))/(2.*Rp(i)))**2.) +
     &            (((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4.
     &             -(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/dRu(i)
     &            +((U(i,kp) +U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz)
     &              )**2.)

      tmpDiv)=(Ru(i)*U(i,k)-Ru(im)*U(im,k))/(Rp(i)*dru(i))
     &              +(      W(i,k) -      W(i,km))/dz

      tmpPk = Pk(i,k) - 2./3.*(rho(i,k)*putink(i,k)+ekmt(i,k)*(div(i,k)))*(div(i,k))

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
      include 'param.txt'
      include 'common.txt'
      real*8 dnew(0:i1,0:k1),dimpl(0:i1,0:k1)
      real*8 Utmp(0:i1,0:k1),Wtmp(0:i1,0:k1),Rtmp(0:i1,0:k1),ftmp(imax,kmax)
      real*8 rho3(0:i1,0:k1)

      real*8     a  (imax)
      real*8     b  (imax)
      real*8     c  (imax)
      real*8     rhs(imax)

      real*8 scl
      real*8 resE

      resE  = 0.0
      dnew  = 0.0; dimpl = 0.0;
      scl=1.0

      call advecc(dnew,dimpl,eNew,utmp,wtmp,Ru,Rp,dru,dz,i1,k1,rank,periodic,.true.)     

      if (turbmod.eq.1) then 
         call prodis_MK(dnew,dimpl,kNew,eNew,Utmp,Wtmp,temp,Rtmp,scl)
      elseif (turbmod.eq.3) then
         call prodis_VF(dnew,dimpl,kNew,eNew,v2New,ftmp,Utmp,Wtmp,temp,Rtmp,scl)  
      endif

      call diffEPS(dnew,eNew,ekm,ekmi,ekmk,ekmt,sigmae,Rtmp,Ru,Rp,dru,dz,rank,modifDiffTerm)

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

         b(1)=b(1)+a(1)
         i=imax
         rhs(i) = dnew(i,k) - c(i)*eNew(i1,k) + (1-alphae)*b(i)*eNew(i,k)

         call matrixIdir(imax,a,b,c,rhs)

         do i=1,imax
            resE = resE + ((eNew(i,k) - rhs(i))/(eNew(i,k)+1.0e-20))**2.0
            eNew(i,k) = max(rhs(i), 1.0e-8)
         enddo
      enddo

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
      subroutine advanceK(resK,Utmp,Wtmp,Rtmp,rho3,ftmp,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      real*8 dnew(0:i1,0:k1),dimpl(0:i1,0:k1)
      real*8 Utmp(0:i1,0:k1),Wtmp(0:i1,0:k1),Rtmp(0:i1,0:k1),ftmp(imax,kmax)
      real*8 rho3(0:i1,0:k1)

      real*8     a  (imax)
      real*8     b  (imax)
      real*8     c  (imax)
      real*8     rhs(imax)

      real*8 scl
      real*8 resE

      resK  = 0.0
      dnew  = 0.0; dimpl = 0.0;
      scl=0.0


      call advecc(dnew,dimpl,kNew,utmp,wtmp,Ru,Rp,dru,dz,i1,k1,rank,periodic,.true.)

      if (turbmod.eq.1) then 
         call prodis_MK(dnew,dimpl,kNew,eNew,Utmp,Wtmp,temp,Rtmp,scl)
      elseif (turbmod.eq.3) then
         call prodis_VF(dnew,dimpl,kNew,eNew,v2New,ftmp,Utmp,Wtmp,temp,Rtmp,scl)
!      elseif (turbmod.eq.5) then
!         call prodis_SST(dnew,dimpl,kNew,eNew,Utmp,Wtmp,temp,Rtmp,scl)   
      endif

      call diffc(dnew,kNew,ekm,ekmi,ekmk,ekmt,sigmak,Rtmp,Ru,Rp,dru,dz,rank,modifDiffTerm)

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
           b(1) = b(1) + a(1)
             
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
           b(1) = b(1) + a(1)
             
           i=imax
           b(i) = b(i) - (c(i) /alphak)
           !b(i) = (rho3(i,k)*(-(a(i)/rho3(i-1,k))-(c(i)/rho3(i+1,k))) - c(i) + dimpl(i,k) )/alphak
           !b(i) = (rho3(i,k)*(-a(i)-c(i)) - rho3(i+1,k)*c(i) + dimpl(i,k) )/alphak
           rhs(i) = dnew(i,k) + (1-alphak)*b(i)*kNew(i,k)

           call matrixIdir(imax,a,b,c,rhs)

           do i=1,imax
              resK = resK + ((kNew(i,k) - rhs(i))/(kNew(i,k)+1.0e-20))**2.0
              kNew(i,k) = max(rhs(i), 1.0e-8)
           end

      endif

      end


      !>********************************************************************
      !! diffusion term for epsilon in the z-direction, set as a source term...
      !!********************************************************************
      subroutine diffEPS(putout,putin,ek,eki,ekk,ekmt,sigma,rho,Ru,Rp,dru,dz,rank1,diffVersion)
      implicit none
      include 'param.txt'
      integer   km,kp, rank1,diffVersion
      real*8     putout(0:i1,0:k1),putin(0:i1,0:k1),
     &      rho(0:i1,0:k1),ek(0:i1,0:k1),eki(0:i1,0:k1),ekk(0:i1,0:k1),
     &      ekmt(0:i1,0:k1),dru(0:i1),dz,Ru(0:i1),Rp(0:i1),sigma, difcp, difcm

         if ((diffVersion == 1) .or. (diffVersion == 2)) then       ! Inverse SLS  and Aupoix
            do k=1,k1-1
               kp=k+1
               km=k-1
               do i=1,i1-1
                  difcp = (ekk(i,k ) + 0.5*(ekmt(i,k)+ekmt(i,kp))/sigma)/sqrt(0.5*(rho(i,k)+rho(i,kp)))
                  difcm = (ekk(i,km) + 0.5*(ekmt(i,k)+ekmt(i,km))/sigma)/sqrt(0.5*(rho(i,k)+rho(i,km)))
                  putout(i,k) = putout(i,k) + 1.0/rho(i,k)/rho(i,k)*(
     3                     (     difcp * ((rho(i,kp)**1.5)*putin(i,kp)-(rho(i,k )**1.5)*putin(i,k ))
     3                          -difcm * ((rho(i,k )**1.5)*putin(i,k )-(rho(i,km)**1.5)*putin(i,km))  )/(dz*dz)   )
               enddo
            enddo
         else                               ! Standard
            do k=1,k1-1
               kp=k+1
               km=k-1
               do i=1,i1-1
                     putout(i,k) = putout(i,k) + 1.0/rho(i,k)*(
     3                    ( (ekk(i,k ) + 0.5*(ekmt(i,k)+ekmt(i,kp))/sigma)*(putin(i,kp)-putin(i,k ))
     3                     -(ekk(i,km) + 0.5*(ekmt(i,k)+ekmt(i,km))/sigma)*(putin(i,k )-putin(i,km))   )/(dz*dz)   )
               enddo
            enddo
         endif


      end



