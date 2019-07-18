!>******************************************************************************************
!>******************************************************************************************
!>******************************************************************************************
!! MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK
!! MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK
!! MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK MK
!!              MK =  parameter (turbmod           = 1) 
!!
!>******************************************************************************************
!!      MK routine to estimate the eddy viscosity
!!******************************************************************************************
      subroutine calculate_mut_MK(U,W,ekmetmp,ekmttmp,ekmtin,step)
      implicit none
      include 'param.txt'
      include 'common.txt'

      integer  im,ip,km,kp,step
      real*8   tauw(0:k1) 
      real*8, dimension(0:i1,0:k1) :: U,W,ekmetmp,ekmttmp,Tt
      real*8, dimension(0:i1) :: ekmtb,ekmtf,ekmtin

      sigmak       = 1.4
      sigmae       = 1.3
      cmu          = 0.09
      ce1          = 1.4
      ce2          = 1.8

      do k=1,kmax
         km=k-1
         kp=k+1

         tauw(k) = ekmi(imax,k)*0.5*(W(imax,km)+W(imax,k))/wallDist(imax)

         do i=1,imax
            im=i-1
            ip=i+1

            yp(i,k)     = sqrt(rNew(i,k))/ekm(i,k)*(wallDist(i))*tauw(k)**0.5           ! ystar
  !          yp(i,k)     = sqrt(rNew(imax,k))/ekm(imax,k)*(wallDist(i))*tauw(k)**0.5    ! yplus
  !          yp(i,:)     = (wallDist(i))*Re*(1/Re*(Win(imax)/(wallDist(imax))))**0.5    ! yplus
            ReTauS(i,k) = 0.5*sqrt(rNew(i,k))/ekm(i,k)*tauw(k)**0.5
            Ret(i,k)    = rNew(i,k)*(kNew(i,k)**2.)/(ekm(i,k)*eNew(i,k))        ! not sure if r2 or r


            fmu(i,k)     = (1.-exp(-yp(i,k)/70.))*(1+3.45/Ret(i,k)**0.5)
            f1(i,k)      = 1.
            f2(i,k)      = (1.-2./9.*exp(-(Ret(i,k)/6.)**2.))*(1.-exp(-yp(i,k)/5.))**2.0
            dterm(i,k)   = 0.
            eterm(i,k)   = 0.
            ekmttmp(i,k) = min(1.,rNew(i,k)*cmu*fmu(i,k)*kNew(i,k)**2./(eNew(i,k)))


         enddo
      enddo

      end


!>******************************************************************************************
!!      MK routine to estimate the equation terms
!!******************************************************************************************
      subroutine prodis_MK(putout,dimpl,putink,putine,U,W,T,rho,scl)
      implicit none
      include 'param.txt'
      include 'common.txt'

      integer im,ip,jm,jp,km,kp,ib,ie,kb,ke !< integers
      real*8, dimension(0:i1,0:k1) :: putout,U,W,T,rho,div,putink,putine,Tt,dimpl
      real*8  scl


      ib = 1
      ie = i1-1

      kb = 1
      ke = k1-1

      do k=kb,ke
         kp=k+1
         km=k-1
         do i=ib,ie
            ip=i+1
            im=i-1
            
            ! Production of turbulent kinetic energy
            Pk(i,k) = ekmt(i,k)*(
     &         2.*(((W(i,k)-W(i,km))/dz)**2. +
     &             ((U(i,k)-U(im,k))/dRu(i))**2. +
     &             ((U(i,k)+U(im,k))/(2.*Rp(i)))**2.) +
     &            (((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4.
     &             -(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/dRu(i)
     &            +((U(i,kp) +U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz)
     &              )**2.)

            div(i,k) =(Ru(i)*U(i,k)-Ru(im)*U(im,k))/(Rp(i)*dru(i))
     &              +(      W(i,k) -      W(i,km))/dz

            Pk(i,k) = Pk(i,k) - 2./3.*(rho(i,k)*putink(i,k)+ekmt(i,k)*(div(i,k)))*(div(i,k))

            ! turbulent time scale
            Tt(i,k)=putink(i,k)/putine(i,k)

            ! Bouyancy prodution
            Gk(i,k)=-ctheta*beta(i,k)*Fr_1*Tt(i,k)
     &           *  (ekmt(i,k)*(((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4.-(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/dRu(i)
     &                         +((U(i,kp)+U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz) )*
     &                                                                              (T(ip,k)-T(im,k))/(dRp(i)+dRp(im))  )
     &           +(2.*ekmt(i,k)*((W(i,k)-W(i,km))/dz-2./3.*(rho(i,k)*putink(i,k)))*(T(i,kp)-T(i,km))/(2.*dz)
     &           )

            Gk(i,k) = Gk(i,k) + ctheta*beta(i,k)*Fr_1*Tt(i,k)*2./3.*ekmt(i,k)*div(i,k)*(T(i,kp)-T(i,km))/(2.*dz)

            if (scl.eq.0) then
               !k equation
               putout(i,k) = putout(i,k) + ( Pk(i,k) + Gk(i,k) )/rho(i,k)
               dimpl(i,k)  = dimpl(i,k) + putine(i,k)/putink(i,k)       ! note, rho*epsilon/(rho*k), set implicit and divided by density

            elseif (scl.eq.1) then
               !epsilon equation
               putout(i,k) = putout(i,k) +(ce1*f1(i,k)*Pk(i,k)/Tt(i,k) +  ce1*f1(i,k)*Gk(i,k)/Tt(i,k) )/rho(i,k)
               dimpl(i,k)  = dimpl(i,k)  + ce2*f2(i,k)/Tt(i,k)              ! note, ce2*f2*rho*epsilon/T/(rho*epsilon), set implicit and divided by density

            endif

         enddo
      enddo

      end



!>******************************************************************************************
!!      MK advancing the turbulence scalars of this model: k and epsilon
!!******************************************************************************************
      subroutine advanceScalar_MK(resK,resE,Utmp,Wtmp,Rtmp,ftmp,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      real*8 Utmp(0:i1,0:k1),Wtmp(0:i1,0:k1),Rtmp(0:i1,0:k1),ftmp(imax,kmax)
      real*8 resK, resE
      real*8 rho3(0:i1,0:k1)

      
      ! modified turb. model
      !    modifDiffTerm = 1, our modification
      !    modifDiffTerm = 2, Aupoix modification
      if ((modifDiffTerm == 1) .or. (modifDiffTerm == 2)) then
         rho3 = Rtmp
      else
         rho3 = 1.0
      endif

      call advanceEpsilon(resE,Utmp,Wtmp,Rtmp,rho3,ftmp,rank)
      call advanceK(resK,Utmp,Wtmp,Rtmp,rho3,ftmp,rank)


      end

