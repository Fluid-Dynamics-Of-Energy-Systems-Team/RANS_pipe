!>******************************************************************************************
!>******************************************************************************************
!>******************************************************************************************
!! V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F 
!! V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F 
!! V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F 
!! V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F V2F 
!!              V2F =  parameter (turbmod           = 3) 
!!
!>******************************************************************************************
!!      V2F routine to estimate the eddy viscosity
!!******************************************************************************************
      subroutine calculate_mut_VF(U,W,ekmetmp,ekmttmp,ekmtin,rank,step)

      implicit none
      include 'param.txt'
      include 'common.txt'
      integer  im,ip,km,kp,step,rank
      real*8   tauwLoc, tauw(0:k1) 
      real*8, dimension(0:i1,0:k1) :: U,W,ekmetmp,ekmttmp,Tt,Srsq
      real*8, dimension(0:i1) :: ekmtb,ekmtf,ekmtin
      real*8  StR



      sigmat = 0.9


      do k=1,kmax
         km=k-1
         kp=k+1

         tauw(k) = ekmi(imax,k)*0.5*(W(imax,km)+W(imax,k))/wallDist(imax)

         do i=1,imax
            im=i-1
            ip=i+1

            yp(i,k)     = sqrt(rNew(i,k))/ekm(i,k)*(wallDist(i))*tauw(k)**0.5           ! ystar
  !          yp(i,k)     = sqrt(rNew(imax,k))/ekm(imax,k)*(wallDist(i))*tauw(k)**0.5    ! yplus
  !          yp(i,:)     = (wallDist(i))*Re*(1/Re*(Win(imax)/(wallDist(imax))))**0.5        ! yplus
            ReTauS(i,k) = 0.5*sqrt(rNew(i,k))/ekm(i,k)*tauw(k)**0.5
            Ret(i,k)    = rNew(i,k)*(kNew(i,k)**2.)/(ekm(i,k)*eNew(i,k))        ! not sure if r2 or r

            sigmak = 1.0
            sigmae = 1.3
            cmu    = 0.22
            ce1    = 1.4
            ce2    = 1.9
           
            StR= (2.*(((W(i,k)-W(i,km))/dz)**2. +
 &                ((U(i,k)-U(im,k))/dru(i))**2. +
 &                ((U(i,k)+U(im,k))/(2.*Rp(i)))**2.) +  
 &                (((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4.
 &                -(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/dru(i)
 &                +((U(i,kp) +U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz)  )**2.)
    
!               Srsq(i,k) = Pk(i,k)*rNew(i,k)/(2.*ekmt(i,k))

            Srsq(i,k) = Str*rNew(i,k)*0.5
            Tt(i,k)   = max(kNew(i,k)/eNew(i,k),6.0*(ekm(i,k)/(rNew(i,k)*eNew(i,k)))**0.5)

            !extras
            !Tt(i,k)   = max(Tt(i,k), 1.0e-8)
            !Tt(i,k)   = min(Tt(i,k),0.6*kNew(i,k)/(3.**0.5*v2New(i,k)*cmu*(2.*Srsq(i,k))**0.5))

            fmu(i,k) = v2New(i,k)*Tt(i,k)/(kNew(i,k)**2./eNew(i,k))
            f1(i,k)  = 1.0 + 0.045*(kNew(i,k)/v2New(i,k))**0.5
            f2(i,k)  = 1.0
            dterm(i,k) = 0.0
            eterm(i,k) = 0.0
            ekmttmp(i,k) = min(1.,rNew(i,k)*cmu*v2New(i,k)*Tt(i,k))
!               ekmttmp(i,k) = min(max(ekmttmp(i,k),1.0e-8), 1.0)
          
         enddo
      enddo

      end

!>******************************************************************************************
!!      V2F prodis subroutine which calculates the production term of the turbulent scalar equation
!>******************************************************************************************
      subroutine prodis_VF(putout,dimpl,putink,putine,putinv2,putinf,U,W,T,rho,scl)
      implicit none
      include 'param.txt'
      include 'common.txt'

      integer im,ip,km,kp,ib,ie,kb,ke !< integers
      real*8, dimension(0:i1,0:k1) :: putout,U,W,T,rho,div,putink,putine,putinv2,Tt,Srsq,dimpl
      real*8, dimension(imax,kmax) :: putinf,putinftmp
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
            call calculate_Pk(tmpPk, tmpDiv, U,W,T,rho,i,im,ip,k,km,kp)
            Pk(i,k) = tmpPk
            div(i,k)= tmpDiv

            ! time scale for v2f model
            StR = (2.*(((W(i,k)-W(i,km))/dz)**2. +
     &                ((U(i,k)-U(im,k))/(Ru(i)-Ru(im)))**2. +
     &                ((U(i,k)+U(im,k))/(2.*Rp(i)))**2.) +
     &                (((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4.
     &                -(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/(Ru(i)-Ru(im))
     &                +((U(i,kp) +U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz)  )**2.)  
    
!               Srsq(i,k) = Pk(i,k)*rho(i,k)/(2.*ekmt(i,k))
            Srsq(i,k) = Str*rho(i,k)*0.5
            Tt(i,k)   = max(putink(i,k)/putine(i,k), 6.0*(ekm(i,k)/(rho(i,k)*putine(i,k)))**0.5)
            Tt(i,k)   = max(Tt(i,k), 1.0e-8)
            Tt(i,k)   = min(Tt(i,k),0.6*putink(i,k)/(3.**0.5*putinv2(i,k)*cmu*(2.*Srsq(i,k))**0.5))
               
               
               ! Bouyancy prodution with a different time scale
            Gk(i,k)=-ctheta*beta(i,k)*Fr_1*Tt(i,k)
     &             *  (ekmt(i,k)*(((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4.-(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/(Ru(i)-Ru(im))
     &                          +((U(i,kp)+U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz) )*
     &                                                                             (T(ip,k)-T(im,k))/(Rp(ip)-Rp(im))  )
     &               +(2.*ekmt(i,k)*((W(i,k)-W(i,km))/dz-2./3.*(rho(i,k)*putink(i,k)))*(T(i,kp)-T(i,km))/(2.*dz)
     &               )


            Gk(i,k) = Gk(i,k) + ctheta*beta(i,k)*Fr_1*Tt(i,k)*2./3.*ekmt(i,k)*div(i,k)*(T(i,kp)-T(i,km))/(2.*dz)


            if (scl.eq.0) then
               !k equation
               putout(i,k) = putout(i,k) + ( Pk(i,k) + Gk(i,k) )/rho(i,k)
               dimpl(i,k)  = dimpl(i,k) + putine(i,k)/putink(i,k)       ! note, rho*epsilon/(rho*k), set implicit and divided by density

            elseif (scl.eq.1) then
               !epsilon equation
               putout(i,k) = putout(i,k) +(ce1*f1(i,k)*Pk(i,k)/Tt(i,k) +  ce1*f1(i,k)*Gk(i,k)/Tt(i,k) )/rho(i,k)
               dimpl(i,k)  = dimpl(i,k)  + ce2*f2(i,k)/Tt(i,k)              ! note, ce2*f2*rho*epsilon/T/(rho*epsilon), set implicit and divided by density

            elseif (scl.eq.2) then
               !v'2 equation
               putout(i,k) = putout(i,k) + putink(i,k)*putinf(i,k)       ! note, source is rho*k*f/rho
               dimpl(i,k)  = dimpl(i,k)  + 6.*putine(i,k)/putink(i,k)    ! note, 6*rho*v'2*epsilon/k/(rho*v'2), set implicit and divided by density

            

            endif
         enddo
      enddo

      end


!>******************************************************************************************
!!      VF advancing the turbulence scalars of this model: k and epsilon and v2
!!******************************************************************************************
      subroutine advanceScalar_VF(resK,resE,resV2,Utmp,Wtmp,Rtmp,ftmp,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      real*8 dnew(0:i1,0:k1),dimpl(0:i1,0:k1)
      real*8 Utmp(0:i1,0:k1),Wtmp(0:i1,0:k1),Rtmp(0:i1,0:k1),ftmp(imax,kmax)
      real*8 resK, resE, resV2
      real*8 rho3(0:i1,0:k1)

      real*8     a  (imax)
      real*8     b  (imax)
      real*8     c  (imax)
      real*8     rhs(imax)

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

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! advanceV2
      resV2 = 0.0
      dnew  = 0.0; dimpl = 0.0;
      scl=2.0

      call advecc(dnew,dimpl,v2New,utmp,wtmp,Ru,Rp,dru,dz,i1,k1,rank,periodic,.true.)
      call prodis_VF(dnew,dimpl,kNew,eNew,v2New,ftmp,Utmp,Wtmp,temp,Rtmp,scl)
      call diffc(dnew,v2New,ekm,ekmi,ekmk,ekmt,sigmak,Rtmp,Ru,Rp,dru,dz,rank,modifDiffTerm)

      if ((modifDiffTerm == 0) .or. (modifDiffTerm == 1)) then
        do k=1,kmax
           do i=1,imax
            
              a(i) = (ekmi(i-1,k)+0.5*(ekmt(i,k)+ekmt(i-1,k))/sigmak)/((0.5*(rho3(i-1,k)+rho3(i,k)))**0.5)
              a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/Rtmp(i,k)/(rho3(i,k)**0.5)
              c(i) = (ekmi(i  ,k)+0.5*(ekmt(i,k)+ekmt(i+1,k))/sigmak)/((0.5*(rho3(i+1,k)+rho3(i,k)))**0.5)
              c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/Rtmp(i,k)/(rho3(i,k)**0.5)

              b(i) = (rho3(i,k)*(-a(i)-c(i)) + dimpl(i,k))/alphav2
              a(i) = a(i)*rho3(i-1,k)
              c(i) = c(i)*rho3(i+1,k)

               rhs(i) =  dnew(i,k) + (1-alphav2)*b(i)*v2New(i,k)
            enddo

            b(1)=b(1)+a(1)
            i=imax
            b(i) = b(i) - (c(i) /alphav2)
            !b(i) = (rho3(i,k)*(-(a(i)/rho3(i-1,k))-(c(i)/rho3(i+1,k)))  - c(i) + dimpl(i,k) )/alphav2
            !b(i) = (rho3(i,k)*(-a(i)-c(i)) - rho3(i+1,k)*c(i) + dimpl(i,k) )/alphav2
            rhs(i) = dnew(i,k) + (1-alphav2)*b(i)*v2New(i,k)


            call matrixIdir(imax,a,b,c,rhs)

            do i=1,imax
               resV2 = resV2 + ((v2New(i,k) - rhs(i))/(v2New(i,k)+1.0e-20))**2.0
               v2New(i,k) = min(2.0/3.0*kNew(i,k), max(rhs(i), 1.0e-8))
            enddo
        enddo

      else if (modifDiffTerm == 2) then
         do k=1,kmax
            do i=1,imax
             
               a(i) = (ekmi(i-1,k)+0.5*(ekmt(i,k)+ekmt(i-1,k))/sigmak)/(0.5*(rho3(i-1,k)+rho3(i,k)))
               a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/Rtmp(i,k)
               c(i) = (ekmi(i  ,k)+0.5*(ekmt(i,k)+ekmt(i+1,k))/sigmak)/(0.5*(rho3(i+1,k)+rho3(i,k)))
               c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/Rtmp(i,k)


               b(i) = (rho3(i,k)*(-a(i)-c(i)) + dimpl(i,k))/alphav2
               a(i) = a(i)*rho3(i-1,k)
               c(i) = c(i)*rho3(i+1,k)

               rhs(i) =  dnew(i,k) + (1-alphav2)*b(i)*v2New(i,k)
             enddo

             b(1)=b(1)+a(1)
             i=imax
             b(i) = b(i) - (c(i) /alphav2)
             !b(i) = (rho3(i,k)*(-(a(i)/rho3(i-1,k))-(c(i)/rho3(i+1,k)))  - c(i) + dimpl(i,k) )/alphav2
             !b(i) = (rho3(i,k)*(-a(i)-c(i)) - rho3(i+1,k)*c(i) + dimpl(i,k) )/alphav2
             rhs(i) = dnew(i,k) + (1-alphav2)*b(i)*v2New(i,k)


             call matrixIdir(imax,a,b,c,rhs)

             do i=1,imax
                resV2 = resV2 + ((v2New(i,k) - rhs(i))/(v2New(i,k)+1.0e-20))**2.0
                v2New(i,k) = min(2.0/3.0*kNew(i,k), max(rhs(i), 1.0e-8))
             enddo
        enddo
      endif


      end
