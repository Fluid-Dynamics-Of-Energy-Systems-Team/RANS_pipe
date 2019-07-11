!>******************************************************************************************
!>******************************************************************************************
!>******************************************************************************************
!! SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST 
!! SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST 
!! SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST 
!! SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST SST 
!!              SST =  parameter (turbmod           = 5) 
!!
!>******************************************************************************************
!!      SST routine to estimate the eddy viscosity
!!******************************************************************************************
      subroutine calculate_mut_SST(U,W,ekmetmp,ekmttmp,ekmtin,step)

      implicit none
      include 'param.txt'
      include 'common.txt'
      integer  im,ip,km,kp,step
      real*8   tauwLoc, tauw(0:k1) 
      real*8, dimension(0:i1,0:k1) :: U,W,ekmetmp,ekmttmp
      real*8, dimension(0:i1) :: ekmtb,ekmtf,ekmtin
      real*8  sigma_om2,betaStar,gradkom,gamma1,gamma2,gamma3,gammaSST,zetaSST,StR, wallD



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

            !constants
            sigma_om2 = 0.856
            betaStar  = 0.09
            wallD     = wallDist(i)

            ! Vorticity rate
            StR = ( ( -( (W(ip,km)+W(ip,k)+W(i,km)+W(i ,k))/4.-(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/dru(i)
     &                +( (U(i,kp) +U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz))**2.)

            StR = StR**0.5
               
            gradkom = ((kNew(ip,k) - kNew(im,k))/(dRp(i)+dRp(im))) * ((omnew(ip,k) - omnew(im,k))/(dRp(i)+dRp(im)))
     &                 +((kNew(i,kp) - kNew(i,km))/(2.0*dz))         * ((omnew(i,kp) - omnew(i,km))/(2.0*dz))

            cdKOM(i,k) = 2.0*sigma_om2*rNew(i,k)/omnew(i,k)*gradkom;
              
            gamma1   = 500.0*ekm(i,k)/(rNew(i,k)*omNew(i,k)*wallD**2.0);
            gamma2   = 4.0*sigma_om2*rNew(i,k)*kNew(i,k)/(wallD*wallD*max(cdKOM(i,k), 1.0e-20));
            gamma3   = (kNew(i,k)**0.5)/(betaStar*omNew(i,k)*wallD)

            gammaSST = min(max(gamma1, gamma3), gamma2)
            bF1(i,k) = tanh(gammaSST**4.0)

            gammaSST = max(2.0*gamma3, gamma1)
            bF2(i,k) = tanh(gammaSST**2.0)

            zetaSST  = max(0.31*omNew(i,k), bF2(i,k)*StR)
            ekmttmp(i,k) = rNew(i,k)*kNew(i,k)/omNew(i,k)


         enddo
      enddo

      ! Boundary condition for bF1 of sst model
      bF1(i1,:) =  bF1(imax,:)
      bF1(0,:)  =  bF1(1,:)

      call shiftf(bF1,ekmtf,rank)
      call shiftb(bF1,ekmtb,rank)
      bF1(:,0)  = ekmtf(:)
      bF1(:,k1) = ekmtb(:)

      if ((periodic.ne.1).and.(rank.eq.0)) then
         bF1(:,0) = bF1(:,1)  ! ATTENTION
      endif

      if ((periodic.ne.1).and.(rank.eq.px-1)) then
         bF1(:,k1) = 2.*bF1(:,kmax)-bF1(:,kmax-1)
      endif


      end
!>******************************************************************************************
!!      SST prodis subroutine which calculates the production term of the turbulent scalar equation
!>******************************************************************************************
      subroutine prodis_SST(putout,dimpl,putink,U,W,T,rho,scl)

      implicit none
      include 'param.txt'
      include 'common.txt'

      integer im,ip,jm,jp,km,kp,ib,ie,jb,je,kb,ke !< integers
      real*8, dimension(0:i1,0:k1) :: putout,U,W,T,rho,div,putink,putine,Tt,dimpl
      real*8  scl
      real*8  StR
      real*8  sigma_om1,sigma_om2,beta_1,beta_2,betaStar,alfa_1,alfa_2,alfaSST,betaSST, GtR

      real*8 tmpPk, tmpDiv ! COMMENT RENE

      sigma_om1 = 0.5
      sigma_om2 = 0.856
      beta_1    = 0.075
      beta_2    = 0.0828
      betaStar  = 0.09
      alfa_1    = beta_1/betaStar - sigma_om1*(0.41**2.0)/(betaStar**0.5)
      alfa_2    = beta_2/betaStar - sigma_om2*(0.41**2.0)/(betaStar**0.5)

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

            ! turbulent time scale
            Tt(i,k)   = 1.0/omNew(i,k)   ! 0.31 cmu/omega
               
            ! Bouyancy prodution with a different time scale
            Gk(i,k)=-ctheta*beta(i,k)*Fr_1*Tt(i,k)
     &             *  (ekmt(i,k)*(((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4.-(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/dRu(i)
     &                          +((U(i,kp)+U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz) )*
     &                                                                             (T(ip,k)-T(im,k))/(dRp(i)+dRp(im))  )
     &               +(2.*ekmt(i,k)*((W(i,k)-W(i,km))/dz-2./3.*(rho(i,k)*putink(i,k)))*(T(i,kp)-T(i,km))/(2.*dz)
     &                )


            Gk(i,k) = Gk(i,k) + ctheta*beta(i,k)*Fr_1*Tt(i,k)*2./3.*ekmt(i,k)*div(i,k)*(T(i,kp)-T(i,km))/(2.*dz)

            if (scl.eq.10) then
               ! k- equation
               putout(i,k) = putout(i,k) + ( Pk(i,k) + Gk(i,k) )/rho(i,k)         
               dimpl(i,k)  = dimpl(i,k)  + 0.09*omNew(i,k)            ! note, betaStar*rho*k*omega/(rho*k), set implicit and divided by density

            elseif (scl.eq.11) then
               ! omega- equation
               alfaSST   = alfa_1*bF1(i,k) + alfa_2*(1-bF1(i,k))
               betaSST   = beta_1*bF1(i,k) + beta_2*(1.0 - bF1(i,k))

               StR = (2.*(((W(i,k)-W(i,km))/dz)**2. +
     &                    ((U(i,k)-U(im,k))/dRu(i))**2. +
     &                    ((U(i,k)+U(im,k))/(2.*Rp(i)))**2.) +
     &                   (((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4.
     &                    -(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/dRu(i)
     &                   +((U(i,kp) +U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz)  )**2.)

               ! Bouyancy prodution divided by mut 
               GtR=-ctheta*beta(i,k)*Fr_1*Tt(i,k)
     &         *  ((((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4.-(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/dRu(i)
     &                      +((U(i,kp)+U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz) )*
     &                                                                           (T(ip,k)-T(im,k))/(dRp(i)+dRp(im))  )
     &         +(2*((W(i,k)-W(i,km))/dz-2./3.*(rho(i,k)*putink(i,k)))*(T(i,kp)-T(i,km))/(2.*dz)
     &         )


               !!! RENE: change to turbulent time scale here!
               GtR = GtR + ctheta*beta(i,k)*Fr_1*Tt(i,k)*2./3.*div(i,k)*(T(i,kp)-T(i,km))/(2.*dz)

               putout(i,k) = putout(i,k) + (alfaSST*StR*rho(i,k) + alfaSST*GtR*rho(i,k) + (1.0-bF1(i,k))*cdKOM(i,k) ) /rho(i,k)
               dimpl(i,k)  = dimpl(i,k)  + betaSST*omNew(i,k) ! note, beta*rho*omega^2/(rho*omega), set implicit and divided by density

            endif

            
         enddo
      enddo

      end

!>******************************************************************************************
!!      SST advancing the turbulence scalars of this model: k and omega
!!******************************************************************************************
      subroutine advanceScalar_SST(resK,resOm,Utmp,Wtmp,Rtmp,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      real*8 dnew(0:i1,0:k1),tempArray(0:i1,0:k1),dimpl(0:i1,0:k1)
      real*8 Utmp(0:i1,0:k1),Wtmp(0:i1,0:k1),Rtmp(0:i1,0:k1),sigmakSST(0:i1,0:k1)
      real*8 rho3(0:i1,0:k1)

      real*8     a  (imax)
      real*8     b  (imax)
      real*8     c  (imax)
      real*8     rhs(imax)

      real*8 scl
      real*8 resK, resOm

      resK  = 0.0
      resOm = 0.0

      ! modified turb. model
      !    modifDiffTerm = 1, our modification
      !    modifDiffTerm = 2, Aupoix modification

      if ((modifDiffTerm == 1) .or. (modifDiffTerm == 2)) then
         rho3 = Rtmp
      else
         rho3 = 1.0
      endif

      ! --------------------------------------------------------------------
      ! --------------------------------------------------------------------
      ! advance K
      dnew=0.0; dimpl = 0.0;
      call advecc(dnew,dimpl,kNew,utmp,wtmp,Ru,Rp,dru,dz,i1,k1,rank,periodic,.true.)
      !if (modifDiffTerm == 1) call advecrho(dnew,kNew,utmp,wtmp,Ru,Rp,dru,dz,i1,k1,rank)
      scl = 10.0
      call prodis_SST(dnew,dimpl,kNew,Utmp,Wtmp,temp,Rtmp,scl)

      ! calculating constant with blending function factor
      sigmakSST = 0.85*bF1 + 1.0*(1.0 - bF1)
      sigmakSST = 1.0/sigmakSST

      call diffcSSTKine(dnew,kNew,ekm,ekmi,ekmk,ekmt,sigmakSST,
     &                 Rtmp,Ru,Rp,dru,dz,rank,modifDiffTerm)


      if ((modifDiffTerm == 0) .or. (modifDiffTerm == 1)) then
          do k=1,kmax
             do i=1,imax
                
                a(i) = (ekmi(i-1,k)+(ekmt(i,k)+ekmt(i-1,k))/(sigmakSST(i,k)+sigmakSST(i-1,k)))/(0.5*(rho3(i-1,k)+rho3(i,k)))**0.5
                a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/Rtmp(i,k)/rho3(i,k)**0.5
                c(i) = (ekmi(i  ,k)+(ekmt(i,k)+ekmt(i+1,k))/(sigmakSST(i,k)+sigmakSST(i+1,k)))/(0.5*(rho3(i+1,k)+rho3(i,k)))**0.5
                c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/Rtmp(i,k)/rho3(i,k)**0.5

                b(i) = (rho3(i,k)*(-a(i)-c(i)) + dimpl(i,k))/alphak
                a(i) = a(i)*rho3(i-1,k)
                c(i) = c(i)*rho3(i+1,k)

                rhs(i) = dnew(i,k)  + (1-alphak)*b(i)*kNew(i,k)
             enddo

             b(1) = b(1)+a(1)
             i=imax
             b(i) = b(i) - (c(i) /alphak)
             !b(i) = (rho3(i,k)*(-(a(i)/rho3(i-1,k))-(c(i)/rho3(i+1,k))) - c(i) + dimpl(i,k) )/alphak
             !b(i) = (rho3(i,k)*(-a(i)-c(i)) - rho3(i+1,k)*c(i) + dimpl(i,k) )/alphak
             rhs(i) = dnew(i,k)  + (1-alphak)*b(i)*kNew(i,k)

             call matrixIdir(imax,a,b,c,rhs)

             do i=1,imax
                resK = resK + ((kNew(i,k) - rhs(i))/(kNew(i,k)+1.0e-20))**2.0
                kNew(i,k) = max(rhs(i), 1.0e-8)
             enddo
          enddo
      else if (modifDiffTerm == 2) then
          do k=1,kmax
             do i=1,imax

                a(i) = (ekmi(i-1,k)+(ekmt(i,k)+ekmt(i-1,k))/(sigmakSST(i,k)+sigmakSST(i-1,k)))/(0.5*(rho3(i-1,k)+rho3(i,k)))
                a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/Rtmp(i,k)
                c(i) = (ekmi(i  ,k)+(ekmt(i,k)+ekmt(i+1,k))/(sigmakSST(i,k)+sigmakSST(i+1,k)))/(0.5*(rho3(i+1,k)+rho3(i,k)))
                c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/Rtmp(i,k)
                
                b(i) = (rho3(i,k)*(-a(i)-c(i)) + dimpl(i,k))/alphak
                a(i) = a(i)*rho3(i-1,k)
                c(i) = c(i)*rho3(i+1,k)

                rhs(i) = dnew(i,k)  + (1-alphak)*b(i)*kNew(i,k)
             enddo

             b(1) = b(1)+a(1)
             i=imax
             b(i) = b(i) - (c(i) /alphak)
             !b(i) = (rho3(i,k)*(-(a(i)/rho3(i-1,k))-(c(i)/rho3(i+1,k))) - c(i) + dimpl(i,k) )/alphak
             !b(i) = (rho3(i,k)*(-a(i)-c(i)) - rho3(i+1,k)*c(i) + dimpl(i,k) )/alphak
             rhs(i) = dnew(i,k)  + (1-alphak)*b(i)*kNew(i,k)

             call matrixIdir(imax,a,b,c,rhs)

             do i=1,imax
                resK = resK + ((kNew(i,k) - rhs(i))/(kNew(i,k)+1.0e-20))**2.0
                kNew(i,k) = max(rhs(i), 1.0e-8)
             enddo
          enddo
      endif

      ! --------------------------------------------------------------------
      ! --------------------------------------------------------------------
      ! advance omega
      dnew=0.0; dimpl = 0.0;
      call advecc(dnew,dimpl,omNew,utmp,wtmp,Ru,Rp,dru,dz,i1,k1,rank,periodic,.true.)
      scl = 11.0
      call prodis_SST(dnew,dimpl,kNew,Utmp,Wtmp,temp,Rtmp,scl)

      ! calculating constant with blending function factor
      sigmakSST = 0.5*bF1 + 0.856*(1.0 - bF1)
      sigmakSST = 1.0/sigmakSST
      call diffcSSTOmega(dnew,omNew,ekm,ekmi,ekmk,ekmt,sigmakSST,Rtmp,Ru,Rp,dru,dz,rank,modifDiffTerm)

      do k=1,kmax
         do i=1,imax
            a(i) = (ekmi(i-1,k)+(ekmt(i,k)+ekmt(i-1,k))/(sigmakSST(i,k)+sigmakSST(i-1,k)))/(0.5*(rho3(i-1,k)+rho3(i,k)))**0.5 
            a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/Rtmp(i,k)/rho3(i,k)**0.5
            c(i) = (ekmi(i  ,k)+(ekmt(i,k)+ekmt(i+1,k))/(sigmakSST(i,k)+sigmakSST(i+1,k)))/(0.5*(rho3(i+1,k)+rho3(i,k)))**0.5 
            c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/Rtmp(i,k)/rho3(i,k)**0.5

            b(i) = ((-a(i)-c(i))*rho3(i,k)**0.5 + dimpl(i,k))/alphae
            a(i) = a(i)*rho3(i-1,k)**0.5
            c(i) = c(i)*rho3(i+1,k)**0.5

            rhs(i) = dnew(i,k)  + (1-alphae)*b(i)*omNew(i,k)
         enddo

         b(1) = b(1) + a(1)
         i = imax
         rhs(i) = dnew(i,k) - c(i)*omNew(i1,k) + (1-alphae)*(b(i)*omNew(i,k))

         call matrixIdir(imax,a,b,c,rhs)

         do i=1,imax
            resOm = resOm + ((omNew(i,k) - rhs(i))/(omNew(i,k)+1.0e-20))**2.0
            omNew(i,k) = max(rhs(i), 1.0e-8)
         enddo
      enddo

      end

      !>********************************************************************
      !! diffusion term for kine of the SST model in the z-direction, set as a source term...
      !!********************************************************************
      subroutine diffcSSTKine(putout,putin,ek,eki,ekk,ekmt,sigma,rho,Ru,Rp,dru,dz,rank1,diffVersion)
      implicit none
      include 'param.txt'
      integer   km,kp,rank1,diffVersion
      real*8     putout(0:i1,0:k1),putin(0:i1,0:k1),
     &     rho(0:i1,0:k1),ek(0:i1,0:k1),eki(0:i1,0:k1),ekk(0:i1,0:k1),
     &     ekmt(0:i1,0:k1),dru(0:i1),dz,Ru(0:i1),Rp(0:i1),sigma(0:i1,0:k1)
c
c
         if (diffVersion == 1) then       ! Inverse SLS
            do k=1,kmax
               kp=k+1
               km=k-1
               do i=1,imax
                     putout(i,k) = putout(i,k) + 1.0/rho(i,k)/sqrt(rho(i,k))*(
     3         ( (ekk(i,k ) + (ekmt(i,k)+ekmt(i,kp))/(sigma(i,k)+sigma(i,kp)))/sqrt(0.5*(rho(i,k)+rho(i,kp)))*(rho(i,kp)*putin(i,kp)-rho(i,k )*putin(i,k ))
     3          -(ekk(i,km) + (ekmt(i,k)+ekmt(i,km))/(sigma(i,k)+sigma(i,km)))/sqrt(0.5*(rho(i,k)+rho(i,km)))*(rho(i,k )*putin(i,k )-rho(i,km)*putin(i,km))
     &                 )/(dz*dz)   )
               enddo
            enddo
         elseif (diffVersion == 2) then   ! Aupoix
            do k=1,kmax
               kp=k+1
               km=k-1
               do i=1,imax
                     putout(i,k) = putout(i,k) + 1.0/rho(i,k)*(
     3         ( (ekk(i,k ) + (ekmt(i,k)+ekmt(i,kp))/(sigma(i,k)+sigma(i,kp)))/(0.5*(rho(i,k)+rho(i,kp)))*(rho(i,kp)*putin(i,kp)-rho(i,k )*putin(i,k ))
     3          -(ekk(i,km) + (ekmt(i,k)+ekmt(i,km))/(sigma(i,k)+sigma(i,km)))/(0.5*(rho(i,k)+rho(i,km)))*(rho(i,k )*putin(i,k )-rho(i,km)*putin(i,km))
     &                 )/(dz*dz)   )
               enddo
            enddo
         else                               ! Standard
            do k=1,kmax
               kp=k+1
               km=k-1
               do i=1,imax
                  putout(i,k) = putout(i,k) + 1.0/rho(i,k)*(
     3                    ( (ekk(i,k ) + (ekmt(i,k)+ekmt(i,kp))/(sigma(i,k)+sigma(i,kp)))*(putin(i,kp)-putin(i,k ))
     3                     -(ekk(i,km) + (ekmt(i,k)+ekmt(i,km))/(sigma(i,k)+sigma(i,km)))*(putin(i,k )-putin(i,km))   )/(dz*dz)   )
               enddo
            enddo
         endif
      end

      !>********************************************************************
      !! diffusion term for omega of the SST model in the z-direction, set as a source term...
      !!********************************************************************
      subroutine diffcSSTOmega(putout,putin,ek,eki,ekk,ekmt,sigma,rho,Ru,Rp,dru,dz,rank1,diffVersion)
      implicit none
      include 'param.txt'
      integer   km,kp,rank1,diffVersion
      real*8     putout(0:i1,0:k1),putin(0:i1,0:k1),
     &     rho(0:i1,0:k1),ek(0:i1,0:k1),eki(0:i1,0:k1),ekk(0:i1,0:k1),
     &     ekmt(0:i1,0:k1),dru(0:i1),dz,Ru(0:i1),Rp(0:i1),sigma(0:i1,0:k1)
c
c
         if ((diffVersion == 1) .or. (diffVersion == 2)) then       ! Inverse SLS & Aupoix
            do k=1,kmax
               kp=k+1
               km=k-1
               do i=1,imax
                     putout(i,k) = putout(i,k) + 1.0/rho(i,k)*(
     3         ( (ekk(i,k ) + (ekmt(i,k)+ekmt(i,kp))/(sigma(i,k)+sigma(i,kp)))/sqrt(0.5*(rho(i,k)+rho(i,kp)))*
     3                                (putin(i,kp)*sqrt(rho(i,kp)) - putin(i,k )*sqrt(rho(i,k )))
     3          -(ekk(i,km) + (ekmt(i,k)+ekmt(i,km))/(sigma(i,k)+sigma(i,km)))/sqrt(0.5*(rho(i,k)+rho(i,km)))*
     3                                (putin(i,k )*sqrt(rho(i,k )) - putin(i,km)*sqrt(rho(i,km)))
     &                 )/(dz*dz)   )
               enddo
            enddo
         else                               ! Standard
            do k=1,kmax
               kp=k+1
               km=k-1
               do i=1,imax
                  putout(i,k) = putout(i,k) + 1.0/rho(i,k)*(
     3                    ( (ekk(i,k ) + (ekmt(i,k)+ekmt(i,kp))/(sigma(i,k)+sigma(i,kp)))*(putin(i,kp)-putin(i,k ))
     3                     -(ekk(i,km) + (ekmt(i,k)+ekmt(i,km))/(sigma(i,k)+sigma(i,km)))*(putin(i,k )-putin(i,km))   )/(dz*dz)   )
               enddo
            enddo
         endif
      end


