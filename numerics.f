!*************************************************************************************
!
!*************************************************************************************
      !> prodis subroutine which calculates the production term of the turbulent scalar equation
      !! 
      subroutine prodis(putout,dimpl,putink,putine,putinv2,nuSAtmp,putinf,U,W,T,rho,scl)
      implicit none
      include 'param.txt'
      include 'common.txt'

C      integer im,ip,jm,jp,km,kp,ib,ie,jb,je,kb,ke !< integers
      real*8, dimension(0:i1,0:k1) :: putout,U,W,T,rho,div,putink,putine,putinv2,dimpl,nuSAtmp !,Tt,Srsq
      real*8, dimension(imax,kmax) :: putinf,putinftmp
      real*8  scl !,mut,a11,a12,a21,a22,a33,A,A2,A3,A2t,epsihh,tscl
C      real*8  cv1_3,cb1,cb2,cb3,cw1,cw2,cw3_6,inv_cb3,kappa_2,chi,fv1SA,fv2SA,r_SA,g_SA,fw_SA,StR,shatSA
C      real*8  sigma_om1,sigma_om2,beta_1,beta_2,betaStar,alfa_1,alfa_2,alfaSST,betaSST, GtR

      if (turbmod.eq.1) then
         call prodis_MK(putout,dimpl,putink,putine,U,W,T,rho,scl)
      elseif (turbmod.eq.3) then
         call prodis_VF(putout,dimpl,putink,putine,putinv2,putinf,U,W,T,rho,scl)
c      elseif (turbmod.eq.4) then
c         call prodis_SA(putout,dimpl,nuSAtmp,U,W,T,rho)
c      elseif (turbmod.eq.5) then
c         call prodis_SST(putout,dimpl,putink,U,W,T,rho,scl)
      endif
c
c      ib = 1
c      ie = i1-1
c
c      kb = 1
c      ke = k1-1
c
c      do k=kb,ke
c         kp=k+1
c         km=k-1
c         do i=ib,ie
c            ip=i+1
c            im=i-1
c
c            ! Production of turbulent kinetic energy
c            Pk(i,k) = ekmt(i,k)*(
c     &         2.*(((W(i,k)-W(i,km))/dz)**2. +
c     &             ((U(i,k)-U(im,k))/(Ru(i)-Ru(im)))**2. +
c     &             ((U(i,k)+U(im,k))/(2.*Rp(i)))**2.) +
c     &            (((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4.
c     &             -(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/dRu(i)
c     &            +((U(i,kp) +U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz)
c     &              )**2.)
c
c            div(i,k)=(Ru(i)*U(i,k)-Ru(im)*U(im,k))/(Rp(i)*dru(i))
c     &              +(      W(i,k) -      W(i,km))/dz
c
c            Pk(i,k) = Pk(i,k) - 2./3.*(rho(i,k)*putink(i,k)+ekmt(i,k)*(div(i,k)))*(div(i,k))
c
c            ! Bouyancy prodution
c            Gk(i,k)=-ctheta*beta(i,k)*Fr_1*putink(i,k)/putine(i,k)
c     &           *  (ekmt(i,k)*(((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4.-(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/dRu(i)
c     &                         +((U(i,kp)+U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz) )*
c     &                                                                              (T(ip,k)-T(im,k))/(dRp(i)+dRp(im))  )
c     &           +(2.*ekmt(i,k)*((W(i,k)-W(i,km))/dz-2./3.*(rho(i,k)*putink(i,k)))*(T(i,kp)-T(i,km))/(2.*dz)
c     &           )
c
c
c!!! RENE: change to turbulent time scale here!
c            Gk(i,k) = Gk(i,k) + ctheta*beta(i,k)*Fr_1*putink(i,k)/putine(i,k)*2./3.*ekmt(i,k)*div(i,k)*(T(i,kp)-T(i,km))/(2.*dz)
c
c
c            Tt(i,k)=putink(i,k)/putine(i,k)
c
c                
c            if (turbmod.eq.3) then
c            ! time scale for v2f model
c                StR = (2.*(((W(i,k)-W(i,km))/dz)**2. +
c     &                ((U(i,k)-U(im,k))/(Ru(i)-Ru(im)))**2. +
c     &                ((U(i,k)+U(im,k))/(2.*Rp(i)))**2.) +
c     &                (((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4.
c     &                -(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/dRu(i)
c     &                +((U(i,kp) +U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz)  )**2.)      
c!               Srsq(i,k) = Pk(i,k)*rho(i,k)/(2.*ekmt(i,k))
c               Srsq(i,k) = Str*rho(i,k)*0.5
c               Tt(i,k)   = max(putink(i,k)/putine(i,k), 6.0*(ekm(i,k)/(rho(i,k)*putine(i,k)))**0.5)
c               Tt(i,k)   = max(Tt(i,k), 1.0e-8)
c               Tt(i,k)   = min(Tt(i,k),0.6*putink(i,k)/(3.**0.5*putinv2(i,k)*cmu*(2.*Srsq(i,k))**0.5))
c               
c               
c               ! Bouyancy prodution with a different time scale
c               Gk(i,k)=-ctheta*beta(i,k)*Fr_1*Tt(i,k)
c     &                *  (ekmt(i,k)*(((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4.-(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/dRu(i)
c     &                             +((U(i,kp)+U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz) )*
c     &                                                                                (T(ip,k)-T(im,k))/(Rp(ip)-Rp(im))  )
c     &                  +(2.*ekmt(i,k)*((W(i,k)-W(i,km))/dz-2./3.*(rho(i,k)*putink(i,k)))*(T(i,kp)-T(i,km))/(2.*dz)
c     &                  )
c
c
c               Gk(i,k) = Gk(i,k) + ctheta*beta(i,k)*Fr_1*Tt(i,k)*2./3.*ekmt(i,k)*div(i,k)*(T(i,kp)-T(i,km))/(2.*dz)
c
c            endif
c
c            if (scl.eq.0) then
c               !k equation
c               putout(i,k) = putout(i,k) + ( Pk(i,k) + Gk(i,k) )/rho(i,k)
c               dimpl(i,k)  = dimpl(i,k) + putine(i,k)/putink(i,k)       ! note, rho*epsilon/(rho*k), set implicit and divided by density
c
c            elseif (scl.eq.1) then
c               !epsilon equation
c               putout(i,k) = putout(i,k) +(ce1*f1(i,k)*Pk(i,k)/Tt(i,k) +  ce1*f1(i,k)*Gk(i,k)/Tt(i,k) )/rho(i,k)
c               dimpl(i,k)  = dimpl(i,k)  + ce2*f2(i,k)/Tt(i,k)              ! note, ce2*f2*rho*epsilon/T/(rho*epsilon), set implicit and divided by density
c
c            elseif (scl.eq.2) then
c               !v'2 equation
c               putout(i,k) = putout(i,k) + putink(i,k)*putinf(i,k)       ! note, source is rho*k*f/rho
c               dimpl(i,k)  = dimpl(i,k)  + 6.*putine(i,k)/putink(i,k)    ! note, 6*rho*v'2*epsilon/k/(rho*v'2), set implicit and divided by density
c
c            elseif (turbmod.eq.4) then
c               !Spalart Allmaras
c               cv1_3     = (7.1)**3.0
c               cb1       = 0.1355
c               cb2       = 0.622
c               cb3       = 2.0/3.0
c               inv_cb3   = 1.0/cb3
c               kappa_2   = (0.41)**2.0   ! von karman constant
c               cw1       = (cb1/kappa_2) + (1.0+cb2)/cb3
c               cw2       = 0.3
c               cw3_6     = (2.0)**6.0
c
c              ! magnitude of rate of rotation: omega=sqrt(2*Wij*Wij), Wrz = 0.5*(dU/dz-dW/dr);  note, utheta=0 d/dtheta=0
c               StR = ( ( -( (W(ip,km)+W(ip,k)+W(i,km)+W(i ,k))/4.-(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/dRu(i)
c     &                   +( (U(i,kp) +U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz)           )**2.)
c
c               StR = StR**0.5
c
c               ! calculating Shat from SA model
c               chi    = nuSAtmp(i,k)/(ekm(i,k)/rho(i,k))
c               fv1SA  = (chi**3.0)/((chi**3.0) + cv1_3);
c               fv2SA  = 1.0 - (chi/(1.0 + (chi*fv1SA)))
c               ShatSA = StR + fv2SA*nuSAtmp(i,k)/(kappa_2*((wallDist(i))**2.0))
c
c               ! production term in SA model
c               Pk(i,k) = cb1*nuSAtmp(i,k)*ShatSA
c
c               ! destruction term in SA model
c               r_SA         = min(nuSAtmp(i,k)/(kappa_2*((wallDist(i))**2.0)*ShatSA), 10.0)
c               g_SA         = r_SA + cw2*((r_SA**6.0) - r_SA)
c               fw_SA        = g_SA*(((1.0 + cw3_6)/(g_SA**6.0 + cw3_6))**(1.0/6.0))
c
c               ! gustavo: i think this is not correct
c               !destrSA(i,k) = cw1/rho(i,k)*fw_SA*nuSAtmp(i,k)/((0.5-rp(i))**2)
c               dimpl(i,k) = dimpl(i,k) + cw1*fw_SA*nuSAtmp(i,k)/((wallDist(i))**2.0)
c
c
c               ! source term
c
c               if ((modifDiffTerm == 1) .or. (modifDiffTerm == 2)) then
c               ! invSLS and Aupoix SA model=  advection + Pk + (1/rho)*cb2/cb3*(d(nuSA*sqrt(rho))/dr)^2 +(d(nuSA*sqrt(rho))/dz)^2
c                   putout(i,k) = putout(i,k) + Pk(i,k) + cb2*inv_cb3/rho(i,k) * (
c     &               (((nuSAtmp(ip,k)*(rho(ip,k)**0.5)) - (nuSAtmp(im,k)*(rho(im,k)**0.5)))/(dRp(i)+dRp(im)))**2.0
c     &             + (((nuSAtmp(i,kp)*(rho(i,kp)**0.5)) - (nuSAtmp(i,km)*(rho(i,km)**0.5)))/(2.0*dz))**2.0  )
c               else
c               ! Conventional SA model=  advection + Pk + cb2/cb3*(dnuSA/dr)^2 +(dnuSA/dz)^2
c                   putout(i,k) = putout(i,k) + Pk(i,k) + cb2*inv_cb3 * (
c     &               ((nuSAtmp(ip,k) - nuSAtmp(im,k))/(dRp(i)+dRp(im)))**2.0 + ((nuSAtmp(i,kp) - nuSAtmp(i,km))/(2.0*dz))**2.0  )
c               endif
c
c
c
c            elseif (turbmod .eq. 5) then
c
c               ! k-omega SST               
c               Tt(i,k)   = 1.0/omNew(i,k)   ! 0.31 cmu/omega
c               
c               ! Bouyancy prodution with a different time scale
c               Gk(i,k)=-ctheta*beta(i,k)*Fr_1*Tt(i,k)
c     &                *  (ekmt(i,k)*(((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4.-(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/dRu(i)
c     &                             +((U(i,kp)+U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz) )*
c     &                                                                                (T(ip,k)-T(im,k))/(dRp(i)+dRp(im))  )
c     &                  +(2.*ekmt(i,k)*((W(i,k)-W(i,km))/dz-2./3.*(rho(i,k)*putink(i,k)))*(T(i,kp)-T(i,km))/(2.*dz)
c     &                  )
c
c
c               Gk(i,k) = Gk(i,k) + ctheta*beta(i,k)*Fr_1*Tt(i,k)*2./3.*ekmt(i,k)*div(i,k)*(T(i,kp)-T(i,km))/(2.*dz)
c
c               if (scl.eq.10) then
c                  ! k- equation of SST model
c                  putout(i,k) = putout(i,k) + ( Pk(i,k) + Gk(i,k) )/rho(i,k)          ! Gk(i,k)   ! Does not take into account the bouyancy term...
c                  dimpl(i,k)  = dimpl(i,k)  + 0.09*omNew(i,k)            ! note, betaStar*rho*k*omega/(rho*k), set implicit and divided by density
c
c               elseif (scl.eq.11) then
c                  ! omega- equation of SST model
c                  sigma_om1 = 0.5
c                  sigma_om2 = 0.856
c                  beta_1    = 0.075
c                  beta_2    = 0.0828
c                  betaStar  = 0.09
c                  alfa_1    = beta_1/betaStar - sigma_om1*(0.41**2.0)/(betaStar**0.5)
c                  alfa_2    = beta_2/betaStar - sigma_om2*(0.41**2.0)/(betaStar**0.5)
c                  alfaSST   = alfa_1*bF1(i,k) + alfa_2*(1-bF1(i,k))
c                  betaSST   = beta_1*bF1(i,k) + beta_2*(1.0 - bF1(i,k))
c
c                  StR = (2.*(((W(i,k)-W(i,km))/dz)**2. +
c     &                       ((U(i,k)-U(im,k))/dRu(i)**2. +
c     &                       ((U(i,k)+U(im,k))/(2.*Rp(i)))**2.) +
c     &                      (((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4.
c     &                       -(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/dRu(i)
c     &                      +((U(i,kp) +U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz)  )**2.)
c
c                  ! Bouyancy prodution divided by mut 
c                  GtR=-ctheta*beta(i,k)*Fr_1*Tt(i,k)
c     &            *  ((((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4.-(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/dRu(i)
c     &                         +((U(i,kp)+U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz) )*
c     &                                                                              (T(ip,k)-T(im,k))/(dRp(i)+dRp(im))  )
c     &            +(2*((W(i,k)-W(i,km))/dz-2./3.*(rho(i,k)*putink(i,k)))*(T(i,kp)-T(i,km))/(2.*dz)
c     &            )
c
c
c                  !!! RENE: change to turbulent time scale here!
c                  GtR = GtR + ctheta*beta(i,k)*Fr_1*Tt(i,k)*2./3.*div(i,k)*(T(i,kp)-T(i,km))/(2.*dz)
c
c
c
c                  putout(i,k) = putout(i,k) + (alfaSST*StR*rho(i,k) + alfaSST*GtR*rho(i,k) + (1.0-bF1(i,k))*cdKOM(i,k) ) /rho(i,k)
c                  dimpl(i,k)  = dimpl(i,k)  + betaSST*omNew(i,k) ! note, beta*rho*omega^2/(rho*omega), set implicit and divided by density
c
c               endif
c
c            endif
c         enddo
c      enddo

      end



!>********************************************************************
!!     Newton solver for wall boundary condition
!!********************************************************************
      subroutine funcNewtonSolve(enth_i1, enth_imax)
      implicit none
      include 'param.txt'
      include 'common.txt'
      real*8 enth_i1, enth_imax

         if (EOSmode.eq.0) call funcNewtonSolveIG(enth_i1, enth_imax)
         if (EOSmode.eq.1) call funcNewtonSolveRG(enth_i1, enth_imax)

      end



!>********************************************************************
!!     Newton solver for wall boundary condition with PIG
!!********************************************************************
      subroutine funcNewtonSolveIG(enth_i1, enth_imax)
      implicit none
      include 'param.txt'
      include 'common.txt'
      real*8 enth_i1, enth_imax, ekh_imax

        ekh_imax = 1./(Re*Pr)
        enth_i1 = enth_imax + (Rp(i1)-Rp(imax))*Qwall/(ekh_imax*Re*Pr) ! RENE not sure

      end



!>********************************************************************
!!     Newton solver for wall boundary condition RG
!!********************************************************************
      subroutine funcNewtonSolveRG(enth_i1, enth_imax)
      implicit none
      integer nIterNewton,success
      real*8 enth_i1, enth_imax, fxValue, fxValue1, ekh_imax
      success = 1
      fxValue = 1000.0
      nIterNewton = 0

      if (enth_imax.gt.2)    enth_imax =  2.0
      if (enth_imax.lt.-0.1) enth_imax = -0.1

      enth_i1 = enth_imax

      do while (abs(fxValue).gt.1.0e-10)
         call funcNewtonBC(enth_i1,        enth_imax, fxValue)
         call funcNewtonBC(enth_i1+1.0e-8, enth_imax, fxValue1)

         enth_i1 = enth_i1 - fxValue/((fxValue1-fxValue)/1.0e-8)

         if (nIterNewton.gt.200) then
            fxValue = 0.0
            success = 0
         endif

         nIterNewton = nIterNewton + 1
!     write (*,*) 'newton iter: ', nIterNewton,enth_i1,enth_imax,fxValue
      enddo

      if (success.eq.0) then
         write (*,*) 'newton didnt converge, enthimax= ',enth_imax,', ', nIterNewton, ', ', enth_i1
         stop
      endif
      end

!>********************************************************************
!!     function for wall boundary condition used by Newton solver
!!********************************************************************
      subroutine funcNewtonBC(enth, enthIMAX, fxValue)
      implicit none
      include 'param.txt'
      include 'common.txt'
      integer tabkhi,tabklo
      real*8 enth,lamOcpinter,enthIMAX,fxValue
      tabkhi = 0
      tabklo = 0
      call splint(enthTab,lamocpTab,lamocp2Tab,nTab,0.5*(enth+enthIMAX),lamOcpinter,tabkhi,tabklo)
      fxValue = enth - enthIMAX - (Rp(i1)-Rp(imax))*Qwall/lamocpinter
      end


!>********************************************************************
!!     PIG equation of state
!!********************************************************************
      subroutine state(enth,rho,mu,lam,tp,be,istap,rank)
      implicit none
      include 'param.txt'
      integer istap
      real*8 enth(0:i1,0:k1)
      real*8 rho(0:i1,0:k1)
      real*8 mu (0:i1,0:k1)
      real*8 lam(0:i1,0:k1),tp(0:i1,0:k1),be(0:i1,0:k1)
      if (EOSmode.eq.0) call stateIG(enth,rho,mu,lam,tp,be,istap,rank)
      if (EOSmode.eq.1) call stateRG(enth,rho,mu,lam,tp,be,istap,rank)
      end

!>********************************************************************
!!     PIG equation of state IG
!!********************************************************************
      subroutine stateIG(enth,rho,mu,lam,tp,be,istap,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      integer istap
      real*8 enth(0:i1,0:k1)
      real*8  rho(0:i1,0:k1)
      real*8  mu (0:i1,0:k1)
      real*8  lam(0:i1,0:k1),tp(0:i1,0:k1),be(0:i1,0:k1)

         mu   = 1.0/Re
         ekmi = 1.0/Re
         ekmk = 1.0/Re

         lam  = 1.0/(Re*Pr)
         ekhi = 1.0/(Re*Pr)
         ekhk = 1.0/(Re*Pr)

         rho = 1.0/(enth+1.0)
         tp  = (enth+1.0)
         be  = 1.0/tp
         Cp  = 1.0
         Cpi = 1.0
         Cpk = 1.0
      end


!>********************************************************************
!!     real gas equation of state RG
!!********************************************************************
      subroutine stateRG(enth,rho,mu,lam,tp,be,istap,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'

      integer tabkhi,tabklo,istap

      real*8 enth(0:i1,0:k1)
      real*8 rho(0:i1,0:k1)
      real*8 mu (0:i1,0:k1)
      real*8 lam(0:i1,0:k1),tp(0:i1,0:k1),be(0:i1,0:k1)
      real*8 con (0:i1,0:k1)
      real*8 cpp (0:i1,0:k1)
      real enthface,cpface,conface,muface,beface

      do k=0,k1
         call funcNewtonSolveRG(enth(i1,k), enth(imax,k))
         if (rank.eq.0.and.k.lt.K_start_heat) enth(i1,k)=enth(imax,k)
         do i=0,i1
            tabkhi = 0
            tabklo = 0
            call splint(enthTab,rhoTab,   rho2Tab,   nTab,enth(i,k),rho(i,k),tabkhi,tabklo)
            call splint(enthTab,muTab,    mu2Tab,    nTab,enth(i,k),mu (i,k), tabkhi,tabklo)
            call splint(enthTab,cpTab,    cp2Tab,    nTab,enth(i,k),Cp(i,k),tabkhi,tabklo)
            call splint(enthTab,lamocpTab,lamocp2Tab,nTab,enth(i,k),lam(i,k),tabkhi,tabklo)
            call splint(enthTab,tempTab,  temp2Tab,  nTab,enth(i,k),tp(i,k),tabkhi,tabklo)
            call splint(enthTab,betaTab,  beta2Tab,  nTab,enth(i,k),be(i,k),tabkhi,tabklo)
            mu(i,k)  = mu(i,k)/Re
            lam(i,k) = lam(i,k)/(Re*Pr)
         enddo
      enddo
      
      do k=0,kmax
         do i=0,imax
            tabkhi = 0
            tabklo = 0
            enthface = 0.5*(enth(i,k)+enth(i+1,k))
            call splint(enthTab,cpTab, cp2Tab, nTab,enthface,cpface,tabkhi,tabklo)
            call splint(enthTab,lamTab,lam2Tab,nTab,enthface,conface,tabkhi,tabklo)
            call splint(enthTab,betaTab,beta2Tab,nTab,enthface,beface,tabkhi,tabklo)
            ekhi(i,k) = conface/cpface/(Re*Pr)
            Cpi(i,k)=cpface
!            betai(i,k)=beface
            call splint(enthTab,muTab, mu2Tab, nTab,enthface,muface, tabkhi,tabklo)
            ekmi(i,k)  = muface/Re

            enthface = 0.5*(enth(i,k)+enth(i,k+1))
            call splint(enthTab,cpTab, cp2Tab, nTab,enthface,cpface,tabkhi,tabklo)
            call splint(enthTab,lamTab,lam2Tab,nTab,enthface,conface,tabkhi,tabklo)
            call splint(enthTab,betaTab,  beta2Tab,nTab,enthface,beface,tabkhi,tabklo)
            ekhk(i,k) = conface/cpface/(Re*Pr)
            Cpk(i,k)=cpface
!            betak(i,k)=beface
            call splint(enthTab,muTab, mu2Tab, nTab,enthface,muface, tabkhi,tabklo)
            ekmk(i,k)  = muface/Re
         enddo
      enddo
      return
      end






!>********************************************************************
!!     poisson solver
!!********************************************************************
      subroutine fillps(rank)
      implicit none
c     
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      integer ierr

      real*8 sumps,sumps_tot
c     
c     *** Fill the right hand for the poisson solver. ***
c     
      sumps = 0.0

      do  k=1,kmax
         do i=1,imax
            p(i,k)  = (
     1     ( Ru(i)*dUdt(i,k) - Ru(i-1)*dUdt(i-1,k) )/( Rp(i)*dru(i))
     3     + (     dWdt(i,k) -         dWdt(i,k-1) )/( dz         ) )/dt
     &           +      (rnew(i,k)-rold(i,k))/(dt*dt)

            qcrit(i,k) = p(i,k)*dt

            sumps = sumps + p(i,k)*dru(i)*dz
         enddo
      enddo


      call mpi_allreduce(sumps,sumps_tot,1,mpi_real8,
     &     mpi_sum,mpi_comm_world,ierr)

      end




!>********************************************************************
!!     helmotz solver
!!********************************************************************
      subroutine fillhm(rank)
      implicit none
c     
      include 'param.txt'
      include 'common.txt'
      real*8   Tt(0:i1,0:k1) ,Srsq(0:i1,0:k1)
      !real*8   Str

c     
c     *** Fill the right hand for the poisson solver. ***
c     


      do  k=1,kmax
         do i=1,imax
         
            ! time scale for v2f model
            Tt(i,k)   = max(knew(i,k)/enew(i,k), 6.0*(ekm(i,k)/(rnew(i,k)*enew(i,k)))**0.5)
             
            !extras
            Tt(i,k)   = max(Tt(i,k), 1.0e-8)
            Tt(i,k)   = min(Tt(i,k),0.6*knew(i,k)/(3.**0.5*v2new(i,k)*cmu*(2.*Srsq(i,k))**0.5))

            ! lenght scale for v2f model
            Lh(i,k)=0.23*max(knew(i,k)**1.5/enew(i,k),70.*((ekm(i,k)/rnew(i,k))**3./enew(i,k))**0.25)
            
            !extras
            Lh(i,k)=min(knew(i,k)**1.5/enew(i,k),knew(i,k)**1.5/(3.**0.5*v2new(i,k)*cmu*(2.*Srsq(i,k))**0.5))
            Lh(i,k)=0.23*max(Lh(i,k),70.*((ekm(i,k)/rnew(i,k))**3./enew(i,k))**0.25)
            
            fv2(i,k)= - (1.4-1.)*(2./3.-v2new(i,k)/knew(i,k))/Tt(i,k)
     &                - 0.3*(Pk(i,k))/(rnew(i,k)*knew(i,k))-5.*v2new(i,k)/(knew(i,k)*Tt(i,k))
            fv2(i,k) = fv2(i,k)/Lh(i,k)**2.0
         enddo
      enddo

      return
      end




!>********************************************************************
!!     correc
!!********************************************************************
      subroutine correc(rank,setold)
      implicit none
      include 'param.txt'
      include 'common.txt'
      integer setold
      real*8 pplus_w(imax)

      do k=1,kmax
         do i=1,imax-1
            dUdt(i,k)=dUdt(i,k)-dt*(p(i+1,k)-p(i,k))/(Rp(i+1)-Rp(i))
         enddo
      enddo

      do k=1,kmax-1
         do i=1,imax
            dWdt(i,k)=dWdt(i,k)-dt*(p(i,k+1)-p(i,k))/dz
         enddo
      enddo

!     if periodic is true, don't overwrite pplus on px-1
      call pshiftb_w(p,pplus_w,rank)
      if ((rank.eq. px-1).and.(periodic.ne.1)) then
         do i=1,imax
            pplus_w(i)=p(i,kmax)
         enddo
      endif
      
      do i=1,imax
         dWdt(i,kmax)=dWdt(i,kmax)-dt*(pplus_w(i)-p(i,kmax))/dz
      enddo
      if (setold.eq.0)then
         do k=0,kmax
            do i=0,imax
               dudt(i,k)=dUdt(i,k)/(0.5*(rnew(i,k)+rnew(i+1,k))) !dudt(i,k)=dUdt(i,k)/(0.5*(drdt(i,k)+drdt(i+1,k)))
               dwdt(i,k)=dWdt(i,k)/(0.5*(rnew(i,k)+rnew(i,k+1))) !dwdt(i,k)=dWdt(i,k)/(0.5*(drdt(i,k)+drdt(i,k+1)))
            enddo
         enddo
      endif
      
      if (setold.eq.1)then
         Uold = Unew
         Wold = Wnew

         do k=0,kmax
            do i=0,imax
               Unew(i,k)=dUdt(i,k)/(0.5*(rnew(i,k)+rnew(i+1,k))) !Unew(i,k)=dUdt(i,k)/(0.5*(drdt(i,k)+drdt(i+1,k)))
               Wnew(i,k)=dWdt(i,k)/(0.5*(rnew(i,k)+rnew(i,k+1))) !Wnew(i,k)=dWdt(i,k)/(0.5*(drdt(i,k)+drdt(i,k+1)))
            enddo
         enddo

      endif


      rold = rnew


      end




!>********************************************************************
!!     chkdt
!!********************************************************************
      subroutine chkdt(rank,istap)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      integer ierr,istap
      real*8  tmp,Courant,dtmp

      dt = 100.0

      Courant = 1

      do k=1,kmax
         do i=1,imax
            tmp = ( abs(Unew(i,k)) / ( Rp(i+1)-Rp(i) ) ) +
     &            ( abs(Wnew(i,k)) /         dz        )
            tmp = Courant/tmp
            dt  = min(dt, tmp)
         enddo
      enddo

      dtmp = dt

      call mpi_allreduce(dtmp,dt,1,mpi_real8,mpi_min,mpi_comm_world,ierr)

      end




!>********************************************************************
!!     chkdiv
!!********************************************************************
      subroutine chkdiv(rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      integer ierr,ll
      real*8   div,divmax,divbar,divmax_tot,divbar_tot,rhoip,rhoim,
     &     rhokp,rhokm
      divbar = 0.0
      divmax = 0.0

      do k=1,kmax
         do i=1,imax
            rhoip = 0.5*(rNew(i,k)+rNew(i+1,k))
            rhoim = 0.5*(rNew(i,k)+rNew(i-1,k))

            rhokp = 0.5*(rNew(i,k)+rNew(i,k+1))
            rhokm = 0.5*(rNew(i,k)+rNew(i,k-1))

            div =
     &           (Ru(i)*Unew(i,k)*rhoip-Ru(i-1)*Unew(i-1,k)*rhoim)*dz/Rp(i)+
     &           (Wnew(i,k)*rhokp-Wnew(i,k-1)*rhokm)*dru(i)+
     &           (rNew(i,k)-rold(i,k))/dt*dru(i)*dz

!     if (abs(div).gt.10e-6) write(6,*) i,k+kmax*rank,div

            divbar = divbar+div
            div    = abs(div)
            divmax = max(divmax,div)
         enddo
      enddo

      call mpi_allreduce(divbar,divbar_tot,1,mpi_real8,mpi_max,mpi_comm_world,ierr)
      call mpi_allreduce(divmax,divmax_tot,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)


!     if (rank.eq.0) write(6,100) divbar_tot,divmax_tot
! 100  format('Mass_loss/gain:Tot= ',e15.4,' Max= ',e15.4)

      end






!>********************************************************************
!!     mkgrid
!!********************************************************************
      subroutine mkgrid(rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      real*8  delta(imax),Yplus,X,rmax,drr,Rei
      real*8  pii,y,y1,y2,fA,fB,fC,fact, const

!     const = 0.25
!      const = 0.65
!     const = 0.50
c******************************************************************
      pi    = 4.0*atan(1.0)
      Rei   = 1.0/Re
      dz    = 1.0*LoD/(kmax*px)
      ru(0) = 0

      fA = 0.12
      fB = 2.4

!      fA = 0.0001
!      fB = 0.0001


      do i = 1,imax
         fact = (i-0.)/(imax-0.)
         ru(i) = (1.-tanh(fB*(fA-fact))/tanh(fA*fB))
         ru(i) = ru(i)/(1.-tanh(fB*(fA-1.))/tanh(fA*fB))
         delta(i)=0.5*(ru(i)-ru(i-1))
      enddo

      do i=0,imax
         ru(i)=ru(i)/(2.*ru(imax))
      enddo

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

      do i = 1,imax
         wallDist(i) = 0.5 - rp(i)
      enddo

      if (rank.eq.0) then
         open(11,file = 'grid.txt')
         write(11,*) Re, imax
         do i=0,imax
            Yplus = (0.5-Rp(i))*Re
            write(11,'(i5,4F12.6)') i,yplus,Ru(i),Rp(i),delta(max(1,i))
         enddo
         close(11)
      endif
      return
      end



!>********************************************************************
!!     readOldGrid
!!********************************************************************
      subroutine readOldGrid(rank)
      implicit none
      integer i,i1,imax,tmp,rank
      real*8 Re,Ru(0:200),Rp(0:200),yplus

      if (rank.eq.0) then
         open(11,file = 'grid.txt')
         read(11,*) Re, i1
         do i=0,i1
            Yplus = (0.5-Rp(i))*Re
            read(11,'(i5,4F12.6)') tmp,yplus,Ru(i),Rp(i)
         enddo
         close(11)
      endif
      stop
      end


!>********************************************************************
!!     diffusion term in the z-direction, set as a source term...
!!********************************************************************
      subroutine diffc(putout,putin,ek,eki,ekk,ekmt,sigma,rho,Ru,Rp,dru,dz,rank1,diffVersion)
      implicit none
      include 'param.txt'
      integer   km,kp,rank1,diffVersion
      real*8     putout(0:i1,0:k1),putin(0:i1,0:k1),
     &     rho(0:i1,0:k1),ek(0:i1,0:k1),eki(0:i1,0:k1),ekk(0:i1,0:k1),
     &     ekmt(0:i1,0:k1),dru(0:i1),dz,Ru(0:i1),Rp(0:i1),sigma
c
c
         if (diffVersion == 1) then       ! Inverse SLS
            do k=1,kmax
               kp=k+1
               km=k-1
               do i=1,imax
                     putout(i,k) = putout(i,k) + 1.0/rho(i,k)/sqrt(rho(i,k))*(
     3         ( (ekk(i,k ) + 0.5*(ekmt(i,k)+ekmt(i,kp))/sigma)/sqrt(0.5*(rho(i,k)+rho(i,kp)))*(rho(i,kp)*putin(i,kp)-rho(i,k )*putin(i,k ))
     3          -(ekk(i,km) + 0.5*(ekmt(i,k)+ekmt(i,km))/sigma)/sqrt(0.5*(rho(i,k)+rho(i,km)))*(rho(i,k )*putin(i,k )-rho(i,km)*putin(i,km))
     &                 )/(dz*dz)   )
               enddo
            enddo
         elseif (diffVersion == 2) then   ! Aupoix
            do k=1,k1-1
               kp=k+1
               km=k-1
               do i=1,i1-1
                     putout(i,k) = putout(i,k) + 1.0/rho(i,k)*(
     3         ( (ekk(i,k ) + 0.5*(ekmt(i,k)+ekmt(i,kp))/sigma)/(0.5*(rho(i,k)+rho(i,kp)))*(rho(i,kp)*putin(i,kp)-rho(i,k )*putin(i,k ))
     3          -(ekk(i,km) + 0.5*(ekmt(i,k)+ekmt(i,km))/sigma)/(0.5*(rho(i,k)+rho(i,km)))*(rho(i,k )*putin(i,k )-rho(i,km)*putin(i,km))
     &                 )/(dz*dz)   )
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




     
!>*****************************************************************
!!     
!!     diffu calculates the diffusion of u-velocity, which is
!!     the velocity in the radial direction.
!!     
!!
!!     In formula:  (4 terms)
!!
!!     1  d                  1  d                     d
!!     - -- (r Sigma(r r)) + - ---- (Sigma(phi r)) + -- (Sigma(z r)) -
!!     r dr                  r dphi                  dz
!!
!!
!!     1
!!     - Sigma(phi phi)
!!     r
!!
!!     r   : direction  ---> explicit (subroutine diffu)
!!     phi : direction  ---> implicit (subroutine predic)
!!     z   : direction  ---> explicit (subroutine diffu)
!!     
!!     on input :
!!     
!!     putout            : advection part
!!     Uvel,Vvel,Wvel    : contain velocities at n-1
!!     ekm               : diffusion coefficients (for velocity) in
!!     center points of the grid cells
!!     dr,dphi,dz        : grid spacing in r, phi and z-direction
!!     i1,j1,k1          : parameters for array-dimensions
!!     ib,ie,jb,je,kb,ke : range of gridpoints for which the
!!     diffusion part has to be calculated
!!     Ru,Rp             : radial positions of the U-velocity
!!     component and the pressure location
!!     respectively
!!     
!!     on output :
!!     
!!     putout            : advection and diffusion part
!!     other parameters  : all unchanged
!!     
!!*****************************************************************
      subroutine diffu (putout,Uvel,Wvel,ekme,Ru,Rp,dru,dz,i1,k1,dif)
      implicit none

      integer  i,k,im,ip,km,kp,i1,k1
      real*8     putout(0:i1,0:k1),Uvel(0:i1,0:k1),Wvel(0:i1,0:k1),
     &     ekme(0:i1,0:k1),dru(0:i1),dz,Ru(0:i1),Rp(0:i1),
     &     epop,epom,dzi,divUim,divUip,divUi,dif
c     
      dzi =1./dz
      do k=1,k1-1
         kp=k+1
         km=k-1
         do i=1,i1-1
            ip=i+1
            im=i-1

            epop = 0.25*(ekme(i,k)+ekme(ip,k) + ekme(ip,kp) + ekme(i,kp))
            epom = 0.25*(ekme(i,k)+ekme(ip,k) + ekme(ip,km) + ekme(i,km))

            divUim = (Ru(i)*Uvel(i,k) - Ru(im)*Uvel(im,k))/(Rp(i)*dru(i))
     &           + (        Wvel(i,k) -        Wvel(i,km))/dz

            divUip = (Ru(ip)*Uvel(ip,k)-Ru(i)*Uvel(i ,k ))/(Rp(ip)*dru(ip))
     &           + (         Wvel(ip,k)-      Wvel(ip,km))/dz

            divUi = ( Rp(ip)*(Uvel(ip,k)+Uvel(i,k)) - Rp(i)*(Uvel(i,k)+Uvel(im,k))
     &                )/(2.*Ru(i)*(Rp(ip)-Rp(i)))
     &           + ((Wvel(ip,k)+Wvel(i,k))-(Wvel(ip,km)+Wvel(i,km)))/(2.*dz)

            putout(i,k) = putout(i,k) +
     1   2.0*( Rp(ip)*ekme(ip,k)*(dif*(Uvel(ip,k)-Uvel(i ,k))/dru(ip) -1./3.*divUip)
     1        -Rp(i )*ekme(i ,k)*(dif*(Uvel(i ,k)-Uvel(im,k))/dru(i ) -1./3.*divUim)
     1            )/(Ru(i)*(Rp(ip)-Rp(i)))
     &           +
     3           ( epop * ( (Uvel(i,kp)-Uvel(i,k))*dzi
     3                    + (Wvel(ip,k)-Wvel(i,k))/(Rp(ip)-Rp(i)) )
     3             -
     3             epom * ( (Uvel(i,k)  -Uvel(i,km))*dzi
     3                    + (Wvel(ip,km)-Wvel(i,km))/(Rp(ip)-Rp(i)) )
     3            )*dzi
     &           -
     4           (ekme(i,k)+ekme(ip,k))/Ru(i)*(Uvel(i,k)/Ru(i)-1./3.*divUi)
         enddo
      enddo
      return
      end




!>*****************************************************************
!!     
!!     diffw calculates the diffusion of w-velocity, which is
!!     the velocity in the axial direction.
!!     
!!     
!!     In formula:  (3 terms)
!!     
!!     1  d                  1  d                     d
!!     - -- (r Sigma(r z)) + - ---- (Sigma(phi z)) + -- (Sigma(z z))
!!     r dr                  r dphi                  dz
!!     
!!     r   : direction  ---> explicit (subroutine diffw)
!!     phi : direction  ---> implicit (subroutine predic)
!!     z   : direction  ---> explicit (subroutine diffw)
!!     
!!     on input :
!!     
!!     putout            : advection part
!!     Uvel,Vvel,Wvel    : contain velocities at n-1
!!     ekm               : diffusion coefficients (for velocity) in
!!     center points of the grid cells
!!     dr,dphi,dz        : grid spacing in r, phi and z-direction
!!     i1,j1,k1          : parameters for array-dimensions
!!     ib,ie,jb,je,kb,ke : range of gridpoints for which the
!!     diffusion part has to be calculated
!!     Ru,Rp             : radial positions of the U-velocity
!!     component and the pressure location
!!     respectively
!!     
!!     on output :
!!     
!!     putout            : advection and diffusion part
!!     other parameters  : all unchanged
!!     
!!*****************************************************************
      subroutine diffw(putout,Uvel,Wvel,ekme,Ru,Rp,dru,dz,i1,k1,dif,rank)
      implicit none
     

      integer  i,k,im,ip,km,kp,i1,k1,rank
      real*8     putout(0:i1,0:k1),Uvel(0:i1,0:k1),Wvel(0:i1,0:k1),
     &     ekme(0:i1,0:k1),dru(0:i1),dz,Ru(0:i1),Rp(0:i1),
     &     epop,emop,divUkm,divUkp,dif
c
      do k=1,k1-1
         kp=k+1
         km=k-1
         do i=1,i1-1
            ip=i+1
            im=i-1

            epop = 0.25*(ekme(i,k)+ekme(i,kp) + ekme(ip,k) + ekme(ip,kp) )
            emop = 0.25*(ekme(i,k)+ekme(i,kp) + ekme(im,k) + ekme(im,kp) )

            divUkm = (Ru(i)*Uvel(i,k) - Ru(im)*Uvel(im,k))/(Rp(i)*dru(i))
     &             + (      Wvel(i,k) -        Wvel(i,km))/dz

            divUkp = (Ru(i)*Uvel(i,kp)- Ru(im)*Uvel(im, kp))/(Rp(i)*dru(i))
     &             + (      Wvel(i,kp)-        Wvel(i,  k ))/dz

            putout(i,k) = putout(i,k) +
     1         (Ru(i )*epop*(  (Uvel(i,kp) - Uvel(i,k))/dz
     1                    +dif*(Wvel(ip,k) - Wvel(i,k))/(Rp(ip)-Rp(i)))
     1         -
     1          Ru(im)*emop*(  (Uvel(im,kp) - Uvel(im,k))/dz
     1                    +dif*(Wvel(i,k)   - Wvel(im,k))/(Rp(i)-Rp(im)))
     1         )/(Rp(i)*dru(i))
     &         +
     3         (2.*ekme(i,kp)*((Wvel(i,kp)-Wvel(i,k ))/dz - 1./3.*divUkp)-
     3          2.*ekme(i,k )*((Wvel(i,k )-Wvel(i,km))/dz - 1./3.*divUkm))/dz
         enddo
      enddo
      return
      end




    
!>********************************************************************
!!     
!!     advecc calculates the advection for a scalar variable, which is
!!     situated in the center point of the grid cell.
!!     
!!     In formula:
!!     
!!     1 d(ruC)     1 d(vC)     d(wC)
!!     - (  - ------  +  - -----  +  -----  )
!!     r   dr       r  dphi      dz
!!     
!!     on input :
!!     
!!     putout            : "empty" (initialised to zero)
!!     putin             : variable for which the advection has
!!     to be calculated
!!     U,V,W             : contain velocities at former timestep
!!     putinn            : contains subgrid energy at oldest timestep
!!     dr,dphi,dz        : grid spacing in r, phi and z-direction
!!     i1,j1,k1          : parameters for array-dimensions
!!     ib,ie,jb,je,kb,ke : range of gridpoints for which the
!!     advection has to be calculated
!!     Ru,Rp             : radial positions of the U-velocity
!!     component and the pressure location
!!     respectively
!!     
!!     on output :
!!     
!!     putout            : advection part
!!     other parameters  : all unchanged
!!     
!!********************************************************************
      subroutine advecc(putout,dimpl,putin,U,W,Ru,Rp,dru,dz,i1,k1,rank,periodic,flagImpl)

      implicit none

      integer  i,k,im,ip,km,kp,i1,k1,ib,ie,kb,ke,rank,periodic
      real*8 putout(0:i1,0:k1),putin(0:i1,0:k1),dimpl(0:i1,0:k1)
      real*8 U(0:i1,0:k1),W(0:i1,0:k1),dru(0:i1),dz,Ru(0:i1),Rp(0:i1)
      real*8 eps,r1,rk1,re1,phi1,phik1,phie1,r2,rk2,re2,phi2,phik2,phie2,r3,rk3,re3,phi3,phik3,phie3,fak
      real*8 rhoip,rhoim,rhokp,rhokm
      real*8 dcubf(0:i1),dcwbf(0:i1)
      real*8 dcubb(0:i1),dcwbb(0:i1)
      real*8 dcu(0:i1,0:k1),dcw(0:i1,0:k1)
      real*8 cu(0:i1,0:k1),cw(0:i1,0:k1)
      logical flagImpl
!
      cu = 0.0
      cw = 0.0

      dcu = 0.0
      dcw = 0.0

      ib = 1
      ie = i1-1

      kb = 1
      ke = k1-1
!
!     compute delta C and distribute cpu boundaries
      do k=0,ke
         do i=0,ie
            dcu(i,k) = putin(i+1,k)-putin(i,k)
            dcw(i,k) = putin(i,k+1)-putin(i,k)
         enddo
      enddo

!     if (periodic.eq.1) then

      call shiftf(dcu,dcubf,rank)
      call shiftf(dcw,dcwbf,rank)

      do i=ib,ie
         dcu(i,0) = dcubf(i)
         dcw(i,0) = dcwbf(i)
      enddo

!     endif

      if ((periodic.ne.1).and.(rank.eq.0)) then
         dcu(:,0) = 0.0
         dcw(:,0) = 0.0
      endif

      eps=1.0e-16
      fak=1.0/3.0

!     calculate face value for C

      do k=kb,ke
         kp=k+1
         km=k-1
         do i=ib,ie
            ip=i+1
            im=i-1
            if (U(i,k).ge.(0.0)) then
               r1=(dcu(i,k)+eps)/(dcu(im,k)+eps)
               phi1=max(0.,min(2.*r1,min(fak*(1.+2.*r1),2.)))
               Cu(i,k)=putin(i,k)+0.5*phi1*(dcu(im,k))
            else
               r1=(dcu(i,k)+eps)/(dcu(ip,k)+eps)
               phi1=max(0.,min(2.*r1,min(fak*(1.+2.*r1),2.)))
               Cu(i,k)=putin(ip,k)+0.5*phi1*(-dcu(ip,k))
            endif
            if (W(i,k).ge.(0.0)) then
               r3=(dcw(i,k)+eps)/(dcw(i,km)+eps)
               phi3=max(0.,min(2.*r3,min(fak*(1.+2.*r3),2.)))
               Cw(i,k)=putin(i,k) + 0.5*phi3*(dcw(i,km))
            else
               r3=(dcw(i,k)+eps)/(dcw(i,kp)+eps)
               phi3=max(0.,min(2.*r3,min(fak*(1.+2.*r3),2.)))
               Cw(i,k)=putin(i,kp)+0.5*phi3*(-dcw(i,kp))
            endif
         enddo
      enddo

!     if (periodic.eq.1) then

      call shiftf(cu,dcubf,rank)
      call shiftf(cw,dcwbf,rank)


      do i=ib,ie
         cu(i,0)  = dcubf(i)
         cw(i,0)  = dcwbf(i)
      enddo

!     endif


      if ((periodic.ne.1).and.(rank.eq.0)) then
         do i=0,ie
            cu(i,0) = putin(i,0)
            cw(i,0) = putin(i,0)
         enddo
      endif

!     adv = u dc/dz = d(u c)/dz -c du/dz
      do k=kb,ke
         km=k-1
         do i=ib,ie
            im=i-1
            putout(i,k) =     - (Ru(i)*U(i,k)*cu(i,k) - Ru(im)*U(im,k)*cu(im,k))/(Rp(i)*dru(i))
     &                        - (      W(i,k)*cw(i,k) -        W(i,km)*cw(i,km))/(dz)
     4         + putin(i,k)*(   (Ru(i)*U(i,k) - Ru(i-1)*U(i-1,k ))/(Rp(i)*dru(i))
     &                        + (      W(i,k) -         W(i  ,km))/dz )
         enddo
      enddo

!     Implicit part
      if (flagImpl .eqv. .true.) then
         do k=kb,ke
            km=k-1
            do i=ib,ie
               im=i-1
                  putout(i,k) = putout(i,k) + W(i,k)*putin(i,k)/dz
                   dimpl(i,k) = dimpl(i,k)  + W(i,k)/dz
            enddo
         enddo
      endif

      end

!>********************************************************************
!!   
!!     advecrho calculate the extra term of the modified turbulence models (our modifications)
!!     coming from the adveccion
!!     
!!     In formula:
!!     
!!       U*C   d(rho)     W*C   d(rho)
!!    (------  ------  + -----  -----  )
!!      2*rho    dr      2*rho    dz
!!     
!!     on input :
!!     
!!     putout            : "empty" (initialised to zero)
!!     Uvel,Vvel,Wvel    : contain velocities at former timestep
!!     Wtmp              : contains velocity at oldest timestep
!!     dr,dphi,dz        : grid spacing in r, phi and z-direction
!!     i1,j1,k1          : parameters for array-dimensions
!!     ib,ie,jb,je,kb,ke : range of gridpoints for which the
!!     advection has to be calculated
!!     Ru,Rp             : radial positions of the U-velocity
!!     component and the pressure location
!!     respectively
!!     
!!     on output :
!!     
!!     putout            : advection part
!!     other parameters  : all unchanged
!!   
!!********************************************************************
      subroutine advecrho(putout,putin,U,W,Ru,Rp,dru,dz,i1,k1,rank)

  
      implicit none
      integer  i,k,im,ip,km,kp,i1,k1,ib,ie,kb,ke,rank
      real*8 putout(0:i1,0:k1),putin(0:i1,0:k1)
      real*8 U(0:i1,0:k1),W(0:i1,0:k1),dru(0:i1),dz,Ru(0:i1),Rp(0:i1), rnew(0:i1,0:k1)
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
            putout(i,k)  = putout(i,k) + 0.5*putin(i,k)*U(i,k)/rnew(i,k)*((rnew(ip,k ) - rnew(im,k ))/(Rp(ip)-Rp(im)))  
     &                                 + 0.5*putin(i,k)*W(i,k)/rnew(i,k)*((rnew(i ,kp) - rnew(i ,km))/(2.0*dz))       
      
         enddo
      enddo
      return
      end

!>********************************************************************
!!     
!!     advecu calculates the advection of the u-velocity, which is
!!     the velocity in the radial direction.
!!     
!!     In formula:
!!     
!!          1 d(ruu)     1 d(uv)     d(uw)     vv
!!     - (  - ------  +  - -----  +  -----  -  --  )
!!          r   dr       r  dphi      dz        r
!!     
!!     on input :
!!     
!!     putout            : "empty" (initialised to zero)
!!     Uvel,Vvel,Wvel    : contain velocities at former timestep
!!     Utmp              : contains velocity at oldest timestep
!!     dr,dphi,dz        : grid spacing in r, phi and z-direction
!!     i1,j1,k1          : parameters for array-dimensions
!!     ib,ie,jb,je,kb,ke : range of gridpoints for which the
!!     advection has to be calculated
!!     Ru,Rp             : radial positions of the U-velocity
!!     component and the pressure location
!!     respectively
!!     
!!     on output :
!!     
!!     putout            : advection part
!!     other parameters  : all unchanged
!!     
!!********************************************************************
      subroutine advecu(putout,Uvel,Wvel,RHO,Ru,Rp,dru,dz,i1,k1)
      implicit none
c     

      integer  i,k,im,ip,km,kp,i1,k1,ib,ie,kb,ke
      real*8     putout(0:i1,0:k1),Uvel(0:i1,0:k1),Wvel(0:i1,0:k1),
     &     dru(0:i1),dz,Ru(0:i1),Rp(0:i1)
      real*8 rho(0:i1,0:k1)
      real*8 rhoip,rhoim,rhokp,rhokm
c     
!     if (adv.eq.1) Uin=Uvel
      ib = 1
      ie = i1-1
      kb = 1
      ke = k1-1
c     
      do k=kb,ke
         kp=k+1
         km=k-1
         do  i=ib,ie
            ip=i+1
            im=i-1

            rhokp=0.25*(rho(i,k)+rho(i,kp)+rho(ip,k)+rho(ip,kp))
            rhokm=0.25*(rho(i,k)+rho(i,km)+rho(ip,k)+rho(ip,km))

            putout(i,k) = 0.0

            putout(i,k) = - 0.25 * (
     1           (Rp(ip)*(Uvel(i,k)+Uvel(ip,k))*(Uvel(i,k)+Uvel(ip,k))
     1           *rho(ip,k) -
     1            Rp(i )*(Uvel(im,k)+Uvel(i,k))*(Uvel(i,k)+Uvel(im,k))
     1           *rho(i ,k)  )
     1           / ( Ru(i) * ( Rp(ip)-Rp(i) ) )
     &           +
     3           ( (Wvel(i,k) +Wvel(ip,k) )*(Uvel(i,k)+Uvel(i,kp))*rhokp -
     3             (Wvel(i,km)+Wvel(ip,km))*(Uvel(i,k)+Uvel(i,km))*rhokm  )
     3           / ( dz )
     &           )
         enddo
      enddo
      return
      end




!>********************************************************************
!!     
!!     advecw calculates the advection of the w-velocity, which is
!!     the velocity in the axial direction.
!!     
!!     In formula:
!!     
!!     1 d(ruw)     1 d(wv)     d(ww)
!!     - (  - ------  +  - -----  +  -----  )
!!     r   dr       r  dphi      dz
!!     
!!     on input :
!!     
!!     putout            : "empty" (initialised to zero)
!!     Uvel,Vvel,Wvel    : contain velocities at former timestep
!!     Wtmp              : contains velocity at oldest timestep
!!     dr,dphi,dz        : grid spacing in r, phi and z-direction
!!     i1,j1,k1          : parameters for array-dimensions
!!     ib,ie,jb,je,kb,ke : range of gridpoints for which the
!!     advection has to be calculated
!!     Ru,Rp             : radial positions of the U-velocity
!!     component and the pressure location
!!     respectively
!!     
!!     on output :
!!     
!!     putout            : advection part
!!     other parameters  : all unchanged
!!     
!!********************************************************************
      subroutine advecw(putout,Uvel,Wvel,RHO,Ru,Rp,dru,dz,ekm,peclet_z)
      implicit none
      include 'param.txt'

      integer   im,ip,km,kp,ib,ie,kb,ke
      real*8     putout(0:i1,0:k1),Uvel(0:i1,0:k1),Wvel(0:i1,0:k1),
     &     dru(0:i1),dz,Ru(0:i1),Rp(0:i1)
      real*8 rho(0:i1,0:k1),ekm(0:i1,0:k1),peclet_z(0:i1,0:k1)
      real*8 rhoip,rhoim,advcecw_w

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

            rhoip=0.25*(rho(i,k)+rho(i,kp)+rho(ip,k)+rho(ip,kp))
            rhoim=0.25*(rho(i,k)+rho(i,kp)+rho(im,k)+rho(im,kp))

            peclet_z(i,k)=Wvel(i,k)*rho(i,k)*dz/(ekm(i,k))
            if (peclet_z(i,k).gt.2.)then
               advcecw_w= 2.0*(rho(i,k)+rho(i,k))*Wvel(i,k)*(Wvel(i,k)-Wvel(i,km))/ dz+
     &              2.0*Wvel(i,k)*((rho(i,k)+rho(i,k))*Wvel(i,k)-(rho(i,km)+rho(i,km))*Wvel(i,km))/ dz
            else
               advcecw_w= ((Wvel(i,k) +Wvel(i,kp) )*(Wvel(i,k)+Wvel(i,kp))*rho(i,kp) -
     &              (Wvel(i,km)+Wvel(i,k)  )*(Wvel(i,k)+Wvel(i,km))*rho(i,k )  ) / dz
            endif

            putout(i,k) = 0.0
            putout(i,k) = - 0.25 * (
     1           (Ru(i )*(Uvel(i,k) +Uvel(i,kp) )*(Wvel(i,k)+Wvel(ip,k))
     1           *rhoip -
     1           Ru(im)*(Uvel(im,k)+Uvel(im,kp))*(Wvel(i,k)+Wvel(im,k))
     1           *rhoim )
     1           / ( Rp(i) * dru(i) ) +advcecw_w)
         enddo
      enddo
      return
      end

