!*************************************************************************************
!
!*************************************************************************************
      !> prodis subroutine which calculates the production term of the turbulent scalar equation
      !! 
      subroutine prodis(putout,dimpl,putink,putine,putinv2,nuSAtmp,putinf,U,W,T,rho,C,scl)
      implicit none
      include 'param.txt'
      include 'common.txt'

      integer im,ip,jm,jp,km,kp,ib,ie,jb,je,kb,ke !< integers
      real*8, dimension(0:i1,0:k1) :: putout,U,W,T,rho,C,div,putink,putine,putinv2,Tt,Srsq,nuSAtmp,dimpl
      real*8, dimension(imax,kmax) :: putinf,putinftmp
      real*8  scl,mut,a11,a12,a21,a22,a33,A,A2,A3,A2t,epsihh,tscl
      real*8  cv1_3,cb1,cb2,cb3,cw1,cw2,cw3_6,inv_cb3,kappa_2,chi,fv1SA,fv2SA,r_SA,g_SA,fw_SA,StR,shatSA
      real*8  sigma_om1,sigma_om2,beta_1,beta_2,betaStar,alfa_1,alfa_2,alfaSST,betaSST, GtR
      real*8  uudtdx,uTdudx,betagT,auT,Ttemp,Tmix,fd1,fd2,feps,d2Tdxdr  !modTemp
      real*8  Q,TauRm,diverg,TauRp,TauZm,TauZp
      character*5 cha
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


            if (rank.eq.0.and.k.lt.K_start_heat) then
               Q=0.0
            else
               Q=Qwall
            endif
            
            !******************************************
            ! Production of turbulent kinetic energy
            Pk(i,k) = ekmt(i,k)*(
     &         2.*(((W(i,k)-W(i,km))/dz)**2. +
     &             ((U(i,k)-U(im,k))/(Ru(i)-Ru(im)))**2. +
     &             ((U(i,k)+U(im,k))/(2.*Rp(i)))**2.) +
     &            (((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4.
     &             -(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/(Ru(i)-Ru(im))
     &            +((U(i,kp) +U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz)
     &              )**2.)

            div(i,k)=(Ru(i)*U(i,k)-Ru(im)*U(im,k))/(Rp(i)*dr(i))
     &              +(      W(i,k) -      W(i,km))/dz

            Pk(i,k) = Pk(i,k) - 2./3.*(rho(i,k)*putink(i,k)+ekmt(i,k)*(div(i,k)))*(div(i,k))

            !******************************************
            !defining the time scale
            Tt(i,k)   = putink(i,k)/putine(i,k) 
            
            if (turbmod.eq.3) then
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
               
               if (modVF.eq.1) then
               ! Modifications: Lien&Kalitzin 2001 "Computations of transonic flow with the v2f turbulence model"
                  Tt(i,k)   = max(Tt(i,k), 1.0e-8)
                  Tt(i,k)   = min(Tt(i,k),0.6*putink(i,k)/(3.**0.5*putinv2(i,k)*cmu*(2.*Srsq(i,k))**0.5))
               endif

            elseif (turbmod .eq. 4) then
            ! time scale for SA  
                         
               Tt(i,k)   = 0.0
               
            elseif (turbmod .eq. 5) then
            ! time scale for k-omega SST               
               Tt(i,k)   = 1.0/omNew(i,k)   ! 0.31 cmu/omega   
                           
            endif

            !******************************************
            ! Bouyancy prodution
            Gk(i,k)=-ctheta*beta(i,k)*Fr_1*Tt(i,k)
     &           *  (ekmt(i,k)*(((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4.-(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/(Ru(i)-Ru(im))
     &                         +((U(i,kp)+U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz) )*
     &                                                                              (T(ip,k)-T(im,k))/(Rp(ip)-Rp(im))  )
     &           +(2.*ekmt(i,k)*((W(i,k)-W(i,km))/dz-2./3.*(rho(i,k)*putink(i,k)))*(T(i,kp)-T(i,km))/(2.*dz)
     &           )

            Gk(i,k) = Gk(i,k) + ctheta*beta(i,k)*Fr_1*Tt(i,k)*2./3.*ekmt(i,k)*div(i,k)*(T(i,kp)-T(i,km))/(2.*dz)
            

            if (tempturbmod.GT.1) then   !modTemp

                ! Bouyancy production with Algebraix Flux model (2012 Zhang et al)
                ! Four terms come into play:
                ! 1) uudtdx= <u'i*u'j>*dT/dxj, 
                ! 2) uTdudx= <u'j*T'>*dui/dxj, 
                ! 3) betagT=beta*g_i*<T'2>, 
                ! 4) auT =aij * <u'j*T'> , where aij=<u'i*u'j>/k - 2/3 rho*delta_ij
                ! where <  > is favre averaged
                if (modGk.GT.0) then
                !1)-----------------
                ! <u'i*u'j>*dT/dxj   OK!!!  
                uudtdx= (ekmt(i,k)*(((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4.-(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/(Ru(i)-Ru(im))
     &                         +((U(i,kp)+U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz) )*
     &                                                                              (T(ip,k)-T(im,k))/(Rp(ip)-Rp(im))  )
     &           +(2.*ekmt(i,k)*((W(i,k)-W(i,km))/dz-2./3.*(rho(i,k)*putink(i,k)))*(T(i,kp)-T(i,km))/(2.*dz)           )
                uudtdx= uudtdx + 2./3.*ekmt(i,k)*div(i,k)*(T(i,kp)-T(i,km))/(2.*dz)

                !2)-----------------
                ! <u'j*T'>*dui/dxj, OK!  rho <u'j T'> [kgK/(m2 s)]=(lambda/cp)/cp dCdx
                uTdudx = ekht(i,k)/cp(i,k)*((C(ip,k)-C(im,k))/(Rp(ip)-Rp(im))
     &          *((U(ip,k)-U(im,k))/(Rp(ip)-Rp(im)) + ((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4.-(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/(Ru(i)-Ru(im))   )
     &                  +(C(i,kp)-C(i,km))/(2.*dz)        
     &          *(((W(i,k)-W(i,km))/dz)             + ((U(i,kp)+U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz)             ))

                !3)-----------------
                ! beta*g_i*<T'2>   OK!!! is the density included in the Fr_1?
                betagT= rho(i,k)*beta(i,k)*Fr_1*ktnew(i,k) 

                !4)-----------------
                ! aij * <u'j*T'> , where aij=<u'i*u'j>/k - 2/3 rho delta_ij
                auT = ekmt(i,k)/putink(i,k)*(((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4.-(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/(Ru(i)-Ru(im))
     &                         +((U(i,kp)+U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz) )
     &                                                                             *ekht(i,k)/cp(i,k)*(C(ip,k)-C(im,k))/(Rp(ip)-Rp(im))  
     &                         +(2.*ekmt(i,k)/putink(i,k)*((W(i,k)-W(i,km)))/dz- 2./3.*ekmt(i,k)/putink(i,k)*div(i,k))
     &                                                                             *ekht(i,k)/cp(i,k)*(C(i,kp)-C(i,km))/(2.*dz) 


                Gk(i,k)=beta(i,k)*Fr_1*(-ctheta0*Tt(i,k)*(ctheta1*uudtdx+(1-ctheta2)*uTdudx+(1-ctheta3)*betagT)+ctheta4*auT)
                
                endif

                !******************************************
                ! Production of temperature fluctuations  Pkt= <u'j T'> dTdxj= (lambda_t/cp)/cp dCdxj dTdxj = = (lambda_t/cp)dTdxj^2
                Pkt(i,k) = ekht(i,k)/cp(i,k)*(((C(i,k)-C(i,km))/dz)            *((T(i,k)-T(i,km))/dz) 
     &                                      + ((C(i,k)-C(im,k))/(Ru(i)-Ru(im)))*((T(i,k)-T(im,k))/(Ru(i)-Ru(im))))
!                Pkt(i,k) = ekht(i,k)*(((T(i,k)-T(i,km))/dz)**2.0+ ((T(i,k)-T(im,k))/(Ru(i)-Ru(im)))**2.0)

!                Checking the value of each term 
!                if (mod(istep,1000 ).eq.0) then 
!                   if (rank.eq.0)  then
!                      if (kmax*0.75.eq.k)  write(*,*) 'Values of B_k', (k+rank*kmax)*dz,ctheta0*Tt(i,k)*(ctheta1*uudtdx),
!     &                     ctheta0*Tt(i,k)*((1-ctheta2)*uTdudx), ctheta0*Tt(i,k)*((1-ctheta3)*betagT),ctheta4*auT,Gk(i,k)/(beta(i,k)*Fr_1),Pkt(i,k)
!                   endif
!                endif


            else
                Pkt(i,k) = 0.0
            endif

            !******************************************
            ! Mass flux contribution term to the turbulent kinetic energy budget: Kawai & Oikawa 2018 ETMM
            if (modMktau.eq.1) then   

c                   kp=k+1,km=k-1, ip=i+1, im=i-1
c                    __________________________________
c                   |  i+1,k-1 |  i+1,k   | i+1,k+1  |
c                   |__________|__________|__________| 
c                   |          |          |          |
c                   |  i,k-1   |   i,k    |   i,k+1  |
c                   |__________|__________|__________| 
c                   |          |          |          |
c                   | i-1,k-1  |  i-1,k   |  i-1,k+1 |
c                   |          |          |          |
c                   |__________|__________|__________| 
c

                ! Mktau=mu_t/(0.15 rho^2) [drho/dr d(tau_rr+tau_rx)/dr +drho/dx d(tau_xx+tau_xr)/dx]

                ! Calculating dTau_ri/dr
                TauRm =  ekmi(im,k)*Ru(im)*(
     &                   (2*((U(i,k)-U(im,k))/(Rp(i)-Rp(im)))**2.)                                                !dudr
     &                  +((W(i,k)-W(im,k))/(Rp(i)-Rp(im))                                                         !dwdr
     &                  +((U(i,kp)+U(i,k)+U(im,kp)+U(im,k))/4.
     &                  - (U(i,k)+U(i,km)+U(im,k)+U(im,km))/4.)/(dz)  )**2.                                       !dudz
     &                                                       )**0.5                      
                diverg=  1/Ru(im)*((Rp(i)*U(i,k)-Rp(im)*U(im,k))/(Rp(i)-Rp(im)))                                  !1/r drudr
     &                   +((W(i,kp)+W(i,k)+W(im,kp)+W(im,k))/4.
     &                   - (W(i,k)+W(i,km)+W(im,k)+W(im,km))/4.)/(dz)                                             !dwdz
                TauRm = TauRm - 2/3 * ekmi(im,k)*Ru(im)*(diverg*diverg)**0.5

                TauRp =  ekmi(i,k)*Ru(i)*(
     &                   (2*((U(ip,k)-U(i,k))/(Rp(ip)-Rp(i)))**2.)                                                !dudr
     &                  +((W(ip,k)-W(i,k))/(Rp(ip)-Rp(i))                                                         !dwdr
     &                  +((U(ip,kp)+U(ip,k)+U(i,kp)+U(i,k))/4.
     &                  - (U(ip,k)+U(ip,km)+U(i,k)+U(i,km))/4.)/(dz)  )**2.                                       !dudz
     &                                                       )**0.5                      
                diverg=  1/Ru(i)*((Rp(ip)*U(ip,k)-Rp(i)*U(i,k))/(Rp(ip)-Rp(i)))                                   !1/r drudr
     &                   +((W(ip,kp)+W(ip,k)+W(i,kp)+W(i,k))/4.
     &                   - (W(ip,k)+W(ip,km)+W(i,k)+W(i,km))/4.)/(dz)                                             !dwdz
                TauRp = TauRp - 2/3 * ekmi(i,k)*Ru(i)*(diverg*diverg)**0.5

                ! Calculating dTau_zi/dz
                TauZm =  ekmk(i,km)*(
     &                   (2.*((W(i,k)-W(i,km))/dz)**2.)
     &                  +(((U(i,k)-U(i,km))/dz)
     &                  +((W(ip,k)+W(ip,km)+W(i,k)+W(i,km))/4.
     &                  -(W(i,k)+W(i,km)+W(im,k)+W(im,km))/4.)/(Ru(i)-Ru(im)))**2.
     &                                                       )**0.5

                diverg=  1/Rp(i)*(
     &                     (Rp(ip)*(U(ip,k)+U(ip,km))+Rp(i) *(U(i,k) +U(i, km)))/4
     &                    -(Rp(i) *(U(i,k) +U(i,km)) +Rp(im)*(U(im,k)+U(im,km)))/4
     &                     )/(Ru(i)-Ru(im))            
     &                   +(W(i,k)-W(i,km))/dz         
                TauZm = TauZm - 2/3 * ekmk(im,km)*(diverg*diverg)**0.5

                TauZp =  ekmk(i,k)*(
     &                   (2.*((W(i,kp)-W(i,k))/dz)**2.)                        
     &                  +(((U(i,kp)-U(i,k))/dz)
     &                  +((W(ip,kp)+W(ip,k)+W(i,kp)+W(i,k))/4.
     &                  -(W(i,kp)+W(i,k)+W(im,kp)+W(im,k))/4.)/(Ru(ip)-Ru(i)))**2.
     &                                                       )**0.5

                diverg=  1/Rp(i)*(
     &                     (Rp(ip)*(U(ip,kp)+U(ip,k))+Rp(i) *(U(i,kp) +U(i, k)))/4
     &                    -(Rp(i) *(U(i,kp) +U(i,k)) +Rp(im)*(U(im,kp)+U(im,k)))/4
     &                     )/(Ru(i)-Ru(im))                                  !1/r drudr
     &                   +(W(i,kp)-W(i,k))/dz       !dwdz  
                TauZm = TauZm - 2/3 * ekmk(im,k)*(diverg*diverg)**0.5
                
                Mktau(i,k)=  ((rho(i,k)-rho(im,k))/(Ru(i)-Ru(im)))*((TauRp-TauRm)/(Ru(i)-Rp(im)))
     &                       +((rho(i,k)-rho(i,km))/dz)*((TauZp-TauZm)/(dz))
                Mktau(i,k)= ekmt(i,k)/(0.15*rho(i,k)**2.0)*Mktau(i,k)
            else
                Mktau(i,k) = 0.0
            endif



            !******************************************
            ! Source term for turbulence models
            
            if (scl.eq.0) then
               !k-equation for MK and V2F
               putout(i,k) = putout(i,k) + ( Pk(i,k) + Gk(i,k) + Mktau(i,k))/rho(i,k)
               dimpl(i,k)  = dimpl(i,k) + putine(i,k)/putink(i,k)         ! note, rho*epsilon/(rho*k), set implicit and divided by density

            elseif (scl.eq.1) then
               !epsilon-equation for MK and V2F
               putout(i,k) = putout(i,k) +(ce1*f1(i,k)*( Pk(i,k) + Gk(i,k) )/Tt(i,k))/rho(i,k)
               dimpl(i,k)  = dimpl(i,k)  + ce2*f2(i,k)/Tt(i,k)            ! note, ce2*f2*rho*epsilon/T/(rho*epsilon), set implicit and divided by density

            elseif (scl.eq.2) then
               !v'2 equation for V2F
               putout(i,k) = putout(i,k) + putink(i,k)*putinf(i,k)        ! note, source is rho*k*f/rho
               dimpl(i,k)  = dimpl(i,k)  + 6.*putine(i,k)/putink(i,k)     ! note, 6*rho*v'2*epsilon/k/(rho*v'2), set implicit and divided by density

            elseif (turbmod.eq.4) then
               !Spalart Allmaras
               cv1_3     = (7.1)**3.0
               cb1       = 0.1355
               cb2       = 0.622
               cb3       = 2.0/3.0
               inv_cb3   = 1.0/cb3
               kappa_2   = (0.41)**2.0   ! von karman constant
               cw1       = (cb1/kappa_2) + (1.0+cb2)/cb3
               cw2       = 0.3
               cw3_6     = (2.0)**6.0

              ! magnitude of rate of rotation: omega=sqrt(2*Wij*Wij), Wrz = 0.5*(dU/dz-dW/dr);  note, utheta=0 d/dtheta=0
               StR = ( ( -( (W(ip,km)+W(ip,k)+W(i,km)+W(i ,k))/4.-(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/(Ru(i)-Ru(im))
     &                   +( (U(i,kp) +U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz)           )**2.)

               StR = StR**0.5

               ! calculating Shat from SA model
               chi    = nuSAtmp(i,k)/(ekm(i,k)/rho(i,k))
               fv1SA  = (chi**3.0)/((chi**3.0) + cv1_3);
               fv2SA  = 1.0 - (chi/(1.0 + (chi*fv1SA)))
               ShatSA = StR + fv2SA*nuSAtmp(i,k)/(kappa_2*((0.5-rp(i))**2.0))

               ! production term in SA model
               Pk(i,k) = cb1*nuSAtmp(i,k)*ShatSA

               ! destruction term in SA model
               r_SA         = min(nuSAtmp(i,k)/(kappa_2*((0.5-rp(i))**2.0)*ShatSA), 10.0)
               g_SA         = r_SA + cw2*((r_SA**6.0) - r_SA)
               fw_SA        = g_SA*(((1.0 + cw3_6)/(g_SA**6.0 + cw3_6))**(1.0/6.0))

               ! gustavo: i think this is not correct
               !destrSA(i,k) = cw1/rho(i,k)*fw_SA*nuSAtmp(i,k)/((0.5-rp(i))**2)
               dimpl(i,k) = dimpl(i,k) + cw1*fw_SA*nuSAtmp(i,k)/((0.5-rp(i))**2.0)


               ! source term

               if ((modifDiffTerm == 1) .or. (modifDiffTerm == 2)) then
               ! invSLS and Aupoix SA model=  advection + Pk + (1/rho)*cb2/cb3*(d(nuSA*sqrt(rho))/dr)^2 +(d(nuSA*sqrt(rho))/dz)^2
                   putout(i,k) = putout(i,k) + Pk(i,k) + cb2*inv_cb3/rho(i,k) * (
     &               (((nuSAtmp(ip,k)*(rho(ip,k)**0.5)) - (nuSAtmp(im,k)*(rho(im,k)**0.5)))/(Rp(ip)-Rp(im)))**2.0
     &             + (((nuSAtmp(i,kp)*(rho(i,kp)**0.5)) - (nuSAtmp(i,km)*(rho(i,km)**0.5)))/(2.0*dz))**2.0  )
               else
               ! Conventional SA model=  advection + Pk + cb2/cb3*(dnuSA/dr)^2 +(dnuSA/dz)^2
                   putout(i,k) = putout(i,k) + Pk(i,k) + cb2*inv_cb3 * (
     &               ((nuSAtmp(ip,k) - nuSAtmp(im,k))/(Rp(ip)-Rp(im)))**2.0 + ((nuSAtmp(i,kp) - nuSAtmp(i,km))/(2.0*dz))**2.0  )
               endif



            elseif (turbmod .eq. 5) then

               if (scl.eq.10) then
                  ! k- equation of SST model
                  putout(i,k) = putout(i,k) + ( Pk(i,k) + Gk(i,k) + Mktau(i,k))/rho(i,k)          ! Gk(i,k)   ! Does not take into account the bouyancy term...
                  dimpl(i,k)  = dimpl(i,k)  + 0.09*omNew(i,k)            ! note, betaStar*rho*k*omega/(rho*k), set implicit and divided by density

               elseif (scl.eq.11) then
                  ! omega- equation of SST model
                  sigma_om1 = 0.5
                  sigma_om2 = 0.856
                  beta_1    = 0.075
                  beta_2    = 0.0828
                  betaStar  = 0.09
                  alfa_1    = beta_1/betaStar - sigma_om1*(0.41**2.0)/(betaStar**0.5)
                  alfa_2    = beta_2/betaStar - sigma_om2*(0.41**2.0)/(betaStar**0.5)
                  alfaSST   = alfa_1*bF1(i,k) + alfa_2*(1-bF1(i,k))
                  betaSST   = beta_1*bF1(i,k) + beta_2*(1.0 - bF1(i,k))

                  StR = (2.*(((W(i,k)-W(i,km))/dz)**2. +
     &                       ((U(i,k)-U(im,k))/(Ru(i)-Ru(im)))**2. +
     &                       ((U(i,k)+U(im,k))/(2.*Rp(i)))**2.) +
     &                      (((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4.
     &                       -(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/(Ru(i)-Ru(im))
     &                      +((U(i,kp) +U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz)  )**2.)

                  ! Bouyancy prodution divided by mut 
                  GtR=-ctheta*beta(i,k)*Fr_1*Tt(i,k)
     &            *  ((((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4.-(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/(Ru(i)-Ru(im))
     &                         +((U(i,kp)+U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz) )*
     &                                                                              (T(ip,k)-T(im,k))/(Rp(ip)-Rp(im))  )
     &            +(2*((W(i,k)-W(i,km))/dz-2./3.*(rho(i,k)*putink(i,k)))*(T(i,kp)-T(i,km))/(2.*dz)
     &            )


                  !!! RENE: change to turbulent time scale here!
                  GtR = GtR + ctheta*beta(i,k)*Fr_1*Tt(i,k)*2./3.*div(i,k)*(T(i,kp)-T(i,km))/(2.*dz)



                  putout(i,k) = putout(i,k) + (alfaSST*StR*rho(i,k) + alfaSST*GtR*rho(i,k) + (1.0-bF1(i,k))*cdKOM(i,k) ) /rho(i,k)
                  dimpl(i,k)  = dimpl(i,k)  + betaSST*omNew(i,k) ! note, beta*rho*omega^2/(rho*omega), set implicit and divided by density

               endif

            endif
            
            !******************************************
            ! Source term for temperature turbulence models
            
            if (tempturbmod.GT.1) then   

               if (scl.eq.20) then
               ! source for the kt-equation (T'2)
                  putout(i,k) = putout(i,k) + 2.0*Pkt(i,k)/rho(i,k) !2.0*Pkt(i,k)/rho(i,k)
                  dimpl(i,k)  = dimpl(i,k)  + 2.0*etnew(i,k)/(ktnew(i,k)+1.0e-20) !2.0*etnew(i,k)/ktnew(i,k)            ! note, rho*epst/(rho*kt), set implicit and divided by density
               
               elseif (scl.eq.21) then
                ! source for the epst-equation (et)
                 
                  if (tempturbmod.eq.2) then
                  !-------------------------
                  ! Nagano and Kim model 1988
                    
                     d2Tdxdr = (1-flambda(i,k))*rho(i,k)*ekht(i,k)*ekh(i,k)
     &                                         *((((T(i,k)-T(i,km))/dz)-((T(im,k)-T(im,km))/dz) )/(Ru(i)-Ru(im)))**2
                     putout(i,k) = putout(i,k) + (1.80*etnew(i,k)/(ktnew(i,k)+1.0e-20)*Pkt(i,k)
     &                                          + 0.72*etnew(i,k)/ knew(i,k)* Pk(i,k) + d2Tdxdr)     /rho(i,k)
                     dimpl(i,k)  = dimpl(i,k)  + 2.20*etnew(i,k)/(ktnew(i,k)+1.0e-20) + 0.8*enew(i,k)/knew(i,k)           
                                            ! note, CD1*rho*et^2/kt/(rho*et)+CD2 rho e et/k/(rho*et), set implicit and divided by density
                  
                  elseif ((tempturbmod.eq.3).or.(tempturbmod.eq.4).or.(tempturbmod.eq.5)) then
                  !-------------------------
                  ! Deng, Wu and Xi model 2001 
                  
                     !time scales: temperature and mix
                     Ttemp  = (ktnew(i,k)+1.0e-20)/(etnew(i,k)+1.0e-20)
                     Tmix   = (Tt(i,k)*Ttemp)**0.5
                     !functions
                     fd1  = 1 - (exp(-Reeps(i,k)/1.7))**2.0
                     feps = 1 - 0.3*exp(-((Ret(i,k)/6.5)**2.0))   
                     fd2  = (1/0.9)*(ce2*feps-1.0)*(1 - (exp(-Reeps(i,k)/5.8))**2.0)  

                     putout(i,k) = putout(i,k) + (2.34*Pkt(i,k)/Tmix)/rho(i,k)
                     
                     ! original CD1=2.0, changed to 1.5   (cdiss1 defined in param.txt)
                     dimpl(i,k)  = dimpl(i,k)  + cdiss1*fd1/Ttemp + 0.9*fd2/Tt(i,k)    
                            ! note, (CD1*fd1/Ttemp + CD2*fd2/Tmomentum)*rho*et/(rho*et) set implicit and divided by density

                  endif
               endif
            
            endif


         enddo
      enddo

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
      real*8 lamcp(0:i1,0:k1) !gustavo
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
            call splint(enthTab,lamocpTab,lamocp2Tab,nTab,enth(i,k),lam(i,k),tabkhi,tabklo)   !! lambda/cp 
            call splint(enthTab,tempTab,  temp2Tab,  nTab,enth(i,k),tp(i,k),tabkhi,tabklo)
            call splint(enthTab,betaTab,  beta2Tab,  nTab,enth(i,k),be(i,k),tabkhi,tabklo)
            mu(i,k)  = mu(i,k)/Re
            lam(i,k) = lam(i,k)/(Re*Pr)                                                       !! (lambda/cp)/(Re*Pr) 
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
            ekhi(i,k) = conface/cpface/(Re*Pr)                                              !! can calculate lambda_i/RePr=ekhi*Cpi
            Cpi(i,k)=cpface
!            betai(i,k)=beface
            call splint(enthTab,muTab, mu2Tab, nTab,enthface,muface, tabkhi,tabklo)
            ekmi(i,k)  = muface/Re

            enthface = 0.5*(enth(i,k)+enth(i,k+1))
            call splint(enthTab,cpTab, cp2Tab, nTab,enthface,cpface,tabkhi,tabklo)
            call splint(enthTab,lamTab,lam2Tab,nTab,enthface,conface,tabkhi,tabklo)
            call splint(enthTab,betaTab,  beta2Tab,nTab,enthface,beface,tabkhi,tabklo)
            ekhk(i,k) = conface/cpface/(Re*Pr)
            Cpk(i,k)=cpface                                                                 !! can calculate lambda_k/RePr=ekhk*Cpk
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
     1           ( Ru(i)*dUdt(i,k) - Ru(i-1)*dUdt(i-1,k) )/( Rp(i)*dr(i))
     3           + (     dWdt(i,k) -         dWdt(i,k-1) )/( dz         ) )/dt
     &           +      (rnew(i,k)-rold(i,k))/(dt*dt)
            qcrit(i,k) = p(i,k)*dt

            sumps = sumps + p(i,k)*dr(i)*dz
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
             
            ! lenght scale for v2f model
            Lh(i,k)=0.23*max(knew(i,k)**1.5/enew(i,k),70.*((ekm(i,k)/rnew(i,k))**3./enew(i,k))**0.25)
            
            if (modVF.eq.1) then
               ! Modifications: Lien&Kalitzin 2001 "Computations of transonic flow with the v2f turbulence model"
               Tt(i,k)   = max(Tt(i,k), 1.0e-8)
               Tt(i,k)   = min(Tt(i,k),0.6*knew(i,k)/(3.**0.5*v2new(i,k)*cmu*(2.*Srsq(i,k))**0.5))

               Lh(i,k)=min(knew(i,k)**1.5/enew(i,k),knew(i,k)**1.5/(3.**0.5*v2new(i,k)*cmu*(2.*Srsq(i,k))**0.5))
               Lh(i,k)=0.23*max(Lh(i,k),70.*((ekm(i,k)/rnew(i,k))**3./enew(i,k))**0.25)
            endif

!            fv2(i,k)= - (1.4-1.)*(2./3.-v2new(i,k)/knew(i,k))/Tt(i,k)
!     &                - 0.3*(Pk(i,k))/(rnew(i,k)*knew(i,k))-5.*v2new(i,k)/(knew(i,k)*Tt(i,k))
            ! f-equation also has Gk: Kenjeres et al 2005 "Contribution to elliptic relaxation modelling of turbulent natural and mixed convection"
            fv2(i,k)= - (1.4-1.)*(2./3.-v2new(i,k)/knew(i,k))/Tt(i,k)
     &                - 0.3*(Pk(i,k)+Gk(i,k))/(rnew(i,k)*knew(i,k))-5.*v2new(i,k)/(knew(i,k)*Tt(i,k))
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
      real*8  dr2,dz2,df,df2,kcoeff,tmp1,tmp2,tmp3,Courant,dtmp
      dt = dtmax
      Courant = 10.5
      do k=1,kmax
         do i=1,imax
            dr2 = dr(i) * dr(i)
            dz2 = dz    * dz
            kcoeff =ekm(i,k)
            tmp1 = ( abs(Unew(i,k)) / ( Rp(i+1)-Rp(i) ) ) +
     &           ( abs(Wnew(i,k)) /         dz        )
            tmp2 = (1.0/dr2 + 1.0/dz2)
            tmp3 = 1.0 / ( 1.0 * tmp2 * kcoeff + tmp1 )
            tmp3 = Courant*tmp3 
            dt = min( dt , tmp3 )
         enddo
      enddo
      dtmp = dt
      call mpi_allreduce(dtmp,dt,1,mpi_real8,mpi_min,mpi_comm_world,ierr)
      return
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
     &           (Wnew(i,k)*rhokp-Wnew(i,k-1)*rhokm)*dr(i)+
     &           (rNew(i,k)-rold(i,k))/dt*dr(i)*dz

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
      const = 0.65
!     const = 0.50
c******************************************************************
      pi    = 4.0*atan(1.0)
      Rei   = 1.0/Re
      dz    = 1.0*LoD/(kmax*px)
      ru(0) = 0

      fA = 0.12
      fB = 2.4

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
         Rp(i) = (Ru(i)+Ru(i-1))/2.0
         dr(i) = (Ru(i)-Ru(i-1))
      enddo

      dr(i1) = dr(imax)
      Ru(i1) = Ru(imax) + dr(i1)
      Rp(i1) = Ru(imax) + dr(i1)/2.0
      dr(0)  = dr(1)
      Rp(0)  = Ru(0) - dr(0)/2.0

      if (rank.eq.0) then
         open(11,file = 'grid.txt')
         write(11,*) Re, imax
         do i=1,imax
            Yplus = (0.5-Rp(i))*Re
            write(11,'(i5,4F12.6)') i,yplus,Ru(i),Rp(i),delta(i)
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
      subroutine diffc(putout,putin,ek,eki,ekk,ekmt,sigma,rho,Ru,Rp,dr,dz,rank1,diffVersion)
      implicit none
      include 'param.txt'
      integer   km,kp,rank1,diffVersion
      real*8     putout(0:i1,0:k1),putin(0:i1,0:k1),
     &     rho(0:i1,0:k1),ek(0:i1,0:k1),eki(0:i1,0:k1),ekk(0:i1,0:k1),
     &     ekmt(0:i1,0:k1),dr(0:i1),dz,Ru(0:i1),Rp(0:i1),sigma
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

      !>********************************************************************
      !!  diffusion term for SA model: in the z-direction as, plus extra for Aupoix modifications...
      !!********************************************************************
      subroutine diffcSA(putout,putin,ek,eki,ekk,ekmt,sigma,rho,Ru,Rp,dr,dz,rank1,diffVersion)
      implicit none
      include 'param.txt'
      integer   km,kp,im,ip,rank1,diffVersion
      real*8     putout(0:i1,0:k1),putin(0:i1,0:k1),
     &     rho(0:i1,0:k1),ek(0:i1,0:k1),eki(0:i1,0:k1),ekk(0:i1,0:k1),
     &     ekmt(0:i1,0:k1),dr(0:i1),dz,Ru(0:i1),Rp(0:i1),sigma
!          Important note: this function takes instead of ek, eki, ekk, ekmt: eknu, eknui, eknuk, nuSANew, respectively.
         ! For, Standard, Inverse SLS and Aupoix
         ! rho=1 for standard
      do k=1,kmax
          kp=k+1
          km=k-1
          do i=1,imax
             putout(i,k) = putout(i,k) + 1.0/rho(i,k)*(
     3          ( (ekk(i,k ) + 0.5*(ekmt(i,k)+ekmt(i,kp))/sigma)*sqrt(0.5*(rho(i,k)+rho(i,kp)))*((rho(i,kp)**0.5)*putin(i,kp)-(rho(i,k )**0.5)*putin(i,k ))
     3           -(ekk(i,km) + 0.5*(ekmt(i,k)+ekmt(i,km))/sigma)*sqrt(0.5*(rho(i,k)+rho(i,km)))*((rho(i,k )**0.5)*putin(i,k )-(rho(i,km)**0.5)*putin(i,km))
     &                 )/(dz*dz)   )
          enddo
      enddo
        
      if (diffVersion == 2) then       ! For Aupoix we need to substract the density gradient diffusion with molecular viscosity
      !-1/rho d/dx[nu*nusa/2 drhodx]
      ! in the z-direction
          do k=1,kmax
              kp=k+1
              km=k-1
              do i=1,imax
                 putout(i,k) = putout(i,k) - 1.0/rho(i,k)*((ekk(i,k )*0.5*(ekmt(i,k)+ekmt(i,kp))/2*(rho(i,kp)-rho(i,k ))
     &                                                     -ekk(i,km)*0.5*(ekmt(i,k)+ekmt(i,km))/2*(rho(i,k )-rho(i,km)))/(dz*dz))
              enddo
          enddo
      ! in the r-direction
          do i=1,imax
              ip=i+1
              im=i-1
              do k=1,kmax
                 putout(i,k) = putout(i,k) - 1.0/rho(i,k)* (  Ru(i )/((Rp(ip)-Rp(i ))*Rp(i)*dr(i))*(eki(i ,k)*0.5*(ekmt(i,k)+ekmt(ip,k))/2*(rho(ip,k)-rho(i ,k)))
     &                                                     -  Ru(im)/((Rp(i )-Rp(im))*Rp(i)*dr(i))*(eki(im,k)*0.5*(ekmt(i,k)+ekmt(im,k))/2*(rho(i ,k)-rho(im,k)))  
     &                                                     )
              enddo
          enddo
      endif




      end

      !>********************************************************************
      !! diffusion term for kine of the SST model in the z-direction, set as a source term...
      !!********************************************************************
      subroutine diffcSSTKine(putout,putin,ek,eki,ekk,ekmt,sigma,rho,Ru,Rp,dr,dz,rank1,diffVersion)
      implicit none
      include 'param.txt'
      integer   km,kp,rank1,diffVersion
      real*8     putout(0:i1,0:k1),putin(0:i1,0:k1),
     &     rho(0:i1,0:k1),ek(0:i1,0:k1),eki(0:i1,0:k1),ekk(0:i1,0:k1),
     &     ekmt(0:i1,0:k1),dr(0:i1),dz,Ru(0:i1),Rp(0:i1),sigma(0:i1,0:k1)
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
      subroutine diffcSSTOmega(putout,putin,ek,eki,ekk,ekmt,sigma,rho,Ru,Rp,dr,dz,rank1,diffVersion)
      implicit none
      include 'param.txt'
      integer   km,kp,rank1,diffVersion
      real*8     putout(0:i1,0:k1),putin(0:i1,0:k1),
     &     rho(0:i1,0:k1),ek(0:i1,0:k1),eki(0:i1,0:k1),ekk(0:i1,0:k1),
     &     ekmt(0:i1,0:k1),dr(0:i1),dz,Ru(0:i1),Rp(0:i1),sigma(0:i1,0:k1)
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


      !>********************************************************************
      !! diffusion term for epsilon in the z-direction, set as a source term...
      !!********************************************************************
      subroutine diffEPS(putout,putin,ek,eki,ekk,ekmt,sigma,rho,Ru,Rp,dr,dz,rank1,diffVersion)
      implicit none
      include 'param.txt'
      integer   km,kp, rank1,diffVersion
      real*8     putout(0:i1,0:k1),putin(0:i1,0:k1),
     &      rho(0:i1,0:k1),ek(0:i1,0:k1),eki(0:i1,0:k1),ekk(0:i1,0:k1),
     &      ekmt(0:i1,0:k1),dr(0:i1),dz,Ru(0:i1),Rp(0:i1),sigma, difcp, difcm

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
      subroutine diffu (putout,Uvel,Wvel,ekme,Ru,Rp,dr,dz,i1,k1,dif)
      implicit none

      integer  i,k,im,ip,km,kp,i1,k1,ib,ie,kb,ke
      real*8     putout(0:i1,0:k1),Uvel(0:i1,0:k1),Wvel(0:i1,0:k1),
     &     ekme(0:i1,0:k1),dr(0:i1),dz,Ru(0:i1),Rp(0:i1),
     &     epop,epom,drp,dzi,divUim,divUip,divUi,dif
c     
      ib = 1
      ie = i1-1

      kb = 1
      ke = k1-1
c     
      dzi =1./dz
      do k=kb,ke
         kp=k+1
         km=k-1
         do i=ib,ie
            ip=i+1
            im=i-1
            drp = Rp(ip)-Rp(i)

            epop = 0.25*(ekme(i,k)+ekme(ip,k) + ekme(ip,kp) + ekme(i,kp))
            epom = 0.25*(ekme(i,k)+ekme(ip,k) + ekme(ip,km) + ekme(i,km))

            divUim = (Ru(i)*Uvel(i,k) - Ru(im)*Uvel(im,k))/(Rp(i)*dr(i))
     &           + (        Wvel(i,k) -        Wvel(i,km))/dz

            divUip = (Ru(ip)*Uvel(ip,k)-Ru(i)*Uvel(i,k))/(Rp(ip)*dr(ip))
     &           + (         Wvel(ip,k) -       Wvel(ip, km))/dz

            divUi = ( Rp(ip)*(Uvel(ip,k)+Uvel(i,k))
     &           -Rp(i )*(Uvel(i ,k)+Uvel(im,k)) )/(2.*Ru(i)*(Rp(ip)-Rp(i)))
     &           +((Wvel(ip,k)+Wvel(i,k))-(Wvel(ip,km)+Wvel(i,km)))/(2.*dz)


            putout(i,k) = putout(i,k) +
     1           ( Rp(ip)*ekme(ip,k)*(dif*(Uvel(ip,k)-Uvel(i,k) )/
     1           dr(ip)-1./3.*divUip)-Rp(i )*ekme(i,k) *(dif*(Uvel(i,k)
     1           -Uvel(im,k))/dr(i) -1./3.*divUim))/(0.5*Ru(i)*(drp))
     &           +
     3           ( epop * (   (Uvel(i,kp)  - Uvel(i,k) ) * dzi
     3           + (Wvel(ip,k)  - Wvel(i,k) ) / (Rp(ip) - Rp(i))
     3           )             -
     3           epom * (   (Uvel(i,k)   - Uvel(i,km)) * dzi
     3           + (Wvel(ip,km) - Wvel(i,km)) / (Rp(ip) - Rp(i))
     3           ) ) * dzi
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
      subroutine diffw(putout,Uvel,Wvel,ekme,Ru,Rp,dr,dz,i1,k1,dif,rank)
      implicit none
     

      integer  i,k,im,ip,km,kp,i1,k1,ib,ie,kb,ke,rank
      real*8     putout(0:i1,0:k1),Uvel(0:i1,0:k1),Wvel(0:i1,0:k1),
     &     ekme(0:i1,0:k1),dr(0:i1),dz,Ru(0:i1),Rp(0:i1),
     &     epop,emop,divUkm,divUkp,dif
c     
      ib = 1
      ie = i1-1

      kb = 1
      ke = k1-1
c     
      do k=kb,ke
         kp=k+1
         km=k-1

         do i=ib,ie
            ip=i+1
            im=i-1
            epop = 0.25*(ekme(i,k)+ekme(i,kp) + ekme(ip,k) + ekme(ip,kp) )
            emop = 0.25*(ekme(i,k)+ekme(i,kp) + ekme(im,k) + ekme(im,kp) )

            divUkm = (Ru(i)*Uvel(i,k) - Ru(im)*Uvel(im,k))/(Rp(i)*dr(i))
     &           + (      Wvel(i,k) -        Wvel(i,km))/dz

            divUkp = (Ru(i)*Uvel(i,kp)- Ru(im)*Uvel(im, kp))/(Rp(i)*dr(i))
     &           + (      Wvel(i,kp)-        Wvel(i,  k ))/dz

!     divUkm = 0.0
!     divUkp = 0.0

            putout(i,k) =  putout(i,k)+
     1           (Ru(i )*epop*( (Uvel(i,kp)  - Uvel(i,k) ) / dz
     1           +dif*(Wvel(ip,k)  - Wvel(i,k) ) / (Rp(ip)-Rp(i)))
     1           -
     1           Ru(im)*emop*( (Uvel(im,kp) - Uvel(im,k)) / dz
     1           +dif*(Wvel(i,k)   - Wvel(im,k)) / (Rp(i)-Rp(im)))
     1           ) / ( Rp(i) * dr(i) )
     &           +
     3           (2.*ekme(i,kp)*((Wvel(i,kp)-Wvel(i,k ))/dz - 1./3.*divUkp)-
     3           2.*ekme(i,k )*((Wvel(i,k )-Wvel(i,km))/dz - 1./3.*divUkm))/dz
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
      subroutine advecc(putout,dimpl,putin,U,W,Ru,Rp,dr,dz,i1,k1,rank,periodic,flagImpl)

      implicit none

      integer  i,k,im,ip,km,kp,i1,k1,ib,ie,kb,ke,rank,periodic
      real*8 putout(0:i1,0:k1),putin(0:i1,0:k1),dimpl(0:i1,0:k1)
      real*8 U(0:i1,0:k1),W(0:i1,0:k1),dr(0:i1),dz,Ru(0:i1),Rp(0:i1)
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
            putout(i,k) =     - (Ru(i)*U(i,k)*cu(i,k) - Ru(im)*U(im,k)*cu(im,k))/(Rp(i)*dr(i))
     &                        - (      W(i,k)*cw(i,k) -        W(i,km)*cw(i,km))/(dz)
     4         + putin(i,k)*(   (Ru(i)*U(i,k) - Ru(i-1)*U(i-1,k ))/(Rp(i)*dr(i))
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
      subroutine advecrho(putout,putin,U,W,Ru,Rp,dr,dz,i1,k1,rank)

  
      implicit none
      integer  i,k,im,ip,km,kp,i1,k1,ib,ie,kb,ke,rank
      real*8 putout(0:i1,0:k1),putin(0:i1,0:k1)
      real*8 U(0:i1,0:k1),W(0:i1,0:k1),dr(0:i1),dz,Ru(0:i1),Rp(0:i1), rnew(0:i1,0:k1)
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
      subroutine advecu(putout,Uvel,Wvel,RHO,Ru,Rp,dr,dz,i1,k1)
      implicit none
c     

      integer  i,k,im,ip,km,kp,i1,k1,ib,ie,kb,ke
      real*8     putout(0:i1,0:k1),Uvel(0:i1,0:k1),Wvel(0:i1,0:k1),
     &     dr(0:i1),dz,Ru(0:i1),Rp(0:i1)
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
     1           Rp(i )*(Uvel(im,k)+Uvel(i,k))*(Uvel(i,k)+Uvel(im,k))
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
      subroutine advecw(putout,Uvel,Wvel,RHO,Ru,Rp,dr,dz,ekm,peclet_z)
      implicit none
      include 'param.txt'

      integer   im,ip,km,kp,ib,ie,kb,ke
      real*8     putout(0:i1,0:k1),Uvel(0:i1,0:k1),Wvel(0:i1,0:k1),
     &     dr(0:i1),dz,Ru(0:i1),Rp(0:i1)
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
     1           / ( Rp(i) * dr(i) ) +advcecw_w)
         enddo
      enddo
      return
      end

