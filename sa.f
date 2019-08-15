!>******************************************************************************************
!>******************************************************************************************
!>******************************************************************************************
!! SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA 
!! SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA 
!! SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA 
!! SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA SA 
!!              SA =  parameter (turbmod           = 1) 
!!
!>******************************************************************************************
!!      SA routine to estimate the eddy viscosity
!!******************************************************************************************
      subroutine calculate_mut_SA(U,W,ekmetmp,ekmttmp,ekmtin,step)

      implicit none
      include 'param.txt'
      include 'common.txt'
      integer  im,ip,km,kp,step
      real*8   tauwLoc, tauw(0:k1) 
      real*8, dimension(0:i1,0:k1) :: U,W,ekmetmp,ekmttmp
      real*8, dimension(0:i1) :: ekmtb,ekmtf,ekmtin
      real*8  cv1_3,chi,fv1SA, wallD

      cv1_3 = (7.1)**3.0

      do k=1,kmax
         km=k-1
         kp=k+1

         tauw(k) = ekmi(imax,k)*0.5*(W(imax,km)+W(imax,k))/wallDist(imax)

         do i=1,imax
            im=i-1
            ip=i+1

            yp(i,k)     = sqrt(rNew(i,k))/ekm(i,k)*(wallDist(i))*tauw(k)**0.5           ! ystar
  !          yp(i,k)     = sqrt(rNew(imax,k))/ekm(imax,k)*wallDist(i)*tauw(k)**0.5    ! yplus
  !          yp(i,:)     = (wallDist(i))*Re*(1/Re*(Win(imax)/(wallDist(imax))))**0.5        ! yplus
            ReTauS(i,k) = 0.5*sqrt(rNew(i,k))/ekm(i,k)*tauw(k)**0.5
            Ret(i,k)    = rNew(i,k)*(kNew(i,k)**2.)/(ekm(i,k)*eNew(i,k))        ! not sure if r2 or r
        
            
            chi          = nuSAnew(i,k)/(ekm(i,k)/rNew(i,k))
            fv1SA        = (chi**3.0)/(chi**3.0 + cv1_3);
            ekmttmp(i,k) = min(100.0,max(0.0,rNew(i,k)*nuSAnew(i,k)*fv1SA))

         enddo
      enddo


      end


!>******************************************************************************************
!!      SA prodis subroutine which calculates the production term of the turbulent scalar equation
!>******************************************************************************************
      subroutine prodisSA(nuSAtmp,U,W,T,rho)
      implicit none
      include 'param.txt'
      include 'common.txt'

      integer im,ip,km,kp,ib,ie,kb,ke !< integers
      real*8, dimension(0:i1,0:k1) :: U,W,T,rho,div,nuSAtmp!,Tt
      real*8  cv1_3,cb1,cb2,cb3,cw1,cw2,cw3_6,inv_cb3,kappa_2,chi,fv1SA,fv2SA,r_SA,g_SA,fw_SA,StR,shatSA

      cv1_3     = (7.1)**3.0
      cb1       = 0.1355
      cb2       = 0.622
      cb3       = 2.0/3.0
      inv_cb3   = 1.0/cb3
      kappa_2   = (0.41)**2.0   ! von karman constant
      cw1       = (cb1/kappa_2) + (1.0+cb2)/cb3
      cw2       = 0.3
      cw3_6     = (2.0)**6.0


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

            ! Bouyancy prodution and time scale, not defined for this model
            Gk(i,k)=0
            Tt(i,k)=1

            ! magnitude of rate of rotation: omega=sqrt(2*Wij*Wij), Wrz = 0.5*(dU/dz-dW/dr);  note, utheta=0 d/dtheta=0
            StR = ( ( -( (W(ip,km)+W(ip,k)+W(i,km)+W(i ,k))/4.-(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/dRu(i)
     &                +( (U(i,kp) +U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz)           )**2.)

            StR = StR**0.5

            ! calculating Shat from SA model
            chi         = nuSAtmp(i,k)/(ekm(i,k)/rho(i,k))
            fv1SA       = (chi**3.0)/((chi**3.0) + cv1_3);
            fv2SA       = 1.0 - (chi/(1.0 + (chi*fv1SA)))
            ShatSA      = StR + fv2SA*nuSAtmp(i,k)/(kappa_2*(wallDist(i)**2.0))

            ! production term in SA model
            Pk(i,k) = cb1*nuSAtmp(i,k)*ShatSA

            div(i,k) =(Ru(i)*U(i,k)-Ru(im)*U(im,k))/(Rp(i)*dru(i))
     &              +(      W(i,k) -      W(i,km))/dz

         enddo
      enddo

      end

!>******************************************************************************************
!!      To calculate the rhs of the nuSA equation
!>******************************************************************************************
      subroutine rhs_SA(putout,dimpl,nuSAtmp,rho)
      implicit none
      include 'param.txt'
      include 'common.txt'

      integer im,ip,km,kp,ib,ie,kb,ke !< integers
      real*8, dimension(0:i1,0:k1) :: putout,rho,div,nuSAtmp,dimpl!,Tt
      real*8  cv1_3,cb1,cb2,cb3,cw1,cw2,cw3_6,inv_cb3,kappa_2,chi,fv1SA,fv2SA,r_SA,g_SA,fw_SA,StR,ShatSA

      cv1_3     = (7.1)**3.0
      cb1       = 0.1355
      cb2       = 0.622
      cb3       = 2.0/3.0
      inv_cb3   = 1.0/cb3
      kappa_2   = (0.41)**2.0   ! von karman constant
      cw1       = (cb1/kappa_2) + (1.0+cb2)/cb3
      cw2       = 0.3
      cw3_6     = (2.0)**6.0


      ib = 1
      ie = i1-1

      kb = 1
      ke = k1-1
      if ((modifDiffTerm == 1) .or. (modifDiffTerm == 2)) then
          do k=kb,ke
             kp=k+1
             km=k-1
             do i=ib,ie
                ip=i+1
                im=i-1
                
                ShatSA = Pk(i,k)/(cb1*nuSAtmp(i,k))
    
                ! destruction term in SA model
                r_SA         = min(nuSAtmp(i,k)/(kappa_2*(wallDist(i)**2.0)*ShatSA), 10.0)
                g_SA         = r_SA + cw2*((r_SA**6.0) - r_SA)
                fw_SA        = g_SA*(((1.0 + cw3_6)/(g_SA**6.0 + cw3_6))**(1.0/6.0))
    
                ! gustavo: i think this is not correct
                !destrSA(i,k) = cw1/rho(i,k)*fw_SA*nuSAtmp(i,k)/(wallDist(i)**2)
                dimpl(i,k) = dimpl(i,k) + cw1*fw_SA*nuSAtmp(i,k)/(wallDist(i)**2.0)

                ! source term
                ! invSLS and Aupoix SA model=  advection + Pk + (1/rho)*cb2/cb3*(d(nuSA*sqrt(rho))/dr)^2 +(d(nuSA*sqrt(rho))/dz)^2
                putout(i,k) = putout(i,k) + Pk(i,k) + cb2*inv_cb3/rho(i,k) * (
     &         (((nuSAtmp(ip,k)*(rho(ip,k)**0.5)) - (nuSAtmp(im,k)*(rho(im,k)**0.5)))/(dRp(i)+dRp(im)))**2.0
     &         +(((nuSAtmp(i,kp)*(rho(i,kp)**0.5)) - (nuSAtmp(i,km)*(rho(i,km)**0.5)))/(2.0*dz))**2.0  )

    
             enddo
          enddo
      else
          do k=kb,ke
             kp=k+1
             km=k-1
             do i=ib,ie
                ip=i+1
                im=i-1

                ShatSA = Pk(i,k)/(cb1*nuSAtmp(i,k))
    
                ! destruction term in SA model
                r_SA         = min(nuSAtmp(i,k)/(kappa_2*(wallDist(i)**2.0)*ShatSA), 10.0)
                g_SA         = r_SA + cw2*((r_SA**6.0) - r_SA)
                fw_SA        = g_SA*(((1.0 + cw3_6)/(g_SA**6.0 + cw3_6))**(1.0/6.0))
    
                ! gustavo: i think this is not correct
                !destrSA(i,k) = cw1/rho(i,k)*fw_SA*nuSAtmp(i,k)/(wallDist(i)**2)
                dimpl(i,k) = dimpl(i,k) + cw1*fw_SA*nuSAtmp(i,k)/(wallDist(i)**2.0)
    
                ! source term
                ! Conventional SA model=  advection + Pk + cb2/cb3*(dnuSA/dr)^2 +(dnuSA/dz)^2
                putout(i,k) = putout(i,k) + Pk(i,k) + cb2*inv_cb3 * (
     &         ((nuSAtmp(ip,k) - nuSAtmp(im,k))/(dRp(i)+dRp(im)))**2.0 + ((nuSAtmp(i,kp) - nuSAtmp(i,km))/(2.0*dz))**2.0  )
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
      subroutine advanceSA(resSA,Utmp,Wtmp,Rtmp,rho3,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      real*8 dnew(0:i1,0:k1),tempArray(0:i1,0:k1),dimpl(0:i1,0:k1)
      real*8 Utmp(0:i1,0:k1),Wtmp(0:i1,0:k1),Rtmp(0:i1,0:k1)
      real*8 rho3(0:i1,0:k1), eknu(0:i1,0:k1),eknui(0:i1,0:k1),eknuk(0:i1,0:k1)
      real*8 cb3

      real*8     a  (imax)
      real*8     b  (imax)
      real*8     c  (imax)
      real*8     rhs(imax)

      real*8 resSA 

      cb3 = 2.0/3.0
      resSA = 0.0
      dnew=0.0; dimpl = 0.0;
      call advecc(dnew,dimpl,nuSANew,utmp,wtmp,Ru,Rp,dru,dz,i1,k1,rank,periodic,.true.)
      call rhs_SA(dnew,dimpl,nuSANew,Rtmp)

      do k=0,kmax+1
         do i=0,imax+1
            eknu (i,k) = ekm(i,k)/Rtmp(i,k)
            if (i.LE.imax) eknui(i,k) = ekmi(i,k)/(0.5*(Rtmp(i,k)+Rtmp(i+1,k))) 
            if (k.LE.kmax) eknuk(i,k) = ekmk(i,k)/(0.5*(Rtmp(i,k)+Rtmp(i,k+1)))
          enddo
      enddo

      tempArray = 0.0
      call diffcSA(tempArray,nuSANew,eknu,eknui,eknuk,nuSANew,1.0,rho3,Ru,Rp,dru,dz,rank,modifDiffTerm)
      dnew = dnew + tempArray/cb3

      !> diffusion term in the r-direction, set implicit!
      !! conv.    d/dy (nu+nuSA) dnuSA/dy
      !! modified (1/rho) d/dy [sqrt(rho) (nu+nuSA) d(sqrt(rho)nuSA)/dy]
      do k=1,kmax
         do i=1,imax

            a(i) = (eknui(i-1,k)+0.5*(nuSANew(i,k)+nuSANew(i-1,k)))/cb3*((0.5*(rho3(i-1,k)+rho3(i,k)))**0.5)
            a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/rho3(i,k)
            c(i) = (eknui(i  ,k)+0.5*(nuSANew(i,k)+nuSANew(i+1,k)))/cb3*((0.5*(rho3(i+1,k)+rho3(i,k)))**0.5)
            c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/rho3(i,k)

            b(i) = ((-a(i)-c(i))*(rho3(i,k)**0.5)  +  dimpl(i,k)  )/alphak

            a(i) = a(i)*(rho3(i-1,k)**0.5)
            c(i) = c(i)*(rho3(i+1,k)**0.5)

            rhs(i) = dnew(i,k) + (1-alphak)*b(i)*nuSANew(i,k)
         enddo

         i=1
         b(i) = b(i)+centerBC*a(i)

         i=imax
         b(i) = b(i) - (c(i) /alphak)
         !b(i) = ((-(a(i)/(rho3(i-1,k)**0.5))-(c(i)/(rho3(i+1,k)**0.5)))*(rho3(i,k)**0.5) - c(i) +  dimpl(i,k)   )/alphak
         !b(i) = ((-a(i)-c(i))*(rho3(i,k)**0.5) - c(i)*rho3(i+1,k)**0.5 +  dimpl(i,k)   )/alphak
         rhs(i) = dnew(i,k) + (1-alphak)*b(i)*nuSANew(i,k)

         call matrixIdir(imax,a,b,c,rhs)

         do i=1,imax
            resSA = resSA + ((nuSANew(i,k) - rhs(i))/(nuSANew(i,k)+1.0e-20))**2.0
            nuSANew(i,k) = max(rhs(i), 1.0e-8)
         enddo
      enddo
      end

!>******************************************************************************************
!!      SA advancing the turbulence scalars of this model: nuSA
!!******************************************************************************************
      subroutine advanceScalar_SA(resSA,Utmp,Wtmp,Rtmp,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'

      real*8 Utmp(0:i1,0:k1),Wtmp(0:i1,0:k1),Rtmp(0:i1,0:k1),ShatSA(0:i1,0:k1)
      real*8 rho3(0:i1,0:k1)
      real*8 resSA 

      ! modified turb. model
      !    modifDiffTerm = 1, our modification
      !    modifDiffTerm = 2, Aupoix modification
      if ((modifDiffTerm == 1) .or. (modifDiffTerm == 2)) then
         rho3 = Rtmp
      else
         rho3 = 1.0
      endif

      call prodisSA(nuSANew,Utmp,Wtmp,temp,Rtmp)
      call advanceSA(resSA,Utmp,Wtmp,Rtmp,rho3,rank)

      
      end

      !>********************************************************************
      !!  diffusion term for SA model: in the z-direction as, plus extra for Aupoix modifications...
      !!********************************************************************
      subroutine diffcSA(putout,putin,ek,eki,ekk,ekmt,sigma,rho,Ru,Rp,dru,dz,rank1,diffVersion)
      implicit none
      include 'param.txt'
      integer   km,kp,im,ip,rank1,diffVersion
      real*8     putout(0:i1,0:k1),putin(0:i1,0:k1),
     &     rho(0:i1,0:k1),ek(0:i1,0:k1),eki(0:i1,0:k1),ekk(0:i1,0:k1),
     &     ekmt(0:i1,0:k1),dru(0:i1),dz,Ru(0:i1),Rp(0:i1),sigma
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
                 putout(i,k) = putout(i,k) - 1.0/rho(i,k)* (  Ru(i )/((Rp(ip)-Rp(i ))*Rp(i)*dru(i))*(eki(i ,k)*0.5*(ekmt(i,k)+ekmt(ip,k))/2*(rho(ip,k)-rho(i ,k)))
     &                                                     -  Ru(im)/((Rp(i )-Rp(im))*Rp(i)*dru(i))*(eki(im,k)*0.5*(ekmt(i,k)+ekmt(im,k))/2*(rho(i ,k)-rho(im,k)))
     &                                                     )
              enddo
          enddo
      endif




      end
