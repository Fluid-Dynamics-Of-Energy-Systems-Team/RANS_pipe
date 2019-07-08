
!>    ****************************************************
!!    Main file of the code
!!    This code simulates Super Critical Fluid in a
!!    heated pipe with constant heat flux
!!    ****************************************************

      implicit none

      include      'param.txt'          !> input file
      include      'common.txt'         !> def. of turb model variables
      include      'mpif.h'             !> mpi stuff

      ! def. of variables
      integer      ploc,ierr,istart,noutput
      integer      iTabFoundHi,iTabFoundLo
      real*8       bulk,stress,stime,time1,time2,timer,time3,dif,adv
      real*8       newTemp,newRho,newMu,newLam,newCP,newEnth,
     &     newTemp2,enth_i1,enth_imax,fxvalue,str,str_tot 
      real*8       resC,resK,resE,resV2,resOm,resSA   ! Residuals for energy, kine, eps, v2, omega, nuSA
      real*8 	   Win(0:i1),kin(0:i1),ein(0:i1),ekmtin(0:i1),v2in(0:i1),omIn(0:i1),nuSAin(0:i1)

      call cpu_time(time1)
      call mpi_init(ierr)
      call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)
      call mpi_comm_size(MPI_COMM_WORLD,ploc,ierr)
      call init_transpose

      dtmax = 1.e-3

      ! generating the table for thermodynamic quantities...
      call readTable(rank)

      call spline(enthTab, rhoTab,    nTab, rho2Tab)
      call spline(enthTab, muTab,     nTab, mu2Tab)
      call spline(enthTab, lamTab,    nTab, lam2Tab)
      call spline(enthTab, cpTab,     nTab, cp2Tab)
      call spline(enthTab, lamocpTab, nTab, lamocp2Tab)
      call spline(enthTab, tempTab,   nTab, temp2Tab)
      call spline(enthTab, betaTab,   nTab, beta2Tab)

      call mkgrid(rank)

      dt = dtmax
       
      istart = 1

      if (select_init < 2) then
         call fkdat(rank)
         istart=1
      else
         call loadRestart(istart,rank)
         istart = istart+1
      endif

      ! periodic = 1, turb flow generator
      ! periodic = 2, heated pipe 
!      if (periodic.ne.1) then
!         if (turbmod.eq.0) open(29,file =  '0/Inflow',form='unformatted')
!         if (turbmod.eq.1) open(29,file = 'MK/Inflow',form='unformatted')
!         if (turbmod.eq.2) open(29,file = 'LS/Inflow',form='unformatted')
!         if (turbmod.eq.3) open(29,file = 'VF/Inflow',form='unformatted')
!         if (turbmod.eq.4) open(29,file = 'SA/Inflow',form='unformatted')
!         if (turbmod.eq.5) open(29,file = 'OM/Inflow',form='unformatted')
!         read(29) Win(:),kin(:),ein(:),v2in(:),omIn(:),nuSAin(:),ekmtin(:),Pk(:,0)
!         close(29)
!      endif

      Win(:) = 1.0
      kin(:) = 0.0
      ein(:) = 0.0
      v2in(:) = 0.0
      omIn(:) = 0.0
      nuSAin(:) = 0.0
      ekmtin(:) = 0.0
      Pk(:,0) = 0.0


      call state(cnew,rnew,ekm,ekh,temp,beta,istart,rank);
      call bound_h(kin,ein,v2in,omIn,nuSAin,rank)
      call state(cnew,rnew,ekm,ekh,temp,beta,istart,rank);   ! necessary to call it twice
      rold = rnew

      call turbprop(Unew,Wnew,ekme,ekmt,ekmtin,rank,istep)
      call bound_v(Unew,Wnew,Win,rank)


      call chkdt(rank,istep)


      ! simulation loop 
      do istep=istart,nstep

         ! calculating turbulent viscosity
         call turbprop(Unew,Wnew,ekme,ekmt,ekmtin,rank,istep)

         if (turbmod.eq.3)  then
            call fillhm(rank)
            call SOLVEhelm(fv2,Ru,Rp,dz,rank,Lh)
         endif
 		 
         call advanceScalar(resC,resK,resE,resV2,resOm,resSA,Unew,Wnew,Rnew,fv2,rank)

         call bound_h(kin,ein,v2in,omIn,nuSAin,rank)
         call state(cnew,rnew,ekm,ekh,temp,beta,istep,rank)
      	

         call advance(rank)
         call bound_m(dUdt,dWdt,wnew,rnew,Win,rank)
         call fillps(rank)
         call SOLVEpois(p,Ru,Rp,dz,rank)
         call correc(rank,1)
         call bound_v(Unew,Wnew,Win,rank)

         if (mod(istep,10) .eq. 0)      call chkdiv(rank)

         call cmpinf(bulk,stress)
         call chkdt(rank,istep)

         if (mod(istep,1000 ).eq.0)     call outputProfile(rank)
         if (mod(istep,1000 ).eq.0)     call outputX_h(rank,istep)
         if (mod(istep,1000 ).eq.0)     call output2d(rank,istep)
         if (mod(istep,5000 ).eq.0)     call saveRestart(rank)

         if ((periodic.eq.1) .and. mod(istep,5000).eq.0) then
            call Inflow_output(rank,istep)   
         endif

         noutput = 100
         if (rank.eq.0) then
            if (istep.eq.istart .or. mod(istep,noutput*20).eq.0) then
                write(6,'(A7,9A14)') 'istep','dt','bulk',
     1           'stress','cResid','kineResid','epsResid','v2Resid','omResid','nuSAresid'
            endif
            if (istep.eq.istart .or. mod(istep,noutput).eq.0) then
                write(6,'(i7,9e14.5)') istep,dt,bulk,stress,resC,resK,resE,resV2,resOm,resSA
            endif
         end if


      enddo

      call outputProfile(rank)
      call output2d(rank,istep)
      call mpi_finalize(ierr)
      stop
      end





!>******************************************************************************************
!!      turbprop routine to estimate the eddy viscosity
!!******************************************************************************************
      subroutine turbprop(U,W,ekmetmp,ekmttmp,ekmtin,rank,step)

      implicit none
      include 'param.txt'
      include 'common.txt'
      integer  im,ip,km,kp,step
      real*8   tauwLoc, tauw(0:k1) 
      real*8, dimension(0:i1,0:k1) :: U,W,ekmetmp,ekmttmp,kd,Tt,Srsq
      real*8, dimension(0:i1) :: ekmtb,ekmtf,ekmtin
      real*8  cv1_3,chi,fv1SA,sigma_om2,betaStar,gradkom,gamma1,gamma2,gamma3,gammaSST,zetaSST,StR, wallD



      sigmat = 0.9


      do k=1,kmax
         km=k-1
         kp=k+1

         tauw(k) = ekmi(imax,k)*0.5*(W(imax,km)+W(imax,k))/wallDist(imax)

         do i=1,imax
            im=i-1
            ip=i+1

            yp(i,k) = sqrt(rNew(i,k))/ekm(i,k)*(wallDist(i))*tauw(k)**0.5           ! ystar
  !          yp(i,k) = sqrt(rNew(imax,k))/ekm(imax,k)*(0.5-Rp(i))*tauw(k)**0.5    ! yplus
  !          yp(i,:)=(0.5-Rp(i))*Re*(1/Re*(Win(imax)/(0.5-Rp(imax))))**0.5        ! yplus
            ReTauS(i,k) = 0.5*sqrt(rNew(i,k))/ekm(i,k)*tauw(k)**0.5
            Ret(i,k)=rNew(i,k)*(kNew(i,k)**2.)/(ekm(i,k)*eNew(i,k))        ! not sure if r2 or r

            if (turbmod.eq.0) then

               ekmttmp(i,k) = 0.

            elseif (turbmod.eq.1) then

               sigmak       = 1.4
               sigmae       = 1.3
               cmu          = 0.09
               ce1          = 1.4
               ce2          = 1.8
               fmu(i,k)     = (1.-exp(-yp(i,k)/70.))*(1+3.45/Ret(i,k)**0.5)
               f1(i,k)      = 1.
               f2(i,k)      = (1.-2./9.*exp(-(Ret(i,k)/6.)**2.))*(1.-exp(-yp(i,k)/5.))**2.0
               dterm(i,k)   = 0.
               eterm(i,k)   = 0.
               ekmttmp(i,k) = min(1.,rNew(i,k)*cmu*fmu(i,k)*kNew(i,k)**2./(eNew(i,k)))

            elseif (turbmod.eq.2) then

               sigmak = 1.0
               sigmae = 1.3
               cmu    = 0.09
               ce1    = 1.44
               ce2    = 1.92
               fmu(i,k)=exp(-3.4/(1.+Ret(i,k)/50.)**2.)
               f1(i,k)=1.
               f2(i,k)=1.-0.3*exp(-(Ret(i,k))**2.)
               ekmttmp(i,k)=min(1.,rNew(i,k)*cmu*fmu(i,k)*kNew(i,k)**2./(eNew(i,k)))
               kd(ip,k)=kNew(ip,k)
               kd(im,k)=kNew(im,k)
               kd(i,kp)=kNew(i,kp)
               kd(i,km)=kNew(i,km)
               if (kd(ip,k).lt.0.0) kd(ip,k)=0.0
               if (kd(im,k).lt.0.0) kd(im,k)=0.0
               if (kd(i,kp).lt.0.0) kd(i,kp)=0.0
               if (kd(i,km).lt.0.0) kd(i,km)=0.0
!
               if (i.eq.imax) then
                  dterm(i,k)=2.*ekm(i,k)/rNew(i,k)*((kd(im,k)**0.5)/dru(im)/2.0+
     &                 (kd(i,kp)**0.5-kd(i,km)**0.5)/(2.*dz))**2.
               else
                  dterm(i,k)=2.*ekm(i,k)/rNew(i,k)*((kd(im,k)**0.5-kd(ip,k)**0.5)/(dRp(i)+dRp(im))+
     &                 (kd(i,kp)**0.5-kd(i,km)**0.5)/(2.*dz))**2.
               endif
                   eterm(i,k)=2.*ekm(i,k)*ekmt(i,k)/rNew(i,k)**2.*(((((U(i,kp)+U(im,kp))/2-(U(i,k)+U(im,k))/2)
     &              /dz-((U(i,k)+U(im,k))/2-(U(i,km)+U(im,km))/2)/dz)/dz)**2
     &              +((((W(im,k)+W(im,km))/2-(W(i,k)+W(i,km))/2)/dRp(im)-((W(i,k)+W(i,km))/2-
     &              (W(ip,k)+W(ip,km))/2)/dRp(i))/(Ru(i)-Ru(im)))**2)

            elseif (turbmod.eq.3) then

               sigmak = 1.0
               sigmae = 1.3
               cmu    = 0.22
               ce1    = 1.4
               ce2    = 1.9
               
               StR= (2.*(((W(i,k)-W(i,km))/dz)**2. +
     &                ((U(i,k)-U(im,k))/(Ru(i)-Ru(im)))**2. +
     &                ((U(i,k)+U(im,k))/(2.*Rp(i)))**2.) +  ! CYLCOORD
     &                (((W(ip,km)+W(ip,k)+W(i,km)+W(i,k))/4.
     &                -(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/(Ru(i)-Ru(im))
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

            elseif (turbmod .eq. 4) then

               cv1_3        = (7.1)**3.0
               chi          = nuSAnew(i,k)/(ekm(i,k)/rNew(i,k))
               fv1SA        = (chi**3.0)/(chi**3.0 + cv1_3);
               ekmttmp(i,k) = min(100.0,max(0.0,rNew(i,k)*nuSAnew(i,k)*fv1SA))

            elseif (turbmod.eq.5) then

           
              !constants
               sigma_om2 = 0.856
               betaStar  = 0.09
               wallD     = wallDist(i)

              ! Vorticity rate
              StR = ( ( -( (W(ip,km)+W(ip,k)+W(i,km)+W(i ,k))/4.-(W(im,km)+W(im,k)+W(i,km)+W(i,k))/4.)/(Ru(i)-Ru(im))
     &                   +( (U(i,kp) +U(im,kp)+U(i,k)+U(im,k))/4.-(U(im,km)+U(i,km)+U(im,k)+U(i,k))/4.)/(dz)           )**2.)

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

            endif

         enddo
      enddo


      ! Boundary condition for bF1 of sst model
      if (turbmod.eq.5) then
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
      endif




      ekmttmp(i1,:) = -ekmttmp(imax,:)
      ekmttmp(0,:)  =  ekmttmp(1,:)

      call shiftf(ekmttmp,ekmtf,rank)
      call shiftb(ekmttmp,ekmtb,rank)
      ekmttmp(:,0)  = ekmtf(:)
      ekmttmp(:,k1) = ekmtb(:)

      if ((periodic.ne.1).and.(rank.eq.0)) then
         ekmttmp(:,0) = ekmtin(:)
      endif

      if ((periodic.ne.1).and.(rank.eq.px-1)) then
         ekmttmp(:,k1) = 2.*ekmt(:,kmax)-ekmt(:,kmax-1)
      endif

      ekmetmp = ekm + ekmttmp

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
      subroutine advanceScalar(resC,resK,resE,resV2,resOm,resSA,Utmp,Wtmp,Rtmp,ftmp,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      real*8 dnew(0:i1,0:k1),tempArray(0:i1,0:k1),dimpl(0:i1,0:k1),tscl
      real*8 Utmp(0:i1,0:k1),Wtmp(0:i1,0:k1),Rtmp(0:i1,0:k1),ftmp(imax,kmax),sigmakSST(0:i1,0:k1)
      real*8 rho2(0:i1,0:k1), rho3(0:i1,0:k1), eknu(0:i1,0:k1),eknui(0:i1,0:k1),eknuk(0:i1,0:k1)
      real*8 cb3,Q
      integer ierr
      real*8     a  (imax)
      real*8     b  (imax)
      real*8     c  (imax)
      real*8     rhs(imax)

      real*8 t1,t2,t3,t4,bc,scl
      real*8 term, adiffm, adiffp
      real*8 resC, resK, resE, resV2, resOm, resSA


! ------------------------------------------------------------------------
!     ENTHAPLY ENTHAPLY ENTHAPLY ENTHAPLY
! ------------------------------------------------------------------------
      resC  = 0.0
      dnew=0.0; dimpl = 0.0;
      call advecc(dnew,dimpl,cnew,Utmp,Wtmp,Ru,Rp,dru,dz,i1,k1,rank,periodic,.true.)
      call diffc(dnew,cnew,ekh,ekhi,ekhk,ekmt,sigmat,Rtmp,Ru,Rp,dru,dz,rank,0)

      do k=1,kmax
         if (rank.eq.0.and.k.lt.K_start_heat) then
            Q=0.0
         else
            Q=Qwall
         endif

         do i=1,imax-1
            a(i) = -Ru(i-1)*(ekhi(i-1,k)+0.5*(ekmt(i,k)+ekmt(i-1,k))/sigmat)/(dRp(i-1)*Rp(i)*dru(i))/Rtmp(i,k)
            c(i) = -Ru(i  )*(ekhi(i  ,k)+0.5*(ekmt(i,k)+ekmt(i+1,k))/sigmat)/(dRp(i  )*Rp(i)*dru(i))/Rtmp(i,k)
            b(i) = (-a(i)-c(i) + dimpl(i,k) )/alphac        ! BUG
            rhs(i) = dnew(i,k) + (1-alphac)*b(i)*cnew(i,k)  ! BUG
         enddo

         b(1)=b(1)+a(1)
         i=imax
            a(i)   = -Ru(i-1)*(ekhi(i-1,k)+0.5*(ekmt(i,k)+ekmt(i-1,k))/sigmat)/(dRp(i-1)*Rp(i)*dru(i))/Rtmp(i,k)
            c(i)   =  0.0
            b(i)   =  (-a(i)-c(i) + dimpl(i,k) )/alphac
            rhs(i) = dnew(i,k) + Ru(i)*Q/(Re*Pr*Rtmp(i,k)*Rp(i)*dru(i)) + (1-alphac)*b(i)*cnew(i,k)

         call matrixIdir(imax,a,b,c,rhs)

         do i=1,imax
            !resC = resC + (cnew(i,k) - rhs(i))**2.0
            resC = resC + ((cnew(i,k) - rhs(i))/(cnew(i,k)+1.0e-20))**2.0
            cnew(i,k) = max(rhs(i), 0.0)
         enddo
      enddo


      if (periodic.eq.1) cnew = 0.0

! ------------------------------------------------------------------------
!     TURBULENCE MODELING
! ------------------------------------------------------------------------
    
      resK  = 0.0
      resE  = 0.0
      resV2 = 0.0
      resOm = 0.0
      resSA = 0.0

      ! modified turb. model
      !    modifDiffTerm = 1, our modification
      !    modifDiffTerm = 2, Aupoix modification
      rho2  = Rtmp
      if ((modifDiffTerm == 1) .or. (modifDiffTerm == 2)) then
         rho3 = Rtmp
      else
         rho3 = 1.0
      endif

      if (turbmod.eq.0) return      !>     in case is laminar this fuction finishes here....


!     0..laminar, 1..MK, 2..LS, 3..V2F
      if (turbmod <= 3) then

! ------------------------------------------------------------------------
! epsilon epsilon epsilon epsilon epsilon epsilon epsilon epsilon
! ------------------------------------------------------------------------
          dnew=0.0; dimpl = 0.0;
          call advecc(dnew,dimpl,eNew,utmp,wtmp,Ru,Rp,dru,dz,i1,k1,rank,periodic,.true.)
          !if (modifDiffTerm == 1) call advecrho(dnew,eNew,utmp,wtmp,Ru,Rp,dru,dz,i1,k1,rank)
          scl=1.0
          call prodis(dnew,dimpl,kNew,eNew,v2New,nuSANew,ftmp,Utmp,Wtmp,temp,Rtmp,scl)
          call diffEPS(dnew,eNew,ekm,ekmi,ekmk,ekmt,sigmae,rho2,Ru,Rp,dru,dz,rank,modifDiffTerm)

          do k=1,kmax
             do i=1,imax
                
                a(i) = (ekmi(i-1,k)+0.5*(ekmt(i,k)+ekmt(i-1,k))/sigmae)/sqrt(0.5*(rho3(i-1,k)+rho3(i,k)))
                a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/rho2(i,k)/rho3(i,k)

                c(i) = (ekmi(i  ,k)+0.5*(ekmt(i,k)+ekmt(i+1,k))/sigmae)/sqrt(0.5*(rho3(i+1,k)+rho3(i,k)))
                c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/rho2(i,k)/rho3(i,k)

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


! ------------------------------------------------------------------------
! TKE  TKE  TKE  TKE  TKE  TKE  TKE  TKE  TKE  TKE  TKE  TKE  TKE  TKE
! ------------------------------------------------------------------------
          dnew=0.0; dimpl = 0.0;
          call advecc(dnew,dimpl,kNew,utmp,wtmp,Ru,Rp,dru,dz,i1,k1,rank,periodic,.true.)
          !if (modifDiffTerm == 1) call advecrho(dnew,kNew,utmp,wtmp,Ru,Rp,dru,dz,i1,k1,rank)
          scl=0.0
          call prodis(dnew,dimpl,kNew,eNew,v2new,nuSAnew,ftmp,Utmp,Wtmp,temp,Rtmp,scl)
          call diffc(dnew,kNew,ekm,ekmi,ekmk,ekmt,sigmak,rho2,Ru,Rp,dru,dz,rank,modifDiffTerm)

          do k=1,kmax
             do i=1,imax
                if ((modifDiffTerm == 0) .or. (modifDiffTerm == 1)) then
                   a(i) = (ekmi(i-1,k)+0.5*(ekmt(i,k)+ekmt(i-1,k))/sigmak)/((0.5*(rho3(i-1,k)+rho3(i,k)))**0.5)
                   a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/rho2(i,k)/(rho3(i,k)**0.5)
                   c(i) = (ekmi(i  ,k)+0.5*(ekmt(i,k)+ekmt(i+1,k))/sigmak)/((0.5*(rho3(i+1,k)+rho3(i,k)))**0.5)
                   c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/rho2(i,k)/(rho3(i,k)**0.5)
                else if (modifDiffTerm == 2) then
                   a(i) = (ekmi(i-1,k)+0.5*(ekmt(i,k)+ekmt(i-1,k))/sigmak)/(0.5*(rho3(i-1,k)+rho3(i,k)))
                   a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/rho2(i,k)
                   c(i) = (ekmi(i  ,k)+0.5*(ekmt(i,k)+ekmt(i+1,k))/sigmak)/(0.5*(rho3(i+1,k)+rho3(i,k)))
                   c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/rho2(i,k)
                endif

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





! ------------------------------------------------------------------------
! v2f  v2f  v2f  v2f  v2f  v2f  v2f  v2f  v2f  v2f  v2f  v2f  v2f
! ------------------------------------------------------------------------
          if (turbmod.eq.3) then
             dnew=0.0; dimpl = 0.0;
             call advecc(dnew,dimpl,v2New,utmp,wtmp,Ru,Rp,dru,dz,i1,k1,rank,periodic,.true.)
             !if (modifDiffTerm == 1) call advecrho(dnew,v2New,utmp,wtmp,Ru,Rp,dru,dz,i1,k1,rank)
             scl=2.0
             call prodis(dnew,dimpl,kNew,eNew,v2New,nuSANew,ftmp,Utmp,Wtmp,temp,Rtmp,scl)
             call diffc(dnew,v2New,ekm,ekmi,ekmk,ekmt,sigmak,rho2,Ru,Rp,dru,dz,rank,modifDiffTerm)

             do k=1,kmax
                do i=1,imax

                   if ((modifDiffTerm == 0) .or. (modifDiffTerm == 1)) then
                      a(i) = (ekmi(i-1,k)+0.5*(ekmt(i,k)+ekmt(i-1,k))/sigmak)/((0.5*(rho3(i-1,k)+rho3(i,k)))**0.5)
                      a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/rho2(i,k)/(rho3(i,k)**0.5)
                      c(i) = (ekmi(i  ,k)+0.5*(ekmt(i,k)+ekmt(i+1,k))/sigmak)/((0.5*(rho3(i+1,k)+rho3(i,k)))**0.5)
                      c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/rho2(i,k)/(rho3(i,k)**0.5)
                   else if (modifDiffTerm == 2) then
                      a(i) = (ekmi(i-1,k)+0.5*(ekmt(i,k)+ekmt(i-1,k))/sigmak)/(0.5*(rho3(i-1,k)+rho3(i,k)))
                      a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/rho2(i,k)
                      c(i) = (ekmi(i  ,k)+0.5*(ekmt(i,k)+ekmt(i+1,k))/sigmak)/(0.5*(rho3(i+1,k)+rho3(i,k)))
                      c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/rho2(i,k)
                   endif

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


! ------------------------------------------------------------------------
! Spalart Allmaras
! ------------------------------------------------------------------------
      elseif (turbmod.eq.4) then

          cb3 = 2.0/3.0

          dnew=0.0; dimpl = 0.0;
          call advecc(dnew,dimpl,nuSANew,utmp,wtmp,Ru,Rp,dru,dz,i1,k1,rank,periodic,.true.)
          !if (modifDiffTerm == 1) call advecrho(dnew,nuSANew,utmp,wtmp,Ru,Rp,dru,dz,i1,k1,rank)
          scl = 4.0
          call prodis(dnew,dimpl,kNew,eNew,v2New,nuSANew,ftmp,Utmp,Wtmp,temp,Rtmp,scl)

          do k=0,kmax+1
             do i=0,imax+1
                eknu (i,k) = ekm(i,k)/Rtmp(i,k)
             enddo
          enddo
          do k=0,kmax+1
             do i=0,imax
                eknui(i,k) = ekmi(i,k)/(0.5*(Rtmp(i,k)+Rtmp(i+1,k)))
             enddo
          enddo
          do k=0,kmax
             do i=0,imax+1
                eknuk(i,k) = ekmk(i,k)/(0.5*(Rtmp(i,k)+Rtmp(i,k+1)))
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

             b(1) = b(1)+a(1)
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


! ------------------------------------------------------------------------
! k-omega SST 
! ------------------------------------------------------------------------
      elseif (turbmod.eq.5) then

          ! --------------------------------------------------------------------
          ! ---------------------------- k equation ----------------------------
          ! --------------------------------------------------------------------
          dnew=0.0; dimpl = 0.0;
          call advecc(dnew,dimpl,kNew,utmp,wtmp,Ru,Rp,dru,dz,i1,k1,rank,periodic,.true.)
          !if (modifDiffTerm == 1) call advecrho(dnew,kNew,utmp,wtmp,Ru,Rp,dru,dz,i1,k1,rank)
          scl = 10.0
          call prodis(dnew,dimpl,kNew,eNew,v2New,nuSANew,ftmp,Utmp,Wtmp,temp,Rtmp,scl)

          ! calculating constant with blending function factor
          sigmakSST = 0.85*bF1 + 1.0*(1.0 - bF1)
          sigmakSST = 1.0/sigmakSST
          call diffcSSTKine(dnew,kNew,ekm,ekmi,ekmk,ekmt,sigmakSST,rho2,Ru,Rp,dru,dz,rank,modifDiffTerm)

          do k=1,kmax
             do i=1,imax

                if ((modifDiffTerm == 0) .or. (modifDiffTerm == 1)) then
                   a(i) = (ekmi(i-1,k)+(ekmt(i,k)+ekmt(i-1,k))/(sigmakSST(i,k)+sigmakSST(i-1,k)))/(0.5*(rho3(i-1,k)+rho3(i,k)))**0.5
                   a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/rho2(i,k)/rho3(i,k)**0.5
                   c(i) = (ekmi(i  ,k)+(ekmt(i,k)+ekmt(i+1,k))/(sigmakSST(i,k)+sigmakSST(i+1,k)))/(0.5*(rho3(i+1,k)+rho3(i,k)))**0.5
                   c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/rho2(i,k)/rho3(i,k)**0.5
                else if (modifDiffTerm == 2) then
                   a(i) = (ekmi(i-1,k)+(ekmt(i,k)+ekmt(i-1,k))/(sigmakSST(i,k)+sigmakSST(i-1,k)))/(0.5*(rho3(i-1,k)+rho3(i,k)))
                   a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/rho2(i,k)
                   c(i) = (ekmi(i  ,k)+(ekmt(i,k)+ekmt(i+1,k))/(sigmakSST(i,k)+sigmakSST(i+1,k)))/(0.5*(rho3(i+1,k)+rho3(i,k)))
                   c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/rho2(i,k)
                endif

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


          ! --------------------------------------------------------------------
          ! -------------------------- omega equation --------------------------
          ! --------------------------------------------------------------------
          dnew=0.0; dimpl = 0.0;
          call advecc(dnew,dimpl,omNew,utmp,wtmp,Ru,Rp,dru,dz,i1,k1,rank,periodic,.true.)
          !if (modifDiffTerm == 1) call advecrho(dnew,omNew,utmp,wtmp,Ru,Rp,dru,dz,i1,k1,rank)
          scl = 11.0
          call prodis(dnew,dimpl,omNew,eNew,v2New,nuSANew,ftmp,Utmp,Wtmp,temp,Rtmp,scl)

          ! calculating constant with blending function factor
          sigmakSST = 0.5*bF1 + 0.856*(1.0 - bF1)
          sigmakSST = 1.0/sigmakSST
          call diffcSSTOmega(dnew,omNew,ekm,ekmi,ekmk,ekmt,sigmakSST,rho2,Ru,Rp,dru,dz,rank,modifDiffTerm)

          do k=1,kmax
             do i=1,imax
                a(i) = (ekmi(i-1,k)+(ekmt(i,k)+ekmt(i-1,k))/(sigmakSST(i,k)+sigmakSST(i-1,k)))/(0.5*(rho3(i-1,k)+rho3(i,k)))**0.5 
                a(i) = -Ru(i-1)*a(i)/(dRp(i-1)*Rp(i)*dru(i))/rho2(i,k)/rho3(i,k)**0.5
                c(i) = (ekmi(i  ,k)+(ekmt(i,k)+ekmt(i+1,k))/(sigmakSST(i,k)+sigmakSST(i+1,k)))/(0.5*(rho3(i+1,k)+rho3(i,k)))**0.5 
                c(i) = -Ru(i  )*c(i)/(dRp(i  )*Rp(i)*dru(i))/rho2(i,k)/rho3(i,k)**0.5

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


      endif


      end










!>*************************************************************************************
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
!!*************************************************************************************
      subroutine advance(rank)
      implicit none
      include 'param.txt'
      include 'common.txt'

      real*8, dimension(imax)   :: a, b, c, rhs
      real*8, dimension(imax-1) :: au, bu, cu, rhsu
      real*8 dnew(0:i1,0:k1)
      real*8 dif,alpha,rhoa,rhob,rhoc

      dif=0.0

!>********************************************************************
!!     CALCULATE advection, diffusion and Force in r-direction
!!     at the old(n-1) and new(n) timelevels
!!********************************************************************
      dnew = 0.0
      call advecu(dnew,Unew,Wnew,Rnew,Ru,Rp,dru,dz,i1,k1)
      call diffu (dnew,Unew,Wnew,ekme,Ru,Rp,dru,dz,i1,k1,dif)

      do k=1,kmax
         do i=1,imax-1
            au(i) = -dt*ekme(i  ,k)*Rp(i  )/(dRp(i)*Ru(i)*dru(i  ))
            cu(i) = -dt*ekme(i+1,k)*Rp(i+1)/(dRp(i)*Ru(i)*dru(i+1))
            bu(i) = -au(i)-cu(i)
            rhoa = 0.5*(rnew(i  ,k)+rnew(i-1,k))
            rhoc = 0.5*(rnew(i+1,k)+rnew(i+2,k))
            rhob = 0.5*(rnew(i+1,k)+rnew(i  ,k))
            au(i) = au(i)/rhoa
            bu(i) = bu(i)/rhob + 1.0
            cu(i) = cu(i)/rhoc
         enddo

         i = imax-1;    cu(i) = 0.0
         i = 1;         bu(i) = bu(i) + au(i)

         do i=1,imax-1
            rhsu(i) = dt*dnew(i,k) + Unew(i,k)*(Rnew(i+1,k)+Rnew(i,k))*0.5
         enddo

         call matrixIdir(imax-1,au,bu,cu,rhsu)
         do i=1,imax-1
            dUdt(i,k)=rhsu(i)
         enddo
      enddo

!********************************************************************
!     CALCULATE advection, diffusion and Force in z-direction
!     at the old(n-1) and new(n) timelevels
!********************************************************************

      dnew = 0.0
      call advecw(dnew,Unew,Wnew,Rnew,Ru,Rp,dru,dz,ekm,peclet)
      call diffw (dnew,Unew,Wnew,ekme,Ru,Rp,dru,dz,i1,k1,dif,rank)

      if (periodic.eq.1) dnew = dnew + dpdz

      if (Qwall.ne.0) then
         dnew(1:imax,1:kmax) = dnew(1:imax,1:kmax) + 0.5*(Rnew(1:imax,1:kmax)+Rnew(1:imax,2:kmax+1))*Fr_1
      endif

      do k=1,kmax
         do i=1,imax
            a(i) = -0.25*dt*(ekme(i,k)+ekme(i,k+1)+ekme(i-1,k)+ekme(i-1,k+1))*Ru(i-1)/(dRp(i-1)*Rp(i)*dru(i))
            c(i) = -0.25*dt*(ekme(i,k)+ekme(i,k+1)+ekme(i+1,k)+ekme(i+1,k+1))*Ru(i  )/(dRp(i  )*Rp(i)*dru(i))
            b(i) = -a(i)-c(i)
            rhoa = 0.5*(rnew(i-1,k+1)+rnew(i-1,k))
            rhoc = 0.5*(rnew(i+1,k+1)+rnew(i+1,k))
            rhob = 0.5*(rnew(i  ,k+1)+rnew(i  ,k))
            a(i) = a(i)/rhoa
            b(i) = b(i)/rhob + 1.0
            c(i) = c(i)/rhoc
         enddo

         i=imax;    b(i) = b(i) - c(i)    ! BC at wall: zero vel: subtract one c
         i = 1;     b(i) = b(i) + a(i)    ! BC at center: Neumann: cancel coeff a

         do i=1,imax
            rhs(i) = dt*dnew(i,k) + Wnew(i,k)*(Rnew(i,k)+Rnew(i,k+1))*0.5
         enddo

         call matrixIdir(imax,a,b,c,rhs)
         do i=1,imax
            dWdt(i,k)=rhs(i)
         enddo
      enddo


      end






!<*************************************************************************************
!!
!!  bound_h equation
!!
!!*************************************************************************************
      subroutine bound_h(kin,ein,v2in,omin,nuSAin,rank)
      implicit none
c     
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'

      real*8 unin,flux,Ub,BCvalue(0:k1)
      real*8 Sk(0:k1),kin(0:i1),ein(0:i1),v2in(0:i1),omin(0:i1),nuSAin(0:i1)
      real*8 tmpShift(0:i1)


!     Radial boundary condition for enthalpy c
      do k=0,k1
         if (rank.eq.0.and.k.lt.K_start_heat) then
            cnew(i1,:) = cnew(imax,:)
         else
            call funcNewtonSolve(cnew(i1,k), cnew(imax,k))
         endif
      enddo


!     Radial boundary condition
      if (turbmod.eq.1) then

         knew(i1,:) = -knew(imax,:)
         BCvalue(:) = 2.0*ekm(imax,:)/rNew(imax,:)*knew(imax,:)/wallDist(imax)**2
         enew(i1,:) = 2.0*BCvalue(:) - enew(imax,:)

      elseif (turbmod.eq.2) then

         knew(i1,:) =  -knew(imax,:)
         enew(i1,:) =  -enew(imax,:)

      elseif (turbmod.eq.3) then

         knew(i1,:)  = -knew(imax,:)
         v2new(i1,:) = -v2new(imax,:)
         BCvalue(:)  = 2.*ekm(imax,:)/rNew(imax,:)*knew(imax,:)/wallDist(imax)**2
         enew(i1,:)  = 2.0*BCvalue(:) - enew(imax,:)

      elseif (turbmod.eq.4) then

         nuSAnew(i1,:) = -nuSAnew(imax,:)

      elseif (turbmod.eq.5) then

         knew(i1,:)  = -knew(imax,:)
         BCvalue(:)  = 60.0/0.075*ekm(imax,:)/rNew(imax,:)/wallDist(imax)**2
!         omNew(i1,:) = 50205072.858428769
         omNew(i1,:) = 2.0*BCvalue(:) - omNew(imax,:)

      endif


!     center line BC
      cnew(0,:)    = cnew(1,:)
      knew(0,:)    = knew(1,:)
      enew(0,:)    = enew(1,:)
      v2new(0,:)   = v2new(1,:)
      omNew(0,:)   = omNew(1,:)
      nuSAnew(0,:) = nuSAnew(1,:)

!     ghost cells between CPU
      call shiftf(cnew,   tmpShift,rank);       cnew(:,0) = tmpShift(:);
      call shiftf(knew,   tmpShift,rank);       knew(:,0) = tmpShift(:);
      call shiftf(enew,   tmpShift,rank);       enew(:,0) = tmpShift(:);
      call shiftf(v2new,  tmpShift,rank);      v2new(:,0) = tmpShift(:);
      call shiftf(omNew,  tmpShift,rank);      omNew(:,0) = tmpShift(:);
      call shiftf(nuSAnew,tmpShift,rank);    nuSAnew(:,0) = tmpShift(:);

      call shiftb(cnew,   tmpShift,rank);      cnew(:,k1) = tmpShift(:);
      call shiftb(knew,   tmpShift,rank);      knew(:,k1) = tmpShift(:);
      call shiftb(enew,   tmpShift,rank);      enew(:,k1) = tmpShift(:);
      call shiftb(v2new,  tmpShift,rank);     v2new(:,k1) = tmpShift(:);
      call shiftb(omNew,  tmpShift,rank);     omNew(:,k1) = tmpShift(:);
      call shiftb(nuSAnew,tmpShift,rank);   nuSAnew(:,k1) = tmpShift(:);
      
      if (periodic.eq.1) return


!     set inlet and outlet BC for developing flow
      if (rank.eq.0) then
          cnew(:,0) = 0.0
          knew(:,0) = kin(:)
          enew(:,0) = ein(:)
         v2new(:,0) = v2in(:)
         omNew(:,0) = omin(:)
       nuSAnew(:,0) = nuSAin(:)
      endif

      if (rank.eq.px-1) then
          cnew(:,k1) = 2.0*   cnew(:,kmax) -    cnew(:,kmax-1)
          knew(:,k1) = 2.0*   knew(:,kmax) -    knew(:,kmax-1)
          enew(:,k1) = 2.0*   enew(:,kmax) -    enew(:,kmax-1)
         v2new(:,k1) = 2.0*  v2new(:,kmax) -   v2new(:,kmax-1)
         omNew(:,k1) = 2.0*  omNew(:,kmax) -   omNew(:,kmax-1)
       nuSAnew(:,k1) = 2.0*nuSAnew(:,kmax) - nuSAnew(:,kmax-1)
      endif

      end



!>*************************************************************************************
!!      bound_v(Ubound,Wbound,Win,rank)
!!
!!
!!*************************************************************************************
      subroutine bound_v(Ubound,Wbound,Win,rank)

      implicit none

      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      character*5 inflow
      integer ierr,tabkhi,tabklo
      real*8  y1,y2,y3,y4
      real*8  Ubound(0:i1,0:k1), Wbound(0:i1,0:k1),flux,Ub,Win(0:i1)
      real*8 Rbb(0:i1)
      real*8 ubb(0:i1)
      real*8 wbb(0:i1)
      real*8 Rbf(0:i1)
      real*8 ubf(0:i1)
      real*8 wbf(0:i1)
      integer ib,ie,kb,ke



!     Radial Boundary condition
      do k=0,k1

         Ubound(0,k)    =   Ubound(1,k)

         Ubound(imax,k) =   0.0
         Ubound(i1,k)   = - Ubound(imax-1,k)

         Wbound(0,k)    =   Wbound(1,k)
         Wbound(i1,k)   = - Wbound(imax,k)

      enddo

      call shiftf(Ubound,ubf,rank);     Ubound(:,0)  = Ubf(:);
      call shiftf(Wbound,wbf,rank);     Wbound(:,0)  = Wbf(:);
      call shiftb(Ubound,ubb,rank);     Ubound(:,k1) = Ubb(:);
      call shiftb(Wbound,wbb,rank);     Wbound(:,k1) = Wbb(:);


!     if periodic is true, no need to overwrite k=0 proofile on rank=0 and to
!     apply advective outflow BC
      if (periodic.eq. 1) return

      if (rank.eq.0) then
         Ubound(:,0) = 0.0
         Wbound(:,0) = Win(:)
      endif

      if (rank.eq.px-1)then
         ubound(:,k1) = 2.*ubound(:,kmax)-ubound(:,kmax-1)
         wbound(:,k1) = 2.*wbound(:,kmax)-wbound(:,kmax-1)
      endif

      end




!>*************************************************************************************
!!
!!     Subroutine bound sets the boundary conditions for all variables,
!!     except for the diffusion coefficients. These are set in submod.
!!     The common boundary conditions for the pressure are set in mkgrid.
!!
!!     Set boundary conditions for j=0 and j=j1. Because of the
!!     fact that the index j denotes the tangential direction,
!!     we have to set the values at j=0 equal to the values at
!!     j=jmax and the values at j=j1 equal to the values at j=1.
!!
!!*************************************************************************************
      subroutine bound_m(Ubound,Wbound,W_out,Rbound,Win,rank)
      implicit none
c     
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      character*5 inflow
c     
      real*8      W_out(0:i1,0:k1)
      integer ierr
      real*8  y1,y2,y3,y4
      real*8  Ubound(0:i1,0:k1),Vbound(0:i1,0:k1)
     1     ,Wbound(0:i1,0:k1),Rbound(0:i1,0:k1),Win(0:i1)
      real*8 Ub,flux,flux_tot,deltaW,rhob,wfunc,wr(1:imax)
      real*8 Rbb(0:i1)
      real*8 ubb(0:i1)
      real*8 vbb(0:i1)
      real*8 wbb(0:i1)
      real*8 Rbf(0:i1)
      real*8 ubf(0:i1)
      real*8 vbf(0:i1)
      real*8 wbf(0:i1)
      integer   ib,ie,kb,ke


      do k=0,k1
         Ubound(0,k)    =   Ubound(1,k)

         Ubound(imax,k) =   0.0
         Ubound(i1,k)   = - Ubound(imax-1,k)

         Wbound(0,k)    =   Wbound(1,k)
         Wbound(i1,k)   = - Wbound(imax,k)
      enddo


      call shiftf(Ubound,ubf,rank)
      call shiftf(Wbound,wbf,rank)
      call shiftb(Ubound,ubb,rank)
      call shiftb(Wbound,wbb,rank)
!!!!!!!!!!!!!!!!!!!!!!!

      do i=0,i1
         Ubound(i,k1) = Ubb(i)
         Wbound(i,k1) = Wbb(i)
         Ubound(i,0)  = Ubf(i)
         Wbound(i,0)  = Wbf(i)
      enddo




!     if periodic is true, no need to overwrite k=0 proofile on rank=0 and to
!     apply advective outflow BC
      if (periodic.eq.1) return



      if (rank.eq.0) then
         Rbound(:,0) = 1.0
         Ubound(:,0) = 0.0
         Wbound(:,0) = Win(:)
      endif


      wr = 0
      Ub = 0.
      flux = 0.0

      if (rank.eq.px-1) then
         Ub = 0.
         do i=1,imax
            wr(i) = W_out(i,kmax)
            Ub = max(Ub,2.0*Wbound(i,kmax)/(Rbound(i,kmax)+Rbound(i,k1)))
         enddo

         do i=0,i1
            Wbound(i,kmax) = 2.0*W_out(i,kmax-1) - W_out(i,kmax-2)
            Wbound(i,kmax) = Wbound(i,kmax)*0.5*(Rbound(i,kmax)+Rbound(i,k1))
         enddo

         Wbound(i1,kmax) = -Wbound(imax,kmax)
         Wbound(0,kmax)  =  Wbound(1,kmax)
      endif

c     compute drho/dt*dvol
      do k=1,kmax
         do i=1,imax
            flux = flux - (rnew(i,k)-rold(i,k))/dt*Rp(i)*dru(i)*dz
         enddo
      enddo

c     compute mf in
      if (rank.eq.0)then
         do i=1,imax
            flux = flux + Wbound(i,0)*dru(i)*rp(i)
         enddo
      endif

      if (rank.eq.px-1)then
         Ub = 0
         wfunc = 0
         do i=1,imax
            flux = flux - Wbound(i,kmax)*dru(i)*rp(i)
            wfunc = wfunc + wr(i)*dru(i)*rp(i) ! based on averaged outflow velocity
         enddo
      endif

!      write(*,*) "delta flux: ", flux

c     
      call mpi_allreduce(flux,flux_tot,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)


      if (rank.eq.px-1)then
         deltaW = (flux_tot)/wfunc
!     write(*,*)'------------------------',deltaW
         do i=1,imax
            Wbound(i,kmax) = Wbound(i,kmax) + deltaW*wr(i) ! based on averaged outflow velocity
         enddo
         Wbound(i1,kmax) = -Wbound(imax,kmax)
         Wbound(0,kmax)  =  Wbound(1,kmax)

         flux = 0
         do i=1,imax
            flux = flux - Wbound(i,kmax)*dru(i)*rp(i)
         enddo
!         write(*,*) "flux out: ", flux

      endif




      return
      end






!>*************************************************************************************
!!
!!           fkdat(rank)
!!
!!*************************************************************************************
      subroutine fkdat(rank)
      implicit none
      include 'param.txt'
      include 'common.txt'

      real*8 yplus,t1,t2,t3,in,chl,ran,Wvel,delta

      character*5 inflow

      delta=0.5
      t1=3.5
      t2=5.
      t3=10.
      Unew =0.
      Uold =0.
      Cnew =0.


      rold =1.
      rnew =1.


        if (select_init.eq.1) then    
			!initialized from inflow
            if (rank.eq.0)  write(*,*) 'Initializing flow with inflow = ', select_init

            do k=0,k1
              if (turbmod.eq.0) open(29,file=  '0/Inflow',form='unformatted')
              if (turbmod.eq.1) open(29,file= 'MK/Inflow',form='unformatted')
              if (turbmod.eq.2) open(29,file= 'LS/Inflow',form='unformatted')
              if (turbmod.eq.3) open(29,file= 'VF/Inflow',form='unformatted')
              if (turbmod.eq.4) open(29,file= 'SA/Inflow',form='unformatted')
              if (turbmod.eq.5) open(29,file= 'OM/Inflow',form='unformatted')
              read(29) Wnew(:,k),knew(:,k),enew(:,k),v2new(:,k),omNew(:,k),nuSAnew(:,k),ekmt(:,k),Pk(:,k)
              close(29)
            enddo
        else
            !initialized from scratch values
            if (rank.eq.0)  write(*,*) 'Initializing flow from scratch = ', select_init

            do i=1,imax
              Wnew(i,:)  = 1.0 !Re/6*3/2.*(1-(rp(i)/0.5)**2)
              knew(i,:)  = 0.1
              enew(i,:)  = 1.0
              omnew(i,:) = 0.001
              v2new(i,:) = 2./3.*knew(i,:)
              nuSAnew(i,:) = 0.001
            enddo
        endif


      end
      


!>*************************************************************************************
!!           cmpinf(Bulk,Stress)
!!
!!*************************************************************************************
      subroutine cmpinf(Bulk,Stress)
      implicit none
     
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      integer ierr
      real*8 waver(imax),waver2(imax)
     
      real*8  Bulk,Stress
     
!     --- Initialization ---
   
      Waver = 0.0
     
!     -------------------------------------------start i,j,k-loop
      do  i=1,imax
         do k=1,kmax
            Waver(i) = Waver(i) + wnew(i,k)
         enddo
      enddo
      call mpi_allreduce(waver,waver2,imax,
     &     mpi_real8,mpi_sum,mpi_comm_world,ierr)
!     -------------------------------------------end i,j,k-loop
      waver = waver2/(kmax*px)

     
!     Stress =  Waver(imax) /(1.0-rp(imax))/Re

!     *** use MIDPOINT INTEGRATION RULE ***
     
      Bulk = 0.0
      do i=1,imax
         Bulk = Bulk + 8.*Waver(i) *Rp(i)* dru(i)
      enddo
      Stress =  Waver(imax)/wallDist(imax)/Re
      return
      end








