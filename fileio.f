
c********************************************************************
c     write inflow file
c********************************************************************
      subroutine Inflow_output(rank,istap) !tubstress,rank,istap)
      implicit none
      include 'param.txt'
      include 'common.txt'
      integer istap
      character*5 inflow

          if (rank.eq.px/2) then
             if (turbmod.eq.0) open(29,file='0/Inflow', form='unformatted')
             if (turbmod.eq.1) open(29,file='MK/Inflow',form='unformatted')
             if (turbmod.eq.2) open(29,file='LS/Inflow',form='unformatted')
             if (turbmod.eq.3) open(29,file='VF/Inflow',form='unformatted')
             if (turbmod.eq.4) open(29,file='SA/Inflow',form='unformatted')
             if (turbmod.eq.5) open(29,file='OM/Inflow',form='unformatted')
             !write(29) Wnew(:,kmax/2),knew(:,kmax/2),enew(:,kmax/2),v2new(:,kmax/2),omNEW(:,kmax/2),nuSAnew(:,kmax/2),ekmt(:,kmax/2),Pk(:,kmax/2)
             write(29) Wnew(:,kmax/2),knew(:,kmax/2),enew(:,kmax/2),v2new(:,kmax/2),omNEW(:,kmax/2),nuSAnew(:,kmax/2),ekmt(:,kmax/2),Pk(:,kmax/2),
     &                ekht(:,kmax/2),Pkt(:,kmax/2),ktnew(:,kmax/2),etnew(:,kmax/2)          ! modTemp
             close(29)
          endif

      end





c********************************************************************
c     read table from file
c********************************************************************
      subroutine readTable(rank)
      implicit none
      include 'param.txt'
      include 'mpif.h'
      if (EOSmode.eq.0) call readTableIG(rank)
      if (EOSmode.eq.1) call readTableRG(rank)
      end



c********************************************************************
c     read table from file
c********************************************************************
      subroutine readTableIG(rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      integer ierr

      if (rank.eq.0) then
!     open(27,file='co2h_table.dat')
         do i=1,nTab
            tempTab(i)   = ((1.0*i-1)/(nTab-1.0)-0.1)*3.0 + 1.0
            rhoTab(i)    = 1.0/tempTab(i)
            betaTab(i)    = 1.0/tempTab(i)
            muTab(i)     = 1.0
            lamTab(i)    = 1.0
            cpTab(i)     = 1.0
            enthTab(i)   = tempTab(i) - 1.0
            lamocpTab(i) = lamTab(i)/cpTab(i)
         enddo
         close(27)
      endif

      call MPI_BCAST(tempTab,   nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rhoTab,    nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(muTab,     nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lamTab,    nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cpTab,     nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(enthTab,   nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lamocpTab, nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)

      end

c********************************************************************
c     read table from file
c********************************************************************
      subroutine readTableRG(rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      integer ierr

      if (rank.eq.0) then
         open(27,file='co2h_table.dat')
         do i=1,nTab
            read (27,*) tempTab(i),rhoTab(i),muTab(i),lamTab(i),cpTab(i),enthTab(i),betaTab(i)

!     rhoTab(i)  = 1.0/tempTab(i)
!     muTab(i)   = 1.0
!     lamTab(i)  = 1.0
!     cpTab(i)   = 1.0
!     enthTab(i) = tempTab(i) - 1.0

            lamocpTab(i) = lamTab(i)/cpTab(i)
         enddo
         close(27)
      endif

      call MPI_BCAST(tempTab,   nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rhoTab,    nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(muTab,     nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lamTab,    nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cpTab,     nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(enthTab,   nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lamocpTab, nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(betaTab,   nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
      end

c********************************************************************
c     dump table to file
c********************************************************************
      subroutine dumpTABLE(rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      character*5 cha
      write(cha,'(I5.5)')rank
      open(207,file='table.'//cha)
      do i=1,nTab
         write (207,'(7E20.12)') tempTab(i),rhoTab(i),muTab(i),lamTab(i),cpTab(i),enthTab(i),lamocpTab(i)
      enddo
      close(207)
      end



c***************************************************************************************
c     read istep from the restart file to identify the corresponding inflow profile
c***************************************************************************************
      subroutine loadRestart(startStep,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      character*5 cha
      integer startStep,ierr

      real*8     UNEWR(0:i1,0:k1old)
      real*8     UOLDR(0:i1,0:k1old)
      real*8     WNEWR(0:i1,0:k1old)
      real*8     WOLDR(0:i1,0:k1old)
      real*8     CNEWR(0:i1,0:k1old)
      real*8     KNEWR(0:i1,0:k1old)
      real*8     ENEWR(0:i1,0:k1old)
      real*8    V2NEWR(0:i1,0:k1old)
      real*8  NUSANEWR(0:i1,0:k1old)
      real*8    omNEWR(0:i1,0:k1old)
	  real*8       PKR(0:i1,0:k1old)
	  real*8    KTNEWR(0:i1,0:k1old)   ! modTemp
      real*8    ETNEWR(0:i1,0:k1old)   ! modTemp
	  real*8      PKTR(0:i1,0:k1old)   ! modTemp

      write(cha,'(I5.5)')rank
      if (turbmod.eq.0) open(19,file= '0/Restart/start_stop.'//cha,form='unformatted')
      if (turbmod.eq.1) open(19,file='MK/Restart/start_stop.'//cha,form='unformatted')
      if (turbmod.eq.2) open(19,file='LS/Restart/start_stop.'//cha,form='unformatted')
      if (turbmod.eq.3) open(19,file='VF/Restart/start_stop.'//cha,form='unformatted')
      if (turbmod.eq.4) open(19,file='SA/Restart/start_stop.'//cha,form='unformatted')
      if (turbmod.eq.5) open(19,file='OM/Restart/start_stop.'//cha,form='unformatted')
         read(19) startStep
         read(19) UNEWR,UOLDR,WNEWR,WOLDR,CNEWR,KNEWR,ENEWR,V2NEWR,omNEWR,nuSANEWR,PKR,KTNEWR,ETNEWR,PKTR   ! modTemp
         close(19)

      do i=0,i1
        do k=0,k1
            UNEW(i,k) =    UNEWR(i,k1old*k/k1)
            UOLD(i,k) =    UOLDR(i,k1old*k/k1)
            WNEW(i,k) =    WNEWR(i,k1old*k/k1)
            WOLD(i,k) =    WOLDR(i,k1old*k/k1)
            CNEW(i,k) =    CNEWR(i,k1old*k/k1)
            KNEW(i,k) =    KNEWR(i,k1old*k/k1)
            ENEW(i,k) =    ENEWR(i,k1old*k/k1)
           V2NEW(i,k) =   V2NEWR(i,k1old*k/k1)
         nuSAnew(i,k) = nuSAnewR(i,k1old*k/k1)
           omNew(i,k) =   omNewR(i,k1old*k/k1)
		      Pk(i,k) =      PKR(i,k1old*k/k1)
		   ktnew(i,k) =   KTNEWR(i,k1old*k/k1)
		   etnew(i,k) =   ETNEWR(i,k1old*k/k1)
		     Pkt(i,k) =     PKTR(i,k1old*k/k1)
        enddo
      enddo


      end

c***************************************************************************************
c     store istep in the restart file to identify the corresponding inflow profile
c***************************************************************************************
      subroutine saveRestart(rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      character*5 cha
      write(cha,'(I5.5)')rank
      if (turbmod.eq.0) open(19,file= '0/Restart/start_stop.'//cha,form='unformatted')
      if (turbmod.eq.1) open(19,file='MK/Restart/start_stop.'//cha,form='unformatted')
      if (turbmod.eq.2) open(19,file='LS/Restart/start_stop.'//cha,form='unformatted')
      if (turbmod.eq.3) open(19,file='VF/Restart/start_stop.'//cha,form='unformatted')
      if (turbmod.eq.4) open(19,file='SA/Restart/start_stop.'//cha,form='unformatted')
      if (turbmod.eq.5) open(19,file='OM/Restart/start_stop.'//cha,form='unformatted')
         write(19) istep  
         write(19) UNEW,UOLD,WNEW,WOLD,CNEW,KNEW,ENEW,V2NEW,omNEW,nuSANEW,Pk,ktNEW,etNEW,Pkt       ! modTemp
         close(19)
      end




c***************************************************************************************
c
c***************************************************************************************
      subroutine outputProfile(rank)
      implicit none
      include 'param.txt'
      include 'common.txt'

    ! radius= 1, velocity=3, temperature=5, density=6, turb. kine=7,
    ! epsilon=8, v2=9, nuSA=10, mut=11
      if (rank == 0) then
        open(19,file='profile')
          k = 1
          do i=1,imax
            write(19,'(15E18.6)') rp(i),unew(i,k),Wnew(i,k),cnew(i,k),
     &             temp(i,k),rnew(i,k),knew(i,k),enew(i,k),v2new(i,k),nuSAnew(i,k),ekmt(i,k),
     &             ekht(i,k),Pkt(i,k),ktnew(i,k),etnew(i,k)          ! modTemp
          enddo
        close(19)
      endif

      end


c***************************************************************************************
c     
c***************************************************************************************
      subroutine outputX_h(rank,istap)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      character*5 cha
      integer ierr,istap,ktabhi,ktablo
      real*8 massflow(kmax),enthflow(kmax),enth_b(kmax),Twall(kmax),Tbulk(kmax),Qp(kmax),
     &     tmp(kmax),massfl,enth,laminter,tempinter,cpinter,Gb(kmax),rhob(kmax),ekmb(kmax),Ub(kmax),
     &     addedHeat,addedHeatTot,subheat,subheatTot,w_c

      twall    = 0.0
      Qp       = 0.0
      massflow = 0.0
      enthflow = 0.0
      w_c=0.
      Gb=0.


      do k=1,kmax

         do i=1,imax
            massfl = 0.5*rnew(i,k)*(Wnew(i,k)+Wnew(i,k-1))*rp(i)*dr(i)
!     Ub(k) = Ub(k)+2.*Wnew(i,k)*Rp(i)*dr(i)

            massflow(k) = massflow(k) + massfl

            enthflow(k) = enthflow(k) + massfl*Cnew(i,k)

         enddo
      enddo
!     Ub=Ub/(ru(imax)**2.)
      enth_b=enthflow/massflow
      do k=1,kmax
         ktabhi = 0
         ktablo = 0
         w_c=(Cnew(i1,k)+Cnew(imax,k))/2.
         call splint(enthTab,rhoTab,rho2Tab, nTab,enth_b(k),rhob(k),ktabhi,ktablo)
         call splint(enthTab,muTab,mu2Tab, nTab,enth_b(k),ekmb(k),ktabhi,ktablo)
         ekmb(k)=ekmb(k)/Re
         call splint(enthTab,tempTab,temp2Tab, nTab,w_c,Twall(k),ktabhi,ktablo) 
         call splint(enthTab,tempTab,temp2Tab, nTab,enth_b(k),Tbulk(k),ktabhi,ktablo)
         Gb(k) =massflow(k)/(4.*atan(1.)*Ru(imax)**2.)
      enddo

      ! bulk stuff
      write(cha,'(I5.5)')rank
      if (turbmod.eq.0) open(9,file='0/profX.'//cha)
      if (turbmod.eq.1) open(9,file='MK/profX.'//cha)
      if (turbmod.eq.2) open(9,file='LS/profX.'//cha)
      if (turbmod.eq.3) open(9,file='VF/profX.'//cha)
      if (turbmod.eq.4) open(9,file='SA/profX.'//cha)
      if (turbmod.eq.5) open(9,file='OM/profX.'//cha)
      do k=1,kmax
         write(9,'(8E18.6)') (k+rank*kmax)*dz, massflow(k),enthflow(k),enth_b(k)
     &        ,Twall(k),Tbulk(k),Gb(k)/Gb(1),Gb(k)/Gb(1)/rhob(k)

      enddo
      close(9)
      return
      end





c***************************************************************************************
c
c***************************************************************************************
      subroutine output2d(rank,istap)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      character*5 cha
      real*8 pecletx,peclety,pecletz
      real*8     THF(0:i1,0:k1),ReSS(0:i1,0:k1)
      real*8 massflow(kmax),enthflow(kmax),enth_b(kmax),Twall(kmax),Tbulk(kmax),
     &     massfl,w_c,ndt
      integer istap,jstart,ktabhi,ktablo
      write(cha,'(I5.5)')rank


      twall    = 0.0
      massflow = 0.0
      enthflow = 0.0
      w_c=0.

      do k=1,kmax

         do i=1,imax
            massfl = 0.5*rnew(i,k)*(Wnew(i,k)+Wnew(i,k-1))*rp(i)*dr(i)

            massflow(k) = massflow(k) + massfl

            enthflow(k) = enthflow(k) + massfl*Cnew(i,k)
            
            

         enddo
      enddo

      enth_b=enthflow/massflow
      do k=1,kmax
         ktabhi = 0
         ktablo = 0
         w_c=(Cnew(i1,k)+Cnew(imax,k))/2.
         call splint(enthTab,tempTab,temp2Tab, nTab,w_c,Twall(k),ktabhi,ktablo)
         call splint(enthTab,tempTab,temp2Tab, nTab,enth_b(k),Tbulk(k),ktabhi,ktablo)
      enddo

      if (turbmod.eq.0) open(15,file='0/tecp.'//cha)
      if (turbmod.eq.1) open(15,file='MK/tecp.'//cha)
      if (turbmod.eq.2) open(15,file='LS/tecp.'//cha)
      if (turbmod.eq.3) open(15,file='VF/tecp.'//cha)
      if (turbmod.eq.4) open(15,file='SA/tecp.'//cha)
      if (turbmod.eq.5) open(15,file='OM/tecp.'//cha)

      if (rank.eq.0) then
         write(15,*) 'VARIABLES ="X","Y","U","W","C","T","k","eps","v2","omega","nuSA","yplus","RHO","Pe","mu",
     &               "mut","lamcp","cp","alphat","kt","epst","Pk","Gk" '        ! modTemp
         write(15,*) 'ZONE I=  ', imax+2,' J=  ',(kmax+2)*px,' F=POINT '
      endif

      do k=0,k1
         do i=0,i1
            write(15,'(23ES24.10E3)')  (k+rank*kmax)*dz, rp(i),unew(i,k), Wnew(i,k), cnew(i,k), temp(i,k),                                              ! modTemp
     &           knew(i,k),enew(i,k),v2new(i,k),omNew(i,k),nuSAnew(i,k),yp(i,k),rnew(i,k),peclet(i,k),ekm(i,k),ekmt(i,k),ekh(i,k),cp(i,k),              ! modTemp
     &           ekht(i,k),ktnew(i,k),etnew(i,k),Pk(i,k),Gk(i,k)                                                                                                            ! modTemp
         enddo
      enddo

      close(15)


!      if (turbmod.eq.0) open(16,file='0/k_eps.'//cha)
!      if (turbmod.eq.1) open(16,file='MK/k_eps.'//cha)
!      if (turbmod.eq.2) open(16,file='LS/k_eps.'//cha)
!      if (turbmod.eq.3) open(16,file='VF/k_eps.'//cha)
!      do i=1,imax
!         k=kmax/2
!         if (rank.eq.3) k=2
!         write(16,'(13E18.5)') rp(i),yp(i,k),ypt(i,k),Wnew(i,k),knew(i,k),enew(i,k),v2new(i,k),fv2(i,k),
!     &        ekmt(i,k),fmu(i,k),rnew(i,k),Pk(i,k)*rnew(i,k),Gk(i,k)*rnew(i,k)
!
!      enddo
!      close(16)
!
!      if (turbmod.eq.0) open(17,file='0/temp.'//cha)
!      if (turbmod.eq.1) open(17,file='MK/temp.'//cha)
!      if (turbmod.eq.2) open(17,file='LS/temp.'//cha)
!      if (turbmod.eq.3) open(17,file='VF/temp.'//cha)
!      do i=1,imax
!         k=kmax/2
!         if (rank.eq.3) k=2
!         ndt=(Twall(k)-temp(i,k))/(Twall(k)-Tbulk(k))
!         write(17,'(7E18.5)') rp(i),yp(i,k),ypt(i,k),cnew(i,k),temp(i,k),ndt,h2new(i,k)
!
!      enddo
!      close(17)

      end
