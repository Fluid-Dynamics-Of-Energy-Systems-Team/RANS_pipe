
      implicit none
      include      'param.txt'
      include      'common.txt'

      integer pxNew, startStep
      real*8, dimension(0:i1,0:kmax*px+1) :: UNEW2,WNEW2,CNEW2,KNEW2,ENEW2,V2NEW2,Pk2,EKMT2,UOLD2,WOLD2,COLD2,KOLD2,EOLD2,V2old2
      character*5 cha


      do i=0,px-1
         call loadRestart(startStep,i)
          UNEW2(0:i1,kmax*i:kmax*(i+1)+1) =  UNEW(0:i1,0:k1)
          WNEW2(0:i1,kmax*i:kmax*(i+1)+1) =  WNEW(0:i1,0:k1)
          CNEW2(0:i1,kmax*i:kmax*(i+1)+1) =  CNEW(0:i1,0:k1)
          KNEW2(0:i1,kmax*i:kmax*(i+1)+1) =  KNEW(0:i1,0:k1)
          ENEW2(0:i1,kmax*i:kmax*(i+1)+1) =  ENEW(0:i1,0:k1)
         V2NEW2(0:i1,kmax*i:kmax*(i+1)+1) = V2NEW(0:i1,0:k1)
            Pk2(0:i1,kmax*i:kmax*(i+1)+1) =    Pk(0:i1,0:k1)
          EKMT2(0:i1,kmax*i:kmax*(i+1)+1) =  EKMT(0:i1,0:k1)
          UOLD2(0:i1,kmax*i:kmax*(i+1)+1) =  UOLD(0:i1,0:k1)
          WOLD2(0:i1,kmax*i:kmax*(i+1)+1) =  WOLD(0:i1,0:k1)
          COLD2(0:i1,kmax*i:kmax*(i+1)+1) =  COLD(0:i1,0:k1)
          KOLD2(0:i1,kmax*i:kmax*(i+1)+1) =  KOLD(0:i1,0:k1)
          EOLD2(0:i1,kmax*i:kmax*(i+1)+1) =  EOLD(0:i1,0:k1)
         V2old2(0:i1,kmax*i:kmax*(i+1)+1) = V2old(0:i1,0:k1)
      enddo

      rank = 0
      write(cha,'(I5.5)')rank
      if (turbmod.eq.0) open(19,file= '0/Restart/start_stop.1core.'//cha,form='unformatted')
      if (turbmod.eq.1) open(19,file='SA/Restart/start_stop.1core.'//cha,form='unformatted')
      if (turbmod.eq.2) open(19,file='MK/Restart/start_stop.1core.'//cha,form='unformatted')
      if (turbmod.eq.3) open(19,file='VF/Restart/start_stop.1core.'//cha,form='unformatted')
      if (turbmod.eq.4) open(19,file='OM/Restart/start_stop.1core.'//cha,form='unformatted')

      write(19) time,istep
      write(19) UNEW2,WNEW2,CNEW2,KNEW2,ENEW2,V2NEW2,Pk2,EKMT2,UOLD2,WOLD2,COLD2,KOLD2,EOLD2,V2old2
      close(19)


      end





c***************************************************************************************
c     read istep from the restart file to identify the corresponding inflow profile
c***************************************************************************************
      subroutine loadRestart(startStep,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      character*5 cha
      integer startStep
      write(cha,'(I5.5)')rank
      if (turbmod.eq.0) open(19,file='0/Restart/start_stop.'//cha,form='unformatted')
      if (turbmod.eq.1) open(19,file='SA/Restart/start_stop.'//cha,form='unformatted')
      if (turbmod.eq.2) open(19,file='MK/Restart/start_stop.'//cha,form='unformatted')
      if (turbmod.eq.3) open(19,file='VF/Restart/start_stop.'//cha,form='unformatted')
      if (turbmod.eq.4) open(19,file='OM/Restart/start_stop.'//cha,form='unformatted')
      read(19) time,startStep
      read(19) UNEW,WNEW,CNEW,KNEW,ENEW,V2NEW,Pk,EKMT,UOLD,WOLD,COLD,KOLD,EOLD,V2old
      close(19)
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
      if (turbmod.eq.0) open(19,file='0/Restart/start_stop.'//cha,form='unformatted')
      if (turbmod.eq.1) open(19,file='SA/Restart/start_stop.'//cha,form='unformatted')
      if (turbmod.eq.2) open(19,file='MK/Restart/start_stop.'//cha,form='unformatted')
      if (turbmod.eq.3) open(19,file='VF/Restart/start_stop.'//cha,form='unformatted')
      if (turbmod.eq.4) open(19,file='OM/Restart/start_stop.'//cha,form='unformatted')

      write(19) time,istep
      write(19) UNEW,WNEW,CNEW,KNEW,ENEW,V2NEW,Pk,EKMT,UOLD,WOLD,COLD,KOLD,EOLD,V2old
      close(19)
      end






