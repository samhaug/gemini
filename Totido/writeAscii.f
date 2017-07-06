      subroutine writeAscii(RecLat,RecLong,nfmax,x,nsamp_out,samprat,tred,nseis)
c
c Writes out the seismograms in column-ASCII-format
c
c
      integer nfmax,i,komp,nseis,ns
      real x(0:nfmax,3,*),SourceLat,SourceLong,SourceDepth,samprat,tred,
     &     RecLat(*),RecLong(*)
      integer nsamp_out
      character date*6,time*10,SourceType*20
      common /outpa/ SourceType,date,time,SourceLat,SourceLong,SourceDepth


c           Output of seismograms
c           ~~~~~~~~~~~~~~~~~~~~~

      open(20,file='seismogr')
      write(20,999) 'SRCE ',SourceType,'S',SourceLat,SourceLong,
     &              SourceDepth,date,time
      write(20,998) 'Samples:',nsamp_out,'Samprat(Hz):',samprat,'Timeshift(sec):',tred
      do ns = 1, nseis
        write(20,'(2(a7,1x,f11.6,1x))')
     &         'RECLat ',RecLat(ns),'RECLong',RecLong(ns)
        write(20,'(4a13)') 'time(s)','Z','NS','EW'
        write(20,'(4e15.6)')
     &      (tred+i/samprat,(x(i,komp,ns),komp=1,3),i=0,nsamp_out-1)
      enddo
      close(20)
      print *,'Wrote ',nsamp_out,' samples to file "seismogr"'
c----------

 998  format(a8,1x,i5,1x,a12,1x,f10.6,1x,a15,1x,f13.6)
 999  format(a5,a20,1x,a1,1x,3f15.6,1x,a6,1x,a10)
      end
