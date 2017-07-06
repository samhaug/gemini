      subroutine writesff(nfmax,x,nseis,RecLat,RecLong,RECname,nsamp_out,
     &                   samprat,tred)
c
c Writes out the seismograms in Sff-format
c
c
      integer nfmax,nfdumm,i
      parameter(nfdumm=8192)
      real x(0:nfmax,3,*),RecLat(*),RecLong(*),fdata(nfdumm),SourceLat,
     &     SourceLong,SourceDepth,samprat,tred,seconds,shisec
      integer ns,nseis,ierr,idata(nfdumm),nsamp_out,minutes,iy,imo,iday,iho,
     &        shimin
      equivalence(fdata,idata)
      character date*6,time*10,SourceType*20,RECname(*)*5
      character*132 wid2line
      common /outpa/ SourceType,date,time,SourceLat,SourceLong,SourceDepth

c Convert seconds of tred into minutes and seconds for wid2line

      shimin = int(tred/60.)
      shisec = mod(tred,60.)
      print *, 'Timeshift in minutes and seconds: ',shimin,shisec
      read(date,'(3i2)') iy,imo,iday
      if(iy.lt.50) then
        iy = 2000 + iy
      else
        iy = 1900 + iy
      endif
      read(time,'(2i2,f6.3)') iho,minutes,seconds

c Add time shift

      seconds = seconds + shisec
      if (seconds.ge.60.) then
        seconds = seconds - 60.
        shimin = shimin + 1
      endif
      minutes = minutes + shimin
      if (minutes.ge.60) then
        minutes = minutes - 60
        iho = iho + 1
        if (iho.ge.24) then
          iho = iho - 24
          iday = iday +1   ! Hopefully we are not at the end of the month !!!
        endif
      endif


c           Output of seismograms
c           ~~~~~~~~~~~~~~~~~~~~~

        call sff_New(10,'plot_z.sff',ierr)
        call sff_WOpenS(10,'plot_z.sff',SourceType,'S',SourceLat,SourceLong,
     &            SourceDepth,date,time,ierr)
        do ns=1,nseis
        call sff_PrepWid2(nsamp_out,samprat,RECname(ns),iy,imo,iday,iho
     &       ,minutes,'Z  ','NSP','NSP',seconds,-1.,-1.,-1.,-1.,wid2line,ierr)
            do i=1,nsamp_out
             fdata(i) = x(i-1,1,ns)
            enddo
            if (ns.lt.nseis) then
               call sff_WTraceI(10,wid2line,nsamp_out,fdata,idata,.false.,'S',
     &        RecLat(ns),RecLong(ns),0.,1,ierr)
            else if(ns.eq.nseis) then
               call sff_WTraceI(10,wid2line,nsamp_out,fdata,idata,.true.,'S',
     &        RecLat(ns),RecLong(ns),0.,1,ierr)
            endif
        enddo

        call sff_New(10,'plot_ns.sff',ierr)
        call sff_WOpenS(10,'plot_ns.sff',SourceType,'S',SourceLat,SourceLong,
     &            SourceDepth,date,time,ierr)
        do ns=1,nseis
        call sff_PrepWid2(nsamp_out,samprat,RECname(ns),iy,imo,iday,iho
     &       ,minutes,'NS ','NSP','NSP',seconds,-1.,-1.,-1.,-1.,wid2line,ierr)
            do i=1,nsamp_out
             fdata(i) = x(i-1,2,ns)
            enddo
            if (ns.lt.nseis) then
               call sff_WTraceI(10,wid2line,nsamp_out,fdata,idata,.false.,'S',
     &        RecLat(ns),RecLong(ns),0.,1,ierr)
            else if(ns.eq.nseis) then
               call sff_WTraceI(10,wid2line,nsamp_out,fdata,idata,.true.,'S',
     &        RecLat(ns),RecLong(ns),0.,1,ierr)
            endif
        enddo

        call sff_New(10,'plot_ew.sff',ierr)
        call sff_WOpenS(10,'plot_ew.sff',SourceType,'S',SourceLat,SourceLong,
     &            SourceDepth,date,time,ierr)
        do ns=1,nseis
        call sff_PrepWid2(nsamp_out,samprat,RECname(ns),iy,imo,iday,iho
     &       ,minutes,'EW ','NSP','NSP',seconds,-1.,-1.,-1.,-1.,wid2line,ierr)
            do i=1,nsamp_out
             fdata(i) = x(i-1,3,ns)
            enddo
            if (ns.lt.nseis) then
               call sff_WTraceI(10,wid2line,nsamp_out,fdata,idata,.false.,'S',
     &        RecLat(ns),RecLong(ns),0.,1,ierr)
            else if(ns.eq.nseis) then
               call sff_WTraceI(10,wid2line,nsamp_out,fdata,idata,.true.,'S',
     &        RecLat(ns),RecLong(ns),0.,1,ierr)
            endif
        enddo

      print *,'Wrote ',nsamp_out,' samples to files "plot_z.sff", "plot_ns.sff" and "plot_ew.sff".'

      end
