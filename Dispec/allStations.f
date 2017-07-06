      subroutine allstations(fsname,recnames,netnames,lat,long,nseimx,n)
c
c Reads a station list in IRIS-Dmc format (download from
c `dmc.iris.washington.edu') and gives back all stations in the file.
c
C lat  = station latitude (degree, -90 -> +90)
C long = station longitude (degree, -180 -> +180)
      integer iu,maxs,nseimx,n,i,clength,ifl,iflp,ib,iborder(6,2)
      parameter(maxs=6000,clength=132)
      character recnames(*)*5,netnames(*)*2
      character totline*(clength)
      character*2 cyear,height*6
      character fsname*(*)
      real lat(maxs),long(maxs)


      iu = 1
      open(iu,file=fsname)
      read(iu,'(1x/1x/1x)')             ! Skip header

c-----------------------------------------------------------------
c Read line by line and parse them.                               |
c-----------------------------------------------------------------

      n = 1
 10   read(iu,'(a)',end=666) totline
        iflp=0
        ib = 0
        do i=1,clength
           if (totline(i:i).eq.' ') then
              ifl = 0
           else
              ifl = 1
           endif
           if (iflp.eq.0.and.ifl.eq.1) then
              ib = ib + 1
              iborder(ib,1)=i
           else if (iflp.eq.1.and.ifl.eq.0) then
              iborder(ib,2)=i-1
           endif
           if (ib.eq.6) goto 100
           iflp=ifl
        enddo
 100    cyear=totline(iborder(1,1):iborder(1,2))
        recnames(n)=totline(iborder(2,1):iborder(2,2))
        netnames(n)=totline(iborder(3,1):iborder(3,2))
        read(totline(iborder(4,1):iborder(4,2)),'(f10.6)') lat(n)
        read(totline(iborder(5,1):iborder(5,2)),'(f11.6)') long(n)
        height=totline(iborder(6,1):iborder(6,2))
        print '(a5,1x,a2,f10.6,f11.6,a61)',
     &           recnames(n),netnames(n),lat(n),long(n)
       n = n+1
       if (n.gt.nseimx) stop 'Array dimension <nseimx> too small!'
      goto 10


 666  close(iu)
      n = n-1
      print *,'Got ',n,' stations.'
      return


      end
