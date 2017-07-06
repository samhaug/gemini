      subroutine getstation(station_file,nsro,cnetwork,lat,long)
c
c Reads a station list in IRIS-Dmc format (download from
c `dmc.iris.washington.edu') and searches for station 'nsro'.
c
C lat  = station latitude (degree, -90 -> +90)
C long = station longitude (degree, -180 -> +180)
      integer iu,ib,iborder(6,2),i,ifl,iflp,clength
      parameter(clength=132)
      character*5 name,nsro,Location*61,totline*(clength)
      character*2 cyear,cnetwork,height*6
      character station_file*(*)
      real lat,long



      iu = 1
      open(iu,file=station_file,status='old',err=2001)
      read(iu,'(1x/1x/1x)')             ! Skip header
      print *,'Searching station...'

c-----------------------------------------------------------------
c Read line by line and parse them.                               |
c-----------------------------------------------------------------

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
        name=totline(iborder(2,1):iborder(2,2))
        cnetwork=totline(iborder(3,1):iborder(3,2))
        read(totline(iborder(4,1):iborder(4,2)),'(f10.6)') lat
        read(totline(iborder(5,1):iborder(5,2)),'(f11.6)') long
        height=totline(iborder(6,1):iborder(6,2))
        location=totline(iborder(6,2)+2:clength)
      if(nsro.ne.name) goto 10                  ! station

      print *
      print '(a6,a5,a2,a61)','Found ',name,': ',Location
      print '(a12,f10.6)', 'Latitude :  ',lat
      print '(a12,f11.6)', 'Longitude:  ',long
      print *
      

      close(iu)
      return

 666  stop 'Station not found, try another one!!'
 2001 stop 'getStation: Did not find station file!'
      end
