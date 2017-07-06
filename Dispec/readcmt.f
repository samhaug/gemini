      Subroutine readCMT (cevent,rmom,clat,clon,cdepth)
c----                                                 ---------------
c Subroutine to read source file containing moment tensor and 
c earthquake parameters in  H a r v a r d  Catalogue Format.
c With the factors 'dyne2Nm' and 'Nm2myUnit' the resulting seimograms
c will have the unit 'nm/s' (velocity) or 'nm' (displacement).
c--------------------------------------------------------------------

      character*(72) cevent
      character*2 year,month,day,hour,min
      character  time*10, date*6
      character*8 code
      character*24 region
      character*3 epsrc
      real lat,lon,mb,ms,ev(3),dyne2Nm,Nm2myUnit,totsec
      real rmom(6),emrr,emss,emee,emrs,emre,emse
      real mwst,mwrec,mwcoff,sec,depth,ahdur,smo,
     1            dt,edt,clat,eclat,clon,eclon,cdepth,ecdepth
      integer bwst,bwrec,bwcoff,expo,strike(2),dip(2),rake(2)
      integer evpl(3),evaz(3),i,inut,ihour,inmin,ios
      parameter(dyne2Nm=1.e-7,Nm2myUnit=1.e-6)
      common / cmt_date / date,time,code

c-----------------------------------------------------------------
c  Read source file in Harvard Catalog Format with moment tensor, |
c  See the file 'CMT_explained' for further information.          |
c-----------------------------------------------------------------

      inut = 1
      OPEN(inut,file=cevent,iostat=ios,status='old',err=999)
      read(inut,1001)code,month,day,year,hour,min,sec,lat,lon,depth,
     1            mb,ms,region
 1001 format(a8,5(1x,a2),1x,f4.1,f7.2,f8.2,f6.1,2f3.1,a24)     
      read(inut,1002)epsrc,bwst,bwrec,bwcoff,mwst,mwrec,mwcoff,
     1            dt,edt,clat,eclat,clon,eclon,cdepth,ecdepth
 1002 format(a3,2(4x,i2,i3,i4),4x,f6.1,f4.1,
     1       f7.2,f5.2,f8.2,f5.2,f6.1,f5.1)
      read(inut,1003)ahdur,expo,rmom(1),emrr,rmom(2),emss,rmom(3),emee,
     1            rmom(4),emrs,rmom(5),emre,rmom(6),emse
 1003 format(4x,f4.1,4x,i2,6(f6.2,f5.2))
      read(inut,1004)ev(1),evpl(1),evaz(1),ev(2),evpl(2),evaz(2),
     2              ev(3),evpl(3),evaz(3),smo,strike(1),dip(1),rake(1),
     4              strike(2),dip(2),rake(2)
 1004 format(3(f7.2,i3,i4),f7.2,2(i4,i3,i5)) 
      print *
      print *, 'Centroid moment tensor and coordinates:'
      print '(a,i3,6f7.3)', 'Moment tensor:', expo, (rmom(i),i=1,6)
      print '(a,3f8.3)', 'Centroid-Lat/Long/Depth:', clat,clon,cdepth

c-----------------------------------------------------------------
c Apply scale factor to the tensor elements and convert from Dyne |
c to Nm.                                                          |
c Correct time for centroid-delay.                                |
c-----------------------------------------------------------------

      do i=1,6
       rmom(i) = rmom(i)*(10.**expo)*dyne2Nm*Nm2myUnit
      enddo
      print *,'PDE time:', hour,min,sec
      read(hour,'(i2)') ihour
      read(min,'(i2)') inmin
      totsec = float(ihour)*3600. + float(inmin)*60. + sec + dt
      ihour = int(totsec/3600.)
      inmin = int(mod(totsec,3600.)/60.)
      sec = mod(totsec,60.)
      print *,'CMT time:', ihour,inmin,sec
      write(time,'(2i2,f6.3)') ihour,inmin,sec
      date = year//month//day
c     print *, 'date,time:',date,time

      close(inut)
      return

999   print *,'File status=',ios
      pause 'Error in readCMT!!'
      end


