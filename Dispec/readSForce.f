	Subroutine ReadSForce (cevent,sforce,clat,clon,cdepth)
c
c Subroutine to read source file containing moment tensor and 
c earthquake parameters in  H a r v a r d  Catalogue Format.
c
c It interprets the first three elements of the moment tensor
c as radial, north-south and east-west components of a single
c force:  M_rr -> f_r   M_tt -> f_t   M_pp -> f_p
c
	character*(50) cevent
      character*2 year,month,day,hour,min
      character  time*10, date*6, SourceType*20
      character*8 code
      character*24 region
      character*3 epsrc
      real lat,lon,mb,ms,ev(3),dyne2Nm,Nm2myUnit
      real rmom(6),emrr,emss,emee,emrs,emre,emse,sforce(3)
      real mwst,mwrec,mwcoff,sec,depth,ahdur,smo,
     1            dt,edt,clat,eclat,clon,eclon,cdepth,ecdepth
      integer bwst,bwrec,bwcoff,expo,strike(2),dip(2),rake(2)
      integer evpl(3),evaz(3),i,inut
      parameter(dyne2Nm=1.e-7,Nm2myUnit=1.e-3)
      common / cmt_date / date,time,SourceType

c  Read source file and convert to moment tensor, see source files
c  for further information

      inut = 1
      OPEN(inut,file=cevent)

c ...Harvard Catalog Format!

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
      read(inut,1004)ev(1),evpl(1),evaz(1),
     1               ev(2),evpl(2),evaz(2),
     2	         ev(3),evpl(3),evaz(3),
     3          smo,strike(1),dip(1),rake(1),
     4              strike(2),dip(2),rake(2)
 1004 format(3(f7.2,i3,i4),f7.2,2(i4,i3,i5)) 

      do i = 1, 3
        sforce(i) = rmom(i)
      enddo

      print *
      print *, 'Centroid single force and coordinates:'
      print '(a,i3,6f7.3)', 'Single Force :', expo, (sforce(i),i=1,3)
      print '(a,f8.3)', 'Latitude      :', clat
      print '(a,f8.3)', 'Longitude     :', clon

c Apply scale factor to the tensor elements and convert from Dyne to Nm

      do i=1,3
       sforce(i) = sforce(i)*(10.**expo)*dyne2Nm*Nm2myUnit
      enddo

      write(time,'(2a2,f6.3)') hour,min,sec
      date = year//month//day
      SourceType=code

      close(inut)

	end


