c
      program              D  I  S  P  E  C
c
c
c  Calculates spectra for displacement in vertical and horicontal
c  direction. Optionally, a section of 'nseis' seismograms between
c  the source and the receiver can be computed. 
c
c  Main program activities:
c     1.  Read in running parameters
c     2.  Get receiver locations
c     3.  Get source parameters
c     4.  Transform from geographical coordinates into source coordinates
c           (epicentral distance and azimuth)
c     5.  Compute spherical harmonics for all receiver source coordinates
c     6.  Read in Green basis solutions calculated by 'gemini' and perform
c           sum over spherical harmonics; transform back to geographical
c           coordinates.
c
c
      complex zi
      real pi,dg2rd
      integer nfr,lmx,mm,ntmx,nseimx
      parameter( nfr  = 3000,  lmx  = 10000,   mm=2,
     &     zi = (0.,1.),  pi = 3.1415926535897932, dg2rd=pi/180.,
     &     ntmx = 200,  nseimx = 500  )
      complex u(0:lmx,-mm:mm),v(0:lmx,-mm:mm),w(0:lmx,-mm:mm),
     &   y(0:lmx,-mm:mm,nseimx),ydt(0:lmx,-mm:mm,nseimx),
     &   dydt(0:lmx,-mm:mm,nseimx),zsrc(6,-mm:mm),dsp(3,nfr,nseimx),
     &   zsw1,zsw2,basout(4,4,0:lmx),bsoutt(2,2,0:lmx)
      real snorm(0:lmx,0:mm),rmom(6),force(3),SourceLat,SourceLong,SourceDepth
     &     ,rho,delta(nseimx),phi(nseimx),epidis,azimuth,bazi
     &     ,RLat(nseimx),RLong(nseimx),csa(nseimx),sna(nseimx)
     &     ,tlen,tau,sigma,dom,fmhz,qm,qk,rsourc,rrecv,om
     &     ,tapend(ntmx),rhor,qmr,qkr,RecLat,RecLong,RecLat2,RecLong2
     &     ,vpverts,vsverts,vphoris,vshoris,ettas
     &     ,vpvertr,vsvertr,vphorir,vshorir,ettar
      integer mmax,nt2,iflso,nb1,nb2,nb3,l1,l2,ns,i,j,iflsor,
     &  l,m,nw,lmin,lmax,lbeg,lend,nbeg,nend,ifrotgeo
     &  ,nli1,nli2,nli3,iks,lm,iu,lstart,ntoro,nseis,lstep
      character*5 RECname(nseimx),netname(nseimx)*2
      character*80 cevent*82,cgrsph,cwind,fsname,model_name*10,SourceMech*1
      character date*6,time*10,SourceType*20,cmotion*7,SourceCode*8
      common/tap/tapend
      common / cmt_date / date,time,SourceCode
      data ifrotgeo /0/

c     call buggy_catch_overflow
c     call buggy_catch_underflow

c-----------------------------------------------------------------
c  Input running parameters.                                      |
c-----------------------------------------------------------------

      iu=5
      print *, 'File with Green-basicsolutions ?'
      read(iu,'(a80)') cgrsph
      print *, 'File with source parameters?'
      read(iu,'(a70)') cevent
      print *, 'Source mechanism: Moment tensor(m) | Single force(f)'
      read(iu,'(a1)') SourceMech
      print *, 'Maximum order in sum over spherical harmonics?'
      read(iu,*) mmax
      print *, 'File with boundaries in om-l-domain ?'
      read(iu,'(a80)') cwind
      print *, 'Taper at the end of l-interval? '
      read(iu,*) nt2
      print *, 'Name of file with stations? '
      read(iu,'(a80)') fsname
      print *,
     &  'Geogr. coordinates: 1, Station abbrev.: 2, File with coord.list: 3'
      read(iu,*) iks
      if (iks.eq.1) then
       read(iu,*) RecLat,RecLong
       RLat(1) = RecLat
       RLong(1) = RecLong
       RECname(1) = 'UNKN '                 ! Station abbreviation UNKNown
       netname(1) = 'UK'
       nseis=1
      elseif (iks.eq.2) then
       read(iu,'(a5)') RECname(1)
       netname(1) = 'UK'
       nseis=1
      elseif (iks.eq.3) then
       read(iu,'(1x)')
      else if (iks.eq.4) then
       read(iu,*) RecLat,RecLong, RecLat2, RecLong2, nseis
       RLat(1) = RecLat
       RLong(1) = RecLong
       RECname(1) = 'UNKN '                 ! Station abbreviation UNKNown
       netname(1) = 'UK'
      endif
      print *,'Rotate to NS-EW-coordinates ( n > 0: yes)?'
      read(iu,*) ifrotgeo

c----------------------------------------------------------
c Read file with source parameters and generate moment     |
c tensor or single force.                                  |
c----------------------------------------------------------

      if (SourceMech.eq.'m') then
        call readCMT(cevent,rmom,SourceLat,SourceLong,SourceDepth)
        mmax = min(mmax,2)
      elseif (SourceMech.eq.'f') then
        call readSForce(cevent,force,SourceLat,SourceLong,SourceDepth)
        mmax = min(mmax,1)
      endif
      SourceType=SourceCode

c-----------------------------------------------------------------
c  Receiver location: coordinates from stdin | file               |
c   make many receivers for section                               |
c-----------------------------------------------------------------

       if (iks.eq.2) then
         call getStation(fsname,RECname(1),netname(1),RecLat,RecLong)
         RLat(1) = RecLat
         RLong(1) = RecLong
       else if (iks.eq.3) then
         call allStations(fsname,RECname,netname,RLat,RLong,nseimx,nseis)
       endif
       if (nseis.gt.nseimx) stop 'Array dimension <nseimx> too small!'

c-----------------------------------------------------------------
c Compute epicentral distance, azimut and backazimut for each     |
c receiver. Then specify angle between epicentral and south       |
c direction for transform to south-north and east-west components.|
c-----------------------------------------------------------------

      if (nseis.eq.1) then
       call epitra(SourceLat,SourceLong,RecLat,RecLong,epidis,azimuth,bazi)
       print 999,'Epicentral distance in degree : ',epidis
       print 999,'Azimuth        : ',azimuth
       print 999,'Back-Azimuth   : ',bazi
       delta(1) = epidis
       phi(1)   = 180.-azimuth
       csa(1) = cos(dg2rd*bazi)
       sna(1) = -sin(dg2rd*bazi)
      elseif (nseis.gt.1.and.iks.eq.4) then       ! section
       call epitra(RecLat,RecLong,RecLat2,RecLong2,epidis,azimuth,bazi)
       call necklace(RecLat,RecLong,epidis,azimuth,nseis,RLat,RLong)
       do ns=1,nseis
         call epitra(SourceLat,SourceLong,RLat(ns),RLong(ns),epidis
     &               ,azimuth,bazi)
         csa(ns)   = cos(dg2rd*bazi)
         sna(ns)   = -sin(dg2rd*bazi)
         delta(ns) = epidis
         phi(ns)   = 180.-azimuth
c        write(77,*) RLat(ns), Rlong(ns), epidis,azimuth,bazi
       enddo 
      elseif (nseis.gt.1.and.iks.eq.3) then       ! Many receivers 
       do ns=1,nseis
           call epitra(SourceLat,SourceLong,RLat(ns),RLong(ns),epidis,
     &                 azimuth,bazi)
           csa(ns) = cos(dg2rd*bazi)
           sna(ns) = -sin(dg2rd*bazi)
           delta(ns)= epidis
           phi(ns)   = 180.-azimuth
           write(*,*) RLat(ns), Rlong(ns), epidis,azimuth,bazi
       enddo 
      endif

c-----------------------------------------------------------------
c Open file with Green - basis solutions.                         |
c-----------------------------------------------------------------

      open(10,file=cgrsph,form='unformatted',status='old')
      read(10) cmotion,model_name,rsourc,rho,
     &  vpverts,vsverts,vphoris,vshoris,ettas,iflso,qm,qk,
     &  rrecv,rhor,vpvertr,vsvertr,vphorir,vshorir,ettar,iflsor,qmr,qkr,nli1,
     &  nli2,nli3,nb1,nb2,nb3,tlen,tau,nbeg,nend,l1,l2,lstep
      if(nend.gt.nfr) then
        print *, 'Frequency array dimension <nfr> not sufficient!'
        print *, 'Need at least',nend,' elements.'
        stop
      endif
      sigma = 1./tau
      dom=2.*Pi/tlen
      if (nli2.eq.3) then
       ntoro=1
      else 
       ntoro = 2
      endif
      lstep = max(1,lstep)      ! l-interval, normally = 1
      call costap(nt2)       !  Taper in om-ell-domain

c-----------------------------------------------------------------
c Control output.                                                 |
c-----------------------------------------------------------------

      print '(2a)','Type of motion:',cmotion
      print '(2a)','Earth model: ',model_name
      print '(a)','Source radius, rho, vpv, vsv, vph, vsh, qm, qk'
      print '(7f6.1,2e10.3)', rsourc,rho,vpverts,vsverts,vphoris,vshoris,
     &       ettas,qm,qk
      print '(a)','Receiver radius, rho, vpv, vsv, vph, vsh, qm, qk'
      print '(7f6.1,2e10.3)', rrecv,rhor,vpvertr,vsvertr,vphorir,vshorir,
     &       ettar,qmr,qkr
      print '(a,1p,g15.5)','Seism. length : ',tlen
      print '(a,1p,g15.5)','Damping time ',tau
      print '(a,3i10)','Degree min/max, delta-ell ',l1,l2,lstep
      
c-----------------------------------------------------------------
c Calculate fully normalized spherical harmonics and derivatives  |
c (for each epicentral distance). Get  Normalization factor for   |
c source terms.                                                   |
c-----------------------------------------------------------------

      lm = min(lmx,l2)
      do ns = 1,nseis
         call sharm(lmx,lm,mm,mmax,delta(ns),phi(ns),y(0,-mm,ns),
     &            ydt(0,-mm,ns),dydt(0,-mm,ns))
      enddo
      call stnorm(lmx,lm,mm,snorm)
c     write(88,'(2i5,6e13.5)') ((l,m,y(l,m,1),dydt(l,m,1),ydt(l,m,1),l=0,lm),m=1,1)


c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c~~~~~~~~~~~~~~~~~~~~~~~~~~ Loop over frequencies  ~~~~~~~~~~~~~~~~
      do 1000 nw=nbeg,nend
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      om=nw*dom
      fmhz = om*1000./(2.*Pi)

c-----------------------------------------------------------------
c Read in Green - basic solutions                                 |
c-----------------------------------------------------------------

      if (cmotion.eq.'P-SV-SH') then
           read(10,err=888,end=888) lmin,lmax,
     &         (((basout(i,j,l),i=nli1,nli2,nli3),
     &       j=nb1,nb2,nb3),((bsoutt(i,j,l),i=1,ntoro),j=nb1,2),
     &        l=lmin,lmax,lstep)
      else if (cmotion.eq.'SH     ') then
           read(10,err=888,end=888) lmin,lmax,
     &     (((bsoutt(i,j,l),i=1,ntoro),j=nb1,2),l=lmin,lmax,lstep)
      else if (cmotion.eq.'P-SV   ') then
           read(10,err=888,end=888) lmin,lmax,
     &            (((basout(i,j,l),i=nli1,nli2,nli3),
     &                   j=nb1,nb2,nb3),l=lmin,lmax,lstep)
      endif

c-----------------------------------------------------------------
c Complex Model parameters in source layer, provided via common   |
c block for 'momfor'.                                             |
c Get lmin and lmax for actual frequency.                         |
c-----------------------------------------------------------------

      call elapa(om,sigma,qm,qk,rho,vpverts,vsverts,vphoris,vshoris,ettas)
      call gtllim(fmhz,lbeg,lend,cwind)
      lmin=max(lmin,lbeg)
      lmax=min(lmax,lend)
      if(lmax.gt.lmx) STOP 'Array size for l needs adjusting  !'

c-----------------------------------------------------------------
c Source terms                                                    |
c-----------------------------------------------------------------

      if ( SourceMech .eq.'m') then
          call momfor(rsourc,mm,rmom,zsrc,'m')
      elseif (SourceMech .eq. 'f' ) then
          call momfor(rsourc,mm,force,zsrc,'f')
      endif

c-----------------------------------------------------------------
c Make expansion coefficients for spherical harmonic series from  |
c basic solutions and source terms.                               |
c-----------------------------------------------------------------

      do 1160 m=-mmax,mmax
        do 1150 l=lmin,lmax,lstep
          u(l,m) = 0.
          v(l,m) = 0.
          do i=nb1,nb2,nb3
             u(l,m)=u(l,m)+basout(1,i,l)*zsrc(i,m)
             v(l,m)=v(l,m)+basout(3,i,l)*zsrc(i,m)
          enddo
          u(l,m)=u(l,m)*snorm(l,abs(m))
          v(l,m)=v(l,m)*snorm(l,abs(m))
          w(l,m) = (bsoutt(1,1,l)*zsrc(5,m)
     &             +bsoutt(1,2,l)*zsrc(6,m))*snorm(l,abs(m))
 1150   continue
 1160 continue

c-----------------------------------------------------------------
c Apply taper to avoid cut-off-effects.                           |
c-----------------------------------------------------------------

      if(nt2.gt.0)then
       lstart = max(lmin,lmax-nt2+1)
       do 10 m=-mmax,mmax
       do 10 l=lstart,lmax,lstep
        u(l,m)=u(l,m)*tapend(l+nt2-lmax)
        v(l,m)=v(l,m)*tapend(l+nt2-lmax)
        w(l,m)=w(l,m)*tapend(l+nt2-lmax)
 10    continue
      endif
      if(mod(nw,10).eq.0) print *, fmhz,lmin,lmax,nt2

c-----------------------------------------------------------------
c Summation over spherical harmonics for all receivers.           |
c r: dsp(1,nw) , theta: dsp(2,nw) , phi: dsp(3,nw)                |
c-----------------------------------------------------------------


      do 1200 ns = 1,nseis

      dsp(1,nw,ns)=0.
      dsp(2,nw,ns)=0.
      dsp(3,nw,ns)=0.
       do 1211 l=lmin,lmax,lstep
        dsp(1,nw,ns) = dsp(1,nw,ns)+u(l,0)*y(l,0,ns)
        dsp(2,nw,ns) = dsp(2,nw,ns)+v(l,0)*dydt(l,0,ns)
        dsp(3,nw,ns) = dsp(3,nw,ns)-w(l,0)*dydt(l,0,ns)
 1211 continue
       if (mmax.ge.1) then
       do 1212 l=lmin,lmax,lstep
        dsp(1,nw,ns) = dsp(1,nw,ns)+u(l,1)*y(l,1,ns)+u(l,-1)*y(l,-1,ns)
        dsp(2,nw,ns) = dsp(2,nw,ns)+
     &                 v(l,1)*dydt(l,1,ns)+w(l,1)*ydt(l,1,ns)+
     &                 v(l,-1)*dydt(l,-1,ns)+w(l,-1)*ydt(l,-1,ns)
        dsp(3,nw,ns) = dsp(3,nw,ns)+v(l,1)*ydt(l,1,ns)-w(l,1)*dydt(l,1,ns)
     &                 +v(l,-1)*ydt(l,-1,ns)-w(l,-1)*dydt(l,-1,ns)
 1212  continue
       endif
       if (mmax.eq.2) then
       do 1213 l=lmin,lmax,lstep
        dsp(1,nw,ns) = dsp(1,nw,ns)+u(l,2)*y(l,2,ns)+u(l,-2)*y(l,-2,ns)
        dsp(2,nw,ns) = dsp(2,nw,ns)+
     &                 v(l,2)*dydt(l,2,ns)+w(l,2)*ydt(l,2,ns)+
     &                 v(l,-2)*dydt(l,-2,ns)+w(l,-2)*ydt(l,-2,ns)
        dsp(3,nw,ns) = dsp(3,nw,ns)+v(l,2)*ydt(l,2,ns)-w(l,2)*dydt(l,2,ns)
     &                 +v(l,-2)*ydt(l,-2,ns)-w(l,-2)*dydt(l,-2,ns)
 1213 continue
      endif

      if (ifrotgeo.ne.0) then 
c-----------------------------------------------------------------
c Transform to north-south and east-west components.              |
c          2 -> North-South  3 -> East-West                       |
c I.e. a positive swing points towards north and east!!           |
c-----------------------------------------------------------------
        zsw1 = dsp(2,nw,ns)
        zsw2 = dsp(3,nw,ns)
        dsp(2,nw,ns) = -csa(ns)*zsw1+sna(ns)*zsw2
        dsp(3,nw,ns) =  sna(ns)*zsw1+csa(ns)*zsw2
      endif

 1200 continue


 1000 continue         !------------------------- End of frequency loop
  889 close(10)

c-----------------------------------------------------------------
c Write spectra for vertical and horicontal components.           |
c-----------------------------------------------------------------

      open(40,file='spec3k',form='unformatted')
       write(40) SourceType,date,time,SourceLat,SourceLong,SourceDepth,
     &    nbeg,nend,tlen,tau,nseis,(RECname(ns)(1:5),netname(ns),RLat(ns)
     &    ,RLong(ns),((dsp(i,nw,ns),i=1,3),nw=nbeg,nend),ns=1,nseis)
       close(40)

      stop
 999  format(a,f8.3)
 888  print *, 'An error occured while reading the basis solutions!'
      print *, 'Re-setting "nend" to the last correct value',nw-1,' and'
      print *, 'saving what we have til now.'
      nend = nw-1
      goto 889
      end





c
c -----   Part for using Legendre polynomes of the second kind (Qlm)
c
clm = min(lmx,l2)
cdo ns = 1,nseis
c        call sharm(lmx,lm,mm,mmax,delta(ns),phi,y(0,-mm,ns),
c    &            ydt(0,-mm,ns),dydt(0,-mm,ns))
c        if (iuseq.eq.1) then
c        call sharmq(lmx,lm,mm,mmax,delta(ns),phi,yq(0,-mm,ns),
c    &            yqdt(0,-mm,ns),dyqdt(0,-mm,ns))
c          do 101 m=-mmax,mmax
c          do 101 l = abs(m),lm
c           y(l,m,ns) = 0.5*(y(l,m,ns)+2.*zi/Pi*yq(l,m,ns))
c           ydt(l,m,ns) = 0.5*(ydt(l,m,ns)+2.*zi/Pi*yqdt(l,m,ns))
c           dydt(l,m,ns) = 0.5*(dydt(l,m,ns)+2.*zi/Pi*dyqdt(l,m,ns))
c101       continue
c        endif
c     enddo
