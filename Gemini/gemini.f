c************************************************************************
      program            G E M I N I
c
c---> This is version 2.2: Improved routine for finding starting radius
c     March 1999
c
c     version 2.1: full implementation of tranversal isotropic
c     February 1998        earth models
c
c     version 2.0: new routines to find starting radius
c     April 1997
c
c CALLS: aasmt, gtmanis, lyrfnd, gtllim, stavani, propag, ixone, stripmn
c
c  Green's function of the Earth by MINor Integration         
c  Program to compute Green's matrix at a given source radius 
c  including damping trick to inhibit time aliasing
c************************************************************************

      integer nlayer,lmx,nvmax,charalength
      parameter(nlayer=20,lmx=10000,nvmax=6,charalength=80)
      external odemin,odepco,oderad,odetor
      double complex ystart(nvmax)
      double complex m1(6),m2(6),bs(4,2),det,zdampm(nlayer),
     &   zdampk(nlayer),zom,zi
      real*8 x1psv,x1psv_old,x1sh,vpv(nlayer,4),vph(nlayer,4)
     &  ,vsv(nlayer,4),vsh(nlayer,4),rho(nlayer,4),eta(nlayer,4),
     &   qm(nlayer),qk(nlayer),om,elp1,rb(0:nlayer),eps,pi,tlendp,taudp,
     &   flow,fhigh,depth,dom,piti2,df,rsourc,sigma,x,
     &   fmhz,vel,x1,h1,x2,fref,omref,phimu,phika
      integer iflso(nlayer),nw,nlay,nl,l,i,j,iunit,iprint
     &   ,nba,nbe,ifgoon,ntop,ifliq,ifcalder,
     &   lmin,lmax,nvar,ll,nb,nlstrt,ifdisp,ifd,ichoose

c Input declaration of header content
      include "rwdata.h"
      
      character*(charalength) cwind,coutfl,cmnotstr
      character csyst*1
      dimension h1(nlayer)
      common /modls/   iflso   /modq/ qm,qk
      common /modanis/ rho,vpv,vph,vsv,vsh,eta
      common/modr/rb/modlay/nlay/lyrdmp/zdampm,zdampk
      common/loop/zom,om,elp1,l
      common/acc/eps /layer/ nl
      common  /calder/ ifcalder  /whichs/ csyst
      data ifdisp / 1 /

ccc      call buggy_catch_overflow
ccc      call buggy_catch_underflow
ccc      call buggy_catch_inexact


c-----------------------------------------------------------------
c Evaluate memo-array for SODE-matrix storage.                    |
c Set some constants (imaginary unit, Pi, 2*Pi).                  |
c No left hand side of SODE required (ifcalder)                   |
c-----------------------------------------------------------------

      call aasmst
      zi = (0.d0,1.d0)
      pi = 4.d0*datan(1.d0)
      piti2 = 2.d0*pi
      ifcalder = 0

c-----------------------------------------------------------------
c       Input.                                                    |
c-----------------------------------------------------------------

      iunit=5
      print *, 'Choose P-SV -> 1      SH -> 2     P-SV + SH -> 3'
      read(iunit,*) ichoose
      print *,'Control-output print level'
      read(iunit,*) iprint
      print *,'Seismogram length in seconds'
      read(iunit,*) tlendp
      print *,'Damping time in seconds'
      read(iunit,*) taudp
      print *,'Lowest frequency in mHz'
      read(iunit,*) flow
      print *,'Highest frequency in mHz'
      read(iunit,*) fhigh
      print *,'Take into account dispersion:  1 -> yes   0 -> no'
      read(iunit,*) ifd
      print *,'Lowest degree of spherical harmonics'
      read(iunit,*) l1
      print *,'Highest degree of spherical harmonics'
      read(iunit,*) l2
      print *,'Step in the degree-domain'
      read(iunit,*) ldel
      print *,'Source depth in km'
      read(iunit,*) depth
      print *,'Accuracy'
      read(iunit,*) eps
      print *,'Earth model '
      read(iunit,'(a)') cmnotstr
      print *,'File name of window in the om-l-plane'
      read(iunit,'(a)') cwind
      print *,'Name of output file'
      read(iunit,'(a)') coutfl

      call stripmn(cmodel,cmnotstr,charalength)

c Transfer to frequency points and adjusting if input seems to
c be senseless.
      df     = 1.d0/tlendp
      dom    = piti2*df
      nfhigh = idnint(fhigh/(1000.d0*df))
      nflow  = idnint(flow/(1000.d0*df))
      nflow  = max(1,nflow)
      nfhigh = max(nflow,nfhigh)
      flow   = nflow*df
      fhigh  = nfhigh*df
      ifdisp = max( 0, min(ifdisp,ifd) )
      sigma  = 1.d0/taudp ! Damping constant in exp(-sigma*t)=exp(-t/tau) 
                          ! for complex frequency.


c Adjusting l-range if need be.
      l1=max(0,l1)
      ldel = max(1,ldel)
      if (l2-l1 .gt. lmx) stop 'l-range too big for arrays, resize them!!'

c Set motion-info-string for basis solution file
      if (ichoose.eq.1) then
        cmotion='P-SV   '
      else if (ichoose.eq.2) then
        cmotion='SH     '
      else if (ichoose.eq.3) then
        cmotion='P-SV-SH'
      endif

c Set what basic solutions are needed 
      nba=1
      nbe=2


c           Echo input and possibility to cancel.
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      print '(2a)','Type of motion:',cmotion
      print '(a,i2)','Control-output print level:',iprint
      print '(a,1p,2g15.5)','Seism. length, df : ',tlendp,df
      print '(a,1p,g15.5)','Damping time ',taudp
      print '(a,1p,2g13.3,2i5)','f1,f1,nf1,nf2 ',flow,fhigh,nflow,nfhigh
      print '(a,i5)','Dispersion ?:',ifdisp
      print '(a,3i10)','Degree min/max, delta-ell ',l1,l2,ldel
      print '(2a)','Earth model: ',cmodel
      print '(a,f10.2)','Source depth: ',depth
      print '(a,1pe10.0e1)','Accuracy: ',eps
      print '(2a)','Window file name: ',cwind
      print '(2a)','Output file name: ',coutfl
      print *, 'Everything allright? Type <1>'
      read(iunit,*) ifgoon
      if(ifgoon.ne.1) stop


c Read model (should be a realistic one. NO OCEAN, please!!).
      CALL gtmanis(cmnotstr,fref)
      omref = piti2*fref   ! Reference frequency of anelasticity (Hertz)



c Set top layer for integration downwards.
      ntop   = nlay
      rsourc = rb(nlay)-depth
      rsrc   = sngl(rsourc)
      rrecv  = sngl(rb(nlay))


c Get model parameters at source and receiver radius for later use
c in program which builds moment tensor and source terms.
      CALL lyrfnd(rsourc,nl,ifliq)
      x   = rsourc/rb(nlay)
      ros  = sngl(rho(nl,1)+x*(rho(nl,2)+x*(rho(nl,3)+x*rho(nl,4))))
      vpverts = sngl(vpv(nl,1)+x*(vpv(nl,2)+x*(vpv(nl,3)+x*vpv(nl,4))))
      vsverts = sngl(vsv(nl,1)+x*(vsv(nl,2)+x*(vsv(nl,3)+x*vsv(nl,4))))
      vphoris = sngl(vph(nl,1)+x*(vph(nl,2)+x*(vph(nl,3)+x*vph(nl,4))))
      vshoris = sngl(vsh(nl,1)+x*(vsh(nl,2)+x*(vsh(nl,3)+x*vsh(nl,4))))
      ettas   = sngl(eta(nl,1)+x*(eta(nl,2)+x*(eta(nl,3)+x*eta(nl,4))))

c x equals 1, therefore:
      ror  = sngl(rho(nlay,1)+rho(nlay,2)+rho(nlay,3)+rho(nlay,4))
      vpvertr = sngl(vpv(nlay,1)+vpv(nlay,2)+vpv(nlay,3)+vpv(nlay,4))
      vsvertr = sngl(vsv(nlay,1)+vsv(nlay,2)+vsv(nlay,3)+vsv(nlay,4))
      vphorir = sngl(vph(nlay,1)+vph(nlay,2)+vph(nlay,3)+vph(nlay,4))
      vshorir = sngl(vsh(nlay,1)+vsh(nlay,2)+vsh(nlay,3)+vsh(nlay,4))
      ettar   = sngl(eta(nlay,1)+eta(nlay,2)+eta(nlay,3)+eta(nlay,4))

      iflqs = iflso(nl)
      iflqr = iflso(nlay)
      qms = sngl(qm(nl))
      qks = sngl(qk(nl))
      qmr = sngl(qm(nlay))
      qkr = sngl(qk(nlay))
      tlen = sngl(tlendp)
      tau  = sngl(taudp)

c Open file for Green's basic solution
      nb1   = 1
      nb2   = 4
      nb3   = 1
      nli1  = 1
      nli2  = 3
      nli3  = 2
      iunit = 50
      open(iunit,file=coutfl,form='unformatted',status='new',err=888)
      write(iunit) cmotion,cmodel,
     &    rsrc,ros,vpverts,vsverts,vphoris,vshoris,ettas,iflqs,qms,qks,
     &    rrecv,ror,vpvertr,vsvertr,vphorir,vshorir,ettar,iflqr,qmr,qkr,
     &    nli1,nli2,nli3,nb1,nb2,nb3,tlen,tau,nflow,nfhigh,l1,l2,ldel


c______________________________________________________________________
c
c  INTEGRATION PART INTEGRATION PART INTEGRATION PART INTEGRATION PART |
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


      do 1000 nw = nflow, nfhigh     ! Loooooooooooooop over frequency.

      
      om   = nw*dom
      fmhz = om*1000.d0/piti2

c Make complex frequency for Laplace transform. Calculate anelasticity
c factor in every layer. Compute initial stepsize for Bulirsch-Stoer
c integration algorithm from wave-length of P-/S-waves at layer bottoms.

      zom = dcmplx(om,-sigma)
      do 1 nl = 1, nlay
       if (ifdisp.eq.1) then         ! Dispersion -> complex moduli
         
         phimu      = atan(qm(nl))
         phika      = atan(qk(nl))
         zdampm(nl) = (zi*zom/omref)**(2.d0*phimu/pi)
         zdampk(nl) = (zi*zom/omref)**(2.d0*phika/pi)

       else if (ifdisp.eq.0) then    ! no Dispersion -> real moduli
         zdampm(nl) = 1.d0
         zdampk(nl) = 1.d0
       endif
       x   = rb(nl-1)/rb(nlay)
       vel = vsv(nl,1)+x*(vsv(nl,2)+x*(vsv(nl,3)+x*vsv(nl,4)))
       if (iflso(nl).eq.1) then
         vel = vpv(nl,1)+x*(vpv(nl,2)+x*(vpv(nl,3)+x*vpv(nl,4)))
       endif
       h1(nl) = piti2/4.d0*vel/om
 1    continue


      CALL gtllim(sngl(fmhz),lmin,lmax,cwind)     !   Set lmin and lmax.
      lmax = min0(lmax,l2)
      lmin = max0(lmin,l1)


       l = lmin
 1100      if (l.le.lmax) then ! Loooooooooooooooooooooooooooooooooop over l.
      elp1 = dble(l*(l+1))
      ll   = l-lmin+1

      CALL ixone(x1psv,nlstrt)       ! Find starting radius for integration
c     CALL wixone(x1psv,nlstrt)   ! Find starting radius-wolle-version
c     print *, 'Startradius:', x1psv
      
      if (ichoose.eq.2) goto 2000        ! SH-motion only


      if (l.eq.0) then

c             ++++++++++++++++++++++++++++++++++++++++++++++
c             + Calculate RADIAL MOTION coefficients (l=0) +
c             ++++++++++++++++++++++++++++++++++++++++++++++

c Get starting values for integration UPWARDS.
      x1     = x1psv                ! Start integration at the lowermost level.
      x1psv_old = x1psv   
      CALL stavani(ystart,x1,nlstrt,'s')  ! Get starting values there.

c Propagate solution from starting radius to source.
      nvar   = 2
      csyst  = 'R'
      x2     = rsourc
      CALL propag(x1,x2,ystart,h1,oderad,nvar)
      m1(1)  = ystart(1)
      m1(2)  = ystart(2)


c  Integration downwards to the source.
      ystart(1) = 1.d0
      ystart(2) = 0.d0
      x1        = rb(nlay)
      x2        = rsourc
      CALL propag(x1,x2,ystart,h1,oderad,nvar)
      m2(1)     = ystart(1)
      m2(2)     = ystart(2)


c Construct Green's basic solution for radial motion
       det            =  m1(1)*m2(2)-m1(2)*m2(1)
       basout(1,1,ll) = cmplx(-m1(2)/det)
       basout(1,2,ll) = cmplx(m1(1)/det)
       basout(1,3,ll) =  0.
       basout(1,4,ll) =  0.
 
c------------(l.eq.0)
      else 
c------------(l.gt.0)

c        ++++++++++++++++++++++++++++++++++++++++++++++++++++++
c        + Calculate SPHEROIDAL MOTION coefficients (l>0)     +
c        ++++++++++++++++++++++++++++++++++++++++++++++++++++++


c                 Integration of the minors upwards.
c                 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      nvar  =  5
      csyst = 'M'
      if (x1psv.ge.rsourc) then
c Leave l-loop, adjusting lmax
        lmax = l-1
        goto 1100
      endif
      x1 = x1psv
      x1psv_old = x1psv
      CALL stavani(ystart,x1,nlstrt,'s')
      x2 = rsourc
      CALL propag(x1,x2,ystart,h1,odemin,nvar)
      do 102 i = 1, nvar
 102  m1(i) =  ystart(i)
      m1(6) = -m1(1)/elp1


c                 Integration of the minors downwards.
c                 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ystart(1) = 0.d0
      ystart(2) = 1.d0
      ystart(3) = 0.d0
      ystart(4) = 0.d0
      ystart(5) = 0.d0
      x1        = rb(nlay)
      x2        = rsourc
      CALL propag(x1,x2,ystart,h1,odemin,nvar)
      do 104 i = 1, nvar
 104  m2(i) =  ystart(i)
      m2(6) = -m2(1)/elp1


c         Calculate two BASIC SOLUTIONS of the Green's function
c         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c         downwards from the earth's surface.
c         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c Initial values of basic solutions with vanishing stress at the 
c free surface.

       bs(1,1) = 1.d0
       bs(2,1) = 0.d0
       bs(3,1) = 0.d0
       bs(4,1) = 0.d0

       bs(1,2) = 0.d0
       bs(2,2) = 0.d0
       bs(3,2) = 1.d0
       bs(4,2) = 0.d0

      nvar  =  4
      csyst = 'C'
      do 3200 nb = nba, nbe   ! Looooooooooooop over basic solutions.
          do 115 i = 1, 4
 115      ystart(i) = bs(i,nb)
          x1 = rb(ntop)
          x2 = rsourc
          CALL propag(x1,x2,ystart,h1,odepco,nvar)
          do 105 i = 1, nvar
 105      bs(i,nb) = ystart(i)
 3200 continue                ! Eeeeend of loop over basic solutions.

c                 Green's basic solutions.
c                 ~~~~~~~~~~~~~~~~~~~~~~~~

      det =  m1(1)*m2(6) + m1(6)*m2(1) - m1(2)*m2(5)
     &     + m1(3)*m2(4) + m1(4)*m2(3) - m1(5)*m2(2)
        basout(1,1,ll)=cmplx(
     &    (m1(4)*bs(4,2)-m1(5)*bs(3,2)+m1(6)*bs(2,2))/det  )
        basout(1,2,ll)=cmplx(
     &    (m1(3)*bs(3,2)-m1(6)*bs(1,2)-m1(2)*bs(4,2))/det  )
        basout(1,3,ll)=cmplx(
     &    (m1(1)*bs(4,2)-m1(3)*bs(2,2)+m1(5)*bs(1,2))/det  )
        basout(1,4,ll)=cmplx(
     &    (m1(2)*bs(2,2)-m1(1)*bs(3,2)-m1(4)*bs(1,2))/det  )
        basout(3,1,ll)=cmplx(
     &        (-m1(4)*bs(4,1)-m1(6)*bs(2,1)+m1(5)*bs(3,1))/det  )
        basout(3,2,ll)=cmplx(
     &        (m1(6)*bs(1,1)-m1(3)*bs(3,1)+m1(2)*bs(4,1))/det  )
        basout(3,3,ll)=cmplx(
     &        (m1(3)*bs(2,1)-m1(5)*bs(1,1)-m1(1)*bs(4,1))/det  )
        basout(3,4,ll)=cmplx(
     &        (m1(4)*bs(1,1)-m1(2)*bs(2,1)+m1(1)*bs(3,1))/det  )

      endif   ! l>0


      if(ichoose.eq.1)  goto 1099     ! P-SV-motion only


 2000 if (l.gt.0) then


c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c          Calculate TOROIDAL MOTION
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c Integr. through mantle (no toroidal motion in liquid outer core!).
      nvar  =  2
      csyst = 'T'
      x1    = x1psv
      CALL stavani(ystart,x1,nlstrt,'t')
      x1sh  = x1
      if (x1sh.ge.rsourc) then
c Leave l-loop, adjusting lmax
        lmax = l-1
        goto 1100
      endif
      x2    = rsourc
      CALL propag(x1,x2,ystart,h1,odetor,nvar)
      m1(1) = ystart(1)
      m1(2) = ystart(2)


c                      Integration downwards.
c                      ~~~~~~~~~~~~~~~~~~~~~~
      ystart(1) = 1.d0
      ystart(2) = 0.d0
      x1        = rb(ntop)
      x2        = rsourc
      CALL propag(x1,x2,ystart,h1,odetor,nvar)
      m2(1)     = ystart(1)
      m2(2)     = ystart(2)


c Construct Green's basic solution for toroidal motion
       det            =  m1(1)*m2(2)-m1(2)*m2(1)
       bsoutt(1,1,ll) = cmplx( -m1(2)/det )
       bsoutt(1,2,ll) = cmplx( m1(1)/det )


c-------------(l.gt.0)
      endif
      
 1099  if (iprint.ge.1) then 
        if (mod(l,iprint).eq.0) write(*,998) fmhz,l,x1psv_old,x1sh
       endif

cEeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeend of loop over l
           l = l + ldel
      goto 1100
      endif

c                 Output
c                 ~~~~~~

      if (ichoose.eq.3) then
        write(iunit) lmin,lmax,(((basout(i,j,l),i=1,3,2),
     &    j=1,4),(bsoutt(1,j,l),j=1,2),l=1,lmax-lmin+1,ldel)
      else if (ichoose.eq.2) then
        write(iunit) lmin,lmax,((bsoutt(1,j,l),j=1,2),l=1,lmax-lmin+1,ldel)
      else if (ichoose.eq.1) then
        write(iunit) lmin,lmax,(((basout(i,j,l),i=1,3,2),
     &                 j=1,4),l=1,lmax-lmin+1,ldel)
      endif
c        write(88,'(d15.5)') (zabs(basout(1,2,l)),l=1,lmax-lmin+1)


       if(iprint.le.0) write(*,999) fmhz,nw,lmin,lmax,x1psv_old,x1sh

 1000 continue   ! Eeeeeeeeeeeeeeeeeeeeeeeeeeend of loop over frequeny.

      stop
 998  format(f7.1,i6,2f8.0)
 999  format(f7.1,3i6,2f8.0)
 888  print *, 'Datei scheint schon zu existieren!'
      end

