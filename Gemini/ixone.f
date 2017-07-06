      subroutine ixone(rstart,nla)
c
ccccalls RTBIS, PHASI
c
c Subroutine to find the starting radius for integration
c
      integer nlayer,layer,j,nlay,iaob,nla,JMAX
      real*8 FAC,epssv
      parameter(nlayer=20,JMAX=50,FAC=1.d-2,epssv=1.d-5)
      real*8 rb(0:nlayer), phasi,tobo,xacc,ss,dr,a,b,eps,
     &   epsnex,rtp,rstart,rtbis,x1,x2,rmin
      integer iflso(nlayer),iswi
      external phasi
      common /modr/ rb /modls/ iflso
      common/modlay/nlay 
      common /whichv/ layer,iaob,iswi /acc/ eps
      save

      rmin = rb(1)*FAC         ! Smallest starting radius
      rtp = rmin
      epsnex=-dlog(epssv)

c            Search for layer with turning point
c            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      layer = 1
      iswi = 1
 10   if (layer.le.nlay) then
        iaob = iflso(layer)
        tobo  = phasi(max(rb(layer-1),rmin)) 
c       print *,'an Schichtanfang:', tobo
        if (tobo.gt.0.d0) then
           tobo  = phasi(rb(layer))
c          print *,'an Schichtende:', tobo
           if (tobo.lt.0.d0)then
c    Search for turning point position inside layer
c    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              x1 = max(rb(layer-1),rmin)
              x2 = rb(layer)
              xacc = (x2-x1)*FAC
              rtp = rtbis(phasi,x1,x2,xacc)
              goto 12
           endif
        else if (tobo.lt.0.d0) then
           if (layer.gt.1) then
             layer=layer-1
           else 
             rstart=rmin
             nla=1
             return
           endif
           iaob = iflso(layer)
           rtp = rb(layer)
           goto 12
        endif
        layer = layer+1
      goto 10
      endif

      layer=nlay
      rtp=rb(layer)


c                  Search starting radius
c                  ~~~~~~~~~~~~~~~~~~~~~~

 12   iswi = 0
c     print *, 'Turning point,layer=',rtp,layer
      a = max(rb(layer-1),rmin)
      b = rtp
      ss = 0.d0
c exponent integration should succeed now
 100  dr=(b-a)/dble(JMAX)
      do 30 j=1,JMAX
          rstart = b - j*dr
          ss = ss + phasi(rstart)*dr
          if(ss.ge.epsnex) then
             nla = layer
c            print *, 'Startradius:',rstart
             return
          endif
  30  continue
c     print *, 'Nochmal runter...'
      if (layer.ge.2) then
         layer = layer - 1
         iaob = iflso(layer)
c     print *, 'Nochmal runter...: layer iaob', layer,iaob
         a = max(rb(layer-1),rmin)
         b = rb(layer)
         goto 100
      else
         rstart=rmin
      endif

      end


      real*8 function phasi(r1)
c
c Calculates the integral of the phase in WKBJ-approximation from
c r1 to r2 analytically by setting velocity constant.
c
      integer nlayer,ly,nlay,iaob,l,iswi
      real*8 pi,rp2,rpi
      parameter(nlayer=20,pi=3.1415926535897932d0
     &  ,rpi=1.d0/Pi,rp2=0.5d0/pi)
      real*8 vpv(nlayer,4),vph(nlayer,4),vsv(nlayer,4),vsh(nlayer,4)
     &   ,eta(nlayer,4),rho(nlayer,4),rb(0:nlayer),qm(nlayer),qk(nlayer)
     & ,r1,x,om,elpl1,vel,e43,vela,velb,lnom
      double complex zom
      common /modanis/ rho,vpv,vph,vsv,vsh,eta
      common /modq/  qm,qk 
      common /modr/ rb
      common/modlay/nlay
      common /whichv/ ly,iaob,iswi
      common/loop/zom,om,elpl1,l
      save


      x = r1/rb(nlay)

c                  Phase velocities.
c                  ~~~~~~~~~~~~~~~~
      lnom=log(om*rp2)*rpi
      if (iaob.eq.1) then
       vela = vpv(ly,1)+x*(vpv(ly,2)+x*(vpv(ly,3)+x*vpv(ly,4)))
       velb = vsv(ly,1)+x*(vsv(ly,2)+x*(vsv(ly,3)+x*vsv(ly,4)))
       e43=4.d0*(velb/vela)**2 /3.d0
       vel=vela*(1.d0+lnom*((1.-e43)*qk(ly) + e43*qm(ly)))
      else if (iaob.eq.0) then
       vel = vsv(ly,1)+x*(vsv(ly,2)+x*(vsv(ly,3)+x*vsv(ly,4)))
       vel=  vel*(1.d0+lnom*qm(ly))
      else
       stop 'PHASE: velocity not determined !!!!'
      endif

      if (iswi.eq.1) then
        phasi = elpl1 - (om*r1/vel)**2
c       print *,  elpl1 ,om,r1,vel
      else if (iswi .eq. 0) then
        phasi = sqrt(elpl1 - (om*r1/vel)**2)/r1
      else
        stop 'PHASE: no phase calculated!!'
      endif

      end


      FUNCTION rtbis(func,x1,x2,xacc)
c Changed, so that approaching is done from the side where func(x)
c remains positive. See in-line comments after changed code.
      INTEGER JMAX
      REAL*8 rtbis,x1,x2,xacc,func
      EXTERNAL func
      PARAMETER (JMAX=40)
      INTEGER j
      REAL*8 dx,f,fmid,xmid
      fmid=func(x2)
      f=func(x1)
      if(f*fmid.ge.0.d0) then
        print *, 'x1,x2,f(x1),f(x2)=',x1,x2,func(x1),func(x2)
        pause 'root must be bracketed in rtbis'
      endif
      if(f.gt.0.d0)then
c...changed .lt. in .gt. 
        rtbis=x1
        dx=x2-x1
      else
        rtbis=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*.5d0
        xmid=rtbis+dx
        fmid=func(xmid)
        if(fmid.ge.0.d0)rtbis=xmid
c...changed .le. in .ge.
        if(abs(dx).lt.xacc .or. fmid.eq.0.d0) return
11    continue
      pause 'too many bisections in rtbis'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software "*3.

