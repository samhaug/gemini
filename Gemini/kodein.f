      subroutine kodein (ystart,nvar,x1,x2,h1,hmin,nok,nbad,derivs)
c
c  Complex Kolmo-version of bulirsch-stoer driver
c
c  Bulirsch-Stoer driver with adaptive step size control. integrate
c  the nvar starting values ystart from x1 to x2,
c  h1 should be set as a guessed first stepsize, hmin as the
c  minimum allowed step size (can be zero). on output nok and
c  nbad are the number of good and bad (but retried and fixed)
c  steps taken, and ystart is replaced by values at the end of
c  the integration interval.
c
c Based on:
c   Press et al., 'NUMERICAL RECIPES, The Art of Scientific Computing',
c   1986, Cambridge University Press.
c--------------------------------------------------------------
	real*8 tiny
	integer maxstp,nvmax,nmax
      parameter(maxstp=20000,nmax=5,nvmax=6,tiny=1.d-10)
c
      external derivs
      double complex ystart(nvar),dyst(nvmax),y(nmax),dydx(nmax)
	real*8 x,x1,x2,h,h1,hmin,yscal,hdid,hnext
	integer nok,nbad,i,nvar,nstp,meseq,mepkt,ifcalder
      common/memopnt/meseq,mepkt
      common  /calder/ ifcalder  /derivatives/ dyst
	save

c Meseq,mepkt denotes the line and column respectively, of the
c memo array which contains the information about points in
c a preceding subdivision where the system matrix of the SODE was
c already calculated; so some unnecessary calculations can be 
c avoided.

      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      do 11 i=1,nvar
 11   y(i)=ystart(i)

      do 16 nstp=1,maxstp
c..........................integration loop

c If step can overshoot end, cut down stepsize.
       if((x+h-x2)*(x+h-x1).gt.0.d0) h=x2-x
       if(dabs(h).lt.tiny) then
        if (ifcalder.eq.1) then      ! Calculate derivatives for strain
            meseq=1
            mepkt=0
            call derivs(x,y,dyst)
        endif
        return
       endif

       meseq=1
       mepkt=0
       call derivs(x,y,dydx)

       call bsstep(y,dydx,nvar,x,h,yscal,hdid,hnext,derivs)
       if(hdid.eq.h) then
        nok=nok+1
       else
        nbad=nbad+1
       endif
c
c Are we done ?
       if((x-x2)*(x2-x1).ge.0.d0) then
          if (ifcalder.eq.1) call derivs(x,y,dyst)
          do 14 i=1,nvar
 14       ystart(i)=y(i)
          return
       endif

       if(dabs(hnext).lt.hmin) then
	  print *,'nvar, h, x = ', nvar,h,x
	  pause ' stepsize smaller than min'
	 endif
       h=hnext
 16   continue

      pause ' too many steps '
      return
      end

