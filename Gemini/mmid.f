      subroutine mmid(y,dydx,nvar,xs,htot,nstep,yout,derivs)
c
c             Modified MIDpoint step
c
c Taken from:
c   Press et al., 'NUMERICAL RECIPES, The Art of Scientific Computing',
c   1986, Cambridge University Press.
c
c
      integer nmax
      parameter (nmax=5)
	external derivs
      double complex y(nvar),dydx(nvar),yout(nvar),ym(nmax),yn(nmax),
     &     swap
      real*8 xs,htot,h,h2,x
      integer meseq,mepkt,nvar,nstep,i,n
      common/memopnt/meseq,mepkt
	save

      h=htot/nstep
      do 11 i=1,nvar
      ym(i)=y(i)
      yn(i)=y(i)+h*dydx(i)
11    continue
      x=xs+h
      mepkt = mepkt + 1
      call derivs(x,yn,yout)
      h2=2.d0*h
      do 13 n=2,nstep
      do 12 i=1,nvar
        swap=ym(i)+h2*yout(i)
        ym(i)=yn(i)
        yn(i)=swap
12    continue
      x=x+h
      mepkt = mepkt + 1
      call derivs(x,yn,yout)
13    continue
      do 14 i=1,nvar
       yout(i)=0.5d0*(ym(i)+yn(i)+h*yout(i))
14    continue

      return
      end
