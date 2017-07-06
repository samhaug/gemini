      subroutine bsstep(y,dydx,nv,x,htry,yscal,hdid,hnext,derivs)
c
c             Bulirsch-Stoer STEP
c
c Taken and slightly modified from:
c   Press et al., 'NUMERICAL RECIPES, The Art of Scientific Computing',
c   1986, Cambridge University Press.
c
	real*8 shrink,grow
 	integer nmax,imax,nuse
      parameter(nmax=5,imax=11,nuse=7,shrink=.95d0,grow=1.2d0)
      external derivs
      double complex y(nv),dydx(nv),yerr(nmax),
     1                   ysav(nmax),dysav(nmax),yseq(nmax)
	real*8 h,htry,x,xsav,eps,yscal,hdid,hnext,xest,errmax
      integer nseq(imax),meseq,mepkt,nv,i,j
      common /memopnt/ meseq,mepkt
      common /acc/ eps
	save
      data nseq /2,4,6,8,12,16,24,32,48,64,96/

      h=htry
      xsav=x
      do 11 i=1,nv
      ysav(i)=y(i)
      dysav(i)=dydx(i)
11    continue

1     do 10 i=1,imax
c.....................loop over subdivisions 

      call mmid(ysav,dysav,nv,xsav,h,nseq(i),yseq,derivs)
      xest=(h/nseq(i))**2
      call rzextr(i,xest,yseq,y,yerr,nv,nuse)

c Guard against spurious early convergence.
      if(i.gt.3)then
         errmax=0.d0
         do 12 j=1,nv
	     yscal=1.d0/zabs(y(j))
           errmax=max(errmax,zabs(yerr(j))*yscal)
12       continue
         errmax=errmax/eps

c Cancel and new step if accuracy is sufficient.
         if(errmax.lt.1.d0) then
           x=x+h
           hdid=h
           if(i.eq.nuse)then
             hnext=h*shrink
           else if(i.eq.nuse-1)then
             hnext=h*grow
           else
             hnext=(h*nseq(nuse-1))/nseq(i)
           endif
           return
         endif

      endif
      meseq = meseq + 1
      mepkt = 0
10    continue

      h=0.25d0*h/2**((imax-nuse)/2)
      write(6,7) ' step size diminished to ',h,' at point ',x
   7  format(a25,f10.2,a10,f10.2)
      if(x+h.eq.x)pause 'step size underflow.'
      meseq = 1
      goto 1
      end
