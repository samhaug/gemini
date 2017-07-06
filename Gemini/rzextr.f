      subroutine rzextr(iest,xest,yest,yz,dy,nv,nuse)
c-----------------------------------------------------------------
c                     R Z E X T R
c
c        Rational function extrapolation 
c
c Taken from:
c   Press et al., 'NUMERICAL RECIPES, The Art of Scientific Computing',
c   1986, Cambridge University Press.
c-----------------------------------------------------------------
	integer nmax,imax,ncol,j,nv,k,m1,iest,nuse
      parameter (imax=11,nmax=5,ncol=7)
      double complex yest(nv),yz(nv),dy(nv),d(nmax,ncol),v,c,b,b1,yy,ddy
      real*8 x(imax),fx(ncol),xest
	save

      x(iest)=xest
      if(iest.eq.1) then
      do 11 j=1,nv
        yz(j)=yest(j)
        d(j,1)=yest(j)
        dy(j)=yest(j)
11    continue
      else
      m1=min(iest,nuse)
      do 12 k=1,m1-1
        fx(k+1)=x(iest-k)/xest
12    continue

      do 14 j=1,nv
c...................extrapolation part
        yy=yest(j)
        v=d(j,1)
        c=yy
        d(j,1)=yy
        do 13 k=2,m1
          b1=fx(k)*v
          b=b1-c
          if(b.ne.0.d0) then
             b=(c-v)/b
             ddy=c*b
             c=b1*b
          else
             ddy=v
          endif
          v=d(j,k)
          d(j,k)=ddy
          yy=yy+ddy
13      continue
        dy(j)=ddy
        yz(j)=yy
c................
14    continue

      endif
      return
      end
