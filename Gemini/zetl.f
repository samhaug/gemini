      Subroutine   Zetl (zx2,l,z)
c
c Computes ratio of (complex) spherical Bessel-functions 
c    z(n):= x*j(n+1)/j(n)
c
c by calculation of the recurrence formula
c    z(n)=x**2/(2*n+3 -z(n+1)).
c (Takeuchi and Saito (1972), Methods in computional Physics, p.243)
c
c The recurrence formula is treated as a continued fraction
c and evaluated with the MODIFIED LENTZ'S METHOD. See Press, W.H. et
c al, 1992, Numerical Recipes, 2nd Edition, Cambridge, Chapter 5.2,
c p. 169 for further information.
c
c Author: J.R. Dalkolmo


	real*8 tiny,eps
	integer maxit,l,i
      parameter(maxit=5000,tiny=1.d-30)
	double complex zx2,z,zd,zc,zdelta
      common /acc/ eps
	save

c First iteration
      z      = tiny
      zd     = dcmplx(2*l+3)
      zc     = zd+zx2/z
      zd     = 1.d0/zd
      zdelta = zc*zd
      z      = z*zdelta

c Remaining iterations
      do 10 i=2,maxit
		zd=dcmplx(2*(l+i)+1)-zx2*zd
		if (zd.eq.0.d0) zd=tiny
		zc=dcmplx(2*(l+i)+1)-zx2/zc
		if (zc.eq.0.d0) zc=tiny
		zd=1.d0/zd
		zdelta=zc*zd
		z=z*zdelta
		if ((zabs(zdelta)-1.d0).lt.eps) then
c	  print *, zsqrt(zx2),l,i
		  return
		endif
 10   continue
      print *, 'zx2,l,eps,z,maxit,tiny=',zx2,l,eps,z,maxit,tiny
      pause '*** Iteration failed to converge in ZETL ! ***'

      end
