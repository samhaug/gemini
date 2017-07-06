	subroutine stnorm(lmx,lm,mm,snorm)
c Calls ----> none.
c
c Makes normalization factors for source terms in  m o m f o r .
c This routine has to be adjusted to the definition of the spherical
c harmonics in use.
c
	real snorm(0:lmx,0:mm),pi
	integer lm,lmx,mm,l
	parameter(pi=3.141592)
c
	do 101 l=0,lm
	  snorm(l,0)=sqrt(float(2*l+1)/(4.*Pi))
	  if (mm.gt.0 .and. l.ge.1) 
     &      snorm(l,1)=snorm(l,0)/sqrt(float(l*(l+1)))
	  if (mm.gt.1 .and. l.ge.2) 
     &      snorm(l,2)=snorm(l,1)*sqrt(float((l+2)*(l-1)))
 101  continue
c
	end
