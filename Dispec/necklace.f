      subroutine necklace(lat1,long1,delta,phi,nfg,rlat,rlong)
c 
c Calculates nfg equal spaced receiver points on the meridian from
c (clat1,long1) with azimuth phi. Angles are give in degree!
c
      integer nfg,i
      real clat1,long1,lat1,rlat(*),rlong(*),delta,phi,ddelta,del,
     &   pi,tete,ctete,stete,cgam,gamma,sgam,d2r,r2d
      parameter(pi=3.1415926535897932,d2r=pi/180.,r2d=1./d2r)


c           Transform to radians 

      clat1 = Pi/2. - lat1*d2r
      delta = delta*d2r
      phi   = phi*d2r

c Subdivision of the arc in equal parts
	ddelta = delta/float(nfg)

c             Calculate coordinates of receivers
c             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	do 100 i=1,nfg

	 del = float(i)*ddelta

c   Cosine theorem on the sphere, applied to the spherical triangle
c               north pole -- Point 1 -- receiver

      ctete  = cos(clat1)*cos(del)-sin(clat1)*sin(del)*cos(phi)

      tete = acos(ctete)
      stete=sqrt((1.-ctete)*(1.+ctete))

c   Latitude of receiver 
	rlat(i) = 90. - tete*r2d

	if (sin(clat1).ne.0..and.stete.ne.0.) then
       cgam=(cos(del)-cos(clat1)*ctete)/(sin(clat1)*stete)
	 sgam = sin(del)*sin(phi)/stete
	else
	 print *, 'Division by Zero in <necklace>!!'
	 stop
	endif
      if(cgam.gt.1.) cgam=1.
      if(cgam.lt.-1.) cgam=-1.
      gamma=acos(cgam)
c
c***Longitude of receiver 
      if(sgam.ge.0.) then
	  rlong(i)=long1+gamma*r2d
	else
        rlong(i)=long1-gamma*r2d
	endif
c
 100  continue
c
      return
      end
