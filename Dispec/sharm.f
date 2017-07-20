      subroutine sharm(lmx,lmax,mm,mmax,delta,phi,y,dydp,dydt)
c.....                                                       ...........
c                                                                      .
c               Computes Spherical Harmonics
c
c                    1       d Ylm         d Ylm
c        Ylm,    ---------*---------- ,  --------                     
c                sin(theta)  d Phi        d Theta
c
c  until lmax,mmax. Note that dYlm/dPhi=i*m*Ylm !                      .
c                                                                      .
c This version in  s i n g l e   p r e c i s i o n !!
c 27.9.99: Changed handling of poles 
c.......................................................................
c
      real snrm,faclm,pi,cdel,sdel,fact,r4pi,rsdel,d2r
      complex y(0:lmx,-mm:mm),dydt(0:lmx,-mm:mm),dydp(0:lmx,-mm:mm),emfi,zi
      real delta,phi
      integer l,m,lmax,mmax,mm,lmx,j,ievod
      logical atpoles
      parameter( Pi=3.1415926535897932, r4pi = 1./(4.*pi), 
     &           zi=(0.,1.), d2r  = Pi/180. )


      cdel = cos(delta*d2r)
      sdel = sqrt((1.-cdel)*(1.+cdel))
      if (sdel.gt.0.) then
         atpoles=.false.
      else 
         print *, '*** We are at one of the poles! ***'
         atpoles=.true.
      endif

c ------------------------------------------------------------------
c                   Compute Legendre-Polynoms
c ------------------------------------------------------------------
      do  m = 0,mmax         ! Loop over order m
c  Initialization 
         Y(m,m)  = 1.0
         fact = 1.0
         do  j = 1, m
             Y(m,m)  = -Y(m,m)*fact*sdel
             fact =  fact + 2.0
         enddo
         Y(m+1,m)    = cdel*(2*m+1)*Y(m,m)
c  Main recurrence
        do l = m+2, lmax
         Y(l,m) = (cdel*float(2*l-1)*Y(l-1,m)-(l+m-1)*Y(l-2,m)) 
     &            /float(l-m)
        enddo
      enddo

c----------------------------------------------------------------------
c        Compute derivatives with respect to theta and phi
c----------------------------------------------------------------------
      if (.not.atpoles) then
        rsdel=1./sdel
        do  m = 0,mmax         ! Loop over order m
         dYdt(m,m)   = m*cdel*Y(m,m)*rsdel
         dYdp(m,m)   = zi*m*Y(m,m)*rsdel
         do l = m+1, lmax
           dYdp(l,m) = zi*m*Y(l,m)*rsdel
           dYdt(l,m) = (l*cdel*Y(l,m)-(l+m)*Y(l-1,m)) *rsdel
         enddo
        enddo
      else
         ievod=sign(1.0,cdel)
         do  l=0,lmax
          dYdt(l,0) = 0.
          dYdt(l,1) = -0.5*float((l*(l+1))*ievod**l)
          dYdt(l,2) = 0.
          dYdp(l,0) = 0.
          dYdp(l,1) = zi*dYdt(l,1)*ievod
          dYdp(l,2) = 0.
         enddo
      endif

c            Normalization and multiplication with zonal functions
c            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      do 200 m = 0,mmax
         faclm = 1.0
         emfi  = cexp(zi*m*phi*d2r)
         do  j = 1, 2*m                  !
           faclm = faclm/float(j)        !  l=m -> Compute 1/(2*mm)!
         enddo                           ! 
         do l = m, lmax
             snrm = sqrt(float(2*l+1)*r4pi*faclm)
             Y(l,m)    = snrm*Y(l,m)*emfi
             dYdt(l,m) = snrm*dYdt(l,m)*emfi
             dYdp(l,m) = snrm*dYdp(l,m)*emfi
             Y(l,-m)    = conjg(Y(l,m))*(-1)**m
             dYdt(l,-m) = conjg(dYdt(l,m))*(-1)**m
             dYdp(l,-m)  = conjg(dYdp(l,m))*(-1)**m
             faclm = faclm*float(l+1-m)/float(l+1+m)
         enddo
 200  continue


      return
      end
