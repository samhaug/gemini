      subroutine momfor(r0,mm,fom,s,cfom)
c
c
c Calulates source terms of a moment tensor or a single force
c for spheroidal motion (indices 1 to 4)
c and toroidal motion (indices 5 and 6)
c      r0 ...      Source radius
c      fom ...      contains moment/force
c      s ...            Source terms
c
c  Used definition of spherical harmonics
c      y(l,m)= ((-1)**m)*Norm(l,m)*p(l,m)*exp(i*m*phi)
c  Normalization factor of source terms already computed in main
c  program DISPEC
c
      integer mm,i,j
      complex s(6,-mm:mm),za,zc,zf,zl,zn,zepc,zepl
      real fom(6),rr2,rr3,rr2h,r0
      character cfom*1
c      ro0 ...      Density at source radius
c      za,zc,zf,
c      zl,zn ...      complex elastic constants there.
      common /sopar/ za,zc,zf,zl,zn
      save
c

      rr2=1./r0**2
      rr3=rr2/r0
      rr2h=0.5*rr2
      do i=-mm,mm
       do j = 1,6
        s(j,i) = 0.
       enddo
      enddo

      if(cfom.eq.'m')then
c *** Source is moment tensor  
        zepc=1./zc
        zepl=1./zl
         s(1,0)  =  rr2*fom(1)*zepc
         s(2,0)  =  rr3*(2.*zf*zepc*fom(1)-fom(2)-fom(3))
         s(4,0)  = -0.5*s(2,0)
         if (mm.ge.1) then
          s(3,-1) = -rr2h*zepl*cmplx(-fom(4),-fom(5))
          s(3,1)  = -rr2h*zepl*cmplx( fom(4),-fom(5))
          s(5,-1) =  rr2h*zepl*cmplx(-fom(5), fom(4))
          s(5,1)  =  rr2h*zepl*cmplx( fom(5), fom(4))
         endif
         if (mm.ge.2) then
          s(4,-2) =  rr3*cmplx(fom(3)-fom(2),-2.*fom(6))*0.25
          s(4,2)  =  conjg(s(4,-2))
          s(6,-2) =  rr3*cmplx(2.*fom(6),fom(3)-fom(2))*0.25
          s(6,2)  =  conjg(s(6,-2))
         endif
      elseif(cfom.eq.'f') then
c *** Source is single force   
         s(2,0)   = -rr2*fom(1)
         if (mm.ge.1) then
              s(4,-1) = -rr2h*cmplx(fom(2),fom(3))
              s(4,1)  = -rr2h*cmplx(-fom(2),fom(3))
              s(6,-1) = -rr2h*cmplx(-fom(3),fom(2))
              s(6,1)  = -rr2h*cmplx(fom(3),fom(2))
         endif
      endif
      end            
