      subroutine odemin(r,m,dmdr)
c.....                           .......................................
c Calls ---> cpanis
c                                                                      .
c   Evaluate ODE-system for minors in a solid region.                  .
c   A minor is a combination of basic solutions, i.e. a                .
c   sub-determinant. In this case the minors are built of              .
c   two basic solutions from spheroidal motion (-->subroutine ODEPCO)  .
c                                                                      .
c          m_ij = y_1i*y_2j - y_1j*y_2i                                .
c                                                                      .
c          m_12 = m(1)        m_23 = m(4)                              .
c          m_13 = m(2)        m_24 = m(5)                              .
c          m_14 = m(3)        m_34 = m(6) = -m(1)/(l*(l+1))            .
c                                                                      .
c   m(i)      .... minor                                               .
c   r         .... radius                                              .
c   dmdr(i)    ... derivative of minor                                 .
c   zepa,zepc, ... elastic parameters in a transversal-isotropic       .
c                  medium (complex because of anelasticity)            .
c   ro        .... density                                             .
c   om         ... frequency                                           .
c   zom      ..... double complex frequency (om-i*sigma)               .
c   dll1      .... abbrev. of l*(l+1)                                  .
c   l          ... degree of Spherical surface harmonics               .
c                                                                      .
c.......................................................................
c
      implicit real*8(a-h,o-y),double complex(z)
	integer imx,ix
      parameter(imx=11,ix=96)
      double complex m(5),dmdr(5),zepa,zepc,zepf,zepl,zepn,zom,
     &       z11(imx,0:ix),z12(imx,0:ix),z13(imx,0:ix),z14(imx,0:ix),
     &       z22(imx,0:ix),z23(imx,0:ix),z24(imx,0:ix),
     &       z31(imx,0:ix),z32(imx,0:ix),z33(imx,0:ix),
     &       z42(imx,0:ix),z51(imx,0:ix),z55(imx,0:ix)
      real*8 ro,om,dll1,rr,rr2
	integer l,memo,meseq,mepkt,ml,mc
      common /paramemo/ memo(imx,0:ix,2)
      common /memopnt/  meseq,mepkt
      common /ep/       zepa,zepc,zepf,zepl,zepn,ro
      common /loop/     zom,om,dll1,l
	save

      if(memo(meseq,mepkt,1).eq.0)then
        ml = meseq
        mc = mepkt
	  rr = 1.d0/r
        rr2 = rr*rr
        call cpanis(r)
        zomro = zom*zom*ro
	  zrpc = 1.d0/zepc
        z1 = zepa - zepf*zepf*zrpc - zepn
        z2 = 2.d0*zepf*zrpc
        z3 = 4.d0*z1*rr2
c Evaluate system matrix of SODE
        z11(ml,mc) = -2.d0*rr
        z12(ml,mc) = -2.d0*dll1*z1*rr2
        z13(ml,mc) =  dll1*rr
        z14(ml,mc) = -dll1*zepf*zrpc*rr
c       z15=0.d0
c       z21=0.d0
        z22(ml,mc) =  (1.d0-z2)*rr
        z23(ml,mc) =  1.d0/zepl
        z24(ml,mc) =  zrpc
c       z25=0.d0
        z31(ml,mc) = -z2*rr
        z32(ml,mc) = -zomro + ( dll1*(z1+zepn) - 2.d0*zepn )*rr2
        z33(ml,mc) = -(3.d0+z2)*rr
c       z34=0.d0
c       z35=z24
c       z41=-z11
        z42(ml,mc) = -zomro+z3
c       z43=0.d0
c       z44=-z22
c       z45=z23
        z51(ml,mc) =  z3
c       z52=0.d0
c       z53=z42
c       z54=z32
        z55(ml,mc) =  (z2-5.d0)*rr
      else
        ml=memo(meseq,mepkt,1)
        mc=memo(meseq,mepkt,2)
      endif
c
      dmdr(1) =  z11(ml,mc)*m(1) + z12(ml,mc)*m(2) + z13(ml,mc)*m(3)
     &           +z14(ml,mc)*m(4)
      dmdr(2) =  z22(ml,mc)*m(2) + z23(ml,mc)*m(3) + z24(ml,mc)*m(4)
      dmdr(3) =  z31(ml,mc)*m(1) + z32(ml,mc)*m(2) + z33(ml,mc)*m(3)
     &           +z24(ml,mc)*m(5)
      dmdr(4) = -z11(ml,mc)*m(1) + z42(ml,mc)*m(2)
     &          -z22(ml,mc)*m(4) + z23(ml,mc)*m(5)
      dmdr(5) =  z51(ml,mc)*m(1) + z42(ml,mc)*m(3)
     &           +z32(ml,mc)*m(4) + z55(ml,mc)*m(5)
c
      end
