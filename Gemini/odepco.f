      subroutine odepco(r,y,dydr)
c.....                           .......................................
c Calls ---> celpar
c                                                                      .
c   Evaluate ODE-system for EXPANSION COEFFICIENTS U, V and R, S       .
c   of displacement and tension (spheroidal motion), respectively.     .
c                                                                      .
c     y(1) --> U --> radial component of displacement-vector           .
c     y(2) --> R --> radial component of tension-vector                .
c     y(3) --> V --> tangential component of displacement-vector       .
c     y(4) --> T --> tangential component of tension-vector            .
c                                                                      .
c   r .... radius                                                      .
c   dydr(i) ... derivative of y(i) with respect to r                   .
c   zepa,zepc,... elastic parameters in a transversal-isotropic        .
c                 medium (complex because of anelasticity).            .
c   ro .... density,   om ... frequency                                .
c   zom ..... complex frequency (om-i*sigma)                           .  
c   l ... degree of Spherical surface harmonics                        .
c   dll1 ....  l*(l+1)                                                 .
c.......................................................................
c
	integer imx,ix
      parameter(imx=11,ix=96)
      double complex y(4),dydr(4),z11(imx,0:ix),z12(imx,0:ix),
     &   z13(imx,0:ix),z21(imx,0:ix),z22(imx,0:ix),z23(imx,0:ix),
     &   z24(imx,0:ix),z31(imx,0:ix),z34(imx,0:ix),z41(imx,0:ix),
     &   z42(imx,0:ix),z43(imx,0:ix),z44(imx,0:ix),zepa,zepc,zepf,
     &   zepl,zepn,zom,z1,z2,zrpc,zomro
	real*8 om,ro,rr,rr2,dll1,r
	integer l,memo,meseq,mepkt,ml,mc
      common/paramemo/memo(imx,0:ix,2)
      common/memopnt/meseq,mepkt
      common/ep/zepa,zepc,zepf,zepl,zepn,ro
      common/loop/zom,om,dll1,l
	save
c
      if(memo(meseq,mepkt,1).eq.0)then
        ml    = meseq
        mc    = mepkt
	  rr    = 1.d0/r
        rr2   = rr*rr
        call cpanis(r)
	  zrpc  = 1.d0/zepc
        z1    = zepa - zepf*zepf*zrpc - zepn
        z2    = zepf*rr*zrpc 
        zomro = zom*zom*ro
c Evaluate system matrix of SODE
        z11(ml,mc) = -2.d0*z2
        z12(ml,mc) =  zrpc
        z13(ml,mc) =  dll1*z2
c       z14=0.d0
        z21(ml,mc) = -zomro + 4.d0*z1*rr2
        z22(ml,mc) = -z11(ml,mc) - 2.d0*rr
        z23(ml,mc) = -2.d0*dll1*z1*rr2
        z24(ml,mc) =  dll1*rr
        z31(ml,mc) = -rr
c       z32=0.d0
c       z33=-z31
        z34(ml,mc) =  1.d0/zepl
        z41(ml,mc) = -2.d0*z1*rr2
        z42(ml,mc) = -z2
        z43(ml,mc) = -zomro + ( dll1*(z1+zepn) - 2.d0*zepn )*rr2
        z44(ml,mc) = -3.d0*rr
      else
        ml=memo(meseq,mepkt,1)
        mc=memo(meseq,mepkt,2)
      endif
c
      dydr(1) = z11(ml,mc)*y(1) + z12(ml,mc)*y(2) + z13(ml,mc)*y(3)
      dydr(2) = z21(ml,mc)*y(1) + z22(ml,mc)*y(2) + z23(ml,mc)*y(3)
     &          +z24(ml,mc)*y(4)
      dydr(3) = z31(ml,mc)*(y(1)-y(3)) + z34(ml,mc)*y(4)
      dydr(4) = z41(ml,mc)*y(1) + z42(ml,mc)*y(2) + z43(ml,mc)*y(3)
     &          +z44(ml,mc)*y(4)
c
      end
