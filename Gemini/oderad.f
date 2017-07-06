      subroutine oderad(r,y,dydr)
c.....                           .......................................
c Calls ---> cpanis
c                                                                      .
c Evaluate ODE-system for radial motion, i.e. degree 'l' of            .
c the sperical harmonics equals zero: l=0                              .
c.......................................................................
      integer imx,ix
      parameter(imx=11,ix=96)
      double complex y(2),dydr(2),zepa,zepc,zepf,zepl,zepn,zom,zrpc,z1,
     &   z11(imx,0:ix),z12(imx,0:ix),z21(imx,0:ix),z22(imx,0:ix)
	real*8 om,dll1,ro,r,rr
	integer memo,meseq,mepkt,ml,mc,l
      common/paramemo/memo(imx,0:ix,2)
      common/memopnt/meseq,mepkt
      common/ep/zepa,zepc,zepf,zepl,zepn,ro
      common/loop/zom,om,dll1,l
	save
c
      if(memo(meseq,mepkt,1).eq.0)then
        ml = meseq
        mc = mepkt
	  rr = 1.d0/r
        call  cpanis(r)
	  zrpc = 1.d0/zepc
        z1   = zepf*zrpc
        z11(ml,mc) = -2.d0*z1*rr
        z12(ml,mc) =  zrpc
        z21(ml,mc) = -zom*zom*ro + 4.d0*( zepa-z1*zepf-zepn )*rr*rr
        z22(ml,mc) =  2.d0*rr*(z1-1.d0)
      else
        ml=memo(meseq,mepkt,1)
        mc=memo(meseq,mepkt,2)
      endif
c
      dydr(1) = z11(ml,mc)*y(1) + z12(ml,mc)*y(2)
      dydr(2) = z21(ml,mc)*y(1) + z22(ml,mc)*y(2)
c
      end
