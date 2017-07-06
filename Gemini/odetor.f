      subroutine odetor(r,y,dydr)
c.....                           .......................................
c Calls ---> cpanis
c                                                                      .
c             Evaluate ODE-system for toroidal motion.                 .
c                                                                      .
c.......................................................................

	integer imx,ix
      parameter(imx=11,ix=96)
      double complex y(2),dydr(2),zepa,zepc,zepf,zepl,zepn,zom,
     &   z11(imx,0:ix),z12(imx,0:ix),z21(imx,0:ix),z22(imx,0:ix)
	real*8 ro,rr,om,dll1,r
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
        call cpanis(r)
        z11(ml,mc) =  rr
        z12(ml,mc) =  1.d0/zepl
        z21(ml,mc) = -zom*zom*ro - zepn*rr*rr*(2.d0-dll1)
        z22(ml,mc) = -3.d0*rr
      else
        ml = memo(meseq,mepkt,1)
        mc = memo(meseq,mepkt,2)
      endif
c
      dydr(1) = z11(ml,mc)*y(1) + z12(ml,mc)*y(2)
      dydr(2) = z21(ml,mc)*y(1) + z22(ml,mc)*y(2)
c
      end
