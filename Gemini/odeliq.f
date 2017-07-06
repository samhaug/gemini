      subroutine odeliq(r,y,dydr)
c.....                           .......................................
c Calls : celpar
c                                                                      .
c             Evaluate ODE-system in a liquid layer.                   .
c             dll1 = l*(l+1)                                           .
c.......................................................................

	integer imx,ix
      parameter(imx=11,ix=96)
      double complex y(2),dydr(2),z11(imx,0:ix),z12(imx,0:ix),
     &  z21(imx,0:ix),zepa,zepc,zepf,zepl,zepn,zom,zomro
	real*8 ro,om,dll1,r
	integer memo,meseq,mepkt,ml,mc,l
      common/paramemo/memo(imx,0:ix,2)
      common/memopnt/meseq,mepkt
      common/ep/zepa,zepc,zepf,zepl,zepn,ro
      common/loop/zom,om,dll1,l
	save
c
      if(memo(meseq,mepkt,1).eq.0)then
        ml=meseq
        mc=mepkt
        call  cpanis(r)
        zomro=zom*zom*ro
        z11(ml,mc) = -2.d0/r
        z12(ml,mc) =  1.d0/zepc - dll1/(zomro*r*r)
        z21(ml,mc) = -zomro
      else
        ml=memo(meseq,mepkt,1)
        mc=memo(meseq,mepkt,2)
      endif
c
      dydr(1) = z11(ml,mc)*y(1) + z12(ml,mc)*y(2)
      dydr(2) = z21(ml,mc)*y(1)
c
      return
      end
