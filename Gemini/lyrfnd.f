      subroutine lyrfnd(r,laysrc,ifliq)
c
c Subroutine to find out layer number of given radius r. rb() are
c the outer border radii of the layers in the used earth model. iflso()
c indicates if the layer is liquid(1) or solid(0). nlay is the 
c total number of layers in the model
c
      integer nlayer
      parameter(nlayer=20)
	real*8 rb(0:nlayer),r
	integer iflso(nlayer),i,nlay,laysrc,ifliq
      common/modr/rb/modls/iflso/modlay/nlay
	save
c 
c If r less then zero, it is depth.
      if(r.lt.0.d0) r=rb(nlay)+r
c
      i=1
 101  if(r.gt.rb(i)) then
        i=i+1
        goto 101
      endif
      laysrc=i
	if (laysrc.gt.nlay) stop 'LYRFND: Radius above surface!! ****'
      ifliq=iflso(laysrc)
c
      end
