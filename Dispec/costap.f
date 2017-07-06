	subroutine costap(ntend)
c
c Calculates a cosine taper of length 'ntend'
c
	integer ntpmx,ntend,l
	real pi
	parameter(pi=3.141592,ntpmx=200)
	real tapend(ntpmx)
	common/tap/tapend
c
      if (ntend.gt.ntpmx) then
          print *, 'Array size not sufficient for taper-length!!'
          print *, '**** Change bounds of array <tapend> !! *****'
          stop ' Stopping in <costap>!!'
      endif
c
      if(ntend.gt.0) then
       do 10 l=1,ntend
 10    tapend(l)=cos(pi/float(2*ntend)*float(l))**2
      else
       print *, 'Applying no taper...'
       do 20 l=1,ntpmx
 20    tapend(l) = 1.
      endif
c
	end
