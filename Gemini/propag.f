      subroutine propag(rbeg,rend,zy,h1,derivs,nv)
c.....                                                    .............
c         Subroutine to propagate function vector through solid       .
c         and liquid layers stacked in arbitrary order.               .
c......................................................................
c
      integer nmax,nlayer
      real*8 hmin
      parameter(nmax=6,nlayer=20,hmin=1.d-3)
      double complex zy(nmax)
      real*8 rb(0:nlayer),h1(nlayer),rbeg,rend,x1,x2
      integer iflso(nlayer),iexit,nlinc,nvar,nok,nbad,nl,nv,ifliq
     &        ,ifcalder
      external odeliq,derivs
      character*1 csy
c       csy='M'        : minor-DGS
c       csy='C' or 'P' : potential coeff.-DGS
      common/layer/nl/modr/rb/modls/iflso 
      common /whichs/ csy
      common /calder/ ifcalder
      save

      iexit=0
      call lyrfnd(rbeg,nl,ifliq)
      x1=rbeg

c Which direction to go?
      nlinc=int(sign(1.d0,rend-rbeg))

      if (nlinc.gt.0) then
         x2=min(rend,rb(nl))
      else
         x2=max(rend,rb(nl-1))
      endif
      if(x2.eq.rend) iexit=1

c Integration loop.
 100  if(iflso(nl).eq.1) then
        nvar=2
        call kodein(zy,nvar,x1,x2,h1(nl),hmin,nok,nbad,odeliq)
      else
        nvar=nv
        call kodein(zy,nvar,x1,x2,h1(nl),hmin,nok,nbad,derivs)
      endif
c Rescale the upward solution (of GEMINI!) to prevent overflow at 
c high frequencies
      if (nlinc.gt.0 .and. ifcalder.eq.0) then
         call rescal(zy,nvar)
      endif  
      if(iexit.eq.1) return
      
c Continuation at layer boundary needs special treatment, if boundary
c separates liquid and solid material.
      if(iflso(nl)+iflso(nl+nlinc).eq.1) then
       if(iflso(nl).lt.iflso(nl+nlinc)) then
c Boundary solid -> liquid
         if(csy.eq.'M') then
          zy(1)=zy(3)
          zy(2)=zy(5)
         endif
       else
c Boundary liquid -> solid
         if(csy.eq.'M') then
          zy(4)=zy(2)
          zy(2)=zy(1)
          zy(1)=0.d0
          zy(3)=0.d0
          zy(5)=0.d0
         else if(csy.eq.'C' .or. csy.eq.'P') then
          zy(3)=1.d0
          zy(4)=0.d0
         endif
       endif
      endif
      nl=nl+nlinc
      
c Determine new integration interval.
      if(nlinc.gt.0) then
       x1=rb(nl-1)
       x2=min(rend,rb(nl))
      else
       x1=rb(nl)
       x2=max(rend,rb(nl-1))
      endif
      if(x2.eq.rend) iexit=1
      goto 100
      
      end
