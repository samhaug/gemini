      subroutine gtllim(f,lmn,lmx,window)
c
c          Computes minimum and maximum degree lmn and lmx
c          for given frequency f (in mHz) by
c          linear interpolation. 
c
c Author:  J.R. Dalkolmo, May 1997     

      integer nst
      parameter(nst=30)
      real om(nst),ellmn(nst),ellmx(nst),f
      integer lmn,lmx,ifirst,n,i,klo,khi,k
      character window*(*)
      data ifirst/1/
      save


c            Open window file and get limits.
c            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if(ifirst.ne.0)then
       ifirst=0
       open(1,file=window,status='old')
       read(1,*) n
       do 10 i=1,n
 10    read(1,*) om(i),ellmn(i),ellmx(i)
       close(1)
      endif

c            Perform bisection to get frequency interval.
c            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      klo=1
      khi=n
 1    if(khi-klo.gt.1)then
        k=(khi+klo)/2
        if(om(k).gt.f)then
         khi=k
        else
         klo=k
        endif
        goto 1
      endif

c            Linear interpolation to get l-limits.
c            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      lmn=int(  (ellmn(khi)*(f-om(klo))-ellmn(klo)*(f-om(khi)))/
     &        (om(khi)-om(klo))  )
      lmx=int(  (ellmx(khi)*(f-om(klo))-ellmx(klo)*(f-om(khi)))/
     &        (om(khi)-om(klo))  )

c                  Adjust if nonsens.
c                  ~~~~~~~~~~~~~~~~~~

      if (lmx.lt.lmn.or.lmx.lt.0) then
        lmn=0
        lmx=-1
        return
      endif

      end
