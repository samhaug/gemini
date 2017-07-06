      subroutine aasmst
c.....                 .................................................
c                                                                      .
c              Address Array for System Matrix STorage                 .
c  Tags the subdivisions of an integration step where the system       .
c  matrix of the SODE will be repeatedly calculated by the             .
c  Bulirsch-Stoer algorithm.                                           .
c.......................................................................
c
	integer imax,i,j,k,l,jfound,jleft,jright
      parameter(imax=11)
      integer nseq(imax),memo(imax,0:96,2)
      common /paramemo/ memo
	save 
      data nseq /2,4,6,8,12,16,24,32,48,64,96/
c
      do 1 i=2,imax
         do 2 j=0,nseq(i)
           k=1
           jfound=0
 3         if (k.lt.i .and. jfound.eq.0) then
             l=0
             jleft=j*nseq(k)
 4           if (l.le.nseq(k) .and. jfound.eq.0) then
               jright=l*nseq(i)
               if(jleft.eq.jright) then
                 memo(i,j,1)=k
                 memo(i,j,2)=l
                 jfound=1
               else
                 memo(i,j,1)=0
               endif
               l=l+1
	         goto 4
             endif
             k=k+1
             goto 3
           endif
2        continue
1     continue
c
      return
      end
