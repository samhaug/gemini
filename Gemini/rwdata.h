c
c Declaration for the write-block
c----------------------------------------------------------------------------
      complex basout(4,4,0:lmx),bsoutt(2,2,0:lmx)
	real  rsrc,rrecv,ros,ror,vpverts,vphoris,vsverts,vshoris,ettas,vpvertr,
     &   vphorir,vsvertr,vshorir,ettar,qms,qmr,qks,qkr,tlen,tau
      integer iflqr,iflqs,nli1,nli2,nli3,nb1,nb2,nb3,ldel,nflow,nfhigh,l1,l2
	character cmodel*10,cmotion*7
c
c Include the following lines to read/write the header
c>    OPEN(iunit,file=.....,form='unformatted',status='old')
c>    read/write(iunit) cmotion,cmodel,
c>   &    rsrc,ros,vpverts,vsverts,vphoris,vshoris,ettas,iflqs,qms,qks,
c>   &    rrecv,ror,vpvertr,vsvertr,vphorir,vshorir,ettar,iflqr,qmr,qkr,
c>   &    nli1,nli2,nli3,nb1,nb2,nb3,tlen,tau,nflow,nfhigh,l1,l2,ldel
      
