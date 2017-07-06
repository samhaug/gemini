      subroutine gtmanis(cmodel,fref)
c------------------------------------------------------
c              G T M A N I S 
c
c Reads tranversal isotropic Earth model parameters described by polynomials
c from file.
c See comments at the read-statements for the meaning of the parameter.
c
c------------------------------------------------------
      integer nlayer
      parameter(nlayer=20)
      real*8 qm(nlayer),qk(nlayer),
     &     rho(nlayer,4),rb(0:nlayer),rboc,fref
      real*8 vpv(nlayer,4),vph(nlayer,4),vsv(nlayer,4),vsh(nlayer,4),
     &       eta(nlayer,4)
      integer iflso(nlayer),nco(nlayer),i,j,nlay,iunit,nloc,ifanis
      character text*80,cnlay*2,form*11,cmodel*(*)
      common /modq/  qm,qk /modr/ rb
      common /modanis/ rho,vpv,vph,vsv,vsh,eta  /tranis/ ifanis
      common /modls/ iflso  /modlay/ nlay /core/ rboc
      save

      do i=1,nlayer
        qm(i)=0.d0
        qk(i)=0.d0
        rb(i)=0.d0
        iflso(i)=0
        nco(i)=0
        do j=1,4
          rho(i,j)=0.d0
          vpv(i,j)=0.d0
          vph(i,j)=0.d0
          vsv(i,j)=0.d0
          vsh(i,j)=0.d0
          eta(i,j)=0.d0
        enddo
      enddo

      iunit = 1
      open(iunit,file=cmodel,status='old')

c          Read header with comments and print them on stdout
c          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  1   read(iunit,'(a72)') text
      if (text(1:1).eq.'#') then
c       print '(a72)', text
        goto 1
      endif
      backspace iunit


      read(iunit,'(i2)') nlay                ! Number of layers

      write(cnlay,'(i2)') nlay               ! 
      form='('//cnlay//'i2)'                 ! Number of polynomal
      read(iunit,form) (nco(i),i=1,nlay)     ! coefficients for each layer 
      
      read(iunit,*) fref               ! reference frequency of Qs in Hertz
      read(iunit,*) ifanis             ! Transversal isotropic? 1=y, else=n
      read(iunit,'(1x/1x/)')


c  Reading radii, density, P-velo-vert, P-velo-hori, S-velo-vert, 
c          S-velo-hori, Qmu, Qkappa, eta

      do i = 1, nlay

        read(iunit,*) rb(i-1),rho(i,1),vpv(i,1),vph(i,1),vsv(i,1),vsh(i,1),
     &                qm(i),qk(i),eta(i,1)
c       write(88,*) rb(i-1),rho(i,1),vpv(i,1),vph(i,1),vsv(i,1),vsh(i,1),
c    &                qm(i),qk(i),eta(i,1)

        qk(i)=1./qk(i)                 ! reciprocal quality factors
        if (qm(i).le.0.) then          !
            qm(i) = 0.                 !
            iflso(i)=1                 !
            nloc = i                   ! layer number of outer core
        else                           ! 
	      qm(i)=1./qm(i)             !
            iflso(i)=0                 !
        endif                          ! 

        do j = 2, nco(i)               ! polynomal coefficients
        read(iunit,*) rho(i,j),vpv(i,j),vph(i,j),vsv(i,j),vsh(i,j),eta(i,j)
c       write(88,*) rho(i,j),vpv(i,j),vph(i,j),vsv(i,j),vsh(i,j),eta(i,j)
        enddo
        read(iunit,'(1x)')
c       write(88,'(1x)')

      enddo


      read(iunit,*) rb(nlay)           ! Get Earth radius

c  Check for ocean.
      if(iflso(nlay).eq.1) stop '~~~~ I don`t like oceans. ~~~~'
      rboc = rb(nloc)                  ! Radius of outer core

      close(iunit)
      return

      end
