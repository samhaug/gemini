      subroutine  cpanis(r)
c.....                             .....................................
c                                                                      .
c   Calculates density and complex elastic parameters in a             .
c   medium for a given radius. In this case parameters for             .
c   an transversal isotropic earth model are calculated.               .
c                                                                      .
c   roh ........... density-array                                      .
c   zpa,zpc, ...... complex transv. isotr. elastic constants           .
c   zdmpm,zdmpk ... damping factor calculated in main program          .
c   nl ............ actual layer number                                .
c The conversion formulae are taken from Dziewonski & Anderson: PREM,  .
c 1980, Phys.Earth.Planet.Int.,25,297-356  (DA1980)                    .
c.......................................................................
c
      integer nlayer
      parameter(nlayer=20)
      double complex zpa,zpc,zpf,zpl,zpn,zdmpm(nlayer),zdmpk(nlayer),
     &    zl2m,zmu,zlm,zdma,zdmc,zdmf,zdml,zdmn
      real*8 ro,roi,rho(nlayer,4), rb(0:nlayer),x,r,vk,vmu,amu,akap,al2m
      real*8 epa,epc,epf,epl,epn,vpvert,vphori,vsvert,vshori,etta
      real*8 vpv(nlayer,4),vph(nlayer,4),vsv(nlayer,4),vsh(nlayer,4),
     &       eta(nlayer,4)
      integer nl,nlay,iflso(nlayer),ifanis
      common/ep/zpa,zpc,zpf,zpl,zpn,ro
      common/epi/zmu,zl2m,roi
      common /modr/ rb
      common /modanis/ rho,vpv,vph,vsv,vsh,eta
      common /lyrdmp/zdmpm,zdmpk
      common /layer/ nl /tranis/ ifanis
      common /modls/ iflso /modlay/ nlay
	save

c Density & velocities
      x = r/rb(nlay)
      ro= rho(nl,1)+x*(rho(nl,2)+x*(rho(nl,3)+x*rho(nl,4)))
	roi = ro
      vpvert = vpv(nl,1)+x*(vpv(nl,2)+x*(vpv(nl,3)+x*vpv(nl,4)))


      if (iflso(nl).eq.0) then                   !   Solid region

        if (ifanis.eq.1) then         ! Transversal isotropic material

          vphori = vph(nl,1)+x*(vph(nl,2)+x*(vph(nl,3)+x*vph(nl,4)))
          vsvert = vsv(nl,1)+x*(vsv(nl,2)+x*(vsv(nl,3)+x*vsv(nl,4)))
          vshori = vsh(nl,1)+x*(vsh(nl,2)+x*(vsh(nl,3)+x*vsh(nl,4)))
          etta   = eta(nl,1)+x*(eta(nl,2)+x*(eta(nl,3)+x*eta(nl,4)))
c Elastic constants A, C, F, L, N
           epa = ro*vphori**2
           epc = ro*vpvert**2
           epn = ro*vshori**2
           epl = ro*vsvert**2
           epf = etta*(epa-2.d0*epl)
c Calculate "equivalent" Voigt bulk and shear modulus (--> DA1980)
           vk  = (4.d0*(epa+epf-epn) + epc)/9.d0
           vmu = (epa+epc-2.d0*epf+5.d0*epn+6.d0*epl)/15.d0
c This is for the starting values:
           zmu = vmu*zdmpm(nl)
           zl2m = vk*zdmpk(nl) + 4.d0*zmu/3.d0
c Damping factors for A,C,...
           zdma = (vk*zdmpk(nl)+4.d0*vmu*zdmpm(nl)/3.d0)/(vk+4.d0*vmu/3.d0)
           zdmc = zdma
           zdmf = (vk*zdmpk(nl)-2.d0*vmu*zdmpm(nl)/3.d0)/(vk-2.d0*vmu/3.d0)
           zdml = zdmpm(nl)
           zdmn = zdmpm(nl)
c Complex moduli, handed over to the differential equations:
           zpa = epa*zdma
           zpc = epc*zdmc
           zpf = epf*zdmf
           zpl = epl*zdml
           zpn = epn*zdmn
           
       else

           vsvert = vsv(nl,1)+x*(vsv(nl,2)+x*(vsv(nl,3)+x*vsv(nl,4)))
           amu  = ro*vsvert**2
           al2m = ro*vpvert**2
           akap = al2m - (4.d0*amu)/3.d0
c Complex moduli.
           zmu = amu*zdmpm(nl)
           zlm = akap*zdmpk(nl)-(2.d0*zmu)/3.d0
           zl2m = zlm + 2.d0*zmu
           zpa = zl2m
           zpc = zpa
           zpf = zlm
           zpl = zmu
           zpn = zpl
           
       endif
       
c----------
      else                                      ! Liquid region
c----------

       zpc  = ro*vpvert*vpvert*zdmpk(nl)      
       zl2m = zpc
       zmu  = 0.d0

      endif

      return
      end
