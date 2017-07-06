	subroutine elapa(om,sigma,qm,qk,rho,vpvert,vsvert,vphori,vshori,etta)
c.....                                              ....................
c                                                                      .
c Calculates complex elastic parameters                                .
c                                                                      .
c.......................................................................
c
	complex zom,zdampm,zdampk,zpa,zpc,zpf,zpl,zpn,zi,zdma,zdmc,
     &        zdmf,zdml,zdmn
	real    om,sigma,omref,rho,vk,vmu,epa,epc,epf,epl,epn,pi,qm,qk,phimu,
     &        phika,vpvert,vsvert,vphori,vshori,etta
	common  /sopar/ zpa,zpc,zpf,zpl,zpn
	parameter(pi=3.141592653,omref=2.*Pi,zi=(0.,1.))

c    Transversal Isotropic elastic constants
c    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      epa = rho*vphori**2
      epc = rho*vpvert**2
      epn = rho*vshori**2
      epl = rho*vsvert**2
      epf = etta*(epa-2.*epl)

c    Dispersion factors
c    ~~~~~~~~~~~~~~~~~~
      phimu  = atan(qm)
      phika  = atan(qk)
      zom    = cmplx(om,-sigma)
      zdampm = (zom/omref)**(2.*phimu/pi)*cexp(zi*phimu)
      zdampk = (zom/omref)**(2.*phika/pi)*cexp(zi*phika)

c Calculate "equivalent" Voigt bulk and shear modulus (--> DA1980)
      vk  = (4.*(epa+epf-epn) + epc)/9.
      vmu = (epa+epc-2.*epf+5.*epn+6.*epl)/15.
c Damping factors for A,C,...
      zdma = (vk*zdampk+4.*vmu*zdampm/3.)/(vk+4.*vmu/3.)
      zdmc = zdma
      zdmf = (vk*zdampk-2.*vmu*zdampm/3.)/(vk-2.*vmu/3.)
      zdml = zdampm
      zdmn = zdampm
c                  Complex moduli in source layer.
c                  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      zpa = epa*zdma
      zpc = epc*zdmc
      zpf = epf*zdmf
      zpl = epl*zdml
      zpn = epn*zdmn
c     print *,zpa,zpc,zpf,zpl,zpn


	end
