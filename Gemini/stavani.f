      subroutine stavani(zyst,xfi,nla,csys)
cccccalls cpanis, zetl
c
* Calculates STARTING VALUES for integration of the spheroidal
* and toroidal motion upwards to the source level.
* The values are derived from a homogeneous isotropic Earth model
* assumed beneath the starting radius r1. The results are
* spherical Bessel functions. To avoid explicit calculation of
* those functions, we use ratios of them. These ratios can 
* be computed from a recurrence relation. See subroutine ZETL().
c Author: J.R. Dalkolmo

      integer nlayer,nla,ifliq,nlay,l,nl
      parameter(nlayer=20)
      double complex zyst(6),zom,zxa2,zeta,zxb2,zetb,zl2m,zmue,za,zb
      real*8 om,elp1,xfi,ro1,el,rsq,rboc
      integer iflso(nlayer)
      character csys*1
      common/loop/zom,om,elp1,l
      common/modlay/nlay /layer/ nl /modls/ iflso
      common/epi/zmue,zl2m,ro1
      common /core/ rboc
      save


c Get model parameters (seismic velocities) at the starting radius
      nl = nla
      ifliq = iflso(nl)
      call cpanis(xfi)
      rsq = xfi*xfi
      el = dble(l)


      if(csys.eq.'s') then            !   Spheroidal motion

       zxa2=ro1*zom*zom*rsq/zl2m
c      print *, 'zxa2=',zxa2
       call zetl(zxa2,l,zeta)   ! Ratio of spherical Bessel functions.
       if(ifliq.eq.0)then
         zxb2=ro1*zom*zom*rsq/zmue
         call zetl(zxb2,l,zetb)
       endif

c Starting values for integration UPWARDS in liquid or solid layer.

      if (l.eq.0) then
        zyst(1)=-zeta
        zyst(2)=-ro1*xfi*zom*zom + 4.d0*zmue/xfi*zeta
        return
      else
        if(ifliq.eq.1)then
          zyst(1) = -(1.d0-zeta/el)/(zom*zom*ro1*xfi)
          zyst(2) = 1.d0/el
          return 
        else
          za=zeta/el
          zb=zetb/(el+1.d0)
          zyst(2)=(-za+zb*(za-1.d0))/el
          zyst(1)=zmue/xfi*(zxb2/el+2.d0*elp1*zyst(2))
          zyst(3)=zmue/xfi*(-2.d0*zyst(2)+zxb2/elp1*(za-1.d0))
          zyst(4)=-zmue/xfi*(zxb2/(el*el)*(1.d0-zb)+4.d0*zyst(2))
          zyst(5)=zmue*zmue/rsq*(-4.d0*(el-1.d0)*(el+2.d0)*zyst(2)+
     1       zxb2/el*(zxb2/elp1-2.d0*(el-1.d0)*(2.d0*el+1.d0)/elp1
     2                       -4.d0/(el+1.d0)*za-2.d0/el*zb))
c          print *,(zyst(i),i=1,5)
        endif
       endif

      else if(csys.eq.'t') then   !    Toroidal motion
      
       if (xfi.gt.rboc) then
          zxb2=ro1*zom*zom*rsq/zmue
          call zetl(zxb2,l,zetb)
          zyst(1) = 1.d0/el
c         print *, zetb
          zyst(2) = zmue/(xfi*el)*(el-1.d0-zetb)
       else         ! Start at core-mantle boundary: adjust xfirst
          xfi=rboc
          zyst(1) = 1.d0/el
          zyst(2) = 0.d0
       endif

      else
       stop 'Error in <stavani>: No motion specified!'
      endif


      return
      end 
