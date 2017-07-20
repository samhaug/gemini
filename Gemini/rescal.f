c *****************************************************************
c ** Skalierung des Funktionenvektors um overflow zu vermeiden
c **
      subroutine rescal(zy,nvar)

      integer nvar,i, nexponent
      double complex zy(nvar)
      double precision maxabs, scalfac

c Suche Maximalbetrag des Funktionenvektors
      maxabs = 0.d0
      do i = 1, nvar
        maxabs = max(zabs(zy(i)),maxabs)
      enddo

c Skaliere die einzelnen Komponenten des Vektors falls maxabs groesser
c als 10 ist.
      if (maxabs.gt.10.d0) then
        nexponent = int(log10(maxabs))
        scalfac = 10.d0**(-nexponent)
        do i = 1, nvar
          zy(i) = zy(i)*scalfac
        enddo
      endif

      return
      end
