C==============================================================================D
C     subroutine DENSITY calculates the density profile of the distribution
C
C     USED BY LIBRARY
C
C==============================================================================D
      SUBROUTINE density()
      INCLUDE 'libdefs.h'

      do ir=1,irmax
         do iv=1,ivmax

            rho(ir,iv) = dL(ir,iv)*ratML(ir,iv)/vol2d(ir)

         enddo
      enddo

      RETURN
      END
