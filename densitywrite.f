C==============================================================================D
C     subroutine DENSITYWRITE writes out the density derived from the
C       original data map
C
C     USED BY LIBRARY
C
C==============================================================================D
      SUBROUTINE densitywrite()
      INCLUDE 'libdefs.h'

      open (unit=74,file='density.out',status='unknown')

741   format (i4,2x,i4,4(2x,e12.6))

      term = totlight/distance/distance/distance*arcsec*arcsec*arcsec
     &  /angrad/angrad/angrad


      do ir=1,irmax
      do iv=1,ivmax

         r = RneeIR(ir)
         v = VneeIV(iv)
         write (74,741)ir,iv,r,v,dL(ir,iv)*term/vol2d(ir),
     &        rho(ir,iv)*term

         write (19,741)ir,iv,r*angrad,v,dL(ir,iv)*term/vol2d(ir),
     &        rho(ir,iv)*term

       enddo
       enddo

       close(74)
       
       RETURN
       END
