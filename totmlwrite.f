C------------------------------------------------------------------------------D
C     subroutine TOTMLWRITE writes out the total mass (sans hole) and
C       light of the galaxy for use by the model program
C
C     USED BY LIBRARY
C
C------------------------------------------------------------------------------D
      SUBROUTINE TOTMLWRITE()
      INCLUDE 'libdefs.h'

      open (unit=70,file='totml.out',status='unknown')
      write (70,*)'totmass  ',totmass
      write (70,*)'totlight ',totlight
      write (70,*)'coremass ',core*totlight
      close(70)

      RETURN
      END
