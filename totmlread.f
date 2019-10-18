C------------------------------------------------------------------------------D
C     subroutine TOTMLREAD reads the total mass and light of the galaxy
C
C     USED BY MODEL
C
C------------------------------------------------------------------------------D
      SUBROUTINE TOTMLREAD()
      INCLUDE 'moddefs.h'
      CHARACTER dummy

      open (unit=70,file='totml.out',status='old')
      read (70,*)dummy,totmass
      read (70,*)dummy,totlight
      read (70,*)dummy,core
      close(70)

      RETURN
      END
