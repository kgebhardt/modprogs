C==============================================================================D
C     subroutine DATAREAD reads in dL.dat and ratML.dat
C
C     USED BY LIBRARY
C
C==============================================================================D
      SUBROUTINE dataread()
      INCLUDE 'libdefs.h'
      real dM(Nrdat,Nvdat)


C---- read in mass-to-light ratios

      open (unit=66,file='ratML.dat',status='old')
      do i=1,ivmax*irmax
        read (66,*)ir,iv,ratML(ir,iv)
      enddo
      close(66)
c      do ir=1,irmax
c         do iv=1,ivmax
c            ratML(ir,iv)=1.
c         enddo
c      enddo

      open (unit=68,file='dL.dat',status='old')
      open (unit=69,file='dM.dat',status='old')
      totlight=0.
      totmass=0.
      do i=1,ivmax*irmax
        read (68,*)ir,iv,dL(ir,iv)
        read (69,*)ir,iv,dM(ir,iv)
        totlight = totlight + dL(ir,iv)
        totmass = totmass + dM(ir,iv)
      enddo
      close(68)

C---- normalize light map to 1

      do ir=1,irmax
      do iv=1,ivmax
        dL(ir,iv) = dL(ir,iv)/totlight
      enddo
      enddo


      RETURN
      END










