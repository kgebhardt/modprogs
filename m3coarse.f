C==============================================================================D
C     subroutine SMCOARSE creates a coarse-binned surface mass map and then
C       smears out the ir=1 bin (i.e. sums over all angles)
C
C     USED BY MODEL
C
C==============================================================================D
      SUBROUTINE m3coarse()
      INCLUDE 'moddefs.h'

C---- read in fine-grid surface brightness map and convert it to surface mass
 
      open (unit=67,file='dM.dat',status='old')
      do i=1,ivmax*irmax
        read (67,*)ir,iv,SBM(ir,iv)
      enddo
      close(67)

      do ivS=1,Nvlib
      do irS=1,Nrlib
        SMctemp(irS,ivS) = 0.0
        do k=1,ivrat
        do j=1,irrat
          irSS = (irS-1)*irrat+j
          ivSS = (ivS-1)*ivrat+k
          SMctemp(irS,ivS) = SMctemp(irS,ivS) + SBM(irSS,ivSS)
        enddo
        enddo
      enddo
      enddo

      do iv=1,Nvlib
        d3c(1) = d3c(1) + SMctemp(1,iv)
      enddo

      total = d3c(1)
      do ibin=2,Nbin
        ir = itor(ibin)
        iv = itov(ibin)
        d3c(ibin) = SMctemp(ir,iv)
        total = total + d3c(ibin)
      enddo

C---- normalize it:
      do ibin=1,Nbin
        d3c(ibin) = d3c(ibin)/total
      enddo
      print *,'total M3 ',total

      close(69)

      RETURN
      END
