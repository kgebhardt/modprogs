C=============================================================================D
C     subroutine COMPARER compares the result of subroutine summ() (i.e.
C       the summed() array) with the data array SMc().  It returns the
C       model fit, and the absolute and relative errors: summed(ibin),
C       errA(ibin), and errR(ibin), respectively.
C
C     USED BY MODEL
C
C=============================================================================D
      SUBROUTINE comparer(errA,errR)
      INCLUDE 'moddefs.h'
      DIMENSION errA(Nbin),errR(Nbin)

      do ibin=1,Nbin
        errA(ibin) = SMc(ibin) - summed(ibin)
        errR(ibin) = summed(ibin)/SMc(ibin) - 1.
c        print *,ibin,errA(ibin),errR(ibin)
      enddo

      RETURN
      END
