C==============================================================================D
C     subroutine QUAL uses the (modified) Press et al. subroutine moment()
C       to calculate statistical quantities of the model fit
C
C     USED BY MODEL
C
C==============================================================================D
      SUBROUTINE qual(errA,errR,aveA,aveR,adevA,adevR,sdevA,sdevR)
      INCLUDE 'moddefs.h'
      DIMENSION errA(Nbin),errR(Nbin)

      call moment(errA,Nbin,aveA,adevA,sdevA,varA)
      call moment(errR,Nbin,aveR,adevR,sdevR,varR)

      RETURN
      END
