C=============================================================================D
C     subroutine MINMAX determines the best and worst bin matches
C
C     USED BY MODEL
C
C=============================================================================D
      SUBROUTINE minmax(errA,errR,errAmin,errRmin,errAmax,errRmax,
     &  ibAmin,ibAmax,ibRmin,ibRmax)
      INCLUDE 'moddefs.h'
      DIMENSION errA(Nbin),errR(Nbin)

      errAmin =  1.e20
      errRmin =  1.e20
      errAmax =  0.
      errRmax =  0.

      do ibin=1,Nbin

        if (abs(errA(ibin)).lt.abs(errAmin)) then
          ibAmin  = ibin
          errAmin = errA(ibin)
        endif

        if (abs(errR(ibin)).lt.abs(errRmin)) then
          ibRmin  = ibin
          errRmin = errR(ibin)
        endif

        if (abs(errA(ibin)).ge.abs(errAmax)) then
          ibAmax  = ibin
          errAmax = errA(ibin)
        endif

        if (abs(errR(ibin)).ge.abs(errRmax)) then
          ibRmax  = ibin
          errRmax = errR(ibin)
        endif

      enddo

      RETURN
      END
