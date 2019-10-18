C==============================================================================D
C     subroutine BINSET sets the binning parameters
C
C     USED BY LIBRARY AND MODEL
C
C==============================================================================D
      SUBROUTINE binset(Nrdat,Nvdat,Nrlib,Nvlib)
      COMMON/rbinc/a,b,rmin
      COMMON/kbin/irmax,ivmax,irrat,ivrat

      COMMON/relzbin/aelz,belz,relzmin,relzmax

      irmax = Nrdat
      ivmax = Nvdat

      irrat = Nrdat/Nrlib
      ivrat = Nvdat/Nvlib

      open(unit=1,file='ab.dat',status='old')
      read(1,*) a,b
      close(1)

      open(unit=1,file='abelz.dat',status='old')
      read(1,*) aelz,belz
      read(1,*) relzmin,relzmax
      close(1)

      rmin = RneeIR(1)

      RETURN
      END
