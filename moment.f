C==============================================================================D
C     subroutine MOMENT is a (modified) Press et al. routine used by
C       qual() to calculate statistical quantities related to the model
C       fit
C
C     USED BY MODEL
C
C==============================================================================D
      SUBROUTINE moment(data,n,ave,adev,sdev,var)
      INTEGER n
      REAL adev,ave,sdev,var,data(n)
      INTEGER j
      REAL p,s,ep
      if(n.le.1)print *,'n must be at least 2 in moment'
      s=0.
      do 11 j=1,n
        s=s+data(j)
11    continue
      ave=s/n
      adev=0.
      var=0.
      ep=0.
      do 12 j=1,n
        s=data(j)-ave
        ep=ep+s
        adev=adev+abs(s)
        p=s*s
        var=var+p
        p=p*s
        p=p*s
12    continue
      adev=adev/n
      var=(var-ep**2/n)/(n-1)
      sdev=sqrt(var)
      RETURN
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0S,)0(5'R3.
