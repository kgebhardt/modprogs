C==============================================================================D
C     function RTBIS is a (modified) Press et al. routine used to locate
C       the zero-velocity surface
C
C     USED BY LIBRARY
C
C==============================================================================D
      FUNCTION rtbis(func,x1,x2,xacc,xopt)
C     FUNCTION rtbis(func,x1,x2,xacc)
      INTEGER JMAX
      REAL rtbis,x1,x2,xacc,func,xopt
C     REAL rtbis,x1,x2,xacc,func
C     EXTERNAL func
      PARAMETER (JMAX=40)
      INTEGER j
      REAL dx,f,fmid,xmid
      fmid=func(x2,xopt)
C     fmid=func(x2)
      f=func(x1,xopt)
C     f=func(x1)
      if(f*fmid.ge.0.) print *,'root must be bracketed in rtbis'
      if(f.lt.0.)then
        rtbis=x1
        dx=x2-x1
      else
        rtbis=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*.5
        xmid=rtbis+dx
        fmid=func(xmid,xopt)
C       fmid=func(xmid)
        if(fmid.le.0.)rtbis=xmid
        if(abs(dx).lt.xacc .or. fmid.eq.0.) return
11    continue
      print *,'too many bisections in rtbis'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0S,)0(5'R3.
