C==============================================================================D
C     subroutine GALAXREAD reads in galaxy.params
C
C     USED BY LIBRARY AND MODEL
C
C==============================================================================D
      SUBROUTINE galreadm(galname)
      INCLUDE 'moddefs.h'

      CHARACTER dummy*20,galname*40

      open (unit=65,file='galaxy.params',status='old')
      open (unit=66,file='gal.dat',status='old')
      read (65,*)galname
      read (66,*)dummy,hole
      read (66,*)dummy,xinclin
      read (65,*)dummy,distance
      read (65,*)dummy,angrad
      read (65,*)dummy,vmin
      read (65,*)dummy,vmax
      read (65,*)dummy,cslope
      read (65,*)dummy,xmtol
      read (65,*)dummy,islitsize
      close(65)
      close(66)
      xmutarget=1./sqrt(xmtol)
      xinclin = xinclin*pi/180.
      sini=sin(xinclin)
      cosi=cos(xinclin)
      distance = distance*1.e6
      fac1 = sqrt(GG*totlight*arcsec/distance/angrad)
      vmult=float(Nvel-1)/(vmax-vmin)*fac1
      vadd=1.-float(Nvel-1)*vmin/(vmax-vmin)

      perimin = xorbitmin/angrad
      apomax = max(xorbitmax/angrad,1.2)

      rpotmin = min(perimin/10.,RneeIR(1)/10.)
      rpotmax = max(apomax*2.,RneeIR(Nrdat)*2.)
      rlgpotmin = log10(rpotmin)
      rlgpotmax = log10(rpotmax)

      RETURN
      END
