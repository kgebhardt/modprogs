C=============================================================================D
C     subroutine STEP calculates the step size for the Runge-Kutta
C       integration routine (RK4)
C
C     USED BY LIBRARY
C
C=============================================================================D
      SUBROUTINE step(pos0,xLz,dt,Nvar,dposdt,epsilon)
      DIMENSION pos0(Nvar),dposdt(Nvar)

      r   = pos0(1)
      v   = sin(pos0(2))
      vc  = cos(pos0(2))
      vr  = pos0(3)
      vth = pos0(4)
      vph = xLz/r/vc
      ar  = dposdt(3)
      ath = dposdt(4)
      aph = ( xLz/(r*r*vc*vc) ) * (vth*v - vr*vc)

      a = sqrt(ar*ar + ath*ath + aph*aph)
      vel = (vr*vr + vth*vth + vph*vph)

      dt = epsilon*sqrt( r / (a + vel/r) )

c      dt = dt*dt*dt*dt/(dt+1.e-6)/(dt+1.e-6)/(dt+1.e-6)+1.e-6

      RETURN
      END
