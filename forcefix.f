C==============================================================================D
C     subroutine FORCEFIX calculates the correction terms for the
C       force in the core and at the outer edges when the mass
C       distribution is not spherical. Also does potential at edge.
C
C     USED BY LIBRARY
C
C==============================================================================D
      SUBROUTINE forcefix()
      INCLUDE 'libdefs.h'

      frLmin = -( core/rmin/rmin + hole/totlight/rmin/rmin )
      frLmax = -( totmass/totlight + hole/totlight )
      VVmax=frLmax
      n=0
      nstep=10
      do i=1,nstep
         v=0.+float(i-1)/float(nstep-1)*1.0
c      do v=0.,1.0001,.1
         n=n+1
         vin=min(1.,v)
         vfix(n)=vin
         call force(1.002*rmin,vin,frL,fvL)
         fmin(n)=frL/frLmin
         call force(0.998,vin,frL,fvL)
         fmax(n)=frL/frLmax
         call potential(0.998,vin,VV)
         pmax(n)=VV/VVmax
      enddo
      return
      end
