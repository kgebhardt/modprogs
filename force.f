C=============================================================================D
C     subroutine FORCE calculates the acceleration of test particle cal-
C       culated by the r and theta derivatives of the Legendre polynomial-
C       expanded potential
C
C     USED BY LIBRARY
C
C=============================================================================D
      SUBROUTINE force(r,v,frL,fvL)
      INCLUDE 'libdefs.h'
      real frLhalo,fvLhalo,frLhalolo,frLhaloup,fvLhalolo,fvLhaloup

      rir=(log10(r)-rlgpotmin)*float(npot-1)/(rlgpotmax-
     &     rlgpotmin)+1.
c      rir=(log10(r)-log10(rpotmin))*float(npot-1)/(log10(rpotmax)-
c     &     log10(rpotmin))+1.
      irlo=max(int(rir),1)
      irup=max(int(rir+1.),irlo+1)

      rlo = RneeIRhalo(irlo)
      rup = RneeIRhalo(irup)

      rdiff = r - rlo
      frLlo=0.
      frLup=0.
      fvLlo=0.
      fvLup=0.

C------ linearly interpolate for the radial sums

      do n=0,nlegup,2
         nleg = n/2+1

         if (irlo.eq.0) then
            term = 0.
         else
            term = Pl(n,v)*tabfr(nleg,irlo)
         endif

         frLlo = frLlo + term
         frLup = frLup + Pl(n,v)*tabfr(nleg,irup)
         
         if(n.gt.0) then
            fvLlo = fvLlo + (v*Pl(n,v) - Pl(n-1,v)) *
     &           float(n) * tabfv(nleg,irlo)
            fvLup = fvLup + (v*Pl(n,v) - Pl(n-1,v)) *
     &           float(n) * tabfv(nleg,irup)
         endif
         
      enddo
      
c --- original version
      frL = frLlo + (frLup-frLlo)/(rup-rlo) * rdiff
      fvL = fvLlo + (fvLup-fvLlo)/(rup-rlo) * rdiff
c --- end


c --- for hernquist potential
c      frLlo = 1./sqrt(frLlo)
c      frLup = 1./sqrt(frLup)
c      frL = frLlo + (frLup-frLlo)/(rup-rlo) * rdiff
c      fvL = fvLlo + (fvLup-fvLlo)/(rup-rlo) * rdiff
c      frL = 1./(frL**2)
c --- end


      frL = -2.*pi*frL - hole/totlight/gdennorm/r/r
      if(abs(v).lt.1.0) then
         fvL = -2.*pi*fvL/sqrt(1.-v*v)
      else
         fvL = 0.
      endif

      RETURN
      END
      
      
