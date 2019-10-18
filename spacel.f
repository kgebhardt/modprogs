C==============================================================================D
C     subroutine SPACEL calculates the spatial arrays vol2d(ir), area2d(ir),
C       area2dC(ir,iv), rect(islit), and determines if the matter has a non-
C       zero quadrupole moment
C
C     USED BY LIBRARY
C
C==============================================================================D
      SUBROUTINE spacel()
      INCLUDE 'libdefs.h'


C---- calculate 2-index 3-space volume and 2-index 2-space sky area (fine grid)

      dv2 = 1./float(Nvdat) / 2.
      do ir=1,Nrdat
        rlo = (RneeIR(ir-1)+RneeIR(ir))/2.
        rup = (RneeIR(ir)+RneeIR(ir+1))/2.
        vol2d(ir) = 4./3.*pi * (rup*rup*rup-rlo*rlo*rlo) / float(ivmax)
        d0 = dL(ir,1)
        do iv=1,Nvdat
          thlo = asin(max(VneeIV(iv)-dv2,0.))
          thup = asin(min(VneeIV(iv)+dv2,1.))
          area2d(ir,iv) = 4.*(rup*rup-rlo*rlo)*(thup-thlo)/2.
        enddo
      enddo


C---- calculate 2-index 2-space sky area (coarse grid)

      do ivS=1,Nvlib
      do irS=1,Nrlib
        area2dC(irS,ivS) = 0.0
        do k=1,ivrat
        do j=1,irrat
          irSS = (irS-1)*irrat+j
          ivSS = (ivS-1)*ivrat+k
          area2dC(irS,ivS) = area2dC(irS,ivS) + area2d(irSS,ivSS)
        enddo
        enddo
      enddo
      enddo


C---- determine if mass distribution has a non-zero quadupole moment

      do ir=1,Nrdat
        d0 = dL(ir,1)
        do iv=1,Nvdat
          if (dL(ir,iv).ne.d0) then
            iquad=1
            goto 123
          endif
        enddo
      enddo
      iquad=0
123   iquad=iquad

      if(qdm.ne.1.0) iquad=1

      if (iquad.eq.0) then
        nlegup=0
      else
        nlegup=nlegmax
      endif

      RETURN
      END
