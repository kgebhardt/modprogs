C==============================================================================D
C     subroutine SPACEM calculates the spatial arrays vol2d(ir), area2d(ir),
C       area2dC(ir,iv), and rect(islit)
C
C     USED BY MODEL
C
C==============================================================================D
      SUBROUTINE spacem()
      INCLUDE 'moddefs.h'


C---- calculate 2-index 3-space volume and 2-index 2-space sky area (fine grid)

      dv2 = 1./float(Nvdat) / 2.
      do ir=1,Nrdat
        rlo = (RneeIR(ir-1)+RneeIR(ir))/2.
        rup = (RneeIR(ir)+RneeIR(ir+1))/2.
        vol2d(ir) = 4./3.*pi * (rup*rup*rup-rlo*rlo*rlo) / float(ivmax)
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
 
      RETURN
      END
