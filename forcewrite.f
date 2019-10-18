C==============================================================================D
C     subroutine FORCEWRITE samples the force components
C
C     USED BY LIBRARY
C
C==============================================================================D
      SUBROUTINE forcewrite()
      INCLUDE 'libdefs.h'

      open (unit=81,file='force.out',status='unknown')

88    format (4(e12.6,2x))

      ddd = log(3.6/b)/999.
      ccc = b/3./exp(ddd)
      do irr=1,1000
        r = ccc*exp(ddd*float(irr))
        call force(r,0.,frL0,fvL)
        call force(r,.5,frL5,fvL)
        call force(r,.99,frL9,fvL)
        write (81,88)r,-frL0,-frL5,-frL9
      enddo

      r=0.1

      do ith=0,899
        th = float(ith)/10.*pi/180.
        v = sin(th)
        call force(0.05,v,frL,fvL05)
        call force(0.5,v,frL,fvL5)
        call force(0.9,v,frL,fvL9)
        write (81,88)180.*th/pi,-fvL05,-fvL5,-fvL9
      enddo

      close(81)

      RETURN
      END
