      SUBROUTINE addtables()
      INCLUDE 'libdefs.h'

      do n=0,nlegup,2
         nleg = n/2+1
         do irtab = 1,npot
            tabv(nleg,irtab)=tabv(nleg,irtab)+htabv(nleg,irtab)
            tabfr(nleg,irtab)=tabfr(nleg,irtab)+htabfr(nleg,irtab)
            tabfv(nleg,irtab)=tabfv(nleg,irtab)+htabfv(nleg,irtab)
         end do
      end do
      RETURN
      END
