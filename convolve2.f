
      subroutine convolve2(sumb,sumbn,v1libt,v1libt2)
      INCLUDE 'moddefs.h'
      real sumb(Nrlib,Nvlib,Nrlib,Nvlib),sumbn(Nrlib,Nvlib,Nrlib,Nvlib)
      real v1libt(Nrlib,Nvlib),v1libt2(Nrlib,Nvlib)

      do iv=1,Nvlib
         do ir=1,Nrlib
            v1libt2(ir,iv)=0.
         enddo
      enddo

      do ivb=1,Nvlib
         do irb=1,Nrlib
            do iv=1,Nvlib
               do ir=1,Nrlib
                  frac1=sumb(irb,ivb,ir,iv)
                  frac2=sumbn(irb,ivb,ir,iv)
                  v1libt2(ir,iv)=v1libt2(ir,iv)+
     $                 frac1*v1libt(irb,ivb)+
     $                 frac2*v1libt(irb,ivb)
               enddo
            enddo
         enddo
      enddo

      return
      end
