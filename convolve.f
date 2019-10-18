
      subroutine convolve(sumb,sumbn,v1libt,v1libt2)
      INCLUDE 'moddefs.h'
      real sumb(Nrlib,Nvlib,Nrlib,Nvlib),sumbn(Nrlib,Nvlib,Nrlib,Nvlib)
      real v1libt(Nvel,Nrlib,Nvlib),v1libt2(Nvel,Nrlib,Nvlib)

      do iv=1,Nvlib
         do ir=1,Nrlib
            do ivel=1,Nvel
               v1libt2(ivel,ir,iv)=0.
            enddo
         enddo
      enddo

      do ivb=1,Nvlib
         do irb=1,Nrlib
            do iv=1,Nvlib
               do ir=1,Nrlib
                  frac1=sumb(irb,ivb,ir,iv)
                  frac2=sumbn(irb,ivb,ir,iv)
                  do ivel=1,Nvel
                     v1libt2(ivel,ir,iv)=v1libt2(ivel,ir,iv)+
     $                    frac1*v1libt(ivel,irb,ivb)+
     $                    frac2*v1libt(Nvel-ivel+1,irb,ivb)
                  enddo
               enddo
            enddo
         enddo
      enddo

      return
      end
