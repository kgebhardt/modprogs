C==============================================================================D
C     functions which determine (ir,iv) from ibin, or ibin from (ir,iv)
C
C     USED BY MODEL
C
C==============================================================================D

C------------------------------------------------------------------------------D
C     function ITOBIN converts (ir,iv) to (ibin)
C------------------------------------------------------------------------------D
      FUNCTION itobin(ir,iv)
      INCLUDE 'moddefs.h'
      if(ir.eq.1) then
         itobin=1
      else
         itobin=1+(ir-2)*Nvlib+iv
      endif
      RETURN
      END


C------------------------------------------------------------------------------D
C     function ITOR converts (ibin) to (ir)
C------------------------------------------------------------------------------D
      FUNCTION itor(ibin)
      INCLUDE 'moddefs.h'
      if(ibin.eq.1) then
         itor=1
      else
         iv=itov(ibin)
         itor=nint(float(ibin-1-iv)/float(Nvlib)+2.)
      endif
      RETURN
      END


C------------------------------------------------------------------------------D
C     function ITOV converts (ibin) to (iv)
C------------------------------------------------------------------------------D
      FUNCTION itov(ibin)
      INCLUDE 'moddefs.h'
      if(ibin.eq.1) then
         itov=1
      else
         x=float(ibin-1)/float(Nvlib)-.00001
         itov=nint((x-int(x))*Nvlib)
      endif
      RETURN
      END
