
      INCLUDE 'moddefs.h'
      real sumb(Nrlib,Nvlib,Nrlib,Nvlib),sumbn(Nrlib,Nvlib,Nrlib,Nvlib)
      real SMcsc(Nbin),SMctemp2(Nrlib,Nvlib)
      character filen(Nveld)*40,galname*40,filesb(9)*40

      filesb(1)='SBsc1.dat'
      filesb(2)='SBsc2.dat'
      filesb(3)='SBsc3.dat'
      filesb(4)='SBsc4.dat'
      filesb(5)='SBsc5.dat'
      filesb(6)='SBsc6.dat'
      filesb(7)='SBsc7.dat'
      filesb(8)='SBsc8.dat'
      filesb(9)='SBsc9.dat'

      call binset(Nrdat,Nvdat,Nrlib,Nvlib)
      call galreadm(galname)
      open(unit=1,file='SB.new',status='old')
      read (1,*) ir,iv,x1
      do iv=1,Nvlib
         SMctemp(1,iv)=x1/float(Nvlib)
      enddo
      do ibin=2,Nbin
         read(1,*) ir,iv,SMctemp(ir,iv)
      enddo
      close(1)
      call velbin(filen)

      open(unit=2,file='see.dat',status='old')
      read(2,*)
      do i=1,Nstot
         read(2,*,end=667) i1,x2
         if(i1.ne.i) print *,'Order see.dat or else....'
         seeb(i1)=x2/angrad
      enddo
 667  continue
      close(2)

      ifile=0
      do is=1,Nstot
         if(seeb(is).gt.0.) then
            call getsum(is,seeb(is),sumb,sumbn)

c - convolve the light
            call convolve2(sumb,sumbn,SMctemp,SMctemp2)

            SMcsc(1)=0.
            do iv=1,Nvlib
               SMcsc(1) = SMcsc(1) + SMctemp2(1,iv)
            enddo
 
            total = SMcsc(1)
            do ibin=2,Nbin
               ir = itor(ibin)
               iv = itov(ibin)
               SMcsc(ibin) = SMctemp2(ir,iv)
               total = total + SMcsc(ibin)
            enddo

C---- normalize it:
            do ibin=1,Nbin
               SMcsc(ibin) = SMcsc(ibin)/total
            enddo

            ifile=ifile+1
            print *,filesb(ifile)
            open (unit=3,file=filesb(ifile),status='unknown')
            write(3,*) SMcsc       
            close(3)

         endif
      enddo
      write(*,*)

      end
