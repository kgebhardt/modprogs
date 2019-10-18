
      subroutine vdataread(filen)
      parameter(nadt=1000)
      INCLUDE 'moddefs.h'
      real frac(Nvelbm),vnorm(Nbin),vnorma(Nbin,Nstot)
      character filen(Nveld)*40,filesb(9)*40

      filesb(1)='SBsc1.dat'
      filesb(2)='SBsc2.dat'
      filesb(3)='SBsc3.dat'
      filesb(4)='SBsc4.dat'
      filesb(5)='SBsc5.dat'
      filesb(6)='SBsc6.dat'
      filesb(7)='SBsc7.dat'
      filesb(8)='SBsc8.dat'
      filesb(9)='SBsc9.dat'

c - get the light which has been convolved for the normalization

      izerofwhm=0
      ns=0
      if(Nvelb.gt.1) then
         do i=1,Nstot
            if(seeb(i).gt.0) then
               ns=ns+1
               open(unit=1,file=filesb(ns),status='old')
               read(1,*) vnorm
               close(1)
               do ib=1,Nbin
                  vnorma(ib,ns)=vnorm(ib)
               enddo
            else
               izerofwhm=1
            endif
         enddo
      endif

c -- sum along model bins according to the individual configurations

      do ivbin=1,Nvelb
         frac(ivbin)=0.
      enddo

      do ivbin=1,Nvelb
         if(seeb(isee(ivbin)).eq.0.or.Nvelb.eq.1) then
c - set equal to the un-convolved light
            do ibinp=1,Nbin
               frac(ivbin)=frac(ivbin)+
     $              SMcorig(ibinp)*areakin(ivbin,ibinp)
            enddo
         else
            iseeb=isee(ivbin)-izerofwhm
            do ibinp=1,Nbin
               frac(ivbin)=frac(ivbin)+
     $              vnorma(ibinp,iseeb)*areakin(ivbin,ibinp)
            enddo
         endif
      enddo

      do ivel=1,Nvel
         velm(ivel)=vmin+(vmax-vmin)*float(ivel-1)/float(Nvel-1)
      enddo

      do ivbin=1,Nvelb
         suma=0.
         open(unit=56,file=filen(ivbin),status='old')
         nad=0
         do i=1,nadt
            read (56,*,end=666) x1,x2,x3,x4
            nad=nad+1
            veld(nad,ivbin)=x1
            ad(nad,ivbin)=x2
            adfer(nad,ivbin)=max(x4-x2,(x4-x3)/2.)
            suma=suma+ad(nad,ivbin)
         enddo
 666     continue
         close(56)
         if(nad.gt.Nveld) print *,'Make Nveld bigger'
         if(frac(ivbin).eq.0) print *,
     $'frac is zero; increase ny in subroutine getbin in velbin.f'

         do i=1,nad
            ad(i,ivbin)=ad(i,ivbin)*frac(ivbin)/suma
            adfer(i,ivbin)=adfer(i,ivbin)*frac(ivbin)/suma
         enddo

      enddo

      return
      end
