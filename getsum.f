
      subroutine getsum(is,see,sumb,sumbn)
      INCLUDE 'moddefs.h'
      parameter(npsfmax=1000)
      real sumb(Nrlib,Nvlib,Nrlib,Nvlib),sumbn(Nrlib,Nvlib,Nrlib,Nvlib)
      real sumbt(Nrlib,Nvlib),sumbnt(Nrlib,Nvlib)
      real psftemp

      real psfin,xin,yin,psf(npsfmax,npsfmax)
      real xpsf(npsfmax,npsfmax),ypsf(npsfmax,npsfmax),rpsf
      integer npsfx,npsfy,ipsfx,ipsfy
      integer np,npx,npy
      CHARACTER psffile*40

      rbig=1000.
      np=100
      npx = np
      npy = np

      do ipsfx=1,npsfmax
         do ipsfy=1,npsfmax
            psf(ipsfx,ipsfy)=0.0
         end do
      end do

      if(see*angrad.eq.9.999) then
         !m31p3 setup!
      else
         if(see*angrad.gt.99.00) then
            !get psf from file!
            npsfx=0
            npsfy=0
            rpsf=0.0
            
            open(unit=43,file='psflist.dat',status='old')
            do ix=1,is-1
               read(43,*)
            end do
            read(43,*) psffile
            close(43)

            open(unit=43,file=psffile,status='old')
 201        continue
            read(43,*,end=202) psfin,xin,yin,ipsfx,ipsfy
            psf(ipsfx,ipsfy)=psfin
            xpsf(ipsfx,ipsfy)=xin/angrad
            ypsf(ipsfy,ipsfy)=yin/angrad
            rpsf=max(rpsf,sqrt(xin*xin+
     &           yin*yin)/angrad)
            npsfx=max(npsfx,ipsfx)
            npsfy=max(npsfy,ipsfy)
            goto 201
 202        continue
            close(43)
            print*,'In getsum: PSF from file',psffile,' (',npsfx,',',
     &           npsfy,',',rpsf,')'
            npx = npsfx
            npy = npsfy
            rmaxs = rpsf
         else
            rmaxs=10.*see/2.35
            den=2.*see*see/2.35/2.35
         end if
      end if
c -  this is for the STIS PSF
c      if(is.eq.1.and.see.gt.0.085/angrad) then
c         den1=2.*(0.0329/angrad)**2
c         den2=2.*(0.0384/angrad)**2
c         r2=0.105/angrad
c         amp2=0.12
c         print *,'In getsum: ',den1,den2,r2,amp2
c      elseif(is.eq.1.and.see.lt.0.085/angrad) then
c         den1=2.*(0.0340/angrad)**2
c         den2=2.*(0.0384/angrad)**2
c         r2=0.105/angrad
c         amp2=0.
c         print *,'In getsum: ',den1,den2,r2,amp2
c      endif

      do ir=1,Nrlib
         r=RneeIR(ir*irrat-1)
c         if(r.gt.rbig/angrad) goto 666
         do iv=1,Nvlib
                  
            do irp=1,Nrlib
               do ivp=1,Nvlib
                  sumbt(irp,ivp)=0.
                  sumbnt(irp,ivp)=0.
               enddo
            enddo

c - get the bin centers

            v=VneeIV(iv*ivrat-1)
            x=r*sqrt(1.-v*v)
            y=r*v
            xmin=x-rmaxs
            xmult=2.*rmaxs/float(np-1)
            ymin=y-rmaxs
            
            sumg=0.
            do ix=1,npx
               if(see*angrad.lt.99.00) then
                  xp=xmin+(ix-1)*xmult
                  xd=x-xp
               end if
               do iy=1,npy
                  if(see*angrad.gt.99.00) then
                     xp = x+xpsf(ix,iy)
                     xd = x-xp
                     yp = y+ypsf(ix,iy)
                     yd = y-yp
                  else
                     yp=ymin+(iy-1)*xmult
                     yd=y-yp
                  end if
                  rp=sqrt(xp*xp+yp*yp)
                  vp=abs(yp/rp)
                  irp=IRCneeR(rp)
                  ivp=IVCneeV(vp)
                  rd=xd*xd+yd*yd
                  if(see*angrad.eq.9.999) then
                                !m31p3 setup!
                     psftemp=(1.+(sqrt(rd)/
     &                    (0.124/angrad))**(1.5))**(-4.4)
                  else
                     if(see*angrad.gt.99.00) then
                                !psf from file psf.dat!
                        psftemp = psf(ix,iy)
                     else
                                !gaussian psf!
                        psftemp = exp(-rd/den)
c -  this is for the STIS PSF
c                        if(is.eq.1.and.see*angrad.gt.0.085) then
c                           den1=2.*(0.0329/angrad)**2
c                           den2=2.*(0.0384/angrad)**2
c                           r2=0.105/angrad
c                           amp2=0.12
c                           psftemp=1.01*exp(-rd/den1)+
c     $                          amp2*exp(-(sqrt(rd)-r2)**2/den2)
c                        endif
                     endif
                  endif
                  sumg=sumg+psftemp
                  if(irp.gt.0.and.irp.le.Nrlib) then
                     if(xp.gt.0.) then
                        sumbt(irp,ivp)=sumbt(irp,ivp)+psftemp
                     else
                        sumbnt(irp,ivp)=sumbnt(irp,ivp)+psftemp
                     endif
                  endif
               enddo
            enddo
            
            do irp=1,Nrlib
               do ivp=1,Nvlib
                  sumb(ir,iv,irp,ivp)=sumbt(irp,ivp)/sumg
                  sumbn(ir,iv,irp,ivp)=sumbnt(irp,ivp)/sumg
                  if(sumb(ir,iv,irp,ivp).ne.0.0) then
                  end if
               enddo
            enddo
         enddo
      enddo
 666  continue

      return
      end
