C=============================================================================D
C     subroutine VELBIN reads velbin.dat which list the filename and then
C       on the next line the x, y, and seeing pointers (given in see.dat)
C       of the LOSVD's positions to determine which vlib bins to use in 
C       the modelling.
C
C     USED BY MODEL
C
C=============================================================================D

      subroutine velbin(filen)
      INCLUDE 'moddefs.h'
      real v(Nvelbm),rval(Nrlib),vval(Nvlib)
      character filen(Nveld)*40

      open(unit=1,file='velbin.dat',status='old')

      open(unit=2,file='bin_r.out',status='old')
      read(2,*)
      read(2,*)
      do i=1,Nrlib
         read(2,*,end=666) i1,x2,x3,x4
         rval(i)=x4
      enddo
 666  continue
      close(2)
      open(unit=2,file='bin_v.out',status='old')
      read(2,*)
      read(2,*)
      do i=1,Nvlib
         read(2,*,end=667) i1,x2,x3,x4
         vval(i)=x4
      enddo
 667  continue
      close(2)
            
      n=0
      do i=1,Nvelbm
         read(1,*,end=668) filen(i),x1,x2,i3,x4,x5,x6,x7
         n=n+1
         ibinset=0
         if(i3.lt.0) then
            i3=-i3
            ibinset=1
         endif
         xbin(n)=x1/angrad
         ybin(n)=x2/angrad
         rbin(n)=sqrt(x1*x1+x2*x2)/angrad
         if(rbin(n).eq.0) then
            v(n)=0.
         else
            v(n)=x2/rbin(n)/angrad
         endif
         if(x1.eq.0) then
            iminor(n)=1
         else
            iminor(n)=0
         endif
         isee(n)=i3
         if(x4.lt.x5.and.ibinset.eq.0) then
            call getbin(i,x4,x5,x6,x7,rval,vval)
         else
            if(x6.le.x7.and.x5.gt.-665.and.ibinset.eq.0) then
               ir=IRCneeR(rbin(n))
               iv=IVCneeV(v(n))
               ibin=itobin(ir,iv)
               do ibinp=1,Nbin
                  areakin(i,ibinp)=0.
               enddo
               areakin(i,ibin)=1.
c            pause "You better check this!"
            else
               do ibinp=1,Nbin
                  areakin(i,ibinp)=0.
               enddo
               ir=IRCneeR(rbin(n))
               if(x6.gt.x7) then
                  ivs=1
                  ive=Nvlib
               else
                  ivs=nint(x6)
                  ive=nint(x7)
               endif
               if(ibinset.eq.0) then
                  do iv=ivs,ive
                     ibin=itobin(ir,iv)
                     areakin(i,ibin)=1.
                  enddo
               else
                  irs=nint(x4)
                  ire=nint(x5)
                  ivs=nint(x6)
                  ive=nint(x7)
                  do ir=irs,ire
                     do iv=ivs,ive
                        ibin=itobin(ir,iv)
                        areakin(i,ibin)=1.
                     enddo
                  enddo
               endif
            endif
         endif
      enddo
 668  continue
      close(1)
      Nvelb=n

      do i=1,Nvelb
         irc=IRCneeR(rbin(i))
         ivc=IVCneeV(v(i))
         iang(i)=ivc
         ivelb(i)=itobin(irc,ivc)
         if(rbin(i).eq.0) ivelb(i)=1
         print *,i,irc,ivc,ivelb(i),isee(i)
      enddo

      return
      end

      subroutine getbin(ibin,xlp,xhp,ylp,yhp,rval,vval)
      INCLUDE 'moddefs.h'
      real rval(Nrlib),vval(Nvlib)
      nx=50
      ny=200
      do i=1,Nrlib
         if(i.eq.1) then
            xl=0.
         else
            xl=rval(i-1)
         endif
         xh=rval(i)
         do j=1,Nvlib
            if(j.eq.1) then
               yl=0.
            else
               yl=vval(j-1)
            endif
            yh=vval(j)
            a=0.
            do ix=1,nx
               xp=xl+(xh-xl)*float(ix-1)/float(nx-1)
               do iy=1,ny
                  yp=yl+(yh-yl)*float(iy-1)/float(ny-1)
                  angle=max(0.,yp*pi/180.)
                  x=xp*cos(angle)
                  y=xp*sin(angle)
                  if(x.ge.xlp.and.x.le.xhp.and.y.ge.ylp.and.y.le.yhp)
     $                 a=a+1.
               enddo
            enddo
            ibinp=itobin(i,j)
            areakin(ibin,ibinp)=a/float(nx*ny)
         enddo
      enddo
      return
      end
