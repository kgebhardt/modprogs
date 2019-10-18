      program plotvoronoi
      INCLUDE 'libdefs.h'
      parameter(nmax=15000)
      integer iseq
      character seqlabel*3
      integer from(nmax),to(nmax),nuc_type(nmax)
      real node_x(nmax),node_y(nmax),nuc_x(nmax),nuc_y(nmax)
      real x(nmax),y(nmax),xl(2),yl(2),yshift,xshift
      real env_x(nmax),env_y(nmax)

      call galaxyread()

      distance=distance/1.e3

      print*,'iseq [9999=integrals]'
      read(*,*) iseq

      if(iseq.lt.10) then
         write(seqlabel,'("00",i1)') iseq
      else
         if(iseq.lt.100) then
            write(seqlabel,'("0",i2)') iseq
         else
            write(seqlabel,'(i3)') iseq
         end if
      end if
      if(iseq.eq.9999) seqlabel="elz"

c --- read edges
      open(unit=12,file=seqlabel//'.vor.edge',status='old')
      iedge=0
      do i=1,nmax
         read(12,*,end=700) i1,i2,i3
         iedge=iedge+1
         from(iedge)=i2
         to(iedge)=i3
      end do
 700  continue
      close(12)
c --- end

c --- read nodes
      open(unit=12,file=seqlabel//'.vor.node',status='old')
      inode=0
      do i=1,nmax
         read(12,*,end=701) i1,x2,x3
         inode=inode+1
         node_x(inode)=x2
         node_y(inode)=x3
      end do
 701  continue
      close(12)
c --- end      

c --- read nuclei and plot info
      open(unit=12,file=seqlabel//'.nuclei',status='old')
      open(unit=13,file=seqlabel//'.vor.plot',status='old')
      inuclei=0
      do i=1,nmax
         read(12,*,end=702) x1,x2,x3,i4
         inuclei=inuclei+1
         nuc_x(inuclei)=x1
         nuc_y(inuclei)=x2
         if(iseq.eq.9999) then
            read(13,*,end=702) x1,x2,i3,x4,x5
            yshift=x4
            xshift=x5
            nuc_type(inuclei)=i3
         else
            read(13,*,end=702) x1,x2,i3,x4
            nuc_type(inuclei)=i3
            yshift=x4
         end if
      end do
 702  continue
      close(12)
      close(13)
c --- end   


c --- read factors
      open(unit=12,file=seqlabel//'.vor.sos',status='old')
      do i=1,nmax
         read(12,*,end=703) x1,x2,x3,i1,x4,xfac,yfac
      end do
 703  continue
      close(12)
c --- end 


c --- non-logarithimcally
      call pgbegin(0,'/cps',1,1)
      call pgsch(1.5)
      call gminmax(inuclei,nuc_x,xmin,xmax)
      call gminmax(inuclei,nuc_y,ymin,ymax)
      if(iseq.eq.9999) then
         xmin=(xmin-0.025)*xfac-abs(xshift)
         xmax=(xmax-0.025)*xfac-abs(xshift)
         ymin=(ymin-0.025)*yfac-1.1*abs(yshift)
         ymax=(ymax-0.025)*yfac-1.1*abs(yshift)
      else
         xmin=(xmin)*xfac
         xmax=(xmax)*xfac
         ymin=(ymin-0.025)*yfac-1.1*abs(yshift)
         ymax=(ymax-0.025)*yfac-1.1*abs(yshift)
      end if
      xmin=xmin*distance/arcsec
      xmax=xmax*distance/arcsec
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      if(iseq.ne.9999) then
         call pglabel('r [kpc]','v\\Dr\\U [km/s]','')
      else
         call pglabel('log(-E)','log(L\\Dz\\U)','')
      end if
      call pgsch(1.0)

      !nuclei
      in=0
      do i=1,inuclei
         if(nuc_type(i).ne.0) then
            in=in+1
            if(iseq.eq.9999) then
               x(in)=(nuc_x(i)-0.025)*xfac-abs(xshift)
               y(in)=(nuc_y(i)-0.025)*yfac-1.1*abs(yshift)
            else
               x(in)=(nuc_x(i))*xfac
               y(in)=(nuc_y(i)-0.025)*yfac-1.1*abs(yshift)
            end if
            x(in)=x(in)*distance/arcsec
         end if
      end do
      call pgsci(2)
      call pgsch(0.7)
      call pgpoint(in,x,y,4)
      call pgsch(1.0)
      call pgsci(1)

      !border
      in=0
      do i=1,inuclei
         if(nuc_type(i).eq.0) then
            in=in+1
            if(iseq.eq.9999) then
               x(in)=(nuc_x(i)-0.025)*xfac-abs(xshift)
               y(in)=(nuc_y(i)-0.025)*yfac-1.1*abs(yshift)
            else
               x(in)=(nuc_x(i))*xfac
               y(in)=(nuc_y(i)-0.025)*yfac-1.1*abs(yshift)
            end if
            x(in)=x(in)*distance/arcsec
         end if
      end do
      call pgsci(4)
      call pgsch(0.6)
      call pgpoint(in,x,y,-32)
      call pgsch(1.0)
      call pgsci(1)

      !nodes
      do i=1,iedge
         if(from(i).ne.-1.and.to(i).ne.-1) then
            if(iseq.eq.9999) then
               xl(1)=(node_x(from(i)+1)-0.025)*xfac-abs(xshift)
               xl(2)=(node_x(to(i)+1)-0.025)*xfac-abs(xshift)
               yl(1)=(node_y(from(i)+1)-0.025)*yfac-1.1*abs(yshift)
               yl(2)=(node_y(to(i)+1)-0.025)*yfac-1.1*abs(yshift)
            else
               xl(1)=(node_x(from(i)+1))*xfac
               xl(2)=(node_x(to(i)+1))*xfac
               yl(1)=(node_y(from(i)+1)-0.025)*yfac-1.1*abs(yshift)
               yl(2)=(node_y(to(i)+1)-0.025)*yfac-1.1*abs(yshift)
            end if
            xl(1)=xl(1)*distance/arcsec
            xl(2)=xl(2)*distance/arcsec
         end if
         iok=1
         if(iok.eq.1) then
            call pgline(2,xl,yl)
         end if
      end do
c --- end 

c --- non-logarithimcally
      call gminmax(inuclei,nuc_x,xmin,xmax)
      call gminmax(inuclei,nuc_y,ymin,ymax)
      if(iseq.eq.9999) then
         xmin=(xmin-0.025)*xfac-abs(xshift)
         xmax=(xmax-0.025)*xfac-abs(xshift)
         ymin=(ymin-0.025)*yfac-1.1*abs(yshift)
         ymax=(ymax-0.025)*yfac-1.1*abs(yshift)
      else
         xmin=(xmin)*xfac
         xmax=(xmax)*xfac
         ymin=(ymin-0.025)*yfac-1.1*abs(yshift)
         ymax=(ymax-0.025)*yfac-1.1*abs(yshift)
      end if
      xmin=xmin*distance/arcsec
      xmax=xmax*distance/arcsec
      call pgenv(xmin,xmax,ymin,ymax,0,0)
      if(iseq.ne.9999) then
         call pglabel('r [kpc]','v\\Dr\\U [km/s]','')
      else
         call pglabel('log(-E)','log(L\\Dz\\U)','')
      end if
      call pgsch(1.0)

      !nuclei
      in=0
      do i=1,inuclei
         if(nuc_type(i).ne.0) then
            in=in+1
            if(iseq.eq.9999) then
               x(in)=(nuc_x(i)-0.025)*xfac-abs(xshift)
               y(in)=(nuc_y(i)-0.025)*yfac-1.1*abs(yshift)
            else
               x(in)=(nuc_x(i))*xfac
               y(in)=(nuc_y(i)-0.025)*yfac-1.1*abs(yshift)
            end if
            x(in)=x(in)*distance/arcsec
         end if
      end do
      call pgsci(2)
      call pgsch(0.7)
      call pgpoint(in,x,y,4)
      call pgsch(1.0)
      call pgsci(1)

      !border
c      in=0
c      do i=1,inuclei
c         if(nuc_type(i).eq.0) then
c            in=in+1
c            if(iseq.eq.9999) then
c               x(in)=(nuc_x(i)-0.025)*xfac-abs(xshift)
c               y(in)=(nuc_y(i)-0.025)*yfac-1.1*abs(yshift)
c            else
c               x(in)=(nuc_x(i))*xfac
c               y(in)=(nuc_y(i)-0.025)*yfac-1.1*abs(yshift)
c            end if
c         end if
c         x(in)=x(in)*distance/arcsec
c      end do
c      call pgsci(4)
c      call pgsch(0.6)
c      call pgpoint(in,x,y,-32)
c      call pgsch(1.0)
c      call pgsci(1)

      !nodes
      do i=1,iedge
         if(from(i).ne.-1.and.to(i).ne.-1) then
            if(iseq.eq.9999) then
               xl(1)=(node_x(from(i)+1)-0.025)*xfac-abs(xshift)
               xl(2)=(node_x(to(i)+1)-0.025)*xfac-abs(xshift)
               yl(1)=(node_y(from(i)+1)-0.025)*yfac-1.1*abs(yshift)
               yl(2)=(node_y(to(i)+1)-0.025)*yfac-1.1*abs(yshift)
            else
               xl(1)=(node_x(from(i)+1))*xfac
               xl(2)=(node_x(to(i)+1))*xfac
               yl(1)=(node_y(from(i)+1)-0.025)*yfac-1.1*abs(yshift)
               yl(2)=(node_y(to(i)+1)-0.025)*yfac-1.1*abs(yshift)
            end if
            xl(1)=xl(1)*distance/arcsec
            xl(2)=xl(2)*distance/arcsec
         end if
         iok=1
         if(iok.eq.1) then
            call pgline(2,xl,yl)
         end if
      end do
c --- end 

      call pgend

      end





