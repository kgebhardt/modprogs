      subroutine gminmax2(n,x1,x2,xmin,xmax)
      real x1(n),x2(n)
      data big/1.e30/
      xmin=big
      xmax=-big
      do i=1,n
         xmin=min(xmin,x1(i),x2(i))
         xmax=max(xmax,x1(i),x2(i))
      enddo
      if(xmax.eq.xmin) then
         xbit=xmax/10.
      else
         xbit=(xmax-xmin)/10.
      end if
      xmin=xmin-xbit
      xmax=xmax+xbit
      return
      end

      subroutine gminmax(n,x1,xmin,xmax)
      real x1(n)
      data big/1.e30/
      xmin=big
      xmax=-big
      do i=1,n
         xmin=min(xmin,x1(i))
         xmax=max(xmax,x1(i))
      end do
      if(abs((xmax-xmin)/xmax).lt.1.e-3) then
         xbit=xmax/10.
      else
         xbit=(xmax-xmin)/10.
      end if
      xmin=xmin-xbit
      xmax=xmax+xbit
      return
      end

      subroutine gminmax3(n,x1,x2,x3,xmin,xmax)
      real x1(n),x2(n),x3(n)
      data big/1.e30/
      xmin=big
      xmax=-big
      do i=1,n
         xmin=min(xmin,x1(i),x2(i),x3(i))
         xmax=max(xmax,x1(i),x2(i),x3(i))
      enddo
      if(abs((xmax-xmin)/xmax).lt.1.e-4) then
         xbit=xmax/10.
      else
         xbit=(xmax-xmin)/10.
      end if
      xmin=xmin-xbit
      xmax=xmax+xbit
      return
      end

      subroutine styleline(n,x,y,il)
      real x(n),y(n)
      call pgsls(il)
      call pgline(n,x,y)
      call pgsls(1)
      return
      end

      subroutine stylepiece(x,y,z,il,xmin,xmax,ymin,ymax)
      real xl(2),yl(2)
      call pgsls(il)
      xl(1)=xmin+x*(xmax-xmin)
      yl(1)=ymin+y*(ymax-ymin)
      yl(2)=yl(1)
      xl(2)=xmin+(x+z)*(xmax-xmin)
      call pgline(2,xl,yl)
      call pgsls(1)
      return
      end 

      subroutine pointline(n,x,y,n2)
      real x(n),y(n)
      call pgline(n,x,y)
      call pgpoint(n,x,y,n2)
      return
      end

      subroutine colorline(n,x,y,il)
      real x(n),y(n)
      call pgsci(il)
      call pgline(n,x,y)
      call pgsci(1)
      return
      end
 
      subroutine puttext(x,y,xrel,st,xmin,xmax,ymin,ymax)
      character*(*) st
      x1=xmin+x*(xmax-xmin)
      x2=ymin+y*(ymax-ymin)
      call pgptxt(x1,x2,0.,xrel,st)
      return
      end
