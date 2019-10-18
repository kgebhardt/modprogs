
      parameter(nmax=10000,nang=5)
      real r(nmax),d(nmax),s(nmax),ang(nang)
      parameter(pi=3.1415926539)
      character filen*40
      external funci
      common/cfunci/axra,axrt

 1    write(*,"('Datafile : '$)")
      read *,filen
      open(unit=1,file=filen,status='old',err=1)
      write(*,"('apparent and true axis ratio : '$)")
      read *,axra,axrt
      fac=axra/axrt

      if(axra.eq.axrt) then
         xinc=pi/2.
      else
         xinc=rtbiso(funci,0.,pi/2.,1.e-5)
      endif
      print *,'inclination = ',xinc*180./pi

      n=0
      do i=1,nmax
         read(1,*,end=666) x1,x2,x3
         n=n+1
         r(n)=x1
         d(n)=x2*fac
         s(n)=x3
      enddo
 666  continue
      close(1)
      
      as=0.
      ae=pi/2.
      do i=1,nang
         ang(i)=as+(ae-as)*float(i-1)/float(nang-1)
      enddo

      open(unit=1,file='gden.out',status='unknown')

      do i=1,n
         do j=1,nang
            x=r(i)*cos(ang(j))
            y=r(i)*sin(ang(j))
            rnew=sqrt(x*x+y*y/axrt/axrt)
            call xlinint(rnew,n,r,d,dn)
            call xlinint(rnew,n,r,s,sn)
            write(1,*) r(i),ang(j),dn,sn
         enddo
      enddo

      close(1)
      end

      subroutine xlinint(xp,n,x,y,yp)
      real x(n),y(n)
      do j=1,n-1
         if(xp.ge.x(j).and.xp.lt.x(j+1)) then
            yp=y(j)+(y(j+1)-y(j))*(xp-x(j))/(x(j+1)-x(j))
            return
         endif
      enddo
      if(xp.le.x(1)) yp=y(1)
      if(xp.ge.x(n)) yp=y(n)
      return
      end

      function funci(x)
      common/cfunci/axra,axrt
      funci=cos(x)**2+axrt*axrt*sin(x)**2-axra*axra
      return
      end

      FUNCTION rtbiso(func,x1,x2,xacc)
      INTEGER JMAX
      REAL rtbiso,x1,x2,xacc,func
      EXTERNAL func
      PARAMETER (JMAX=40)
      INTEGER j
      REAL dx,f,fmid,xmid
      fmid=func(x2)
      f=func(x1)
      if(f*fmid.ge.0.) print *,'root must be bracketed in rtbis'
      if(f.lt.0.)then
        rtbiso=x1
        dx=x2-x1
      else
        rtbiso=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*.5
        xmid=rtbiso+dx
        fmid=func(xmid)
        if(fmid.le.0.)rtbiso=xmid
        if(abs(dx).lt.xacc .or. fmid.eq.0.) return
11    continue
      print *,'too many bisections in rtbis'
      END
