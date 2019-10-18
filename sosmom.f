C==============================================================================D
C     subroutine SOSMOM determines whether a surface of section is a
C       single-island or triple-island
C
C     USED BY LIBRARY
C
C==============================================================================D
      SUBROUTINE sosmom(n,xm,x,y,isos)
      DIMENSION x(n),y(n)


      if (isos.eq.0) then
        xm=.999
        goto 999
      endif

      xmin =  1e10
      xmax = -1e10
      ymax = -1e10

      do i=1,isos
        if (x(i).lt.xmin.and.y(i).ne.0.0) then
          xmin = x(i)
          imin = i
        endif
        if (x(i).gt.xmax.and.y(i).ne.0.0) then
          xmax = x(i)
          imax = i
        endif
        if (y(i).gt.ymax) ymax=y(i)
      enddo

      xm2=0.
      xm4=0.

      do i=imin,imax-1
        xUP = (x(i+1)-xmin)/(xmax-xmin)
        xLO = (x(i)-xmin)/(xmax-xmin)
        yUP = y(i+1)/ymax
        yLO = y(i)/ymax
        xm2 = xm2 + (yUP*xUP*xUP + yLO*xLO*xLO ) * (xUP-xLO)/2.
        xm4 = xm4 + (yUP*xUP*xUP*xUP*xUP + yLO*xLO*xLO*xLO*xLO)*
     &    (xUP-xLO)/2.
      enddo

      xm = xm4/xm2

999   RETURN
      END
