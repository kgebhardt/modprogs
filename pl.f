C==============================================================================D
C     function PL calculates the Legendre polynomial of order n
C
C     USED BY LIBRARY
C
C==============================================================================D
      FUNCTION Pl(n,x)
      if(n.eq.0) then
         pl=1.
         return
      elseif(n.eq.1) then
         pl=x
         return
      elseif(n.eq.2) then
         pl=(3.*x*x-1.)/2.
         return
      elseif(n.eq.3) then
         pl=(5.*x*x*x-3.*x)/2.
         return
      elseif(n.eq.4) then
         pl=(35.*x*x*x*x-30.*x*x+3.)/8.
         return
      elseif(n.eq.5) then
         pl=(63./8.*x*x*x*x*x-35./4.*x*x*x+15./8.*x)
         return
      elseif(n.eq.6) then
         pl=231./16.*x*x*x*x*x*x-315./16.*x*x*x*x+105./16.*x*x-5./16.
         return
      elseif(n.eq.7) then
         pl=429./16.*x*x*x*x*x*x*x-693./16.*x*x*x*x*x+315./16.*x*x*x-
     &        35./16.*x
         return
      elseif(n.eq.8) then
         pl=6435./128.*x*x*x*x*x*x*x*x-3003./32.*x*x*x*x*x*x+3465./
     &        64.*x*x*x*x-315./32.*x*x+35./128.
         return
      elseif(n.eq.9) then
         pl=12155./128.*x**9-6435./32.*x**7+9009./64.*x**5-1155./
     &        32.*x**3+315./128.*x
         return
      elseif(n.eq.10) then
         pl=46189./256.*x**10-109395/256.*x**8+45045./128.*x**6-
     &        15015./128.*x**4+3465./256.*x**2-63./256.
         return
      elseif(n.eq.11) then
         pl=88179./256.*x**11-230945/256.*x**9+109395./128.*x**7-
     &        45045./128.*x**5+15015./256.*x**3-693./256.*x
         return
      elseif(n.eq.12) then
         pl=676039./1024.*x**12-969969./512.*x**10+2078505./1024.*x**8-
     &        255255./256.*x**6+225225./1024.*x**4-9009./512.*x**2+
     &        231./1024.
         return
      else
         print *,'not: -1 < n < 13 in pl!'
      endif
      return
      end
