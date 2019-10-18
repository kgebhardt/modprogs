C=============================================================================D
C     Press et. al. routine -- modified
C
C     USED BY LIBRARY
C
C=============================================================================D
      SUBROUTINE RK4(N,YSTART,DYDX,H,Y,DERIVS)
      DIMENSION YSTART(N),DYDX(N),Y(N)
      DIMENSION YT(4),DYT(4),DYM(4)
C     DIMENSION YT(N),DYT(N),DYM(N)
      HH=H*0.5
      H6=H/6.
      DO 11 I=1,N
        YT(I)=YSTART(I)+HH*DYDX(I)
11    CONTINUE
      CALL DERIVS(YT,DYT)
      DO 12 I=1,N
        YT(I)=YSTART(I)+HH*DYT(I)
12    CONTINUE
      CALL DERIVS(YT,DYM)
      DO 13 I=1,N
        YT(I)=YSTART(I)+H*DYM(I)
        DYM(I)=DYT(I)+DYM(I)
13    CONTINUE
      CALL DERIVS(YT,DYT)
      DO 14 I=1,N
        Y(I)=YSTART(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
14    CONTINUE
      RETURN
      END
