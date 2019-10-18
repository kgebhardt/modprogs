C==============================================================================D
C     subroutine ITERWRITE writes out minimization iteration information
C
C     USED BY MODEL
C
C==============================================================================D
      SUBROUTINE iterwrite(iter,Smax,time,avt)
      INCLUDE 'moddefs.h'
      DOUBLE PRECISION Smax

      open (unit=69,file='iter.out',status='unknown')

      write (69,*)' '
C      write (69,*)' # of iterations performed = ',iter-1

C888   format (' maximum value of entropy achieved = ',e12.6)
C      write (69,888)Smax

C      write (69,*)''
C      write (69,*)'  number of comparison bins:'
C      write (69,*)'      spatial bins = ',Nbin
C      write (69,*)'     velocity bins = ',Nslit
C      write (69,*)''

900   format ('            total job time = ',f12.6,$)
901   format ('average time per iteration = ',f12.6,$)
 
      if (time.lt.60.) then
        write ( 6,900)time
        write (69,900)time
        write ( 6,*)' seconds'
        write (69,*)' seconds'
      else if (time.lt.3600.) then
        write ( 6,900)time/60.
        write (69,900)time/60.
        write ( 6,*)' minutes'
        write (69,*)' minutes'
      else
        write ( 6,900)time/3660.
        write (69,900)time/3660.
        write ( 6,*)' hours'
        write (69,*)' hours'
      endif

      if (avt.lt.60.) then
        write ( 6,901)avt
        write (69,901)avt
        write ( 6,*)' seconds'
        write (69,*)' seconds'
      else if (avt.lt.3600.) then
        write ( 6,901)avt/60.
        write (69,901)avt/60.
        write ( 6,*)' minutes'
        write (69,*)' minutes'
      else
        write ( 6,901)avt/3660.
        write (69,901)avt/3660.
        write ( 6,*)' hours'
        write (69,*)' hours'
      endif
 
      close(69)

      RETURN
      END
