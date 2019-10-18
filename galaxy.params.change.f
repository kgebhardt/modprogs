      program galaxy_params_change

      character(len=14) junk
      character(len=40) galname
      real junk2,faktor
      integer ijunk

      open(unit=11,file='gden.norm')
      read(11,*) faktor
      close(11)

      open(unit=2,file='galaxy.params',status='old')
      open(unit=3,file='galaxy.params.new',status='unknown')
      read(2,*) galname
      write(3,'(a)') galname
      read(2,*) junk,junk2
      write(3,'(a,f6.2)') junk,junk2
      read(2,*) junk,junk2
      write(3,'(a,f7.2)') junk,junk2
      read(2,*) junk,junk2
      junk2=junk2/sqrt(faktor)
      write(3,'(a,f11.2)') junk,junk2
      read(2,*) junk,junk2
      junk2=junk2/sqrt(faktor)
      write(3,'(a,f10.2)') junk,junk2
      read(2,*) junk,junk2
      write(3,'(a,f6.2)') junk,junk2
      read(2,*) junk,junk2
      write(3,'(a,f7.2)') junk,junk2
      read(2,*) junk,ijunk
      write(3,'(a,i1)') junk,ijunk

      close(2)
      close(3)

      end program
