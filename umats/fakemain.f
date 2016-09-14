c     Fake main program to test/debug umat
c     This function will be combined with an external umat file called in test.py

      program main
      write(*,*)'Hi, someone testing a umat'
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1     DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2     STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3     PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4     JSTEP(4)

      CMNAME='FAKE_UMAT_TEST'

      NTENS=6
      
      call umat() ## umat will be included..
      call xit()
      end program main
      
      subroutine xit
      write(*,*) 'Hello, Subroutine xit was called.'
      return
      end subroutine xit

      subroutine xit
      write(*,*) 'Hello, Subroutine xit was called.'
      return
      end subroutine xit      
