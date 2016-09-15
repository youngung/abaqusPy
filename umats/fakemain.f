c-----------------------------------------------------------------------
c     fakemain.f developed to help abaqus UMAT test/debug.
c     Fake main program to test/debug umat
c     This function will be combined with an external
c     umat file called in test.py
c------------------------------------------------------------------------
      program main
      CHARACTER*80 CMNAME
      integer i,j, kinc, npt, kspt
      parameter(NTENS=6)
      parameter(nprops=2)
      parameter(nstatv=2)
      parameter(ndi=3)
      parameter(nshr=3)
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1     DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2     STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3     PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4     JSTEP(4)

      real*8 celent, dtime,temp,dtemp
      integer ndi,nshr,ntens

c     Giving a fake name
c      CMNAME(1:15)='FAKE_UMAT_TEST'
c     Testing for elements/conditions suitable for 6d stress states

c     Initial values
      do i=1,ntens
         stress(i)=0.d0
         stran(i)=0.d0
         dstran(i)=0.d0
         ddsddt(i)=0.d0
         drplde(i)=0.d0
         do j=1,ntens
            ddsdde(i,j)=0.d0
         enddo
      enddo

      predef(1)=0.
      dpred(1)=0.
      do i=1,3
         coords(i)=0.
         do j=1,3
            drot(i,j) = 0.    ! rotation increment
            dfgrd0(i,j) = 0.  ! deformation gradient at the beinning of increment
            dfgrd1(i,j) = 0.  ! deformation gradient at the end of increment
         enddo
      enddo

      jstep(1) = 1 !step number
      jstep(2) = 2 !procedure key
      jstep(3) = 0 ! 1 if NLGEOM=yes 0: otherwise
      jstep(4) = 1 ! 1 if current step is a linear perturbation procedure; 0 otherwise.


      do i=1,nprops
         props(i) = 0.d0
      enddo

      kinc=1 ! increment number
      npt =1 ! integration point number
      noel=1 ! element number
      layer=1 ! layer number
      kspt=1 ! section point number iwthin the current layer
      celent = 1.d-3
      dtime=1.d-3
      temp=273.
      dtemp=0.
      time(1)=0.
      time(2)=dtime

c     Calling umat
      call umat(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1     RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,
     2     DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3     NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4     CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
c     Calling xit
c     call xit()
      end program main

c-----------------------------------------------------------------------
c     Fake xit function
c-----------------------------------------------------------------------
      subroutine xit
      write(*,*) 'Hello, Subroutine xit was called.'
      return
      end subroutine xit

c-----------------------------------------------------------------------
c     Fake rotsig function
c-----------------------------------------------------------------------
      subroutine rotsig(s,r,sprime,lstr,ndi,nshr)
      dimension s(ndi+nshr),r(3,3),sprime(ndi+nshr)
      integer lstr,ndi,nshr

      return
      end subroutine rotsig
