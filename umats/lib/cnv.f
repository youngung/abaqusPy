c-----------------------------------------------------------------------
c     various converting routines
c-----------------------------------------------------------------------
c     33 bases "stress" tensor to 6 bases stress tensor
c-----------------------------------------------------------------------
      subroutine tens33_trans(a33,b33)
      implicit none
      dimension a33(3,3),b33(3,3)
      real*8 a33,b33
      integer i,j
      b33(:,:) = 0d0
      do 10 i=1,3
      do 10 j=1,3
         b33(i,j) = a33(j,i)
 10   continue
      end subroutine
c-----------------------------------------------------------------------
c     33 bases "stress" tensor to 6 bases "stress" tensor
c-----------------------------------------------------------------------
      subroutine voigt1(a33,a6)
      implicit none
      dimension a33(3,3),a6(6)
      real*8 a33, a6
      integer i
cf2py intent(in) a33
cf2py intent(out) a6
      a6(:) = 0.d0
      do 5 i=1,3
         a6(i) = a33(i,i)
 5    continue
      a6(4) = a33(2,3)
      a6(5) = a33(1,3)
      a6(6) = a33(1,2)
      return
      end subroutine voigt1
c-----------------------------------------------------------------------
c     6 bases "stress" tensor to 33 bases stress tensor
c-----------------------------------------------------------------------
      subroutine voigt2(a6,a33)
      implicit none
      dimension a33(3,3),a6(6)
      real*8 a33, a6
      integer i,j,k,ijv(6,2),n,m
cf2py intent(in)  a6
cf2py intent(out) a33
      data ((ijv(n,m),m=1,2),n=1,6)/1,1,2,2,3,3,2,3,1,3,1,2/
      a33(:,:) = 0.d0
      do 5 k=1,6
         i=ijv(k,1)
         j=ijv(k,2)
         a33(i,j) = a6(k)
         a33(j,i) = a6(k)
 5    continue
      return
      end subroutine voigt2
c-----------------------------------------------------------------------
c     33 bases "strain" tensor to 6 bases strain tensor
c-----------------------------------------------------------------------
c     !! Caution, strain tensor is converted such that
c     the shear strains are multiplied by 2 to follow Abaqus convention
      subroutine voigt3(a33,a6)
      implicit none
      dimension a33(3,3),a6(6)
      real*8 a33, a6
      integer i
cf2py intent(in)  a33
cf2py intent(out) a6
      a6(:) = 0.d0
      do 5 i=1,3
         a6(i) = a33(i,i)
 5    continue
      a6(4) = a33(2,3)*2.
      a6(5) = a33(1,3)*2.
      a6(6) = a33(1,2)*2.
      return
      end subroutine voigt3
c-----------------------------------------------------------------------
c     6 bases "strain" tensor to 33 bases strain tensor
c-----------------------------------------------------------------------
c     !! Caution, strain tensor is converted such that
c     the shear strains are multiplied by 1/2. to follow Abaqus convention
      subroutine voigt4(a6,a33)
      implicit none
      dimension a33(3,3), a6(6)
      real*8 a33, a6, fact
      integer i,j,k,ijv(6,2),n,m
cf2py intent(in)  a6
cf2py intent(out) a33
      data ((ijv(n,m),m=1,2),n=1,6)/1,1,2,2,3,3,2,3,1,3,1,2/
      a33(:,:) = 0.d0
      do 5 k=1,6
         i=ijv(k,1)
         j=ijv(k,2)
         if (i.ne.j) fact=0.5d0
         if (i.eq.j) fact=1.0d0
         a33(i,j) = a6(k)*fact
         a33(j,i) = a6(k)*fact
 5    continue
      return
      end subroutine voigt4
c-----------------------------------------------------------------------
      subroutine reduce_6to3(a6,a3)
      dimension a6(6),a3(3)
      real*8 a6,a3
cf2py intent(in) a6
cf2py intent(out) a3
      a3(1) = a6(1)
      a3(2) = a6(2)
      a3(3) = a6(6)
      return
      end subroutine reduce_6to3
c-----------------------------------------------------------------------
      subroutine reduce_3to6(a3,a6)
      dimension a6(6),a3(3)
      real*8 a6,a3
cf2py intent(in) a3
cf2py intent(out) a6
      a6(1) = a3(1)
      a6(2) = a3(2)
      a6(6) = a3(3)
      return
      end subroutine reduce_3to6
c-----------------------------------------------------------------------
c     Convert 4 dimensional plane-stress tensor
c     - convention:
c      vec  ij component
c        1:      11
c        2:      22
c        3:      33
c        4:      12
      subroutine reduce_6to4(a6,a4)
      dimension a6(6),a4(4)
      real*8 a6,a4
      a4(1) = a6(1)
      a4(2) = a6(2)
      a4(3) = 0d0
      a4(4) = a6(6)
      return
      end subroutine reduce_6to4
c-----------------------------------------------------------------------
      subroutine reduce_4to6(a4,a6)
      dimension a6(6),a4(4)
      real*8 a6,a4
      a6(1) = a4(1)
      a6(2) = a4(2)
      a6(3) = 0d0
      a6(4) = 0d0
      a6(5) = 0d0
      a6(6) = a4(4)
      return
      end subroutine reduce_4to6
c-----------------------------------------------------------------------
