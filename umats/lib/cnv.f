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
      do 10 i=1,3
      do 10 j=1,3
         b(i,j) = a(j,i)
 10   continue
      end subroutine
c-----------------------------------------------------------------------
c     33 bases "stress" tensor to 6 bases "stress" tensor
c-----------------------------------------------------------------------
      subroutine voigt1(a33,a6)
      implicit none
      real*8 a33(3,3), a6(6)
      integer i
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
      real*8 a33(3,3), a6(6)
      integer i,j,k,ijv(6,2),n,m
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
      real*8 a33(3,3), a6(6)
      integer i
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
      real*8 a33(3,3), a6(6), fact
      integer i,j,k,ijv(6,2),n,m
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
