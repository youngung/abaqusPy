c     various converting routines
c-----------------------------------------------------------------------
c     a33 to a6
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
c     a6 to a33
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
