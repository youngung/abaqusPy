c-----------------------------------------------------------------------
!     various testing subroutines
      logical function is_inf(val)
c     val: value subjected to test
      use, intrinsic :: ieee_arithmetic
      implicit none
      real*8 val
      is_inf = .not.(ieee_is_finite(val))
      return
      end function is_inf
c-----------------------------------------------------------------------
!     various testing subroutines
      logical function is_finite(val)
c     val: value subjected to test
      use, intrinsic :: ieee_arithmetic
      implicit none
      real*8 val
      is_finite = ieee_is_finite(val)
      return
      end function is_finite
c$$$      program main
c$$$      implicit none
c$$$      real*8 val
c$$$      logical is_inf
c$$$      val = 1.0
c$$$      write(*,*) 'result:',is_inf(val)
c$$$      end program
