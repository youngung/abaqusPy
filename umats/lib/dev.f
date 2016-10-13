c     Given stress s, calculate deviator and hydrostatic pressure
c-----------------------------------------------------------------------
      subroutine deviat(cauchy,ntens,sdev)
      implicit none
      integer ntens
      dimension cauchy(ntens),sdev(ntens)
      real*8 cauchy,sdev,p
      if (ntens.eq.3) then
         call deviat3(cauchy,sdev,p)
      endif
      end subroutine deviat
c-----------------------------------------------------------------------
      subroutine deviat6(s,sd,p)
      implicit none
      real*8 s(6),sd(6),p
      p = (s(1)+s(2)+s(3))/3.
      sd(1) = s(1)-p
      sd(2) = s(2)-p
      sd(3) = s(3)-p
      sd(4) = s(4)
      sd(5) = s(5)
      sd(6) = s(6)
      return
      end subroutine deviat6
c-----------------------------------------------------------------------
c     Given stress s3; plane stress condition with s(3) = 0.
      subroutine deviat3(s,sd,p)
      implicit none
      real*8 s(3),sd(3),p
      p = (s(1)+s(2))/3.
      sd(1) = s(1)-p
      sd(2) = s(2)-p
      sd(3) = s(3)              !! shear component s12
      return
      end subroutine deviat3
c-----------------------------------------------------------------------
c     Given array a33, find deviatoric part and its trace
      subroutine deviat33(a,ad,amean)
      implicit none
      real*8 a(3,3),ad6(6),ad(3,3),amean,a6(6)
      call voigt1(a,a6)
      call deviat6(a6,ad6,amean)
      call voigt2(ad6,ad)
      return
      end subroutine
c     on Pal
c$$$  include "/home/younguj/repo/abaqusPy/umats/lib/cnv.f"
c     on Mac
c      include "/Users/yj/repo/abaqusPy/umats/lib/cnv.f"
