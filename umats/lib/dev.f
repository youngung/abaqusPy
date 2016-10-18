c     Given stress s, calculate deviator and hydrostatic pressure
c-----------------------------------------------------------------------
      subroutine deviat_inv(ntens,cauchy,sdev,p)
c     Arguments
c     ntens  : Len of cauchy, sdev
c     cauchy : cauchy stress
c     sdev   : deviatoric stress
c     p      : hydrostatic pressure
      implicit none
      integer ntens
      dimension cauchy(ntens),sdev(ntens)
      real*8 cauchy,sdev,p
cf2py intent(in) ntens, sdev,p
cf2py intent(out) cauchy
cf2py depend(ntens) sdev,cauchy
      if (ntens.eq.3) then
         cauchy(1) = sdev(1) + p
         cauchy(2) = sdev(2) + p
         cauchy(3) = sdev(3)
      else
         write(*,*) 'deviat_invert is not ready for this dimension'
         call exit(-1)
      endif
      return
      end subroutine deviat_inv
c-----------------------------------------------------------------------
c     Calculate deviatoric part of cauchy tensor
      subroutine deviat(ntens,cauchy,sdev,p)
c     Arguments
c     ntens   : Len of <cauchy> and <sdev>
c     cauchy  : cauchy stress (can be other tensors of interest)
c     sdev    : deviatoric stress
c     p       : hydrostatic pressure
      implicit none
      integer, intent(in):: ntens
      dimension cauchy(ntens),sdev(ntens)
      real*8, intent(in)  :: cauchy
      real*8, intent(out) :: sdev,p
      if (ntens.eq.3) then
         call deviat3(cauchy,sdev,p)
      elseif(ntens.eq.6) then
         call deviat6(cauchy,sdev,p)
      else
         write(*,*) 'unexpected case given to deviat'
         call exit(-1)
      endif
      return
      end subroutine deviat
c-----------------------------------------------------------------------
c     Calculate deviator and hydrostatic pressure
      subroutine deviat6(s,sd,p)
c     Arguments
c     s   : cauchy stress under plane stress space (s11,s22,s12 with s33=0)
c     sd  : stress deviator
c     p   : hydrostatic pressure
      implicit none
      real*8, intent(in) ::s(6)
      real*8, intent(out)::sd(6),p
      p = (s(1)+s(2)+s(3))/3d0
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
c     Arguments
c     s   : cauchy stress under plane stress space (s11,s22,s12 with s33=0)
c     sd  : stress deviator
c     p   : hydrostatic pressure
      implicit none
      real*8, intent(in) ::s(3)
      real*8, intent(out)::sd(3),p
      p = (s(1)+s(2))/3d0
      sd(1) = s(1)-p
      sd(2) = s(2)-p
      sd(3) = s(3)              !! shear component s12
      return
      end subroutine deviat3
c-----------------------------------------------------------------------
c     Given array a33, find deviatoric part and its trace
      subroutine deviat33(a,ad,amean)
c     a    : cauchy stress under plane stress space (s11,s22,s12 with s33=0)
c     ad   : stress deviator
c     amean: hydrostatic pressure
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
