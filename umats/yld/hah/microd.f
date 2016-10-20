c-----------------------------------------------------------------------
c     Microstructure deviator

c     Reference
c     [1] Manual of abaqusPy

c     Youngung Jeong
c     youngung.jeong@gmail.com
c-----------------------------------------------------------------------
      subroutine micro_dev(emic,target)
c     Arguments
c     emic  : microstructure deviator
c     target: the target with which emic tries to realign
      implicit none
      dimension emic(6),target(6)
      real*8, intent(in) ::  emic, target
      real*8 dh, cos_chi2

      return
      end subroutine
c-----------------------------------------------------------------------
      subroutine micro_dev_derv(ntens,emic,emod,dphi)
c     Arguments
c     emic  : microstructure deviator
c     emod  : elastic modulus
c     dphi  : derivative of yield surface
      implicit none
      integer, intent(in) :: ntens
      dimension emic(ntens), emod(ntens,ntens), dphi(ntens)
      real*8 emic, emod, dphi
      end subroutine
