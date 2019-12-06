c$$$      subroutine deviat(cauchy,ntens,sdev)
c$$$      implicit none
c$$$      integer ntens
c$$$      dimension cauchy(ntens), sdev(ntens)
c$$$      real*8 cauchy, sdev
c$$$      if (ntens.eq.6) then
c$$$         sdev(1) = cauchy(1)
c$$$      endif
c$$$      end subroutine deviat
