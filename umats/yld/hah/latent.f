c------------------------------------------------------------------------
c     Latent hardening subroutines to be wrapped by f2py.
c     This module should be also used in Abaqus UMAT (under development...)
c     Find more information of Abaqus UMAT in
c       www.github.com/youngung/abaqusPy
c
c     General references
c     [1] Barlat et al. IJP 58, 2014 p201-218
c     [2] Jeong et al., IJP, 2016 (in press)

c     Youngung Jeong
c     youngung.jeong@gmail.com
c------------------------------------------------------------------------
c     calculate incremental gL (dgL)
      subroutine latent(eL,gL,ekL,ys_iso,ys_hah,emic,target,ntens,
     $     debar,dgL)
c     Arguments
c     eL     : L parameter
c     gL     : gL parameter at step n
c     ekL    : kL constant
c     ys_iso : \bar{\sigma}(0) - the initial yield surface (without HAH effect)
c     ys_hah : \bar{\sigma}(\bar{\varepsilon}) - current yield surface
c     emic   : microstructure deviator
c     target : the target directino towards which the microstructure
c               deviator aims to realign.
c     debar  : incremental equivalent strain
c               (used as multiplier when updating the state variables...)
c     dgL    : incremental gL
      implicit none
c     Arguments passed in
      integer ntens
      dimension emic(ntens)
      dimension target(ntens)
      real*8 eL,gL,ekL,ys_iso,ys_hah,emic,target,debar,dgL
c     local
      real*8 cos2chi,term

cf2py intent(in) eL,gL,ekL,ys_iso,ys_hah,emic,target,debar
cf2py intent(out) dgL

      call calc_cos2chi(emic,target,ntens,cos2chi)

c     Eq 16 --
      term = dsqrt(eL * (1d0-cos2chi) + cos2chi)-1d0
      dgL = ekL *( (ys_hah-ys_iso) / ys_hah * term  + 1d0 - gL )
      dgL = dgL * debar

      return
      end subroutine latent
