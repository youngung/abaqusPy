c------------------------------------------------------------------------
c     References
c     [1] Barlat et al. IJP 58, 2014 p201-218
c------------------------------------------------------------------------
      subroutine cross_hardening(ntens,emic,target,c_kS,ss,
     $     gS,dgS_deps)
c     Arguments
c     ntens      : Len of tensor
c     emic       : microstructure deviator
c     target     : reference tensor (does not have to be a hat tensor)
c     c_kS       : ks parameter
c     ss         : S parameter
c     gS         : gS parameter
c     dgS_deps   : dgs/deps
      integer, intent(in) :: ntens
      dimension emic(ntens),target(ntens)
      real*8, intent(in)  :: c_kS,ss,gS,emic
      real*8, intent(out) :: dgS_deps
c     local
      real*8 coschi
c     intent(in) c_kS,ss,gS,emic,ntens,target
c     intent(out) dgS_deps

      call calc_coschi(ntens,target,emic,coschi)
c     Equation 29 in Ref [1]
      dgS_deps = c_kS * (1d0 + (ss-1d0) * (coschi*coschi) - gS)
      return
      end subroutine cross_hardening
c------------------------------------------------------------------------
