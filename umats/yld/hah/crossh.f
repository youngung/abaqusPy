c------------------------------------------------------------------------
      subroutine cross_hardening(kS,ss,gS,emic,ntens,tensor_ref,dgS)
c     Arguments
c     kS         : ks parameter
c     ss         : S parameter
c     gS         : gS parameter
c     emic       : microstructure deviator
c     ntens      : Len of tensor
c     tensor_ref : reference tensor
c     dgS        : incremental gS
      integer ntens
      dimension emic(ntens),tensor_ref(ntens)
      real*8 ks,ss,gS,emic
      real*8 dgS
c     local
      real*8 cos2chi
c     intent(in) kS,ss,gS,emic,ntens,tensor_ref
c     intent(out) dgS

      call calc_cos2chi(tensor_ref,mic,ntens,cos2chi)
      dgS = kS * (1d0 + (ss-1d0) * cos2chi - gS)

      return
      end subroutine cross_hardening
c------------------------------------------------------------------------
