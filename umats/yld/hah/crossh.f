c------------------------------------------------------------------------
      subroutine cross_hardening(ntens,ndi,nshr,emic,tensor_ref,kS,ss,
     $     gS,dgS)
c     Arguments
c     kS         : ks parameter
c     ss         : S parameter
c     gS         : gS parameter
c     emic       : microstructure deviator
c     ntens      : Len of tensor
c     ndi        : Number of normal components
c     nshr       : Number of shear components
c     tensor_ref : reference tensor
c     dgS        : incremental gS
      integer, intent(in) :: ntens,ndi,nshr
      dimension emic(ntens),tensor_ref(ntens)
      real*8, intent(in)  :: kS,ss,gS,emic
      real*8, intent(out) :: dgS
c     local
      real*8 coschi
c     intent(in) kS,ss,gS,emic,ntens,tensor_ref
c     intent(out) dgS

      call calc_coschi(ntens,ndi,nshr,tensor_ref,emic,coschi)
      dgS = kS * (1d0 + (ss-1d0) * (coschi*coschi) - gS)
      return
      end subroutine cross_hardening
c------------------------------------------------------------------------
