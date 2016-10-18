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
      real*8 ks,ss,gS,emic
      real*8 dgS
c     local
      real*8 cos2chi
c     intent(in) kS,ss,gS,emic,ntens,tensor_ref
c     intent(out) dgS

      call calc_cos2chi(ntens,ndi,nshr,tensor_ref,mic,cos2chi)
      dgS = kS * (1d0 + (ss-1d0) * cos2chi - gS)

      return
      end subroutine cross_hardening
c------------------------------------------------------------------------
      subroutine crossh(gS,so,phi,phi_chi,phi_x)
c     Arguments
c     kS         : ks parameter
c     ss         : S parameter
c     gS         : gS parameter
c     emic       : microstructure deviator
c     ntens      : Len of tensor
c     tensor_ref : reference tensor
c     dgS        : incremental gS
      implicit none
      integer ntens
      dimension emic(ntens),tensor_ref(ntens)
      real*8 ks,ss,gS,emic
      real*8 dgS
c     local
      real*8 cos2chi


      return
      end subroutine crossh
