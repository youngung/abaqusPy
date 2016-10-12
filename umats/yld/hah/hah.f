c-----------------------------------------------------------------------
c     Homogeneous Anisotropic Hardening

c     subroutine hah has several dependents
c     read_alpha: yld2000_2d.f
c     yld2000_2d: yld2000_2d.f
c     hah_yieldsurface: yld_yieldsurface.f
c-----------------------------------------------------------------------
      subroutine hah(iyld_choice,stress,phi,dphi,d2phi,
     $     yldc,yldp,nyldc,nyldp,ntens)
c     Arguments
c     iyld_choice - isotropic yield surface kernel
c     stress      - stress tensor
c     phi         - HAH yield surface
c     dphi        - HAH yield surface 1st derivatives
c     d2phi       - HAH yield surface 2nd derivatives
c     yldc        - isotropic yield surface constants
c     yldp        - HAH yield surface contants/state variables
c     nyldc       - Len of yldc
c     nyldp       - Len of yldp
c     ntens       - Len of stress tensor
      implicit none
      integer ntens,nyldc,nyldp
      dimension yldc(nyldc),yldp(nyldp),stress(ntens)
      dimension dphi(ntens),d2phi(ntens,ntens)
      real*8 yldc,yldp,stress,phi,dphi,d2phi
      dimension dphi_chi(ntens),d2phi_chi(ntens,ntens)
      real*8 phi_chi,dphi_chi,d2phi_chi
      integer iyld_choice,i,j
c**   iyld_choice=2             ! yld2000-2d
c     dummy alphas for testing purpose

      dphi_chi(:)=0d0
      d2phi_chi(:,:)=0d0

      call read_alpha(
     $     '/home/younguj/repo/abaqusPy/umats/yld/alfas.txt',yldc)

c$$$      if (iyld_choice.eq.2) then
c$$$c**      phi_chi, dphi_chi, d2phi_chi
c$$$         call yld2000_2d(stress,phi_chi,dphi_chi,d2phi_chi,yldc)
c$$$      else
c$$$         write(*,*) 'unexpected iyld_choice'
c$$$         stop -1
c$$$      endif
c$$$
c$$$      call hah_yieldsurface(iyld_choice,yldc,nyldc,yldp,nyldp,stress,
c$$$     $     phi_chi,dphi_chi,d2phi_chi,ntens,phi,dphi,d2phi)

      return
      end subroutine hah
c      include '/home/younguj/repo/abaqusPy/umats/yld/yld.f'
c      include '/home/younguj/repo/abaqusPy/umats/yld/hah/hah_lib.f'
