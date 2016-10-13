c-----------------------------------------------------------------------
c     Homogeneous Anisotropic Hardening

c     subroutine hah has several dependents
c     read_alpha: yld2000_2d.f
c     yld2000_2d: yld2000_2d.f
c     hah_yieldsurface: yld_yieldsurface.f
c-----------------------------------------------------------------------
      subroutine hah(iyld_choice,cauchy,phi,dphi,d2phi,
     $     yldc,yldp,nyldc,nyldp,ntens)
c     Arguments
c     iyld_choice - isotropic yield surface kernel
c     cauchy      - cauchy stress tensor
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
      dimension yldc(nyldc),yldp(nyldp),cauchy(ntens)
      dimension dphi(ntens),d2phi(ntens,ntens)
      real*8 yldc,yldp,cauchy,phi,dphi,d2phi
      dimension dphi_chi(ntens),d2phi_chi(ntens,ntens)
      real*8 phi_chi,dphi_chi,d2phi_chi
      integer iyld_choice,i,j,imsg
      imsg=0

      dphi_chi(:)=0d0
      d2phi_chi(:,:)=0d0

      call w_chr(imsg,'cauchy stress passed to subroutine hah')
      call w_dim(imsg,cauchy,ntens,1d0,.false.)
c      call exit(-1)

c$$$      call read_alpha(
c$$$     $     '/home/younguj/repo/abaqusPy/umats/yld/alfas.txt',yldc)
      call w_chr(imsg,'In hah.f')
      call w_ival(imsg,'iyld_choice:',iyld_choice)

      if (iyld_choice.eq.2) then ! plane-stress condition
         if (ntens.ne.3) then
            write(*,*)'When iyld_choice.eq.2, it should be',
     $           ' plane stress condition with ntens=3'
            call exit(-1)
         endif
c**      phi_chi, dphi_chi, d2phi_chi
         call yld2000_2d(cauchy,phi_chi,dphi_chi,d2phi_chi,yldc)
      else
         write(*,*) 'unexpected iyld_choice'
         stop -1
      endif

      call w_chr(imsg,'hah.f')
      call w_chr(imsg,'cauchy')
      call w_dim(imsg,cauchy,ntens,1.,.false.)
      call w_chr(imsg,'Just before hah_yieldsurface')
      call hah_yieldsurface(iyld_choice,yldc,nyldc,yldp,nyldp,cauchy,
     $     phi_chi,dphi_chi,d2phi_chi,ntens,phi,dphi,d2phi)
      call w_chr(imsg,'right after hah_yieldsurface')
c      call exit(-1)
      return
      end subroutine hah
c      include '/home/younguj/repo/abaqusPy/umats/yld/yld.f'
c      include '/home/younguj/repo/abaqusPy/umats/yld/hah/hah_lib.f'
