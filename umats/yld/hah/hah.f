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
      dimension yldc(nyldc),yldp(nyldp),cauchy(ntens),cauchy_ref(ntens),
     $     cauchy_test(ntens),sdev(ntens),sdev_test(ntens),
     $     sdev_ref(ntens)
      dimension dphi(ntens),d2phi(ntens,ntens)
      real*8 yldc,yldp,cauchy,phi,dphi,d2phi,cauchy_ref,cauchy_test,
     $     sdev,sdev_test,sdev_ref
      dimension dphi_chi(ntens),d2phi_chi(ntens,ntens)
      real*8 phi_chi,dphi_chi,d2phi_chi
      integer iyld_choice,i,j

      real*8 ref,q,hydro
c     local controls
      integer imsg
      logical idiaw
      idiaw=.false.
      imsg=0

c     HAH yield surface depends only on the deviatoric stress
      call deviat(ntens,cauchy,sdev,hydro)

      dphi_chi(:)=0d0
      d2phi_chi(:,:)=0d0

      if (idiaw) then
         call w_chr(imsg,'cauchy stress passed to subroutine hah')
         call w_dim(imsg,cauchy,ntens,1d0,.false.)
c        call exit(-1)
      endif
c$$$      call read_alpha(
c$$$     $     '/home/younguj/repo/abaqusPy/umats/yld/alfas.txt',yldc)
      if (idiaw) then
         call w_chr(imsg,'In hah.f')
         call w_ival(imsg,'iyld_choice:',iyld_choice)
      endif

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
         call exit(-1)
      endif

      if (idiaw) then
         call w_chr(imsg,'hah.f')
         call w_chr(imsg,'cauchy')
         call w_dim(imsg,cauchy,ntens,1.,.false.)
         call w_chr(imsg,'Just before hah_yieldsurface')
      endif

      call latent(iyld_choice,ntens,nyldp,nyldc,cauchy,yldp,yldc,phi)
      call w_val(imsg,'phi:',phi)
      call exit(-1)



c$$$      !call deviat(cauchy_test,ntens,sdev_test)
c$$$      call deviat(ntens,cauchy_test,ntens,sdev_test)
c$$$
c$$$c     idiaw=.true.
c$$$
c$$$      call hah_yieldsurface(iyld_choice,yldc,nyldc,yldp,nyldp,
c$$$     $     sdev_test,phi_chi,dphi_chi,d2phi_chi,ntens,phi,dphi,d2phi)
c$$$
c$$$      q   = yldp(9)
c$$$      ref = (1d0/phi)**(1d0/q)
c$$$      if (idiaw) then
c$$$         call w_val(imsg,'ref:',ref)
c$$$      endif
c$$$      do i=1,ntens
c$$$         sdev_ref(i) = sdev(i)/ref
c$$$      enddo
c$$$      call hah_yieldsurface(iyld_choice,yldc,nyldc,yldp,nyldp,
c$$$     $     sdev_ref,phi_chi,dphi_chi,d2phi_chi,ntens,phi,dphi,d2phi)
c$$$      if (idiaw) then
c$$$         call w_val(imsg,'phi:',phi)
c$$$         call w_chr(imsg,'right after hah_yieldsurface')
c$$$      endif

c      call exit(-1)
      return
      end subroutine hah
c      include '/home/younguj/repo/abaqusPy/umats/yld/yld.f'
c      include '/home/younguj/repo/abaqusPy/umats/yld/hah/hah_lib.f'
