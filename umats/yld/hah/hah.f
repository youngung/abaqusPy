c-----------------------------------------------------------------------
c     Homogeneous Anisotropic Hardening

c     General references
c     [1] Barlat et al. IJP 58, 2014 p201-218
c     [2] Jeong et al., IJP, 2016 (in press)

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
c     Arguments
      integer, intent(in):: iyld_choice,ntens,nyldc,nyldp
      dimension yldc(nyldc),yldp(nyldp),cauchy(ntens)
      dimension dphi(ntens),d2phi(ntens,ntens),sdev(ntens)
      real*8, intent(in)::yldc
      real*8, intent(out):: phi,dphi,d2phi
      real*8 yldp,cauchy,sdev

      dimension dphi_chi(ntens),d2phi_chi(ntens,ntens)
      real*8 phi_chi,dphi_chi,d2phi_chi
      integer i

      real*8 hydro
c     local controls
      integer imsg
      logical idiaw
      idiaw=.false.
c      idiaw=.true.
      imsg=0

      yldp(:)=yldp(:)
      d2phi(:,:)=d2phi(:,:)*1d0

c     HAH yield surface depends on the deviatoric stress

      call deviat(ntens,cauchy,sdev,hydro)

      dphi_chi(:)   =0d0
      d2phi_chi(:,:)=0d0

      if (idiaw) then
         call w_ival(imsg,'ntens:',ntens)
         call w_chr(imsg,'cauchy stress passed to subroutine HAH')
         call w_dim(imsg,cauchy,ntens,1d0,.false.)
         call w_chr(imsg,'In hah.f')
         call w_ival(imsg,'iyld_choice:',iyld_choice)
c         call exit(-1)
      endif

c     Calculate yield surface and its derivatives as a function of
c     only yield function kernel without HAH distortion.
c     These are saved to phi_chi, dphi_chi, d2phi_chi
      if (iyld_choice.eq.2) then ! yld2000-2d; plane-stress condition
         if (ntens.ne.3) then
            write(*,*)'When iyld_choice.eq.2, it should be',
     $           ' plane stress condition with ntens=3'
            call exit(-1)
         endif
         if (idiaw) then
            call w_chr(imsg,'cauchy before yld2000_2d')
            call w_dim(imsg,cauchy,ntens,1.d0,.false.)
         endif

         call yld2000_2d(cauchy,phi_chi,dphi_chi,d2phi_chi,yldc)

         if (idiaw) then
            call w_chr(imsg,'cauchy before yld2000_2d')
            call w_dim(imsg,cauchy,ntens,1.d0,.false.)
            call w_val(imsg,'phi_chi:',phi_chi)
            call w_chr(imsg,'dphi_chi:')
            call w_dim(imsg,dphi_chi,ntens,1d0,.false.)
            call w_chr(imsg,'d2phi_chi:')
            call w_mdim(imsg,d2phi_chi,ntens,1d0,.false.)
         endif
      else
         write(*,*) 'unexpected iyld_choice'
         call exit(-1)
      endif

      if (idiaw) then
         call w_chr(imsg,'hah.f')
         call w_chr(imsg,'cauchy')
         call w_dim(imsg,cauchy,ntens,1.d0,.false.)
         call w_chr(imsg,'Just before hah_yieldsurface')
c         call exit(-1)
      endif

c$$$c     test
      phi=phi_chi
      dphi(:)=dphi_chi(:)
      d2phi(:,:)=d2phi_chi(:,:)

      call hah_calc_ref(ntens,nyldp,nyldc,yldp,yldc,iyld_choice)

      call hah_yieldsurface(iyld_choice,yldc,nyldc,yldp,nyldp,
     $     cauchy,phi_chi,dphi_chi,d2phi_chi,ntens,phi,dphi,d2phi)

      return
      end subroutine hah

c$$$      call exit(-1)
c$$$      do 10 i=1,ntens
c$$$      do 10 j=1,ntens
c$$$         write(*,*)i,j
c$$$         d2phi(i,j)=d2phi_chi(i,j)
c$$$ 10   continue
c$$$      call exit(-1)

c     call exit(-1)

c$$$c**   saves ref to yldp
c$$$      call hah_calc_ref(ntens,nyldp,nyldc,yldp,yldc,iyld_choice)
c$$$
c$$$      call hah_yieldsurface(iyld_choice,yldc,nyldc,yldp,nyldp,
c$$$     $     cauchy,phi_chi,dphi_chi,d2phi_chi,ntens,phi,dphi,d2phi)
c$$$      if (idiaw) then
c$$$         call fill_line(imsg,'*',72)
c$$$         call w_chr(imsg,'Exiting subroutine hah')
c$$$         call fill_line(imsg,'*',72)
c$$$      endif




c      call exit(-1)

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



c      include '/home/younguj/repo/abaqusPy/umats/yld/yld.f'
c      include '/home/younguj/repo/abaqusPy/umats/yld/hah/hah_lib.f'
