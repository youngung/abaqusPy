c-----------------------------------------------------------------------
c     Homogeneous Anisotropic Hardening

c     General references
c     [1] Barlat et al. IJP 58, 2014 p201-218
c     [2] Jeong et al., IJP, 2016 (in press)

c     subroutine hah has several dependents
c     read_alpha: yld2000_2d.f
c     yld2000_2d: yld2000_2d.f
c     hah_yieldsurface: hah_yieldsurface.f
c-----------------------------------------------------------------------
      subroutine hah(iyld_choice,ntens,ndi,nshr,nyldc,nyldp,
     $     cauchy,yldc,yldp,phi,dphi,d2phi)
c     Arguments
c     iyld_choice : isotropic yield surface kernel
c     cauchy      : cauchy stress tensor
c     phi         : HAH yield surface
c     dphi        : HAH yield surface 1st derivatives
c     d2phi       : HAH yield surface 2nd derivatives
c     yldc        : isotropic yield surface constants
c     yldp        : HAH yield surface contants/state variables
c     nyldc       : Len of yldc
c     nyldp       : Len of yldp
c     ntens       : Len of stress tensor
c     ndi         : Number of normal components
c     nshr        : Number of shear components
      implicit none
c     Arguments
      integer, intent(in):: iyld_choice,ntens,ndi,nshr,nyldc,nyldp
      dimension yldc(nyldc),yldp(nyldp),cauchy(ntens)
      dimension dphi(ntens),d2phi(ntens,ntens),sdev(ntens)
      real*8, intent(in)::yldc
      real*8, intent(out):: phi,dphi,d2phi
      real*8 yldp,cauchy,sdev

      dimension dphi_chi(ntens),d2phi_chi(ntens,ntens)
      real*8 phi_chi,dphi_chi,d2phi_chi,hydro
c     local vars
      dimension sdev6(6)
      real*8 sdev6
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
         call w_empty_lines(imsg,1)
         call fill_line(imsg,'-',57)
         call w_chr(imsg,'Beginning of subroutine HAH')
         call w_ival(imsg,'** ntens:',ntens)
         call w_chr(imsg,'** cauchy stress')
         call w_dim(imsg,cauchy,ntens,1d0,.true.)
         call w_ival(imsg,'** iyld_choice:',iyld_choice)
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
            call w_chr(imsg,'** cauchy before yld2000_2d')
            call w_dim(imsg,cauchy,ntens,1.d0,.true.)
         endif

         call yld2000_2d(cauchy,phi_chi,dphi_chi,d2phi_chi,yldc)

         if (idiaw) then
            call w_chr(imsg,'** cauchy after yld2000_2d')
            call w_dim(imsg,cauchy,ntens,1.d0,.true.)
            call w_val(imsg,'phi_chi:',phi_chi)
            call w_chr(imsg,'dphi_chi:')
            call w_dim(imsg,dphi_chi,ntens,1d0,.true.)
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

c     calling hah_calc_ref stores (sqrt(phi(so)**2+phi(sdp)**2))**q to ref
      call hah_calc_ref(ntens,ndi,nshr,nyldp,nyldc,yldp,yldc,
     $     iyld_choice)
      call hah_yieldsurface(ntens,ndi,nshr,iyld_choice,nyldc,
     $     nyldp,yldc,yldp,cauchy,phi_chi,dphi_chi,d2phi_chi,
     $     phi,dphi,d2phi,.false.)

c      dphi(:)=dphi_chi(:)
c      d2phi(:,:)=d2phi_chi(:,:)

      if (ntens.eq.3) then
         sdev6(1:2)=sdev(1:2)
         sdev6(3)  = -sdev(1)-sdev(2)
         sdev6(4:5)=0d0
         sdev6(6) = sdev(3)
      elseif (ntens.eq.6) then
         sdev6(:)=sdev(:)
      else
         call w_chr(imsg,'Unexpected ntens given in hah.f')
         call exit(-1)
      endif

c      call hah_deriv(nyldp,sdev6,yldp


      return
      end subroutine hah
