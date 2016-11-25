c-----------------------------------------------------------------------
c     Module used to calculate yield surface distorted by HAH approach

c     General references
c     [1] Barlat et al. IJP 58, 2014 p201-218
c     [2] Jeong et al., IJP, 2016 (in press)

c     Dependents
c     hah_io in hah_lib.f
c     hah_decompose in hah_lib.f
c     yld in  ../yld.f
c-----------------------------------------------------------------------
c     hah yield function that is composed two components: phi, phi_b.
c     yield surface =[phi^q+phi_b^q]^(1/q).
c     phi is a homogeneous function of first degree.

      subroutine hah_yieldsurface(ntens,ndi,nshr,iyld_choice,nyldc,
     $     nyldp,yldc,yldp,cauchy,phi_chi,dphi_chi,d2phi_chi,phi,dphi,
     $     d2phi,idiaw)
c     Arguments
c     iyld_choice: choice of yield surface kernel
c     yldc       : yield surface constants
c     nyldc      : Len of yldc
c     yldp       : yield surface parameters
c     nyldp      : Len of yldp
c     cauchy     : cauchy stress tensor
c     phi_chi    : isotropic yield surface
c     dphi_chi   : isotropic yield surface 1st derivative
c     d2phi_chi  : isotropic yield surface 2nd derivative
c     ntens      : Len of deviatoric stress tensor
c     ndi        : Number of normal components
c     nshr       : Number of shear components
c     phi        : HAH yield surface
c     dphi       : HAH yield surface 1st derivative
c     d2phi      : HAH yield surface 2nd derivative
      implicit none
      integer,intent(in)::iyld_choice,ntens,ndi,nshr,nyldp,nyldc
      dimension yldc(nyldc),yldp(nyldp),sdev(ntens),cauchy(ntens)
      real*8, intent(in)::cauchy
      real*8 yldc,yldp,sdev
c     isotropic yield surface
      dimension dphi_chi(ntens),d2phi_chi(ntens,ntens)
      real*8 phi_chi,dphi_chi,d2phi_chi
c     HAH yield surface
      dimension dphi(ntens),d2phi(ntens,ntens)
      real*8,intent(out):: phi,dphi,d2phi
c     local - microstructure deviator
      dimension emic(ntens),demic(ntens),krs(4),target(ntens)
      real*8 emic,demic,krs,target,dgr
c     local - Bauschinger parameters
      dimension gk(4),e_ks(5),f_ks(2),phi_bs(2)
      real*8 gk,e_ks,f_ks,eeq,phi_bs
c     local - Latent hardening parameters
      real*8 gL,ekL,eL
c     local - cross hardening parameters
      real*8 gS,c_ks,ss
c     local
      dimension sc(ntens),so(ntens),spp(ntens),sp(ntens),dphi_hah(6),
     $     dphi_hah_fin(6)
      real*8 sc,so,spp,sp,ref0,ref1,dphi_hah,dphi_hah_fin

      dimension cauchy6(ntens),dpsi_dspp6(ntens),dpsi_dsp6(ntens)
      real*8 cauchy6,dpsi_dspp6,dpsi_dsp6
c     local-latent
      dimension dphi_h(ntens),d2phi_h(ntens,ntens),dpsi_dspp(ntens),
     $     dpsi_dsp(ntens)
      real*8 phi_h,dphi_h,d2phi_h,dpsi_dspp,dpsi_dsp,
     $     psi_spp,psi_sp
c     local-cross
      dimension dphi_x(ntens),d2phi_x(ntens,ntens)
      real*8 phi_x,dphi_x,d2phi_x,phi_omega
c     local-bau
      dimension phibs(2)
      real*8 phibs
c     local-control
      integer imsg
      logical, intent(in)::idiaw
      real*8 phi_isoh, hydro, fht, rah
cf2py intent(in) iyld_choice,yldc,nyldc,yldp,nylpd
cf2py intent(in) sdev,phi_chi,dphi_chi,d2phi_chi,ntens
cf2py intent(out) phi,dphi,d2phi

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      imsg = 0
c     Fooling the compiler - by-pass warning of unused arguments.
      d2phi(:,:)     = d2phi(:,:)
      d2phi_chi(:,:) = d2phi_chi(:,:)
      d2phi_h(:,:) = d2phi_h(:,:)
      d2phi_x(:,:)   = d2phi_x(:,:)
      dphi(:)        = dphi(:)
      dphi_chi(:)    = dphi_chi(:)
      dphi_h(:)    = dphi_h(:)
      dphi_x(:)      = dphi_x(:)
      phi_omega      = phi_omega
      phi_x          = phi_x
      phibs(:)       = phibs(:)
      spp(:)         = spp(:)
      sp(:)          = sp(:)
      target(:)      = target(:)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c$$$c     obtain deviatoric stress
      call deviat(ntens,cauchy,sdev,hydro)
      if (idiaw) then
         call fill_line(imsg,'#',72)
         call w_chr(imsg,'Enter HAH_YIELDSURFACE')
         call fill_line(imsg,'#',72)
         call w_chr(imsg,'deviatoric stress')
         call w_dim(imsg,sdev,ntens,1d0,.true.)
         call fill_line(imsg,'-',23)
c     call exit(-1)
      endif
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Restore yldp into state variables/parameters
      call hah_io(0,nyldp,ntens,yldp,emic,demic,dgr,gk,e_ks,f_ks,eeq,
     $     ref0,ref1,gL,ekL,eL,gS,c_ks,ss,krs,target)
c     calculate yield surface
c     decompose deviatoric stress
      if (idiaw) then
         call w_chr(imsg,'calling to hah_decompose')
      endif
c      call hah_decompose(sdev,ntens,emic,sc,so)
      call hah_decompose(ntens,ndi,nshr,sdev,emic,sc,so)
      if (idiaw) then
         call w_chr(imsg,'deviatoric stress decomposition')
         call w_dim(imsg,sdev,ntens,1d0,.true.)
         call w_chr(imsg,'emic')
         call w_dim(imsg,emic,ntens,1d0,.true.)
         call w_chr(imsg,'sc')
         call w_dim(imsg,sc,  ntens,1d0,.true.)
         call w_chr(imsg,'so')
         call w_dim(imsg,so,  ntens,1d0,.true.)
         call fill_line(imsg,'-',23)
c         call exit(-1)
      endif

c---  anisotropic yield surface with isotropic hardening / phi_isoh
c     phi_chi, dphi_chi, d2phi_chi
      phi_isoh  = phi_chi/ref1
      if (idiaw) then
         call w_val(imsg,'ref1        :',ref1)
         call w_val(imsg,'phi_chi     :',phi_chi)
         call w_val(imsg,'phi_isoh    :',phi_isoh)
         call fill_line(imsg,'-',23)
c         call exit(-1)
      endif
c      call exit(-1)

c---  anisotropic yield surface with latent hardening + cross hardening
c     phi_h = (phi(sp)**2 + phi(spp)**2)^{1/2}
c     phi_h being the homogeneous function of degree 1.
      call latent(iyld_choice,ntens,ndi,nshr,nyldp,nyldc,
     $     cauchy,yldp,yldc,dpsi_dspp,dpsi_dsp,psi_spp,psi_sp,
     $     phi_h)
      if (idiaw) call w_val(imsg,'phi_h     :',phi_h)
      phi_h = phi_h**yldc(9)
      if (idiaw) then
         call w_val(imsg,'phi_h_norm:',phi_h)
         call fill_line(imsg,'-',23)
c         call exit(-1)
      endif

c---  anisotropic yield surface with full HAH
      call bauschinger(ntens,ndi,nshr,emic,sdev,f_ks,yldc(9),phi_bs(1),
     $     phi_bs(2))
      fht = phi_h+phi_bs(1)+phi_bs(2)
c      phi = (1.d0/fht)**(1d0/yldc(9))
      phi = (ref1/fht)**(1d0/yldc(9))
      if (idiaw) then
         call w_chr(imsg,'phi_bs (Bauschinger)')
         call w_dim(imsg,phi_bs,2,1d0,.true.)
         call w_val(imsg,'fht          :',fht)
         call w_val(imsg,'ref1         :',ref1)
         call w_val(imsg,'rah          :',rah)
         call w_val(imsg,'phi          :',phi)
         call w_chr(imsg,'cauchy stress:')
         call w_dim(imsg,cauchy,ntens,1d0,.true.)
         call w_empty_lines(imsg,2)
c         call exit(-1)
      endif

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Calculate the derivatives of the yield surface
c     with respect to the cauchy stress.
      !! dphi_hah is the guess.
      !! dphi_hah_fin is the newly found dphi_hah_fin

c     hah_deriv is compatible with full-dimensional stress space.

      if (ntens.eq.3) then
         call reduce_3to6(cauchy   ,cauchy6)
         call reduce_3to6(dpsi_dspp,dpsi_dspp6)
         call reduce_3to6(dpsi_dsp ,dpsi_dsp6)
      elseif(ntens.eq.6) then
         cauchy6(:)   =cauchy(:)
         dpsi_dspp6(:)=dpsi_dspp(:)
         dpsi_dsp6(:) =dpsi_dsp(:)
      endif

c     !! use derivatives of isotropic yield surface kernel
      if (ntens.eq.3) then
         call reduce_3to6(dphi_chi,dphi_hah)
      elseif(ntens.eq.6) then
         dphi_hah(:)=dphi_chi(:)
      endif

      call hah_deriv(nyldp,cauchy6,yldp,yldc(9),phi_h,
     $     dphi_h,psi_spp,dpsi_dspp,psi_sp,dpsi_dsp,phi,
     $     dphi_hah,dphi_hah_fin)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      dphi(1:2) = dphi_hah_fin(1:2)
      dphi(3)   = dphi_hah_fin(3)

      return
      end subroutine hah_yieldsurface
c-----------------------------------------------------------------------
