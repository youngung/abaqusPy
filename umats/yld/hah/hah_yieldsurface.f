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
      subroutine hah_yieldsurface(iyld_choice,yldc,nyldc,yldp,nyldp,
     $     cauchy,phi_chi,dphi_chi,d2phi_chi,ntens,phi,dphi,d2phi,idiaw)
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
c     phi        : HAH yield surface
c     dphi       : HAH yield surface 1st derivative
c     d2phi      : HAH yield surface 2nd derivative
      implicit none
      integer,intent(in)::iyld_choice,ntens,nyldp,nyldc
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
      dimension emic(6)
      real*8 emic
c     local - Bauschinger parameters
      dimension gk(4),e_ks(5),f_ks(2),target(ntens),phi_bs(2)
      real*8 gk,e_ks,f_ks,eeq,target,phi_bs
c     local - Latent hardening parameters
      real*8 gL,ekL,eL
c     local - cross hardening parameters
      real*8 gS,c_ks,ss
c     local
      dimension sc(ntens),so(ntens),sdp(ntens),sp(ntens) ! stress double prime
      real*8 sc,so,sdp,sp,ref
c     local-latent
      dimension dphi_lat(ntens),d2phi_lat(ntens,ntens)
      real*8 phi_lat,dphi_lat,d2phi_lat,phi_lat_norm
c     local-cross
      dimension dphi_x(ntens),d2phi_x(ntens,ntens)
      real*8 phi_x,dphi_x,d2phi_x,phi_omega
c     local-bau
      dimension phibs(2)
      real*8 phibs
c     local-control
      integer imsg
      logical idiaw
      real*8 phi_chi_ref, phi_isoh, hydro, fht, rah
cf2py intent(in) iyld_choice,yldc,nyldc,yldp,nylpd
cf2py intent(in) sdev,phi_chi,dphi_chi,d2phi_chi,ntens
cf2py intent(out) phi,dphi,d2phi

c-----------------------------------------------------------------------
      imsg = 0
c      idiaw=.true.
c      idiaw=.false.

c     Fooling the compiler - by-pass warning of unused arguments.
      d2phi(:,:)=d2phi(:,:)
      d2phi_chi(:,:)=d2phi_chi(:,:)
      d2phi_lat(:,:)=d2phi_lat(:,:)
      d2phi_x(:,:)=d2phi_x(:,:)
      dphi(:)=dphi(:)
      dphi_chi(:)=dphi_chi(:)
      dphi_lat(:)=dphi_lat(:)
      dphi_x(:)=dphi_x(:)
      phi_omega=phi_omega
      phi_x=phi_x
      phibs(:)=phibs(:)
      sdp(:)=sdp(:)
      sp(:)=sp(:)
      target(:)=target(:)
c -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -
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

c-----------------------------------------------------------------------
c     Restore yldp into state variables/parameters
c      write(*,*)'ref:',ref
      call hah_io(0,nyldp,ntens,yldp,emic,gk,e_ks,f_ks,eeq,ref,
     $     gL,ekL,eL,gS,c_ks,ss)
c      call exit(-1)

c     calculate yield surface
c     decompose deviatoric stress
      if (idiaw) then
         call w_chr(imsg,'calling to hah_decompose')
      endif
      call hah_decompose(sdev,ntens,emic,sc,so)
      if (idiaw) then
         call w_chr(imsg,'deviatoric Stress decomposition')
         call w_dim(imsg,sdev,ntens,1d0,.true.)
         call w_chr(imsg,'emic')
         call w_dim(imsg,emic,ntens,1d0,.true.)
         call w_chr(imsg,'sc')
         call w_dim(imsg,sc,ntens,1d0,.true.)
         call w_chr(imsg,'so')
         call w_dim(imsg,so,ntens,1d0,.true.)
         call fill_line(imsg,'-',23)
c         call exit(-1)
      endif

c---  anisotropic yield surface with isotropic hardening / phi_isoh
c     phi_chi, dphi_chi, d2phi_chi
      phi_isoh  = phi_chi/ref
      if (idiaw) then
         call w_val(imsg,'ref (unix)  :',ref)
         call w_val(imsg,'phi_chi     :',phi_chi)
         call w_val(imsg,'phi_isoh    :',phi_isoh)
         call fill_line(imsg,'-',23)
c         call exit(-1)
      endif
c      call exit(-1)

c---  anisotropic yield surface with latent hardening + cross hardening
c     (sqrt(phi(sp)**2 + phi(sdp)**2)) ** q
      call latent(iyld_choice,ntens,nyldp,nyldc,cauchy,
     $     yldp,yldc,phi_lat)
      phi_lat_norm = phi_lat**(1d0/yldp(9)) !! homogeneous function of degree 1.
      if (idiaw) then
         call w_val(imsg,'phi_lat     :',phi_lat)
         call w_val(imsg,'phi_lat_norm:',phi_lat_norm)
         call fill_line(imsg,'-',23)
c         call exit(-1)
      endif

c---  anisotropic yield surface with full HAH
      call bauschinger(f_ks,yldp(9),emic,sdev,ntens,phi_bs(1),phi_bs(2))
      fht = (phi_lat_norm/ref)**yldp(9)+
     $     phi_bs(1)**yldp(9) + phi_bs(2)**yldp(9)
      phi = (fht)**(1d0/yldp(9))
      if (idiaw) then
         call w_chr(imsg,'phi_bs (Bauschinger)')
         call w_dim(imsg,phi_bs,2,1d0,.true.)
         call w_val(imsg,'fht          :',fht)
         call w_val(imsg,'rah          :',rah)
         call w_val(imsg,'phi          :',phi)
         call w_chr(imsg,'cauchy stress:')
         call w_dim(imsg,cauchy,ntens,1d0,.true.)
c         call exit(-1)
      endif

      return
      end subroutine hah_yieldsurface
c-----------------------------------------------------------------------
