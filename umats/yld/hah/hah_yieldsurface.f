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
     $     cauchy,phi_chi,dphi_chi,d2phi_chi,ntens,phi,dphi,d2phi)
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
      idiaw=.true.
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
c         call exit(-1)
      endif

c-----------------------------------------------------------------------
c     Restore yldp into state variables/parameters
      call hah_io(0,nyldp,ntens,yldp,emic,gk,e_ks,f_ks,eeq,ref,
     $     gL,ekL,eL,gS,c_ks,ss)

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
c         call exit(-1)
      endif

c---  anisotropic yield surface with isotropic hardening / phi_isoh
c     phi_chi, dphi_chi, d2phi_chi
      phi_chi_ref = phi_chi**yldc(9)
      phi_isoh    = (ref/phi_chi_ref)**(1d0/yldc(9))
      if (idiaw) then
         call w_val(imsg,'ref',ref)
         call w_val(imsg,'phi_chi',phi_chi)
         call w_val(imsg,'phi_chi_ref',phi_chi_ref)
         call w_val(imsg,'phi_isoh',phi_isoh)
         call exit(-1)
      endif
c      call exit(-1)

c---  anisotropic yield surface with latent hardening
      call latent(iyld_choice,ntens,nyldp,nyldc,cauchy,
     $     yldp,yldc,phi_lat)

      phi_lat_norm = phi_lat**yldp(9)
      phi_lat_norm = (ref/phi_lat_norm)**(1d0/yldp(9))

      if (idiaw) then
         call w_val(imsg,'phi_lat',phi_lat)
         call w_val(imsg,'phi_lat_norm',phi_lat_norm)
      endif
c      call exit(-1)

c---  anisotropic yield surface with full HAH
      call bauschinger(f_ks,yldp(9),emic,sdev,ntens,phi_bs(1),phi_bs(2))
      fht = phi_lat**yldp(9) + phi_bs(1)**yldp(9) + phi_bs(2)**yldp(9)
      rah = (ref/fht) ** (1/yldp(9))
      phi = rah*1d0
      if (idiaw) then
         call w_dim(imsg,phi_bs,2,1d0,.false.)
         call w_val(imsg,'rah',rah)
      endif

c      call exit(-1)

c$$$      sdp(:) = sc(:) + so(:)/gL
c$$$      if (idiaw) then
c$$$         call w_val(imsg,'gL:',gL)
c$$$         call w_chr(imsg,'sdp')
c$$$         call w_dim(imsg,sdp,ntens,1d0,.false.)
c$$$      endif
c$$$
c$$$c------------------------------
c$$$c     Latent extension
c$$$c------------------------------
c$$$c***  Target direction
c$$$      target(:) = sdev(:)
c$$$c***  stress double prime following eq 25 in Ref [1]
c$$$      if (gL.eq.0) then
c$$$         call w_empty_lines(imsg,2)
c$$$         call fill_line(imsg,'*',72)
c$$$         call w_chr(imsg,'**** Error gL is zero ****')
c$$$         call fill_line(imsg,'*',72)
c$$$         call exit(-1)
c$$$      endif
c$$$c      call exit(-1)
c$$$c------------------------------
c$$$      if (.false.) then
c$$$         call w_chr(imsg,'** calling yld for phi_lat **')
c$$$c     call exit(-1)
c$$$         call yld(iyld_choice,yldp,yldc,nyldp,nyldc,sdp,phi_lat,
c$$$     $        dphi_lat,d2phi_lat,ntens)
c$$$c     call exit(-1)
c$$$      else
c$$$         phi_lat=1d0
c$$$      endif
c$$$
c$$$c------------------------------
c$$$c     Cross load hardening
c$$$c------------------------------
c$$$      sp(:) = 4d0*(1d0-gS)*so(:)
c$$$      if (idiaw) then
c$$$         call w_chr(imsg,'** calling yld for phi_x **')
c$$$      endif
c$$$      call yld(iyld_choice,yldp,yldc,nyldp,nyldc,sp,phi_x,dphi_x,
c$$$     $     d2phi_x,ntens)
c$$$c      call exit(-1)
c$$$
c$$$      phi_omega = (dsqrt(phi_lat**2+phi_x**2))**8d0
c$$$
c$$$c------------------------------
c$$$c     Bauschinger
c$$$c------------------------------
c$$$      if (idiaw) then
c$$$         call w_chr(imsg,'** calling Bauschinger for phibs **')
c$$$         call w_chr(imsg,'f_ks passed to Bauschinger')
c$$$         call w_dim(imsg,f_ks,2,1d0,.false.)
c$$$      endif
c$$$      call bauschinger(f_ks,yldc(9),emic,sdev,ntens,phibs(1),phibs(2))
c$$$      if (idiaw) then
c$$$         call w_chr(imsg,'phib1, phib2')
c$$$         call w_dim(imsg,phibs,2,1d0,.false.)
c$$$         call w_chr(imsg,'**')
c$$$      endif
c$$$
c$$$c--------------------------------------------------
c$$$c     HAH Yield surface
c$$$      phi = (phi_omega+phibs(1)+phibs(2))**(1d0/8d0)
c$$$c--------------------------------------------------
c$$$
c$$$      if (idiaw)  then
c$$$         call w_empty_lines(imsg,2)
c$$$         call fill_line(imsg,'*',72)
c$$$         call w_val(imsg,'phi_x    :',phi_x)
c$$$         call w_val(imsg,'phi_chi  :',phi_chi)
c$$$         call w_val(imsg,'phi_lat  :',phi_lat)
c$$$         call w_val(imsg,'phi_omega:',phi_omega)
c$$$         call w_val(imsg,'phi      :',phi)
c$$$         call w_chr(imsg,'Exit HAH_YIELDSURFACE')
c$$$         call fill_line(imsg,'*',72)
c$$$         call w_empty_lines(imsg,2)
c$$$      endif
c$$$
c$$$c      call exit(-1)

      return
      end subroutine hah_yieldsurface

c-----------------------------------------------------------------------
c$$$      subroutine hah_ys_ref(iyld,cauchy,ntens,yldc,nyldc)
c$$$      implicit none
c$$$      integer iyld,ntens,nyldc
c$$$      dimension cauchy(ntens),yldc(nyldc)
c$$$      real*8 cauchy,yldc
c$$$      dimension dphi(ntens),d2phi(ntens,ntens)
c$$$      real*8 ref,dphi,d2phi
c$$$
c$$$
c$$$
c$$$      call yld2000_2d(cauchy,ref,dphi,d2phi,yldc)
c$$$
c$$$      return
c$$$  end subroutine hah_ys_ref
