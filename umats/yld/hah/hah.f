c-----------------------------------------------------------------------
c     Homogeneous Anisotropic Hardening
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
      integer iyld_choice
c**   iyld_choice=2             ! yld2000-2d
c     dummy alphas for testing purpose
      call read_alpha(
     $     '/home/younguj/repo/abaqusPy/umats/yld/alfas.txt',yldc)

      if (iyld_choice.eq.2) then
c**      phi_chi, dphi_chi, d2phi_chi
         call yld2000_2d(stress,phi_chi,dphi_chi,d2phi_chi,yldc)
      else
         write(*,*) 'unexpected iyld_choice'
         stop -1
      endif

      call hah_yieldsurface(iyld_choice,yldc,nyldc,yldp,nyldp,stress,
     $     phi_chi,dphi_chi,d2phi_chi,ntens,phi,dphi,d2phi)

      return
      end subroutine hah
c-----------------------------------------------------------------------
      subroutine hah_yieldsurface(iyld_choice,yldc,nyldc,
     $     yldp,nyldp,stress,
     $     phi_chi,dphi_chi,d2phi_chi,ntens,
     $     phi_n,dphi_n,d2phi_n)
c     Arguments
c     iyld_choice: choice of yield surface kernel
c     yldc      : yield surface constants
c     nyldc     : Len of yldc
c     yldp      : yield surface parameters
c     nyldp     : Len of yldp
c     stress    : stress tensor
c     phi_chi   : isotropic yield surface
c     dphi_chi  : isotropic yield surface 1st derivative
c     d2phi_chi : isotropic yield surface 2nd derivative
c     ntens     : Len of stress tensor
c     phi_n     : HAH yield surface at step n
c     dphi_n    : HAH yield surface 1st derivative at step n
c     d2phi_n   : HAH yield surface 2nd derivative at step n
      implicit none
      integer iyld_choice,ntens,nyldp,nyldc
      dimension yldc(nyldc),yldp(nyldp),stress(ntens)

      real*8 yldc,yldp,stress

c     isotropic yield surface
      dimension dphi_chi(ntens),d2phi_chi(ntens,ntens)
      real*8 phi_chi,dphi_chi,d2phi_chi

c     HAH yield surface
      dimension dphi_n(ntens),d2phi_n(ntens,ntens)
      real*8 phi_n,dphi_n,d2phi_n
      dimension phi_ns(0:1),dphi_ns(0:1,ntens),
     $     d2phi_ns(0:1,ntens,ntens)
      real*8 phi_ns,dphi_ns,d2phi_ns

c     local - microstructure deviator
      dimension emic(6)
      real*8 emic
c     local - Bauschinger parameters
      dimension gk(4)
      dimension e_ks(5)
      dimension f_ks(2)
      dimension target(ntens)
      real*8 gk,e_ks,f_ks,eeq,target
c     local - Latent hardening parameters
      real*8 gL,ekL,eL
c     local - cross hardening parameters
      real*8 gS,c_ks,ss
c     local
      dimension sc(ntens),so(ntens),sdp(ntens),sp(ntens) ! stress double prime
      real*8 sc,so,sdp,sp
c     local-latent
      dimension dphi_lat(ntens),d2phi_lat(ntens,ntens)
      real*8 phi_lat,dphi_lat,d2phi_lat
c     local-cross
      dimension dphi_x(ntens),d2phi_x(ntens,ntens)
      real*8 phi_x,dphi_x,d2phi_x,phi_omega


cf2py intent(in) yldp,nylpd,stress,phi_chi,dphi_chi,d2phi_chi,ntens
cf2py intent(out) phi,dphi,d2phi

c-----------------------------------------------------------------------

      phi_ns(0)       = phi_n
      dphi_ns(0,:)    = dphi_n(:)
      d2phi_ns(0,:,:) = d2phi_n(:,:)

c     Restore yldp into state variables/parameters
      call hah_io(yldp,nyldp,eeq,emic,gk,e_ks,f_ks,
     $     gL,ekL,eL,gS,c_ks,ss,0)

c     calculate yield surface

c     decompose stress
      call hah_decomposition(stress,ntens,emic,sc,so)

c------------------------------
c     Latent extension
c------------------------------
c***  Target direction
      target(:) = stress(:)
c***  stress double prime following eq 25 in Ref [1]
      sdp(:) = sc(:) + 1d0/gL

c------------------------------
c     Not sure if below would be okay since it seems like a recursive call
c     if that's the case, use yld2000_2d or vm_shell, hill48_shell directly
c------------------------------
      call yld(iyld_choice,yldp,yldc,nyldp,nyldc,sdp,
     $     phi_lat,dphi_lat,d2phi_lat,ntens)


c------------------------------
c     Cross load hardening
c------------------------------
      sp(:) = 4d0*(1d0-gS)*so(:)
      call yld(iyld_choice,yldp,yldc,nyldp,nyldc,sp,
     $     phi_x,dphi_x,d2phi_x,ntens)

      phi_omega = (phi_chi**2+phi_x**2)**(0.5d0)

      return
      end subroutine hah_yieldsurface

c-----------------------------------------------------------------------
c     subroutine that converts back and forth
      subroutine hah_io(yldp,nyldp,
     $     eeq,

     $     emic,

     $     gk,e_ks,f_ks,

     $     gL,ekL,eL,

     $     gS,c_ks,ss,

     $     iopt)
c     Arguments
c     yldp  : yield surface state variables
c     nyldp : Len of yldp
c     eeq   : equivalent plastic strain
c     emic   : microstructure deviator
c     gk    : gk parameters
c     gL    : gL
c     ekL   : kL
c     eL    : L parameters
c     gS    : gS parameter
c     c_ks  : ks parameter
c     ss    : S parameter
c     iopt  :  if 0, state variables <- yldp
c              if 1, state variables -> yldp
      implicit none
      integer iopt,nyldp
      dimension yldp(nyldp)
      real*8 yldp

c     local - microstructure deviator
      dimension emic(6)
      real*8 emic
c     local - Bauschinger parameters
      dimension gk(4)
      dimension e_ks(5)
      dimension f_ks(2)
      real*8 gk,e_ks,f_ks,eeq
c     local - Latent hardening parameters
      real*8 gL,ekL,eL
c     local - cross hardening parameters
      real*8 gS,c_ks,ss

c     local
      integer i

      if (iopt.eq.0) then       ! state variables <- yldp

c     HAH yield surface/state variables (and a few constants if any)
c***  equivalent plastic strain (cumulative)
         eeq     = yldp(1)
c***  microstructure deviator
         emic(:) = yldp(2:7)
c***  Bauschinger effect
         gk(:)   = yldp(8:11)   ! state variables
         e_ks(:) = yldp(12:16)  ! k1,k2,k3,k4,k5 constants
         do i=1,2
            f_ks(i) = yldp(i+16) ! f_k state parameters
         enddo
c***  Latent hardening
         gL      = yldp(19)
         ekL     = yldp(20)
         eL      = yldp(21)
c***  cross hardening
         gS      = yldp(22)
         c_ks    = yldp(23)
         ss      = yldp(24)
      elseif (iopt.eq.1) then   ! state variables -> yldp

c     HAH yield surface/state variables (and a few constants if any)
c***  equivalent plastic strain (cumulative)
         yldp(1)   = eeq
c***  microstructure deviator
         yldp(2:7) = e mic(:)
c***  Bauschinger effect
         yldp(8:11) = gk(:)      ! state variables
         yldp(12:16) = e_ks(:)  ! k1,k2,k3,k4,k5 constants
         do i=1,2
            yldp(i+16) = f_ks(i) ! f_k state parameters
         enddo
c***  Latent hardening
         yldp(19) = gL
         yldp(20) = ekL
         yldp(21) = eL
c***  cross hardening
         yldp(22) = gS
         yldp(23) = c_ks
         yldp(24) = ss
      endif

c     yield surface constants

      return
      end subroutine hah_io
c-----------------------------------------------------------------------
      include '/home/younguj/repo/abaqusPy/umats/yld/yld.f'
      include '/home/younguj/repo/abaqusPy/umats/yld/hah_lib.f'
