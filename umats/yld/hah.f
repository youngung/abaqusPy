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

      if (iyld_choice.eq.0) then
c**      phi_chi, dphi_chi, d2phi_chi
         call yld2000_2d(stress,phi_chi,dphi_chi,d2phi_chi,yldc)
      else
         write(*,*) 'unexpected iyld_choice'
         stop -1
      endif


      call hah_yieldsurface(yldp,nyldp,stress,
     $     phi_chi,dphi_chi,d2phi_chi,ntens,phi,dphi,d2phi)


      return
      end subroutine hah
c-----------------------------------------------------------------------
      subroutine hah_yieldsurface(yldp,nyldp,stress,
     $     phi_chi,dphi_chi,d2phi_chi,ntens,
     $     phi,dphi,d2phi)
c     Arguments
c     yldp      : yield surface parameters
c     nyldp     : Len of yldp
c     stress    : stress tensor
c     phi_chi   : isotropic yield surface
c     dphi_chi  : isotropic yield surface 1st derivative
c     d2phi_chi : isotropic yield surface 2nd derivative
c     ntens     : Len of stress tensor
      implicit none
      integer ntens,nyldp
      dimension yldp(nyldp),stress(ntens)

      real*8 yldp,stress

c     isotropic yield surface
      dimension dphi_chi(ntens),d2phi_chi(ntens,ntens)
      real*8 phi_chi,dphi_chi,d2phi_chi

c     HAH yield surface
      dimension dphi(ntens),d2phi(ntens,ntens)
      real*8 phi,dphi,d2phi



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



c     Restore yldp into state variables/parameters
      call hah_io(yldp,nyldp,eeq,emic,gk,e_ks,f_ks,
     $     gL,ekL,eL,0)





      return
      end subroutine hah_yieldsurface


c-----------------------------------------------------------------------
c     subroutine that converts back and forth
      subroutine hah_io(yldp,nyldp,
     $     eeq,

     $     emic,

     $     gk,e_ks,f_ks,

     $     gL,ekL,eL,

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

      elseif (iopt.eq.1) then   ! state variables -> yldp

c     HAH yield surface/state variables (and a few constants if any)
c***  equivalent plastic strain (cumulative)
         yldp(1)   = eeq
c***  microstructure deviator
         yldp(2:7) = e mic(:)
c***  Bauschinger effect
         yldp(8:11) =gk(:)      ! state variables
         yldp(12:16) = e_ks(:)  ! k1,k2,k3,k4,k5 constants
         do i=1,2
            yldp(i+16) = f_ks(i) ! f_k state parameters
         enddo
c***  Latent hardening
         yldp(19) = gL
         yldp(20) = ekL
         yldp(21)= eL
      endif

c     yield surface constants

      return
      end subroutine hah_io
c-----------------------------------------------------------------------
      include '/home/younguj/repo/abaqusPy/umats/yld/yld.f'
