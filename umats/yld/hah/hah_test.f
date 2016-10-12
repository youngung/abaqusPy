c     read_alpha is located in
      program test
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
      parameter(ntens=6,nyldc=20,nyldp=20)
      dimension yldc(nyldc),yldp(nyldp),stress(ntens)
      dimension dphi(ntens),d2phi(ntens,ntens)
      real*8 yldc,yldp,stress,phi,dphi,d2phi
      dimension dphi_chi(ntens),d2phi_chi(ntens,ntens)
      real*8 phi_chi,dphi_chi,d2phi_chi
      integer iyld_choice
c**   iyld_choice=2             ! yld2000-2d

      call read_alpha(
     $     '/home/younguj/repo/abaqusPy/umats/yld/alfas.txt',yldc)

      call hah(iyld_choice,stress,phi,dphi,d2phi,yldc,yldp,nyldc,
     $     nyldp,ntens)


c**   to suppress -wnused-dummy-argument
      write(*,*) iyld_choice,yldp,d2phi_chi,
     $     dphi_chi,phi_chi,d2phi_chi

      end program test
