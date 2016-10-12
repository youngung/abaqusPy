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
      parameter(ntens=6,nyldc=9,nyldp=30)
      dimension yldc(nyldc),yldp(nyldp),stress(ntens)
      dimension dphi(ntens),d2phi(ntens,ntens)
      real*8 yldc,yldp,stress,phi,dphi,d2phi
      dimension dphi_chi(ntens),d2phi_chi(ntens,ntens)
      real*8 phi_chi,dphi_chi,d2phi_chi
      integer iyld_choice,i,imsg
      iyld_choice=2             ! yld2000-2d

      imsg=0

      call fill_line(imsg,'*',72)

      call read_alpha(
     $     '/home/younguj/repo/abaqusPy/umats/yld/alfas.txt',yldc)
      stress(:)=0d0
      stress(1)=1d0
      call w_chr(imsg,'stress')
      call w_dim(imsg,stress,ntens,1d0,.true.)
      call w_chr(imsg,'yldc')
      call w_dim(imsg,yldc,nyldc,1d0,.true.)
      call fill_line(imsg,'*',72)
      
      call hah(iyld_choice,stress,phi,dphi,d2phi,yldc,yldp,nyldc,
     $     nyldp,ntens)

      call w_ival(imsg,'iyld_choice:',iyld_choice)
      call w_val( imsg,'phi_chi    :',phi_chi)
      call w_val( imsg,'phi        :',phi)

      end program test
