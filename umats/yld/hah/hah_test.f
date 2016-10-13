c     read_alpha is located in
      program test
c     Arguments
c     iyld_choice - isotropic yield surface kernel
c     stress      - cauchy stress tensor
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
      parameter(ntens=3,nyldc=9,nyldp=30)
      dimension yldc(nyldc),yldp(nyldp),stress(ntens)
      dimension dphi(ntens),d2phi(ntens,ntens)
      real*8 yldc,yldp,stress,phi,dphi,d2phi
      dimension dphi_chi(ntens),d2phi_chi(ntens,ntens)
      real*8 phi_chi,dphi_chi,d2phi_chi
      integer iyld_choice,i,imsg

c     local - microstructure deviator
      dimension emic(ntens)
      real*8 emic
c     local - Bauschinger parameters
      dimension gk(4)
      dimension e_ks(5),aux_ten(ntens)
      dimension f_ks(2)
      real*8 gk,e_ks,f_ks,eeq,aux_ten
c     local - Latent hardening parameters
      real*8 gL,ekL,eL
c     local - cross hardening parameters
      real*8 gS,c_ks,ss

      aux_ten(:)=0d0
      aux_ten(1)=1d0
      gL = 1d0

      call deviat(aux_ten,ntens,emic)

      call hah_io(yldp,nyldp,eeq,ntens,emic,gk,e_ks,f_ks,
     $     gL,ekL,eL,gS,c_ks,ss,1)

      iyld_choice=2             ! yld2000-2d

      imsg=0

      call fill_line(imsg,'*',72)

c      call read_alpha(
c     $     '/home/younguj/repo/abaqusPy/umats/yld/alfas.txt',yldc)
      yldc(:8) = 1d0
      yldc(9) = 8d0

      stress(:)=0d0
      stress(1)=1d0
      call w_chr(imsg,'cauchy stress')
      call w_dim(imsg,stress,ntens,1d0,.true.)
      call w_chr(imsg,'yldc')
      call w_dim(imsg,yldc,nyldc,1d0,.true.)
      call fill_line(imsg,'*',72)

      call w_chr(imsg,'just before entering hah')
      call w_chr(imsg,'cauchy stress')
      call w_dim(imsg,stress,ntens,1d0,.true.)
      call hah(iyld_choice,stress,phi,dphi,d2phi,yldc,yldp,nyldc,
     $     nyldp,ntens)
      call w_chr(imsg,'right after exit hah')

c$$$      call w_ival(imsg,'iyld_choice:',iyld_choice)
c$$$      call w_val( imsg,'phi_chi    :',phi_chi)
c$$$      call w_val( imsg,'phi        :',phi)

      end program test
