c-----------------------------------------------------------------------
c     Microstructure deviator

c     Reference
c     [1] Manual of abaqusPy

c     Youngung Jeong
c     youngung.jeong@gmail.com
c-----------------------------------------------------------------------
c     Calculates 1) \frac{\partial \hat{h}}{\partial\bar{\varepsilon}}
c                2) dgr = ekrs(3)*(ekrs(4)*(1-coschi*coschi)-ekrs(5))
      subroutine micro_dev_deriv(ntens,ndi,nshr,nyldp,sdev,yldp)
c     Arguments
c     ntens : Len of tensor
c     ndi   : Len of normal components
c     nshr  : Len of shear components
c     nyldp : Len of yldp
c     sdev  : deviatoric stress
c     yldp  : yldp
      implicit none
      integer, intent(in) :: ntens,ndi,nshr,nyldp
      dimension yldp(nyldp),sdev(ntens)
      real*8, intent(in) :: sdev
      real*8, intent(inout) :: yldp
c     Additional local variables to transfer from yldp
      dimension emic(ntens),demic(ntens),target(ntens),emod(ntens,ntens)
     $     ,ds_dcauchy(6,6),aux(ntens),aux33(3,3),aux33_inv(3,3),aux6(6)
     $     ,gk(4),e_ks(5),f_ks(2),ekrs(5)
      real*8 emic,demic,dgr,emod,ds_dcauchy,aux,aux6,aux33,aux33_inv,
     $     target
      dimension
      real*8 gk,e_ks,f_ks,eeq,ref,gL,ekL,eL,gS,c_ks,ss,ekrs,coschi,H
      integer i,j

      H=8d0/3d0
c     restore variables from yldp
      call hah_io(0,nyldp,ntens,yldp,emic,demic,dgr,gk,e_ks,f_ks,eeq,
     $     ref,gL,ekL,eL,gS,c_ks,ss,ekrs,target)

c     Calculate coschi using target deviatoric tensor <target>
c     and current microstructure deviator <sdev>
      call calc_coschi(ntens,sdev,target,coschi) ! located in <hah_lib.f>
c     Calculate dh/debar - see Eqs 2&3 in Ref [1]
c     demic: \frac{\partial \hat{h}}{\partial\bar{\varepsilon}}
      demic(:) = ekrs(1) * sign(1d0,coschi)*
     $     (dabs(coschi/H)**(1d0/ekrs(2))+ekrs(5))*
     $     (target(:)-coschi*emic(:))
c     Calculate dgR:kr3 * [kr4 * ( 1-cos2chi) - kr5]
      dgr = ekrs(3)*(ekrs(4)*(1-coschi*coschi)-ekrs(5))
c     Save state variables to yldp
      call hah_io(1,nyldp,ntens,yldp,emic,demic,dgr,gk,e_ks,f_ks,eeq,
     $     ref,gL,ekL,eL,gS,c_ks,ss,ekrs,target)
      return
      end subroutine micro_dev_deriv
c-----------------------------------------------------------------------
c$$$      program test
c$$$      implicit none
c$$$      integer n
c$$$      parameter(n=3)
c$$$      dimension d(n,n)
c$$$      character(len=82) fmt
c$$$      integer i,j
c$$$      real*8 d
c$$$      if (n.eq.3) then
c$$$         fmt='(3f7.3)'
c$$$      elseif (n.eq.6) then
c$$$         fmt='(6f7.3)'
c$$$      endif
c$$$      call calc_ds_dcauchy(n,d)
c$$$      write(*,*)'d:'
c$$$      do i=1,n
c$$$         write(*,fmt,advance='no') (d(i,j),j=1,n)
c$$$         write(*,*)
c$$$      enddo
c$$$      end program
