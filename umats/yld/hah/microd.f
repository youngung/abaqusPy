c-----------------------------------------------------------------------
c     Microstructure deviator

c     Reference
c     [1] Manual of abaqusPy

c     Youngung Jeong
c     youngung.jeong@gmail.com
c$$$c-----------------------------------------------------------------------
c$$$      subroutine micro_dev(emic,target,nyldp,yldp)
c$$$c     Arguments
c$$$c     emic  : microstructure deviator
c$$$c     target: the target with which emic tries to realign
c$$$      implicit none
c$$$      integer, intent(in):: nyldp
c$$$      dimension emic(6),target(6),dh_deps(6),yldp(nyldp)
c$$$      real*8, intent(in) :: emic,target,yldp
c$$$      real*8 dh, cos_chi2, dh_deps
c$$$
c$$$
c$$$
c$$$      return
c$$$      end subroutine
c$$$c-----------------------------------------------------------------------
c     Obtain derivatives: Eqs. 25 and 26 in Ref [1]
      subroutine micro_dev(ntens,ndi,nshr,sdev,nyldp,yldp)
c     Arguments
c     ntens : Len of tensor
c     ndi   : Len of normal components
c     nshr  : Len of shear components
c     emic  : microstructure deviator
c     emod  : elastic modulus
c     dphi  : derivative of yield surface
      implicit none
      integer, intent(in) :: nyldp,ntens,ndi,nshr
      dimension yldp(nyldp),sdev(ntens)
      real*8 yldp,sdev
      dimension emic(ntens),demic(ntens),target(ntens),emod(ntens,ntens)
     $     ,dphi(ntens),ds_dcauchy(6,6),aux(ntens),aux33(3,3),
     $     aux33_inv(3,3),aux6(6)
      real*8 emic,demic,dgr,emod,dphi,ds_dcauchy,aux,aux6,aux33,
     $     aux33_inv,target
      integer i,j
c**   Arguments passed
      dimension gk(4),e_ks(5),f_ks(2),krs(5)
      real*8 gk,e_ks,f_ks,eeq,ref,gL,ekL,eL,gS,c_ks,ss,krs,coschi,H
      H=8d0/3d0
c     restore variables from yldp
      call hah_io(0,nyldp,ntens,yldp,emic,demic,dgr,gk,e_ks,f_ks,eeq,
     $     ref,gL,ekL,eL,gS,c_ks,ss,krs,target)

c     calculate {-C : dphi}^-1
      aux(:) = 0
      do 5 i=1,ntens
      do 5 j=1,ntens
         aux(i) = aux(i) - emod(i,j)*dphi(j)
 5    continue

c     Convert aux to its second order tensor
      if (ntens.eq.3) then
         aux6(:)   = 0d0
         aux6(1:2) = aux(1:2)
         aux6(6)   = aux(3)
         call voigt2(aux6,aux33)
      elseif (ntens.eq.6) then
         call voigt2(aux,aux33)
      endif

      aux33_inv(:,:)=aux33(:,:)
      call lu_inverse(aux33_inv,3)
      call calc_ds_dcauchy(6,ds_dcauchy)

c     Calculate dh
      call calc_coschi(ntens,ndi,nshr,sdev,target,coschi)
      demic(:) = krs(1) * sign(1d0,coschi)*
     $     (dabs(coschi/H)**(1d0/krs(2))+krs(3))*
     $     (target(:)-coschi*emic(:))
      dgr = krs(4)*(krs(5)*(1-coschi*coschi)-krs(3))
      call hah_io(1,nyldp,ntens,yldp,emic,demic,dgr,gk,e_ks,f_ks,eeq,
     $     ref,gL,ekL,eL,gS,c_ks,ss,krs,target)
      return
      end subroutine micro_dev
c-----------------------------------------------------------------------
      subroutine calc_ds_dcauchy(ntens,d)
      implicit none
      integer, intent(in) :: ntens
      dimension d(ntens,ntens)
      real*8, intent(out) :: d
      real*8 onethird
      parameter (onethird=1d0/3d0)
      integer i,j,ndi,nshr
cf2py intent(in) ntens
cf2py intent(out) d
cf2py depends(ntens) d
      d(:,:)=0d0
      if (ntens.eq.3) then
         ndi=2
         nshr=1
      elseif (ntens.eq.6) then
         ndi=3
         nshr=3
      else
         write(*,*)'ntens should be either 3 or 6 in',
     $        ' microd/calc_ds_dcauchy'
         call exit(-1)
      endif
      do 1 i=1,ntens
         d(i,i) = 1d0
 1    continue
      do 5 j=1,ndi
      do 5 i=1,ndi
         d(i,j) = d(i,j) -  onethird
 5    continue
      return
      end subroutine calc_ds_dcauchy

c$$$c     -------------
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
