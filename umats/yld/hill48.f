c-----------------------------------------------------------------------
c     Hill quadratic yield surface

c     Youngung Jeong, Clemson University
c     youngung.jeong@gmail.com
c-----------------------------------------------------------------------
      subroutine hill48_shell(cauchy,phi,dphi,d2phi,yldc4)
c-----------------------------------------------------------------------
c     Arguments
c     cauchy : cauchy stress
c     phi    : yield surface
c     dphi   : yield surface 1st derivative
c     d2phi  : yield surface 2nd derivative
c     yldc4  : Hill parameters
c-----------------------------------------------------------------------
      implicit none
      dimension cauchy(3),yldc4(4),dphi(3),d2phi(3,3),cauchy6(6),
     $     dphi6(6),d2phi66(6,6),yldc6(6)
      real*8 cauchy, phi,dphi,d2phi,yldc4,dphi6,d2phi66,cauchy6,yldc6
      integer i,j,ii,jj
c-----------------------------------------------------------------------
!     Inflate cauchy3 to cauchy6
      cauchy6(:) = 0.d0
      cauchy6(1) = cauchy(1)
      cauchy6(2) = cauchy(2)
      cauchy6(6) = cauchy(3)
!     Inflate yldc4 to yldc6
      yldc6(1) = yldc4(1)
      yldc6(2) = yldc4(2)
      yldc6(3) = yldc4(3)
      yldc6(4) = 1.5d0
      yldc6(5) = 1.5d0
      yldc6(6) = yldc4(4)
      call hill48_gen(cauchy6,phi,dphi6,d2phi66,yldc6)
      dphi(:)=0.d0
      dphi(1)=dphi6(1)
      dphi(2)=dphi6(2)
      dphi(3)=dphi6(6)
      d2phi(:,:)=0.d0
      do 10 i=1,3
      do 10 j=1,3
         if (i.eq.3) ii=6
         if (i.lt.3) ii=i
         if (j.eq.3) jj=6
         if (j.lt.3) jj=j
         d2phi(i,j) = d2phi66(ii,jj)
 10   continue
      return
      end subroutine hill48_shell
c-----------------------------------------------------------------------
      subroutine hill48_gen(cauchy,phi,dphi,d2phi,yldc)
c-----------------------------------------------------------------------
c     Arguments
c     cauchy : cauchy stress
c     phi    : yield surface
c     dphi   : yield surface 1st derivative
c     d2phi  : yield surface 2nd derivative
c     yldc   : Hill surface parameters
c-----------------------------------------------------------------------
      implicit none
      dimension cauchy(6),dphi(6),d2phi(6,6),dh(6),d2h(6,6),
     $     yldc(6),s(6)
      real*8 cauchy,s,dphi,d2phi,phi,dh,d2h,psi,yldc,dff
      real*8 hh,hf,hg,hl,hm,hn
      integer i,j
c     Local Hill parameters
      hh=yldc(1)
      hf=yldc(2)
      hg=yldc(3)
      hl=yldc(4)
      hm=yldc(5)
      hn=yldc(6)
c     psi: homogeneous function
c-----------------------------------------------------------------------
      psi =  hH * (cauchy(1) - cauchy(2))**2 +
     $       hF * (cauchy(2) - cauchy(3))**2 +
     $       hG * (cauchy(3) - cauchy(1))**2 +
     $   2d0*hL * cauchy(4)**2 +
     $   2d0*hM * cauchy(5)**2 +
     $   2d0*hN * cauchy(6)**2
      phi = psi**5d-1
      s(:) = cauchy(:)

      dh(1) = 2d0*hG*s(1) - 2d0*hG*s(3) + 2d0*hH*s(1) - 2d0*hH*s(2)
      dh(2) = 2d0*hF*s(2) - 2d0*hF*s(3) + 2d0*hH*s(2) - 2d0*hH*s(1)
      dh(3) = 2d0*hF*s(3) - 2d0*hF*s(2) + 2d0*hG*s(3) - 2d0*hG*s(1)
      dh(4) = 4d0*hL *  s(4)
      dh(5) = 4d0*hM *  s(5)
      dh(6) = 4d0*hN *  s(6)

      d2h(:,:) = 0d0
      d2h(1,1) = 2d0*hh + 2d0*hg
      d2h(1,2) =    -hh
      d2h(1,3) =    -hg

      d2h(2,1) =    -hh
      d2h(2,2) = 2d0*hh + 2d0*hf
      d2h(2,3) =    -hf

      d2h(3,1) =    -hg
      d2h(3,2) =    -hf
      d2h(3,3) = 2d0*hf + 2d0*hg

      d2h(4,4) = 4d0*hl
      d2h(5,5) = 4d0*hm
      d2h(6,6) = 4d0*hn

c     1st derivatives
      dff     = 1d0 / (2d0*phi)
      dphi(:) = 0d0
      do 10 i=1,6
         dphi(i) =  dff * dh(i)
10    continue

c     2nd derivatives
      do 30 i=1,6
      do 30 j=1,6
         d2phi(i,j) = d2phi(i,j)
     $        - 1d0/4d0 * psi**(-3d0/2d0) * dh(i) * dh(j)
     $        + 1d0/2d0 * psi**(-1d0/2d0) * d2h(i,j)
 30   continue
      return
      end subroutine hill48_gen
c----------------------------------------------------------------------
c     Calculate Hill48 parameters using three r-values
c     Eq 3 in Dasappa et al. IJSS, vol 49, (2012)
      subroutine tuneH48(rvs,hp)
c     Arguments
c     rvs
c     hp
c----------------------------------------------------------------------
      implicit none
      dimension rvs(3),hp(4)
      real*8 rvs,hp,r0,r45,r90
      real*8 hh,hg,hf,hn
      r0 =rvs(1)
      r45=rvs(2)
      r90=rvs(3)
c-----------------------------------------------------------------------
      hh = r0/(r0+1d0)
      hg = 1d0-hh
      hf = hg * r0/r90
      hn = (r45+0.5d0)*(r0/r90+1d0)*hg
c-----------------------------------------------------------------------
      hp(1)=hh
      hp(2)=hg
      hp(3)=hf
      hp(4)=hn
      return
      end subroutine tuneH48
c-----------------------------------------------------------------------
c$$$      program test
c$$$      dimension rvs(3),hp(4)
c$$$      real*8 rvs,hp,r0,r45,r90
c$$$      real*8 hh,hg,hf,hn
c$$$      rvs(1)=0.9
c$$$      rvs(2)=1.5
c$$$      rvs(3)=0.7
c$$$      call tuneH48(rvs,hp)
c$$$      write(*,*)'hp:',hp
c$$$      return
c$$$      end program test
