      subroutine
      end subroutine









c-----------------------------------------------------------------------
c     Derivatives associated with Homogeneous Anisotropic Hardening
c     model

c     General references
c     [1] Lee et al., IJP 29, (2012) p 13-41
c     [2] Lee et al., Comput. Methods. Appl. Mech. Engrg. 247 (2012) p73-92
c     [3] Lee et al., Comput. Methods. Appl. Mech. Engrg. 286 (2015) p63-86

c-----------------------------------------------------------------------
      subroutine main_deriv(nyldp,ntens,yldp,sdev)
      implicit none
      integer, intent(in) nyldp,ntens
      dimension yldp(nyldp),sdev(ntens)
      real*8, intent(in) yldp,sdev

      return
      end subroutine main_deriv
c-----------------------------------------------------------------------
c     Calculate derivatives between various stress components
c     dsc_ds, dso_ds,ds_dcauchy,dsp_dcauchy
      subroutine deriv_stress_components(iopt,emic,gS,f1,f2,
     $     dsc_ds,dso_ds,ds_dcauchy,dsp_dcauchy)
      implicit none
      integer, intent(in) :: iopt
      dimension emic(3,3),dsc_ds(3,3,3,3),dso_ds(3,3,3,3),
     $     ds_dcauchy(3,3,3,3),kro(3,3)
      real*8, intent(in) :: emic,gS,f1,f2
      real*8, intent(out) :: dsc_ds,dso_ds,ds_dcauchy,dsp_dcauchy
      real*8 temp
      integer kro,i,j,krok,onethird
      parameter(onethird=1d0/3d0)
      real H
      H=8d0/3d0
      kro(:,:)=0
c     Kronecker
      kro(1,1) = 1
      kro(2,2) = 1
      kro(3,3) = 1



c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Eq 48-a in Ref [3]
c     dsc_ds: \partial{s_c} / \partial{sdev}
      temp = H*0.5d0
      do 10 i=1,3
      do 10 j=1,3
      do 10 k=1,3
      do 10 l=1,3
         dsc_ds(i,j,k,l) = 0d0
      do 10 m=1,3
      do 10 n=1,3
         krok = (kro(m,k) * kro(n,l) + kro(m,l) * kro(n,k))
         dsc_ds(i,j,k,l) = dsc_ds(i,j,k,l) + temp * (
     $        krok * emic(m,n) * emic(i,j))
 10   continue
c     Eq 48-b
      temp = 0.5d0
      do 20 i=1,3
      do 20 j=1,3
         dso_ds(i,j,:,:)=0d0
      do 20 k=1,3
      do 20 l=1,3
         dso_ds(i,j,k,l,) = dso_ds(i,j,k,l) +
     $        (kro(i,k)*kro(j,l)+kro(i,l)*kro(j,k))*temp-
     $        dsc_ds(i,j,k,l)
 20   continue
c     Eq 49-a
      do 30 i=1,3
      do 30 j=1,3
         ds_dcauchy(i,j,:,:)=0d0
      do 30 k=1,3
      do 30 l=1,3
         ds_dcauchy(i,j,k,l)=ds_dcauchy(i,j,k,l) +
     $        kro(i,k)*kro(j,l)-onethird*kro(i,j)*kro(k,l)
 30   continue
c     Eq 49-b
      temp = 4d0*(1d0-gS)
      do 40 i=1,3
      do 40 j=1,3
      do 40 k=1,3
      do 40 l=1,3
         dsp_dcauchy(i,j,k,l) = 0d0
      do 40 m=1,3
      do 40 n=1,3
         dsp_dcauchy(i,j,k,l) = dsp_dcauchy(i,j,k,l) +
     $        temp * dso_ds(i,j,m,n) * ds_dcauchy(m,n,k,l)
 40   continue
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Eq 47-a in Ref [3]

      return
      end subroutine
