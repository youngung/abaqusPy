c-----------------------------------------------------------------------
c     Derivatives associated with Homogeneous Anisotropic Hardening
c     model

c     General references
c     [1] Lee et al., IJP 29, (2012) p 13-41
c     [3] Lee et al., Comput. Methods. Appl. Mech. Engrg. 286 (2015) p63-86
      subroutine deriv_generals_exp(emic,gS,f1,f2)
      integer ntens
      dimension emic(3,3),phi_c(3,3,3,3),phi_o(3,3,3,3),phi_s(3,3,3,3),
     $     kro(3,3)
      real*8, intent(in) :: emic,gS,f1,f2
      real*8 phi_c,phi_o,phi_s,phi_p
      integer kro,i,j,krok
      real H
      H=8d0/3d0
      kro(:,:)=0
      do 5 i=1,3
         kro(i,i) = 1
 5    continue
c-----------------------------------------------------------------------
c     For hhat:s<0
c     Eq 48-a in Ref [3]
c     phi_c: \partial{s_c} / \partial{sdev}
      do 10 i=1,3
      do 10 j=1,3
      do 10 k=1,3
      do 10 l=1,3
         phi_c(i,j,k,l) = 0d0
      do 10 m=1,3
      do 10 n=1,3
         krok = (kro(m,k) * kro(n,l) + kro(m,l) * kro(n,k))
         phi_c(i,j,k,l) = phi_c(i,j,k,l) + H/2d0 * (
     $        krok * emic(m,n) * emic(i,j))
 10   continue
c     Eq 48-b
      do 20 i=1,3
      do 20 j=1,3
      do 20 k=1,3
      do 20 l=1,3
         phi_o(i,j,k,l,) = (kro(i,k)*kro(j,l)+kro(i,l)*kro(j,k))/2d0-
     $        phi_c(i,j,k,l)
 20   continue
c     Eq 49-a
      do 30 i=1,3
      do 30 j=1,3
      do 30 k=1,3
      do 30 l=1,3
         phi_s(i,j,k,l)=kro(i,k)*kro(j,l)-1d0/2d0*kro(i,j)*kro(k,l)
 30   continue
c     Eq 49-b
      do 40 i=1,3
      do 40 j=1,3
      do 40 k=1,3
      do 40 l=1,3
         phi_p(i,j,k,l) = 0d0
      do 40 m=1,3
      do 40 n=1,3
         phi_p(i,j,k,l) = phi_p(i,j,k,l) + 4d0*(1d0-gS)
     $        * phi_o(i,j,m,n) * phi_s(m,n,k,l)
 40   continue
c-----------------------------------------------------------------------
      return
      end subroutine
