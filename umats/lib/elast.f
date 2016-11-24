c----------------------------------------------------------------------c
c     Calculate complete 6x6 ISOTROPIC elastic modulus (c) using
c     Young's modulus (e) and poisson ratio (nu)
c----------------------------------------------------------------------c
c$$$  3D shell with S11,S22 and S12
      subroutine emod_iso_shell(e,enu,G,ekappa,cel)
c     intent(in) e, nu
c     intent(out) G,ekappa,c
c      e    : Young's modulus
c      enu  : Poisson ratio
c      G    : Shear modulus
c     ekappa: bulk modulus
c      cel  : elastic constants
      implicit none
      real*8 cel(3,3)
      real*8 enu, e, x, G, ekappa
      integer i,j
      cel(:,:)=0d0
c     Multiplier
      x = e/(1d0+enu)/(1d0-2d0*enu)
      do 10 i=1,2
         do 5 j=1,2
            cel(i,j) = x*enu
 5       continue
         cel(i,i) = x*(1d0-enu)
 10   continue
      cel(3,3) = x* (1d0-2d0*enu)/2d0
      ekappa = e /3d0 / (1d0-2d0*enu)
      G=E/2d0/(1d0+enu)
c      call w_mdim(0,cel,3,1d0)
      return
      end subroutine
c     -------------------------------------
      subroutine emod_iso(e,nu,c,ndi,nshr)
      implicit none
c     intent(in): e,nu,ndi,nshr
c     intent(out): c
      integer, intent(in) :: ndi,nshr
      dimension c(ndi+nshr,ndi+nshr)
      real*8, intent(in) ::  e, nu
      real*8, intent(out) :: c
      real*8 x
      integer i,j,imsg,ntens
      ntens=ndi+nshr
      imsg=7

c     initialization
      c(:,:)=0d0
c$$$      do 20 i=1,ntens
c$$$      do 20 j=1,ntens
c$$$         c(i,j) = 0d0
c$$$ 20   continue
c      write(imsg,*) 'after initialization c matrix'
c
c     construct elastic tensor (6x6) with assuming that
c     \gamma_ij = 2\varepsilon_ij is the engineering shear strain
c
c     Multiplier
      x = e/(1d0+nu)/(1d0-2d0*nu)
c     off-diagonal terms
      do 50 i=1,3
      do 50 j=1,3
         c(i,i) = nu * x
 50   continue
      do 100 i=1,3
         c(i,i)     = (1d0-nu)*x    !! overwrite the diganogal term
         c(i+3,i+3) = (1d0-2d0*nu)/2d0 * x
 100  continue
c      write(imsg,*) 'just before returning'
      return
      end subroutine emod_iso
c----------------------------------------------------------------------c
      subroutine el_iso_jacobian(e,nu,flow,syield,predictor,
     $     ndi,ntens,ddsdde)
c     intent(in):: e,nu,flow,syield,ndi,ntens
c     intent(out):: ddsdde
      parameter(one=1.d0,two=2.d0,three=3.d0)
      dimension flow(ntens), ddsdde(ntens,ntens), c(ntens,ntens)
      real*8 e, nu, syield, predictor,ddsdde
      real*8 mu, k              !! elastic constants (shear modulus and bulk modulus)
      integer i,j,ndi,ntens,nshr
      real*8 mustar, lamstar, fact1,fact2,fact3,c
      nshr = ntens - ndi
      call emod_iso(e,nu,c,ndi,nshr)
      k=e/(three*(one-two*nu))
      mu=e/(  two*(one+    nu))
      mustar = mu * syield / predictor
      lamstar = k - two/three*mustar
      fact1 = one*lamstar
      fact2 = two*mustar
      fact3 = (h/(one+h/three/mu)-three*mustar)
      do 10 i=1,ndi
      do 10 j=1,ndi
         ddsdde(i,j) = fact1
 10   continue
      do 20 i=1,ntens
      do 20 j=1,ntens
         ddsdde(i,j) = ddsdde(i,j) + fact2
 20   continue
      do 30 i=1,ntens
      do 30 j=1,ntens
         ddsdde(i,j) = ddsdde(i,j) + fact3 * flow(i) * flow(j)
 30   continue
      end subroutine el_iso_jacobian
c----------------------------------------------------------------------c

c     Place holder for anisotropic elasticity
      subroutine emod_aniso
      return
      end subroutine emod_aniso
