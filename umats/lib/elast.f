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
c----------------------------------------------------------------------c
      subroutine emod_iso(e,enu,cel,ndi,nshr)
c     Arguments
c     ---------
c     e   : Young's modulus
c     enu : Poisson ratio
c     Cel : Elastic modulus
c     ndi : Number of direct components
c     nshr: Number of shear components
c     Obtain elastic modulus tensor using
c       1. Young's modulus (e)
c       2. poisson ratio   (enu)
c     sigma_i = Cel_ij g_j
c     where e_j is the engineering strain for which the shear strain components
c     are such that g_j = 2 * e_j (here j index loops over only shear components)
      implicit none
      integer, intent(in) :: ndi,nshr
      dimension cel(ndi+nshr,ndi+nshr)
      real*8, intent(in) :: e,enu
      real*8, intent(out) :: cel
      real*8 x
      integer i,j,imsg,ntens
c     intent(in): e,enu,ndi,nshr
c     intent(out): cel
      ntens=ndi+nshr
      imsg=7

c     initialization
      cel(:,:)=0d0
c      write(imsg,*) 'after initialization c matrix'
c
c     construct elastic tensor (6x6) with assuming that
c     \gamma_ij = 2\varepsilon_ij is the engineering shear strain
c
c     Multiplier
      x = e/(1d0+enu)/(1d0-2d0*enu)
c     off-diagonal terms for 'direct' components
      do 50 i=1,ndi
      do 50 j=1,ndi
         cel(i,j) = enu * x
 50   continue
      do 60 i=1,ndi
         cel(i,i) = (1d0-enu)*x !! overwrite the diganogal term
 60   continue
      do 70 i=ndi+1,ndi+nshr
         cel(i,i) = (1d0-2d0*enu)/2d0 * x
 70   continue
      call fill_line(0,'-',40)
      call w_val(0,'mod:',e)
      call w_val(0,'nu:',enu)
      call w_chr(0,'Cel:')
      call w_mdim(0,cel,ntens,1d0)
      call fill_line(0,'-',40)
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
