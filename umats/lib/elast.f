c----------------------------------------------------------------------c
c     Calculate complete 6x6 ISOTROPIC elastic modulus (c) using
c     Young's modulus (e) and poisson ratio (nu)
c----------------------------------------------------------------------c

c$$$  3D shell with S11,S22 and S12
      subroutine emod_iso_shell(e,nu,c)
      real*8 c(3,3)
      real*8 nu, e
      integer i,j
      c(:,:)=0.d0
c     Multiplier
      x = e/(1.+nu)/(1.-2.*nu)
      do i=1,2
         do j=1,2
            c(i,j) = x*nu
         enddo
         c(i,i) = x*(1.-nu)
      enddo
      c(3,3) = x* (1.-2.*nu)/2.
      return
      end subroutine

c$$$  3D shell with S11,S22 and S12
      subroutine emod_iso(e,nu,c,ndi,nshr)
c     intent(in): e,nu,ndi,nshr
c     intent(out): c
      integer ndi,nshr
      real*8 e, nu, c(ndi+nshr,ndi+nshr), x
      integer i,j,imsg

      imsg=7

c     initialization
      do i=1,ntens
         do j=1,ntens
            c(i,j) = 0.d0
         enddo
      enddo

      write(imsg,*) 'after initialization c matrix'

c
c     construct elastic tensor (6x6) with assuming that
c     \gamma_ij = 2\varepsilon_ij is the engineering shear strain
c

c     Multiplier
      x = e/(1.+nu)/(1.-2.*nu)

c     off-diagonal terms
      do i=1,3
         do j=1,3
            c(i,i) = nu * x
         enddo
      enddo

      do i=1,3
         c(i,i) = (1.-nu)*x     !! overwrite the diganogal term
         c(i+3,i+3) = (1.-2.*nu)/2. * x
      enddo

      write(imsg,*) 'just before returning'

      return
      end subroutine emod_iso
c----------------------------------------------------------------------c
      subroutine el_iso_jacobian(e,nu,hard,flow,syield,predictor,
     $     ndi,ntens,ddsdde)
c     intent(in): e,nu,hard,flow,syield,ndi,ntens
c     intent(out): ddsdde
      parameter(one=1.d0,two=2.d0,three=3.d0)
      dimension flow(ntens), ddsdde(ntens,ntens), c(ntens,ntens)
      real*8 e, nu, hard, syield, predictor
      real*8 mu, k              !! elastic constants (shear modulus and bulk modulus)
      integer i,j,ndi,ntens
      real*8 mustar, lamstar, fact1,fact2,fact3
      call emod_iso(e,nu,c)
      k=e/(three*(one-two*nu))
      mu=e/(  two*(one+    nu))

      mustar = mu * syield / predictor
      lamstar = k - two/three*mustar
      fact1 = one*lamstar
      fact2 = two*mustar
      fact3 = (h/(one+h/three/mu)-three*mustar)

      do i=1,ndi
         do j=1,ndi
            ddsdde(i,j) = fact1
         enddo
      enddo

      do i=1,ntens
         do j=1,ntens
            ddsdde(i,j) = ddsdde(i,j) + fact2
         enddo
      enddo

      do i=1,ntens
         do j=1,ntens
            ddsdde(i,j) = ddsdde(i,j) + fact3 * flow(i) * flow(j)
         enddo
      enddo



c$$$
c$$$
c$$$      ebulk3  = e / (one - two * nu)  !  e /  (1-2nu) : bulk mod x 3
c$$$      eg2     = e/(one+nu)            !  e /  (1+ nu) :  shr mod x 2
c$$$      eg      = eg2/two               !  e /{2(1+ nu)}:  shr mod
c$$$      eg3     = three *eg             ! 3e /{2(1+ nu)}:  shr mod x 3
c$$$      elam    = (ebulk3-eg2)/three    !  Lame constant lambda
c$$$
c$$$      !      = (3 mu x h) / (3 mu + h) -
c$$$      effhrd = eg3 * hard / (eg3+hard) - effg3
c$$$
c$$$      do i=1,ndi
c$$$         do j=1,ndi
c$$$            ddsdde(j,i) = efflam
c$$$         enddo
c$$$         ddsdde(i,i) = effg2+efflam
c$$$      enddo
c$$$      do i=ndi+1,ntens
c$$$         ddsdde(i,i) = effg
c$$$      enddo
c$$$      do i=1,ntens
c$$$         do j=1, ntens
c$$$            ddsdde(j,i) = ddsdde(j,i) +  effhrd*flow(j)*flow(i)
c$$$         enddo
c$$$      enddo
c$$$

      end subroutine el_iso_jacobian
c----------------------------------------------------------------------c

c     Place holder for anisotropic elasticity
      subroutine emod_aniso
      return
      end subroutine emod_aniso

