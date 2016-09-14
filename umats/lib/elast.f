c----------------------------------------------------------------------c
      subroutine emod_iso(e,nu,c)
c     intent(in): e,nu
c     intent(out): c
      parameter(ndi=3,ntens=6)
      real*8 e, nu, c(ntens,ntens), x
      integer i,j

c     initialization
      do i=1,ntens
         do j=1,ntens
            c(i,j) = 0.d0
         enddo
      enddo
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

      return
      end subroutine emod_iso
c----------------------------------------------------------------------c

c     Place holder for anisotropic elasticity
      subroutine emod_aniso
      return
      end subroutine emod_aniso

