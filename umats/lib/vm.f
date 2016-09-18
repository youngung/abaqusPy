c     Von Mises
      subroutine vm(s,phi)
      implicit none
      real*8 s(3),h,phi,dff,dphi(3),d2phi(3,3),d2h(3,3),d2ff
Cf2py intent(in,out) s
Cf2py intent(out) phi,dphi,d2phi
      h        = 5d-1*((s(1)-s(2))**2+s(1)**2+s(2)**2+6*s(3)**2)
      phi      = h**5d-1        ! yield surface
      s(:)     = s(:)/phi       ! stress on the yield locus
!     analytic solutions of first and second derivatives of the yield surface (plastic potential)
!     dff = 1./(2*(h**0.5))
      dff      = 1d0/(2*phi)
      dphi(1)  = dff*(2*s(1)-s(2))
      dphi(2)  = dff*(2*s(2)-s(1))
      dphi(3)  = dff*6*s(3)

      d2h(1,1) = 2d0
      d2h(1,2) =-1d0
      d2h(2,2) = 2d0
      d2h(3,3) = 6d0

      d2ff     = -(phi**(-3d0))/4
      d2phi(1,1) = d2ff*dphi(1)*dphi(1) + dff*d2h(1,1)
      d2phi(1,2) = d2ff*dphi(1)*dphi(2) + dff*d2h(1,2)
      d2phi(1,3) = 0d0
      d2phi(2,1) = d2phi(1,2)
      d2phi(2,2) = d2ff*dphi(2)*dphi(2) + dff*d2h(2,2)
      d2phi(2,3) = 0d0
      d2phi(3,1) = 0d0
      d2phi(3,2) = 0d0
      d2phi(3,3) = d2ff*dphi(3)*dphi(3) + dff*d2h(3,3)
      dphi(3) = dphi(3)/2
      return                    !! returns phi, dphi, d2phi
      end subroutine vm

c     Calculate Von Mises deviator and flow direction
      subroutine vm_devi_flow(stress,devi,shydro,flow,ntens,ndi)
      integer ntens,ndi
      dimension stress(ntens),devi(ntens),flow(ntens)
      real*8 shydro, smises
      call vm(stress,smises)
      call deviat(stress,devi,shydro)
      do i=1,ndi
         flow(i) = (stress(i)-shydro)/smises
      enddo
      do i=ndi+1,ntens
         flow(i) = stress(i) / smises
      enddo
      return
      end subroutine vm_devi_flow

      include "/home/younguj/repo/abaqusPy/umats/lib/dev.f"
