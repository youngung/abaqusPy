c-----------------------------------------------------------------------
      subroutine vm_shell(cauchy,phi,dphi,d2phi)
      implicit none
      real*8 cauchy(3), phi,dphi(3),d2phi(3,3),
     $     dphi6(6),d2phi66(6,6),cauchy6(6)
      integer i,j

      cauchy6(:) = 0.d0
      cauchy6(1) = cauchy(1)
      cauchy6(2) = cauchy(2)
      cauchy6(6) = cauchy(3)

      call vm_gen(cauchy6,phi,dphi6,d2phi66)
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
      end subroutine vm_shell
c-----------------------------------------------------------------------
c     General Von Mises plastic potential (yield surface)
      subroutine vm_gen(cauchy,phi,dphi,d2phi)
c     cauchy (cauchy stress of 6 Dimension)
c     phi    (the potential)
c     dphi  (1st stress derivative of the potential)
c     d2phi (2nd stress derivative of the potential)
      implicit none
      real*8 cauchy(6), phi, dphi(6), d2phi(6,6), h, dh(6),d2h(6,6),
     $     dff
      integer i,j

      dh(:)      = 0.d0
      d2h(:,:) = 0.d0
      dphi(:)    = 0.d0
      d2phi(:,:) = 0.d0

      h =  (cauchy(1)-cauchy(2))**2+
     $     (cauchy(2)-cauchy(3))**2+
     $     (cauchy(3)-cauchy(1))**2+
     $     6.d0*(cauchy(4)**2+cauchy(5)**2+cauchy(6)**2)
      h = h / 2.d0
      phi = dsqrt(h)

c     round(H)/round(s_i)
      dh(1) = 2.d0*cauchy(1)-cauchy(2)-cauchy(3)
      dh(2) = 2.d0*cauchy(2)-cauchy(1)-cauchy(3)
      dh(3) = 2.d0*cauchy(3)-cauchy(1)-cauchy(2)
      dh(4) = 6.d0*cauchy(4)
      dh(5) = 6.d0*cauchy(5)
      dh(6) = 6.d0*cauchy(6)

c     Only non-zero components
      d2h(1,1) =  2.d0
      d2h(1,2) = -1.d0
      d2h(1,3) = -1.d0
      d2h(2,1) = -1.d0
      d2h(2,2) =  2.d0
      d2h(2,3) = -1.d0
      d2h(3,1) = -1.d0
      d2h(3,2) = -1.d0
      d2h(3,3) =  2.d0
      d2h(4,4) =  6.d0
      d2h(5,5) =  6.d0
      d2h(6,6) =  6.d0
c
c     1st derivatives
c     dphi_i = round(phi) / round(sig_i)
c     dphi_i = round(phi)/ round (H) x  round(H)/round(s_i)
      dff     = 1.d0 / (2.d0*phi)
      dphi(1) =  dff * dh(1)
      dphi(2) =  dff * dh(2)
      dphi(3) =  dff * dh(3)
      dphi(4) =  dff * dh(4)
      dphi(5) =  dff * dh(5)
      dphi(6) =  dff * dh(6)

c     2nd derivatives
c     dphi_i/round(s_j) = round(round(phi)/ round (H) x  round(H)/round(s_i))/round(s_j)
c     dphi_i/round(s_j) = rr(phi)/{round(H)*round(s_j)} x round(H)/round(s_i) +
c                         round(round(phi)/ round (H) x
      do 10 i=1,6
      do 10 j=1,6
         d2phi(i,j) = d2phi(i,j) + 1.d0/4.d0 * phi**(-3.d0)*
     $        dh(j) * dh(i) + dff * d2h(i,j)
 10   continue
      return
      end subroutine vm_gen

c     Calculate Von Mises deviator and flow direction
      subroutine vm_devi_flow(stress,devi,shydro,flow,ntens,ndi)
c     stress: stress tensor
c     devi  : stress deviator
c     shydro: hydrostatic pressure
c     flow  : flow direction by associated flow rule
c     ntens
c     ndi
      implicit none
      integer ntens,ndi
      dimension stress(ntens),devi(ntens),flow(ntens),s6,
     $     dphi(ntens),d2phi(ntens,ntens)
      real*8 shydro,smises,flow,stress,devi
      integer i

      if (ntens.eq.6 .and. ndi.eq.3) then
         call vm_gen(stress,smises,dphi,d2phi)
         call deviat6(stress,devi,shydro)
      elseif (ntens.eq.3.and.ndi.eq.2) then
c        Assuming plane-stress of (s11,s22,s12)
         call vm_shell(stress,phi,dphi,d2phi)
         call deviat3(stress,devi,shydro)
      else
         write(*,*) 'Unexpected case in VM'
         stop -1
      endif

      do i=1,ndi
         flow(i) = (stress(i)-shydro)/smises
      enddo
      do i=ndi+1,ntens
         flow(i) = stress(i) / smises
      enddo

c$$$      write(*,'(a)',advance='no')'----'
c$$$      write(*,'(e13.3)', advance='no') smises
c$$$      write(*,'(3e13.3)',advance='no') (stress(i),i=1,ntens)
c$$$      write(*,'(3e13.3)',advance='no') (flow(i),i=1,ntens)
c$$$      write(*,'(a)',advance='no')'----'
c$$$      write(*,*)

      return
      end subroutine vm_devi_flow
c-----------------------------------------------------------------------
      program test_vm
      implicit none
      integer nth,i,j
      parameter(nth=100)
      real*8 th(nth)



      return
      end program
c-----------------------------------------------------------------------




      include "/home/younguj/repo/abaqusPy/umats/lib/dev.f"
