      real*8 function vm(stress,ndi,ntens)
      integer ndi,ntens,i, real
      dimension stress(ntens)
      parameter(six=6.d0,two=2.d0)

      vm = (stress(1)-stress(2))**2 + (stress(2)-stress(3))**2
     $     + (stress(3)-stress(1))**2

      do i=1,ndi+1,ntens
         vm = vm + six*stress(i)**2
      enddo
      vm = dsqrt(vm/two)
      return
      end function vm

c     calculate Von Mises deviator and flow direction
      subroutine vm_devi_flow(stress,devi,hydro,flow,ntens)
      dimension stress(ntens),devi(ntens),flow(ntens)
      real*8 hydro
      enddo
