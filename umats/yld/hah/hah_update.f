c-----------------------------------------------------------------------
c     subroutine used to update hah parameter updates.
c-----------------------------------------------------------------------
      subroutine hah_update(iyld_choice,yldp_ns,nyldp,deeq)
c     Arguments
c     iyld_choice   : choice of yield surface kernel
c     yldp_ns       : yield surface parameters including HAH parameters
c     nyldp         : len of 2nd axis of yldp_ns
c     deeq          : incremental value of equivalent plastic strain

      implicit none
      integer nyldp
      dimension yldp_ns(0:1,nyldp)
      real*8 yldp_ns,deeq
      integer iyld_choice



c**   to suppress -wnused-dummy-argument
      write(*,*) deeq,iyld_choice,yldp_ns


      return
      end subroutine hah_update
