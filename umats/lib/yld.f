c-----------------------------------------------------------------------
      subroutine yld(iyld_law,yldp,yldc,nyldp,nyldc,
     $     phi,dphi,d2phi)
c-----------------------------------------------------------------------
c***  Arguments
c     iyld_law  : choice of yield function
c     yldp      : state variables associated with yield function
c     yldc      : constants for yield function
c     nyldp     : Len of yldp
c     nyldc     : Len of yldc
c     phi       : yield surface
c     dphi      : 1st derivative of yield surface w.r.t. stress
c     d2phi     : 2nd derivative of yield surface w.r.t. stress
c-----------------------------------------------------------------------
c     intent(in) iyld_law,yldp,yldc,nyldp,nyldc
c     intent(out) phi,dphi,d2phi
c-----------------------------------------------------------------------
      implicit none
      integer iyld_law,nyldp,nyldc,ntens
      dimension yldp(nyldp),yldc(nyldc),dphi(ntens),d2phi(ntens,ntens)
      real*8 yldp,yldc,dphi,d2phi

c***  Local variables for better readibility
      real*8 stress(ntens),strain(ntens)


c***  Restore local variables from yldp
      call restore_yld_statev(
     $     iyld_law,yldp,nyldp,stress,strain,ntens,1)


c***  Define phi,dphi,d2phi
      if (iyld_law.eq.1) then
         call vm_shell(stress,phi,dphi,d2phi)
      else
         write(*,*)'unexpected iyld_law given'
         stop -1
      endif

      end subroutine yld
c-----------------------------------------------------------------------
      subroutine restore_yld_statev(iyld_law,yldp,nyldp,
     $     stress,strain,ntens,iopt)
c***  Arguments
c     iopt: behavior (0: save to yldp; 1: read from yldp)
      implicit none
      integer iyld_law,nyldp,ntens
      dimension yldp(nyldp),stress(ntens),strain(ntens)
      real*8 yldp,stress,strain
      integer i

      if (iyld_law.eq.1) then ! von mises
         if (iopt.eq.0) then
            do 10 i=1,ntens
               stress(i) = yldp(i)
 10         continue
         else
            do 10 i=1,ntens
               yldp(i) = stress(i)
 10         continue
         endif
      else
         write(*,*)'unexpected iyld_law given to restore_yld_statev'
      endif

      end subroutine restore_yld_statev
c-----------------------------------------------------------------------
c     Von Mises
      include "/home/younguj/repo/abaqusPy/umats/lib/vm.f"
