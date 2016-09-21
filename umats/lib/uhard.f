c-----------------------------------------------------------------------
c     custom hardening law manager...
      subroutine uhard(ihard_law,hrdp,nhrdp,hrdc,nhrdc,flow_stress,
     $     dflow_stress,fact)
c     intent(in) ihard_law,hrdp,nhrdp,hrdc,nhrdc, fact
c     intent(out) flow_stress,dflow_stress

c     ihard_law: hardening law
c     hrdp     : state parameter for hardening
c     hrdc     : hardening constants (invariant)
c      flow_stress: flow stress
c     dflow_stress: slope of flow stress (df/de^eq)
c     fact        : multiplicative factor applied to flow/dflow
      implicit none
      integer ihard_law,nhrdp,nhrdc
      dimension hrdp(nhrdp),hrdc(nhrdc)
      real*8 flow_stress,dflow_stress,hrdp,hrdc,fact
c-----------------------------------------------------------------------
      if (ihard_law.eq.1) then  ! voce
c        hrdp(1) : equivalent plastic strain
         call voce(hrdp(1),hrdc(1),hrdc(2),hrdc(3),hrdc(4),flow_stress,
     $        dflow_stress)
         flow_stress  =  flow_stress*fact
         dflow_stress = dflow_stress*fact
      else
         write(*,*) 'Err: Unexpected hardening law given in uhard.f'
         stop -1
      endif

      return
      end subroutine
c-----------------------------------------------------------------------
c     Voce - rate independent
      subroutine voce(e,a,b0,c,b1,sig,dsig)
      implicit none
      real*8 e,a,b0,c,b1,sig,dsig
!     Voce
      sig  = a - b0 * exp(-c*e) + b1*e
      dsig = c*b0*exp(-c*e) + b1
      return
      end subroutine voce
c-----------------------------------------------------------------------
