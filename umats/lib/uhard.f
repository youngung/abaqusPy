c     user subrutine to define yield surface size
c     and hardening parameters for isotropic plasticity or
c     combined hardening models
      SUBROUTINE UHARD(SYIELD,HARD,EQPLAS,EQPLASRT,TIME,DTIME,TEMP,
     1     DTEMP,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,CMNAME,NSTATV,
     2     STATEV,NUMFIELDV,PREDEF,DPRED,NUMPROPS,PROPS)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION HARD(3),STATEV(NSTATV),TIME(*),
     $     PREDEF(NUMFIELDV),DPRED(*),PROPS(*)
      real*8 a,b0,c,b1

c     user coding to define SYIELD,HARD(1),HARD(2),HARD(3)

c     hard(1): \partial(syield)/\partial(equivalent plastic strain)
c     hard(2): \partial(syield)/\partial(equivalent plastic strain rate)
c     hard(3): \partial(syield)/\partial(temperature)

c     hard(1)=0 !! default


c     IF steel: 303.23, 273.35, 11.95, 191.29
      a = 403.23d6
      b0= 273.35d6
      c = 11.95d6
      b1=191.29d6
      call voce(eqplas,a,b0,c,b1, syield, hard(1))
      RETURN
      END SUBROUTINE UHARD

c     Voce - rate independent
      subroutine voce(e,a,b0,c,b1,sig,dsig)
      implicit none
      real*8 e,a,b0,c,b1,sig,dsig
!     Voce
      sig  = a - b0 * dexp(-c*e) + b1*e
      dsig = c * b0 * dexp(-c*e) + b1
      return
      end subroutine voce
