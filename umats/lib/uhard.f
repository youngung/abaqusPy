c     user subrutine to define yield surface size
c     and hardening parameters for isotropic plasticity or
c     combined hardening models
      SUBROUTINE UHARD(SYIELD,HARD,EQPLAS,EQPLASRT,TIME,DTIME,TEMP,
     1     DTEMP,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,CMNAME,NSTATV,
     2     STATEV,NUMFIELDV,PREDEF,DPRED,NUMPROPS,PROPS)
C
      implicit none
c      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      DIMENSION HARD(3),STATEV(NSTATV),TIME(*),
     $     PREDEF(NUMFIELDV),DPRED(*),PROPS(*)
      real*8 SYIELD,HARD,EQPLAS,EQPLASRT,TIME,DTIME,TEMP,DTEMP,STATEV,
     $     PREDEF,DPRED,PROPS
      integer NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NSTATV,NUMFIELDV,NUMPROPS

c     local variables
      real*8 a,b0,c,b1,empa,dsig
      integer imsg
c$$$      integer imsg,nodel,npt,layer,kstp,kstep,kinc,nstatv
c$$$      integer numfieldv,numprops

      imsg = 7
      empa = 1.e6

c     user coding to define SYIELD,HARD(1),HARD(2),HARD(3)
      hard(:)=0.d0

c     hard(1): \partial(syield)/\partial(equivalent plastic strain)
c     hard(2): \partial(syield)/\partial(equivalent plastic strain rate)
c     hard(3): \partial(syield)/\partial(temperature)

c     hard(1)=0 !! default

c     IF steel: 303.23, 273.35, 11.95, 191.29
      a=479.00408591
      b0=339.71480479
      c=7.68395984
      b1=70.86783572

c$$$      a  = 403.23d0
c$$$      b0 = 273.35d0
c$$$      c  =  11.95d0
c$$$      b1 = 191.29d0

c$$$      write(imsg,*) 'Voce params:'
c$$$      write(imsg,*) a, b0, c, b1
      dsig=0.d0
      call voce(eqplas,a,b0,c,b1,syield,dsig)
      syield = syield * empa
c$$$      write(7,*)'dsig:::::  ',dsig
      hard(1) = dsig*1.d0 * empa
c$$$      write(imsg,*) 'syield :',syield
c$$$      write(imsg,*) 'hard(1):',hard(1)

      RETURN
      END SUBROUTINE UHARD

c     Voce - rate independent
      subroutine voce(e,a,b0,c,b1,sig,dsig)
      implicit none
      real*8 e,a,b0,c,b1,sig,dsig
!     Voce
      sig  = a - b0 * exp(-c*e) + b1*e
      dsig = c*b0*exp(-c*e) + b1
c$$$      write(7,*)'strain:',e
c$$$      write(7,*)'in voce'
c$$$      write(7,*)'c,b0,c,b1'
c$$$      write(7,*)c,b0,c,b1
c$$$      write(7,*)'dsig:'
c$$$      write(7,*) c*b0*exp(-c*e) + b1
c$$$      write(7,'(a)') '---===---'
      return
      end subroutine voce
