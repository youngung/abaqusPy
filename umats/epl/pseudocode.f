c$$$  ABAQUS UMAT Interface.
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,DRPLDE
     $     ,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
     $     CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     $     CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
c      implicit none  !! to test the namespace in UMAT
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     $     DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     $     TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),DROT(3,3),
     $     DFGRD0(3,3),DFGRD1(3,3),JSTEP(4)
      integer NDI,NSHR,NTENS,NSTATV,NPROPS,NOEL,NPT,LAYER,KSPT,JSTEP,
     $     KINC,kstep
      real*8 STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,DRPLDE,DRPLDT,
     $     STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,PROPS,COORDS,
     $     DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,syield
c-----------------------------------------------------------------------
c$$$  User coding to define DDSDDE, STRESS, STATEV, SSE, SPD, SCD
c$$$  and, if necessary, RPL, DDSDDT, DRPLDE, DRPLDT, PNEWDT
c$$$  End of ABAQUS UMAT Interface
c$$$  implicit none - variables passed through umat


c$$$  local arrays
c     Moduluar pseudo code for stress integration

c     i.   Check the current (given) variables

c     ii.  Trial stress calculation

c     iii. See if the trial stress calculation is in the plastic or elastic regime

c     if elastic,
c              1. Save jacobian as elastic moduli
c              2. Update stress,
c              3. Update strain.
c              4. Update other state varaiables (if necessary)
c     elif plastic, go to iv

c     iv. return mapping (loop over k)
c         1. Find normal of current predictor stress (s_(n+1)^k)
c             save the normal to m_(n+alpha)
c         2. Configure NR condition
c             f   = yield - hardening             (objective function)
c             fp  = r(s^eq_(n+1)^k)/r(s_(n+1)^k) (-C^el : r(s^eq_(n+1)^k / r(s_(n+1)^k))
c         3.  Update the multiplier^(k+1)  (dlamb)
c             dlamb^(k+1) = dlamb^k - f/fp
c             find the new predictor stress (using  dE = dE^(el)^(k+1) + dlamb^(k+1))
c             s_(n+1)^(k+1) = C^e dE^(el)

c     v. Exit from iv. means
c       s_(n+1) is obtained.
c       dlamb   is obtained
c       de^(el)_(n+1) is obtained

c       update plastic strain = depl_ij = depl_ij + n_ij dlamb

c     vi. Caculate jacobian (ddsdde)
c       C^el - [ C^el:m_(n+1) cross C^el:m_(n+1) ]   /  [ m_(n+1):C^el:m_(n+1) + h(depl_ij_(n+1)) ]
      return
      end subroutine umat


c     iso elastic
      include "/home/younguj/repo/abaqusPy/umats/lib/elast.f"
c     UHARD using Voce
      include "/home/younguj/repo/abaqusPy/umats/lib/uhard.f"
c     Von Mises
      include "/home/younguj/repo/abaqusPy/umats/lib/vm.f"
c     lib.f
      include "/home/younguj/repo/abaqusPy/umats/lib/lib.f"
c     diag.f
      include "/home/younguj/repo/abaqusPy/umats/lib/diag.f"
