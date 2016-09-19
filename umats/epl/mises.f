c$$$  ABAQUS UMAT Interface.
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1     RPL,DDSDDT,DRPLDE,DRPLDT,
     2     STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3     NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4     CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
c      implicit none  !! to test the namespace in UMAT
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1     DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2     STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3     PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4     JSTEP(4)

c$$$  user coding to define DDSDDE, STRESS, STATEV, SSE, SPD, SCD
c$$$  and, if necessary, RPL, DDSDDT, DRPLDE, DRPLDT, PNEWDT
c$$$  End of ABAQUS UMAT Interface


c$$$  implicit none - variables passed through umat
c-----------------------------------------------------------------------
      integer NDI,NSHR,NTENS,NSTATV,NPROPS,NOEL,NPT,LAYER,KSPT,JSTEP,
     $     KINC,kstep
      real*8 STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,DRPLDE,DRPLDT,
     $     STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,PROPS,COORDS,
     $     DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,syield

c$$$  local arrays
      dimension eelas(ntens),eplas(ntens),flow(ntens),sdeviator(ntens),
     $     hard(3),stress_old(ntens),stress_new(ntens)
      character*255 fndia,cwd,fnstr
      real*8 zero,one,two,three,six,enumax,toler
      real*8 eqplas,syiel0,hard,sdeviator,empa,gpa,emod,enu,f,fp,eg3,
     $     deqpl,smises,shydro,flow,G,ekappa,emus,hs,elabs,eelas,
     $     eplas,eqplasrt,dphi(ntens),d2phi(ntens,ntens)
      integer imsg,idia,istr,k,i,j,NUMFIELDV,numprops,kewton,newton
      logical isnan,idiaw,istrw,isnan_in_marr
      parameter(zero=0.d0,one=1d0,two=2d0,three=3d0,six=6d0,
     $     enumax=0.4999d0,newton=10,toler=1d-6)
      empa=1.e6
      gpa=1.e9

c     Moduluar pseudo code for stress integration

c     i.   Check the current (given) variables

c     ii.  Trial stress calculation

c     iii. See if the trial stress calculation is in the plastic or elastic regime
c          if plastic go to iv
c          else elastic,
c              1. Save jacobian as elastic moduli
c              2. Update stress,
c              3. Update strain.
c              4. Update other state varaiables (if necessary)

c     iv. return mapping (look over k)
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
