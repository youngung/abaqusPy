c$$$  ABAQUS UMAT Interface.
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1     RPL,DDSDDT,DRPLDE,DRPLDT,
     2     STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3     NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4     CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
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
      
      RETURN
      END

c$$$  CMNAME: User material name - useful in case there are multiple types of user materials discussed
c$$   in the umat subroutine

c$$$  STRESS: stress tensor array (vectorized)
c$$$  This array is paassed as the stress tensor at the beginning of the increment
c$$$  This array should be updated to be the one at the end of increment.
c$$$  i.e., STRESS(:) = dSTRESS(:)*dtime


c$$$  NDI: Number of direct (normal) stress components
c$$$  NSHR: Number of shear stress component
c$$$  NTENS = NDI + NSHR

c$$$  NPROPS: User-defined number of material constants associated with this user material
c$$$  PROPS(NPROPS): Array of material constants


c$$$  DDSDDE: Jacobian matrix of the constitutive model: round incr(sigma) / round incr(epsilon)
c$$$  DDSDDE(ntens,ntens)
      DDSDDE(i,j) defines the change in ith stress component at the end of the time increment caused by an infinitesimal

c$$$  STATEV(NSTATV) an array containing the solution-dependent state variables.

c$$$  STRANS(NSTENS): Array of total strains at the beinning of the increment.
c$$$  DSTRANS(NSTENS): Array of strain increments.

c$$$  DTIME: Time increment
c$$$  TIME(1): value of step time at the beginning of the current increment or frequency
c$$$  TIME(2): value of total time at the beinning of the current increment

c$$$  PNEWDT: ratio of suggested new time increment tothe current time increment passed (i.e., DTIME)
c$$$  if PNEWDT<1, PNEWDT*DTIME will be the new time increment

c$$$  NOEL: Element number

c$$$  NPT : Integration point number

c$$$  JSTEP(1): Step Number

c$$$  JSTEP(2): Procedure type key
