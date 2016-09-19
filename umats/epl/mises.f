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


c$$$  local arrays
      dimension Cel(ntens,ntens),dphi_n(ntens),d2phi_n(ntens,ntens),
     $     stran_el(ntens),stran_pl(ntens)
      real*8 e,enu,G,eK,Cel
c     predictor stress
      real*8 spr(ntens)
      real*8 voce_params(4)

      real*8 eeq_n,yld_h,dh
!     eeq_n: eq plastic strain at step n
!     yld_h: flow stress at hardening curve
!     dh: dH/dE^eq

      real*8 phi_n,dphi_n,d2phi_n

      real*8 toler_yield,stran_el,stran_pl
      parameter(toler_yield=1e-6)

      voce_params(1) = 479.0
      voce_params(2) = 339.7
      voce_params(3) = 7.784
      voce_params(4) = 70.87

      eeq_n = statev(1) ! previous equivalent plastic strain (at step n)
      do i=1,ntens
         stran_el(i) = statev(i+1)
         stran_pl(i) = statev(i+1+ntens)
      enddo

c     Moduluar pseudo code for stress integration

c     i.   Check the current (given) variables
      call emod_iso_shell(e,enu,G,eK,Cel)
c     ii.  Trial stress calculation
      call mult_array(e,stress,ntens,nshr,spr)
c     iii. See if the pr stress calculation is in the plastic or elastic regime
      call voce(eeq_n, voce_params(1),voce_params(2),voce_params(3),
     $     voce_params(1),yld_h,dh)
      call vm_shell(stress,phi_n,dphi_n,d2phi_n)

c     if elastic,
      if (phi_n<yld_h) then
c              1. Save jacobian as elastic moduli
         ddsdde(:,:) = Cel(:,:)
c              2. Update strain.
         call add_array(stran,dstran,ntens)
c              3. Update stress,
         call multi_array(ddsdde,strans)
c              4. Update other state varaiables (if necessary)
         return
c     elif plastic, go to iv
      else
         call return_mapping(Cel,spr,phi_n,eeq_n,dphi_n,voce_params,
     $        dstran,stran,stran_el,ntens)

      endif

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
c     return_mapping.f
      include "/home/younguj/repo/abaqusPy/umats/lib/return_mapping.f"
