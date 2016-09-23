c-----------------------------------------------------------------------
c     General USER MAT subroutine that is highly modulized for each
c     individual constitutitve components to allow easy modifications
c     using new material models.
c
c
c     Youngung Jeong, Clemson University
c     youngung.jeong@gmail.com
c-----------------------------------------------------------------------
c$$$  ABAQUS UMAT Interface.
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,DRPLDE
     $     ,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,
     $     CMNAME,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     $     CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
C      INCLUDE 'ABA_PARAM.INC'
      implicit none  ! to test the namespace in UMAT
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     $     DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     $     TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),DROT(3,3),
     $     DFGRD0(3,3),DFGRD1(3,3),JSTEP(4)
      INTEGER NDI,NSHR,NTENS,NSTATV,NPROPS,NOEL,NPT,LAYER,KSPT,JSTEP,
     $     KINC,KSTEP
      REAL*8 STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,DRPLDE,DRPLDT,
     $     STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,PROPS,COORDS,
     $     DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,SYIELD
c-----------------------------------------------------------------------
c$$$  User coding to define DDSDDE, STRESS, STATEV, SSE, SPD, SCD
c$$$  and, if necessary, RPL, DDSDDT, DRPLDE, DRPLDT, PNEWDT
c$$$  End of ABAQUS UMAT Interface

c-----------------------------------------------------------------------
c$$$  local arrays
      character*255 fndia,fnstr
      integer nhrdc,nhrdp
      parameter(nhrdc=4,nhrdp=1)
      dimension Cel(ntens,ntens),phi_ns(0:1),dphi_n(ntens),
     $     d2phi_n(ntens,ntens),stran_el_ns(0:1,ntens),
     $     stran_pl_ns(0:1,ntens),dstran_el(ntens),dstran_pl(ntens),
     $     spr(ntens),epr(ntens),hrdc(nhrdc),hrdp(nhrdp),
     $     stress_ns(0:1,ntens)
      real*8 e,enu,G,eK,Cel
c     predictor stress
      real*8 spr,epr,stress_ns
      real*8 hrdc,hrdp
      dimension eeq_ns(0:1)
      real*8 eeq_ns,flow_stress,dflow_stress
!     eeq_ns: eq plastic strain at steps n,n+1

      real*8 phi_ns,dphi_n,d2phi_n
      real*8 stran_el_ns,stran_pl_ns,dstran_el,dstran_pl
      integer imsg,idia,i,istv,ihrd_law,iyld_law
      logical idiaw,failnr
      real*8 empa,gpa
      parameter(empa=1d6,gpa=1d9)

c$$$  yld function parameters
      integer nyldp,nyldc
      parameter(nyldp=1,nyldc=4) ! this depends on iyld_law...
      dimension yldp_ns(0:1,nyldp),yldc(nyldc)
      real*8 yldp_ns,yldc,deq_pr,deq


c$$$  Material data
      dimension rvs(4)
      real*8 rvs

c-----------------------------------------------------------------------
c***  hardwired material parameters
c-----------------------------------------------------------------------
c**   material constitutive laws
      ihrd_law=0                ! Voce isotropic hardening

c      iyld_law=0                ! Generic von Mises (shell)
      iyld_law=1                ! Generic Hill48 (shell)

      rvs(1)=1.0d0
      rvs(2)=2.2d0
      rvs(3)=1.4d0
      call tuneH48(rvs,yldc)
c-----------------------------------------------------------------------
c**   hardening parameters
      hrdc(1) = 479.0d0
      hrdc(2) = 339.7d0
      hrdc(3) = 7.784d0
      hrdc(4) = 70.87d0
c**   elastic constants
      e = 210d0*gpa
      enu=0.3d0
c-----------------------------------------------------------------------
      imsg=7
      idia=315 ! 0 (in case stdo is preferred)
      istv=425 ! reserved for state variable output file

      stress_ns(0,:) = stress(:)

      idiaw=.false.
!cc   if (kspt.eq.1 .and. noel.eq.1 .and. npt.eq.1.) then
      if (idiaw) then
         fndia='/home/younguj/repo/abaqusPy/examples/one/diagnose.txt'
         if (idia.ne.0) then
            open(idia,position='append',file=fndia)
         endif
         call print_head(idia)
      endif
c     print head

c     restore stran_el, stran_pl, and yldp
      call restore_statev(statev,nstatv,eeq_ns(0),stran_el_ns(0,:),
     $     stran_pl_ns(0,:),ntens,yldp_ns(0,:),nyldp,0,.false.,idia,
     $     .false.,kinc,noel,npt,time(0),stress)

c     yld parameters pertaining to step n - need to find
c     the updated yld parameters for step n+1 later...

      if (idiaw) then
         call w_ival(idia,'* noel :', noel)
         call w_ival(idia,'* kinc :', kinc)
         call w_ival(idia,'* npt  :', npt)
         call w_ival(idia,'* kspt :', kspt)
         call w_ival(idia,'* layer:', layer)
         call w_val(idia,'* dtime:',dtime)
         call w_val(idia,'* time1:',time(1))
         call w_val(idia,'* time2:',time(2))
         call w_chr(idia,'* JSTEP:')
         call w_dim(idia,jstep,4,1d0,.true.)
         call w_chr(idia,'** state variable restored **')
         call w_chr(idia,'* stran at step n')
         call w_dim(idia,stran,ntens,1d0,.true.)
         call w_chr(idia,'* stran_el at step n')
         call w_dim(idia,stran_el_ns(0,:),ntens,1d0,.true.)
         call w_chr(idia,'* stran_pl at step n')
         call w_dim(idia,stran_pl_ns(0,:),ntens,1d0,.true.)
         call w_chr(idia,'* dstran')
         call w_dim(idia,dstran,ntens,1d0,.true.)
         call w_chr(idia,'* DROT')
         call w_mdim(idia,drot,3,1d0)
      endif

      !call stop_debug(0)        ! debug

c     Moduluar pseudo code for stress integration
c-----------------------------------------------------------------------
c     i.   Check the current (given) variables
      call emod_iso_shell(e,enu,G,eK,Cel)
c-----------------------------------------------------------------------
c     ii.  Trial stress calculation
c     Calculate predict elastic strain - assuming dstran = dstran_el
      call add_array2(stran_el_ns(0,:),dstran,epr,ntens)
      call mult_array(Cel,epr,ntens,stress_ns(1,:))
      spr(:) = stress_ns(1,:)

      if (idiaw) then
         call w_chr( idia,'* C^el [GPa]')
         call w_mdim(idia,Cel,ntens,1d0/gpa)
         call w_chr( idia,'* Trial elastic strain ')
         call w_dim( idia,epr,ntens,1d0,.true.)
         call w_chr( idia,'* Trial stress [MPa]')
         call w_dim( idia,stress_ns(1,:),ntens,1d0/empa,.true.)
      endif
c-----------------------------------------------------------------------
c     iii. See if the pr stress (stress_ns(1,:)) calculation is in the plastic or
c          elastic regime
      deq_pr = 0.               !assuming no plastic increment
      hrdp(1) = eeq_ns(0)+deq_pr
      call uhard(ihrd_law,hrdp,nhrdp,hrdc,nhrdc,
     $     flow_stress,dflow_stress,empa)
      if (idiaw) then
         call w_val(idia,'* E^eq pl          :',eeq_ns(0))
         call w_val(idia,'* flow_stress [MPa]:',flow_stress/empa)
         call w_val(idia,'* dflow       [MPa];',dflow_stress/empa)
         call w_val(idia,'yldp_ns(0,1):',yldp_ns(0,1))
      endif

      call update_yldp(iyld_law,yldp_ns,nyldp,deq_pr)
      if (idiaw) call w_val(idia,'yldp_ns(1):',yldp_ns(1,1))
c     predicting yield surface
      call yld(iyld_law,yldp_ns(1,:),yldc,nyldp,nyldc,stress_ns(1,:),
     $     phi_ns(0),dphi_n,d2phi_n,ntens)
      if (idiaw) then
         call w_val(idia,'*         phi step n [MPa]:',
     $        phi_ns(0)/empa)
         call w_val(idia,'* flow stress step n [MPa]:',
     $        flow_stress/empa)
      endif

c-----------------------------------------------------------------------
      if (phi_ns(0).lt.flow_stress) then ! elastic
         call update_elastic(idia,idiaw,iyld_law,ntens,nyldp,nstatv,
     $        ddsdde,cel,stran,stran_el_ns,stran_pl_ns,dstran,stress,
     $        eeq_ns,deq,yldp_ns,statev,kinc)
      else !! plastic
         if (idiaw) then
            call w_chr(idia,'PLASTIC')
            call print_foot(idia)
         endif
c-----------------------------------------------------------------------
c     vi. Return mapping
c        Return mapping subroutine updates stress/statev
         if (idiaw) then
            call w_chr(idia,'** Inquiry for state variable before RM')
            call restore_statev(statev,nstatv,eeq_ns(1),
     $           stran_el_ns(1,:),stran_pl_ns(1,:),ntens,yldp_ns(1,:),
     $           nyldp,0,.true.,idia,
     $           .false.,kinc,noel,npt,time(0),stress)
         endif
         if (idiaw) call w_chr(idia,'** Begin return-mapping **')
         call return_mapping(Cel,stress_ns(1,:),phi_ns(0),eeq_ns(0),
     $        dphi_n,dstran,stran_el_ns(0,:),stran_pl_ns(0,:),
     $        ntens,idiaw,idia,hrdp,nhrdp,hrdc,nhrdc,ihrd_law,
     $        iyld_law,yldc,nyldc,yldp_ns,nyldp,
     $        spd,stress,statev,nstatv,ddsdde,failnr,kinc)
         if (idiaw) call w_chr(idia,'** Exit return-mapping **')
         if (failnr) then
            ! reduce time step?
            pnewdt = 0.5d0
            if (idiaw) then
               call w_val(idia,'** dtime :',dtime)
               call w_val(idia,'** pnewdt:',pnewdt)
               call fill_line(idia,'*',72)
            endif
            return ! out of umat
         endif

         stress_ns(1,:) = stress(:)

         if (idiaw) then
            call w_chr(idia,'** inquiry for state variable after RM')
            call restore_statev(statev,nstatv,eeq_ns(1),
     $           stran_el_ns(1,:),stran_pl_ns(1,:),ntens,yldp_ns(1,:),
     $           nyldp,0,.true.,idia,
     $           .false.,kinc,noel,npt,time(0),stress)
            call w_chr( idia,'* stress at n   [MPa]')
            call w_dim( idia,stress_ns(0,:),ntens,1d0/empa,.true.)
            call w_chr( idia,'* spr           [MPa]')
            call w_dim( idia,spr,ntens,1d0/empa,.true.)
            call w_chr( idia,'* stress at n+1 [MPa]')
            call w_dim( idia,stress_ns(1,:),ntens,1d0/empa,.true.)
            call w_chr( idia,'* stress        [MPa]')
            call w_dim( idia,stress(:),ntens,1d0/empa,.true.)

            call w_chr( idia,'* stran_pl at n    :')
            call w_dim( idia,stran_pl_ns(0,:),ntens,1d0,.true.)
            call w_chr( idia,'* stran_pl at n+1  :')
            call w_dim( idia,stran_pl_ns(1,:),ntens,1d0,.true.)
            call w_chr( idia,'* stran n+1        :')
            call w_dim( idia,stran,ntens,1d0,.true.)
            call w_chr( idia,'* ddsdde     [GPa] :')
            call w_mdim(idia,ddsdde,ntens,1d0/gpa)
            call w_chr( idia,'* stran_el_ns(n+1) :')
            call w_dim( idia,stran_el_ns(1,:),ntens,1d0,.true.)
         endif

         !call stop_debug(0)     ! debug

c-----------------------------------------------------------------------
c     v. Exit from iv. means
c        s   _(n+1) is obtained.
c        e^pl_(n+1) is obtained.
c        e^el_(n+1) is obtained.
c       de^pl_(n+1) is obtained.
c       de^el_(n+1) is obtained.
c
c       Update all state variables. (yldp included.)
c-----------------------------------------------------------------------
c     vi. Caculate jacobian (ddsdde)
c       C^el - [ C^el:m_(n+1) cross C^el:m_(n+1) ]   /  [ m_(n+1):C^el:m_(n+1) + h(depl_ij_(n+1)) ]
      endif

c     Write statev
      if (noel.eq.1 .and. dtime.gt.0 .and. npt.eq.1) then
         call restore_statev(statev,nstatv,eeq_ns(1),
     $        stran_el_ns(1,:),stran_pl_ns(1,:),ntens,yldp_ns(1,:),
     $        nyldp,0,.false.,idia,
     $        .true.,kinc,noel,npt,time(1),stress)
      endif

      if (idia.ne.0.and.idiaw) close(idia)
      if (idiaw) call print_foot(idia)
      return
      end subroutine umat
c-----------------------------------------------------------------------
c     iso elastic
      include "/home/younguj/repo/abaqusPy/umats/lib/elast.f"
c     elastic update
      include "/home/younguj/repo/abaqusPy/umats/lib/update_elastic.f"
c     UHARD using Voce
      include "/home/younguj/repo/abaqusPy/umats/lib/uhard.f"
c     lib.f
      include "/home/younguj/repo/abaqusPy/umats/lib/lib.f"
c     lib_w.f
      include "/home/younguj/repo/abaqusPy/umats/lib/lib_write.f"
c     return_mapping.f
      include "/home/younguj/repo/abaqusPy/umats/lib/return_mapping.f"
c     is.f - testing subroutines/functions (e.g., is_inf)
      include "/home/younguj/repo/abaqusPy/umats/lib/is.f"
c     yld.f - yield function associated calculations
      include "/home/younguj/repo/abaqusPy/umats/yld/yld.f"
c     restore_state.f - save/read from statev
      include "/home/younguj/repo/abaqusPy/umats/lib/restore_statev.f"
