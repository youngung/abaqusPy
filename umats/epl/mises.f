c-----------------------------------------------------------------------
c     General USER MAT subroutine that is highly modulized for each
c     individual constitutitve components to allow easy modifications
c     using new material models.

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
      integer imsg,idia,i,istr,ihrd_law,iyld_law
      logical idiaw,failnr
      real*8 empa,gpa
      parameter(empa=1d6,gpa=1d9)

c$$$  yld function parameters
      integer nyldp,nyldc
      parameter(nyldp=1,nyldc=1) ! this depends on iyld_law...
      dimension yldp_ns(0:1,nyldp),yldc(nyldc)
      real*8 yldp_ns,yldc,deq_pr,deq

c-----------------------------------------------------------------------
c***  hardwired material parameters
c-----------------------------------------------------------------------
c**   material constitutive laws
      ihrd_law=0                ! Voce isotropic hardening
      iyld_law=0                ! von Mises (shell)
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
      idia=315
      istr=425

      stress_ns(0,:) = stress(:)


      idiaw=.true.
!cc   if (kspt.eq.1 .and. noel.eq.1 .and. npt.eq.1.) then
      if (idiaw) then
         fndia='/home/younguj/repo/abaqusPy/examples/one/diagnose.txt'
         fnstr='/home/younguj/repo/abaqusPy/examples/one/strstr.txt'
         if (idia.ne.0) then
            open(idia,position='append',file=fndia)
         endif
         if (istr.ne.0) then
            open(istr,position='append',file=fnstr)
         endif
         call print_head(idia)
      endif
c     print head

c     restore stran_el, stran_pl, and yldp
      call restore_statev(statev,nstatv,eeq_ns(0),stran_el_ns(0,:),
     $     stran_pl_ns(0,:),ntens,yldp_ns(0,:),nyldp,0,.false.,idia)

c     yld parameters pertaining to step n - need to find
c     the updated yld parameters for step n+1 later...

      if (idiaw) then
         call w_val(idia,'* noel :', noel)
         call w_val(idia,'* kinc :', kinc)
         call w_val(idia,'* npt  :', npt)
         call w_val(idia,'* kspt :', kspt)
         call w_val(idia,'* layer:', layer)
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
      endif

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
         call w_val(idia,'% E^eq pl          :',eeq_ns(0))
         call w_val(idia,'% flow_stress [MPa]:',flow_stress/empa)
         call w_val(idia,'% dflow       [MPa];',dflow_stress/empa)
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
     $        eeq_ns,deq,yldp_ns,statev)
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
     $           nyldp,0,.true.,idia)
         endif
         if (idiaw) call w_chr(imsg,'** Begin return-mapping **')
         call return_mapping(Cel,stress_ns(1,:),phi_ns(0),eeq_ns(0),
     $        dphi_n,dstran,stran_el_ns(0,:),stran_pl_ns(0,:),
     $        ntens,idiaw,idia,hrdp,nhrdp,hrdc,nhrdc,ihrd_law,
     $        iyld_law,yldc,nyldc,yldp_ns,nyldp,
     $        spd,stress,statev,nstatv,ddsdde,failnr)
         if (idiaw) call w_chr(imsg,'** Exit return-mapping **')
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
     $           nyldp,0,.true.,idia)
            call w_chr( idia,'* stress at n')
            call w_dim( idia,stress_ns(0,:),ntens,1d0/empa,.true.)
            call w_chr( idia,'* spr')
            call w_dim( idia,spr,ntens,1d0/empa,.true.)
            call w_chr( idia,'* stress at n+1')
            call w_dim( idia,stress_ns(1,:),ntens,1d0/empa,.true.)
            call w_chr( idia,'* stress')
            call w_dim( idia,stress(:),ntens,1d0/empa,.true.)

            call w_chr( idia,'* stran_pl at n')
            call w_dim( idia,stran_pl_ns(0,:),ntens,1d0/empa,.true.)
            call w_chr( idia,'* stran_pl at n+1')
            call w_dim( idia,stran_pl_ns(1,:),ntens,1d0/empa,.true.)
            call w_chr( idia,'* stran n+1')
            call w_dim( idia,stran,ntens,1d0/empa,.true.)
            call w_chr( idia,'* ddsdde')
            call w_mdim(idia,ddsdde,  ntens,1d0/gpa)
            call w_chr( idia,'* stran_el_ns(n+1)')
            call w_dim( idia,stran_el_ns(1,:),ntens,1d0,.true.)
         endif

         call stop_debug(0) ! debug

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
      if (idia.ne.0.and.idiaw) close(idia)
      close(istr)

      if (idiaw) then
         call print_foot(idia)
      endif
      return
      end subroutine umat
c-----------------------------------------------------------------------
      subroutine update_elastic(
     $     idia,idiaw,iyld_law,
     $     ntens,nyldp,nstatv,
     $     ddsdde,cel,
     $     stran,stran_el_ns,stran_pl_ns,dstran,stress,eeq_ns,deq,
     $     yldp_ns,
     $     statev)
      implicit none
      integer idia,iyld_law,nstatv,ntens,nyldp
      logical idiaw
      dimension ddsdde(ntens,ntens),cel(ntens,ntens),stran(ntens),
     $     stran_el_ns(0:1,ntens),stran_pl_ns(0:1,ntens),dstran(ntens),
     $     stress(ntens),yldp_ns(0:1,nyldp),statev(nstatv),eeq_ns(0:1)
      real*8 ddsdde,cel,stran,stran_el_ns,stran_pl_ns,dstran,
     $     stress,yldp_ns,statev,eeq_ns,deq
      real*8 empa,gpa
      empa = 1d6
      gpa  = 1d9
      deq=0d0
      if (idiaw) then
         call w_chr(idia,'* stran at step n')
         call w_dim(idia,stran,ntens,1d0/empa,.true.)
         call w_chr(idia,'* stress at step n')
         call w_dim(idia,stress,ntens,1d0/empa,.true.)
      endif

c$$$  1. Save jacobian as elastic moduli
      ddsdde(:,:) = Cel(:,:)
c$$$  2. Update strain.
      stran_el_ns(1,:) = stran_el_ns(0,:) + dstran(:)
c     call add_array(stran,dstran,ntens)
c$$$  3. Updates stress
      call mult_array(ddsdde,stran_el_ns(1,:),ntens,stress)

      if (idiaw) then
         call w_chr( idia,'* stress at n+1')
         call w_dim( idia,stress,ntens,1d0/empa,.true.)
c     call w_chr( idia,'* stran n+1')
c     call w_dim( idia,stran,ntens,1d0/empa,.true.)
         call w_chr( idia,'* ddsdde')
         call w_mdim(idia,ddsdde,  ntens,1d0/gpa)
         call w_chr( idia,'* stran_el_ns(n+1)')
         call w_dim( idia,stran_el_ns(1,:),ntens,1d0,.true.)
      endif

      if (idiaw) call w_chr(idia,'5')
c$$$  4. Update all other state varaiables

      eeq_ns(1)        = eeq_ns(0) + deq
      stran_pl_ns(1,:) = stran_pl_ns(0,:)
      call update_yldp(iyld_law,yldp_ns,nyldp,deq)

c$$   5. Store updated state variables to statev
      call restore_statev(statev,nstatv,eeq_ns(1),stran_el_ns(1,:),
     $     stran_pl_ns(1,:),ntens,yldp_ns(1,:),nyldp,1,.false.,idia)

      if (idiaw) call w_chr(idia,'6')
      return
      end subroutine update_elastic
c-----------------------------------------------------------------------
      subroutine restore_statev(statev,nstatv,eqpl,stran_el,stran_pl,
     $     ntens,yldp,nyldp,iopt,verbose,iunit)
c-----------------------------------------------------------------------
c     Arguments
c     statev  : state variable array
c     nstatv  : len of statev
c     eqpl    : equivalent plastic strain
c     stran_el: cumulative elastic strain
c     stran_pl: cumulative plastic strain
c     ntens   : len of stran_el and stran_pl
c     iopt    : option to define the behavior of restore_statev
c                  0: read from statev
c                  1: save to statev
c     verbose
c     iunit   : file unit (if not zero) or std to which
c               inquiry stream will be written
      implicit none
      integer nstatv,ntens
      dimension statev(nstatv),stran_el(ntens),stran_pl(ntens),
     $     yldp(nyldp)
      real*8 statev,eqpl,stran_el,stran_pl,yldp
      integer iopt,i,nyldp,iunit
      logical verbose


      if (iopt.eq.0) then
         ! read from statev
         eqpl = statev(1)
         do i=1,ntens
            stran_el(i) = statev(i+1)
            stran_pl(i) = statev(i+1+ntens)
         enddo
         do i=1,nyldp
            yldp(i)     = statev(i+1+ntens*2)
         enddo
      elseif (iopt.eq.1) then
         ! save to statev
         ! read from statev
         statev(1) = eqpl
         do i=1,ntens
            statev(i+1)       = stran_el(i)
            statev(i+1+ntens) = stran_pl(i)
         enddo
         do i=1,nyldp
            statev(i+1+ntens*2) = yldp(i)
         enddo
      else
         write(*,*) 'Unexpected iopt given'
         stop
      endif

      if (verbose) then
         call w_empty_lines(iunit,2)
         call fill_line(iunit,'*',72)
         call w_chr(iunit,'inquiry request on state variables')
         call w_val(iunit,' eqpl   ',eqpl)
         call w_chr(iunit,'stran_el ')
         call w_dim(iunit,stran_el,ntens,1.d0,.true.)
         call w_chr(iunit,'stran_pl ')
         call w_dim(iunit,stran_pl,ntens,1.d0,.true.)
         call w_chr(iunit,'yldp')
         call w_dim(iunit,yldp,nyldp,1.d0,.true.)
         call fill_line(iunit,'*',72)
         call w_empty_lines(iunit,2)
      endif

      end subroutine
c-----------------------------------------------------------------------
c     iso elastic
      include "/home/younguj/repo/abaqusPy/umats/lib/elast.f"
c     UHARD using Voce
      include "/home/younguj/repo/abaqusPy/umats/lib/uhard.f"
c     lib.f
      include "/home/younguj/repo/abaqusPy/umats/lib/lib.f"
c     lib_w.f
      include "/home/younguj/repo/abaqusPy/umats/lib/lib_write.f"
c     diag.f
      include "/home/younguj/repo/abaqusPy/umats/lib/diag.f"
c     return_mapping.f
      include "/home/younguj/repo/abaqusPy/umats/lib/return_mapping.f"
c     is.f - testing subroutines/functions (e.g., is_inf)
      include "/home/younguj/repo/abaqusPy/umats/lib/is.f"
c     yld.f - yield function associated calculations
      include "/home/younguj/repo/abaqusPy/umats/lib/yld.f"
