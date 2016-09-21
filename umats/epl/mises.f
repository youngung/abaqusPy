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


c$$$  local arrays
      character*255 fndia,fnstr
      dimension Cel(ntens,ntens),dphi_n(ntens),d2phi_n(ntens,ntens),
     $     stran_el(ntens),stran_pl(ntens),spr(ntens),epr(ntens),
     $     voce_params(4)
      real*8 e,enu,G,eK,Cel
c     predictor stress
      real*8 spr,epr
      real*8 voce_params
      real*8 eeq_n,yld_h,dh
!     eeq_n: eq plastic strain at step n
!     yld_h: flow stress at hardening curve
!     dh: dH/dE^eq

      real*8 phi_n,dphi_n,d2phi_n
      real*8 toler_yield,stran_el,stran_pl
      integer imsg,idia,i,istr
      logical idiaw
      real*8 empa,gpa
      parameter(toler_yield=1d-6,empa=1d6,gpa=1d9)

      voce_params(1) = 479.0d0
      voce_params(2) = 339.7d0
      voce_params(3) = 7.784d0
      voce_params(4) = 70.87d0
      e = 210d0*gpa
      enu=0.3d0

      imsg=7
      idia=315
      istr=425
      if (kspt.eq.1 .and. noel.eq.1 .and. npt.eq.1) then
         idiaw=.true.
         fndia='/home/younguj/repo/abaqusPy/examples/one/diagnose.txt'
         fnstr='/home/younguj/repo/abaqusPy/examples/one/strstr.txt'
         open(idia,position='append',file=fndia)
         open(istr,position='append',file=fnstr)
      endif
      write(*,*)'** File opened **'

c     print head
      call print_head(0)
      call print_head(imsg)

      call restore_statev(statev,nstatv,eeq_n,stran_el,stran_pl,
     $     ntens,0)
      write(*,*)'** state variable stored **'
      if (idiaw) then
         write(*,*)'* stran'
c         call w_dim(idia,stran,ntens,1d0,.true.)
         call w_dim(0,stran,ntens,1d0,.true.)
         write(*,*)'* stran_el'
c         call w_dim(idia,stran_el,ntens,1d0,.true.)
         call w_dim(0,stran_el,ntens,1d0,.true.)
         write(*,*)'* stran_pl'
c         call w_dim(idia,stran_pl,ntens,1d0,.true.)
         call w_dim(0,stran_pl,ntens,1d0,.true.)
         write(*,*)'* dstran'
c         call w_dim(idia,dstran,ntens,1d0,.true.)
         call w_dim(0,dstran,ntens,1d0,.true.)
      endif

c     Moduluar pseudo code for stress integration

c     i.   Check the current (given) variables
      call emod_iso_shell(e,enu,G,eK,Cel)
c$$$      write(*,*)    '** EMOD_ISO_SHELL **'
c$$$      write(imsg,*) '** EMOD_ISO_SHELL **'
      write(*,*)'* C^el [GPa]'
c      call w_mdim(idia,Cel,ntens,1d0/gpa)
      call w_mdim(   0,Cel,ntens,1d0/gpa)
c     ii.  Trial stress calculation
c     Calculate predict elastic strain
      call add_array2(stran_el,dstran,epr,ntens)
      call mult_array(Cel,epr,ntens,spr)
c$$$      write(*,*)    '** MULT_ARRAY **'
c$$$      write(imsg,*) '** MULT_ARRAY **'
      write(*,*)'* trial stress [MPa]'
c      call w_dim(idia,spr,ntens,1d0/empa,.true.)
      call w_dim(0,spr,ntens,1d0/empa,.true.)

c     iii. See if the pr stress (spr) calculation is in the plastic or elastic regime
      call voce(eeq_n, voce_params(1),voce_params(2),voce_params(3),
     $     voce_params(1),yld_h,dh)
      yld_h = yld_h * empa
      dh    = dh    * empa
c$$$      write(*,*)    '** VOCE **'
c$$$      write(imsg,*) '** VOCE **'
      call vm_shell(spr,phi_n,dphi_n,d2phi_n)
c$$$      write(*,*)    '** VM_SHELL **'
c$$$      write(imsg,*) '** VM_SHELL **'

      write(*,*)'* phi_n [MPa]', phi_n/empa
      write(*,*)'* yld_h [MPa]', yld_h/empa

      if (phi_n.lt.yld_h) then
c$$$         write(*,*)    'ELASTIC'
c$$$         write(imsg,*) 'ELASTIC'
c              1. Save jacobian as elastic moduli
         ddsdde(:,:) = Cel(:,:)
c              2. Update strain.
         call add_array(stran,dstran,ntens)
c$$$         write(*,*)    'ADARRAY'
c$$$         write(imsg,*) 'ADARRAY'
c              3. Updates
c                 3.1 update stress
         call mult_array(ddsdde,stran,ntens,stress)
c$$$         write(*,*)    'MULT_array'
c$$$         write(imsg,*) 'MULT_array'
c                 3.2 update strain (incremental strain belongs to elastic)
         do i=1,ntens
            stran_el(i) = stran_el(i) + dstran(i)
         enddo
c              4. Update all other state varaiables
         call restore_statev(statev,nstatv,eeq_n,stran_el,stran_pl,
     $        ntens,1)
c$$$         write(*,*)
c$$$         write(*,'(a7,i3)',   advance='no')'kinc:',kinc
c$$$         write(imsg,*)
c$$$         write(imsg,'(a7,i3)',advance='no')'kinc:',kinc
c$$$         call w_dim(imsg,stran,ntens,1.d0,.true.)
c$$$         call w_dim(imsg,stran_el,ntens,1.d0,.true.)
c$$$         call w_dim(imsg,stran_pl,ntens,1.d0,.true.)
c$$$         call w_dim(imsg,stress,ntens,1.d0,.true.)
c$$$         call print_foot(0)
c$$$         call print_foot(imsg)
         return
      else
         write(*,*)   'PLASTIC'
         write(imsg,*)'PLASTIC'

         call print_foot(0)
         call print_foot(imsg)
c     vi. Return mapping
         call return_mapping(Cel,spr,phi_n,eeq_n,dphi_n,voce_params,
     $        dstran,stran,stran_el,stran_pl,ntens,idiaw)
         stop -1
         write(imsg,*)'return-mapping'
c     v. Exit from iv. means
c       s_(n+1) is obtained.
c       dlamb   is obtained
c       de^(el)_(n+1) is obtained

c       update plastic strain = depl_ij = depl_ij + n_ij dlamb

c     vi. Caculate jacobian (ddsdde)
c       C^el - [ C^el:m_(n+1) cross C^el:m_(n+1) ]   /  [ m_(n+1):C^el:m_(n+1) + h(depl_ij_(n+1)) ]
      endif
      close(idia)
      close(istr)
      call print_foot(0)
      call print_foot(imsg)
      return
      end subroutine umat
c-----------------------------------------------------------------------
      subroutine restore_statev(statev,nstatv,eqpl,stran_el,stran_pl,
     $     ntens,iopt)
c     statev  : state variable array
c     nstatv  : len of statev
c     eqpl    : equivalent plastic strain
c     stran_el: cumulative elastic strain
c     stran_pl: cumulative plastic strain
c     ntens   : len of stran_el and stran_pl
c     iopt    : option to define the behavior of restore_statev
c                  0: read from statev
c                  1: save to statev
      implicit none
      integer nstatv,ntens
      dimension statev(nstatv),stran_el(ntens),stran_pl(ntens)
      real*8 statev,eqpl,stran_el,stran_pl
      integer iopt, i
      if (iopt.eq.0) then
         ! read from statev
         eqpl = statev(1)
         do i=1,ntens
            stran_el(i) = statev(i+1)
            stran_pl(i) = statev(i+1+ntens)
         enddo
      elseif (iopt.eq.1) then
         ! save to statev
         ! read from statev
         statev(1) = eqpl
         do i=1,ntens
            statev(i+1)       = stran_el(i)
            statev(i+1+ntens) = stran_pl(i)
         enddo
      else
         write(*,*) 'Unexpected iopt given'
         stop
      endif
      end subroutine
c-----------------------------------------------------------------------
c     print header
      subroutine print_head(i)
c     i: file unit (use std if 0, use imsg elsewhere)
      integer i
      call w_empty_lines(i,2)
      if (i.eq.0) then ! USE std
         write(*,*) '*------------------*'
         write(*,*) '|       UMAT       |'
         write(*,*) '*------------------*'
      else
         write(i,*) '*------------------*'
         write(i,*) '|       UMAT       |'
         write(i,*) '*------------------*'
      endif
      call w_empty_lines(i,1)
      return
      end subroutine print_head
c-----------------------------------------------------------------------
c     print footer
      subroutine print_foot(i)
c     i: file unit (use std if 0, use imsg elsewhere)
      integer i
      call w_empty_lines(i,2)
      if (i.eq.0) then ! USE std
         write(*,*) '*------------------*'
         write(*,*) '|       END        |'
         write(*,*) '*------------------*'
      else
         write(i,*) '*------------------*'
         write(i,*) '|       END        |'
         write(i,*) '*------------------*'
      endif
      call w_empty_lines(i,1)
      return
      end subroutine print_foot
c-----------------------------------------------------------------------
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
