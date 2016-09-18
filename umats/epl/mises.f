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
     $     DROT,PNEWDT,CELENT,DFGRD0,DFGRD1

c$$$  local arrays
      dimension eelas(ntens),eplas(ntens),flow(ntens),sdeviator(ntens),
     $     hard(3)
      real*8 zero,one,two,three,six,enumax,toler
      real*8 eqplas,syiel0,hard,sdeviator,mpa,gpa,emod,enu,f,fp,eg3,
     $     deqpl,smises,shydro,flow,G,kappa,mus,lams,hs,labs
      integer imsg,idia,k,i,j,NUMFIELDV,numprops,kewton,newton
      logical isnan,idiaw
      parameter(zero=0.d0,one=1d0,two=2d0,three=3d0,six=6d0,
     $     enumax=0.4999d0,newton=10,toler=1d-6)
      mpa=1.e6
      gpa=1.e9
      if (dtime.eq.0) then
         write(imsg,'(a)') 'DTIME .eq. 0'
         return
      endif

      imsg=7
      idia=315
      idiaw=.false.
      !! specify when to write the diagnose report
      if (kspt.eq.1 .and. noel.eq.1 .and. npt.eq.1) then
         idiaw=.true.
      endif

      if (idiaw)
     $     open(idia,position='append',
!     $     file='/home/younguj/repo/abaqusPy/umats/epl/diagnose.txt')
     $     file='/home/younguj/repo/abaqusPy/examples/one/diagnose.txt')

      call w_empty_lines(imsg,3)
      write(imsg,'(a)') "Beginning of the UMAT"
      write(imsg,'(a,i5)') 'NOEL:',NOEL
      write(imsg,'(a,i5)') 'NPT:',NPT
      write(imsg,'(3a13)')'time1','time2', 'dtime'
      write(imsg,'(3e13.3)') (time(i),i=1,2), dtime
      write(imsg,*)
      write(imsg,'(a10,i5)')"kinc :",kinc
      write(imsg,'(a10,4i5)')"jstep:",jstep
      write(imsg,'(a10,i5)')"kspt :",kspt

c$$$  Elastic properties
c$$$      props(1) = 200.d9
c$$$      props(2) = 0.3
      emod = props(1)
      enu  = props(2)
      ! call emod_iso(emod,enu,ddsdde) ! isotropic elastic modulus
      call emod_iso_shell(emod,enu,G,kappa,ddsdde)
      write(imsg,'(a,f7.1)')'G     [GPa]:',G/gpa
      write(imsg,'(a,f7.1)')'kappa [GPa]:',kappa/gpa
      eg3 = emod/(one+enu)/two*three

      write(imsg,'(a)') "After emod_iso"
      write(imsg,'(a,i)') 'ntens:', ntens
      write(imsg,'(a,i)') 'nshr:', nshr
      write(imsg,'(a4)') 'DROT'
      call w33f(imsg,drot)
      write(imsg,'(a)') '-----------------'

c     printout elastic modulus
      if (kinc.eq.1 .and. noel.eq.1) then
         write(imsg,'(a)')'DDSDDE'
         do i=1,ntens
            write(imsg,'(6e13.2)') (ddsdde(i,j),j=1,ntens)
         enddo
         write(imsg,'(a)')'----------------------------------------------------
     $--------------------'
      endif
c$$$
c$$$  Recover elastic and plastic strains and rotate forward
c$$$
      write(imsg,'(a)') 'dstrans'
      call w_dim(imsg,dstrans,ntens,1.,.true.)
      write(imsg,'(a)') "before rotsig"
      write(imsg,'(a)') 'eelas'
      call w_dim(imsg,eelas,ntens,1.,.true.)
      write(imsg,'(a)') 'eplas'
      call w_dim(imsg,eplas,ntens,1.,.true.)
      if (idiaw) then
         write(idia,'(i3)',advance='no')kinc
         call w_dim(idia,statev(1:ntens),ntens,1.,.false.)
         write(idia,'(x,a1,x)',advance='no') '|'
         call w_dim(idia,statev(ntens+1:2*ntens),ntens,1.,.false.)
         write(idia,'(f10.3)', advance='no') statev(1+ntens*2)
         write(idia,'(x,a1,x)',advance='no') '|'
      endif

      call rotsig(statev(      1),drot,eelas,2,ndi,nshr)
      call rotsig(statev(ntens+1),drot,eplas,2,ndi,nshr)
      eqplas=statev(1+2*ntens)  ! equivalent plastic strain
      write(imsg,'(a)') "After rotsig"
      write(imsg,'(a)') 'eelas'
      call w_dim(imsg,eelas,ntens,1.,.true.)
      write(imsg,'(a)') 'eplas'
      call w_dim(imsg,eplas,ntens,1.,.true.)

      if (idiaw) then
         !write(idia,'(i3)',advance='no')kinc
         call w_dim(idia,eelas,ntens,1.,.false.)
         write(idia,'(x,a1,x)',advance='no') '|'
         call w_dim(idia,eplas,ntens,1.,.false.)
         write(idia,'(e10.3)',advance='no') eqplas
         write(idia,'(x,a1,x)',advance='no') '|'
      endif

c$$$
c$$$  Caclulate predictor stress and elastic strain
c$$$
      do 10 i=1,ntens
         do 5 j=1,ntens
c     !multiply by total strain increment
            stress(j)=stress(j) + ddsdde(j,i)*dstran(i)
 5       continue
         eelas(i) = eelas(i) + dstran(i)
 10   continue
      write(imsg,'(a)') 'stress [MPa]'
      call w_dim(imsg,stress,ntens,1./mpa,.true.)
c$$$
c$$$  Calculate equivalent Von Mises stress
c$$$
      call vm(stress,smises)

      call UHARD(SYIEL0,HARD,EQPLAS,EQPLASRT,TIME,DTIME,TEMP,
     1     DTEMP,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,CMNAME,NSTATV,
     2     STATEV,NUMFIELDV,PREDEF,DPRED,NUMPROPS,PROPS)
      write(imsg,'(a)') "After uhard"
      write(imsg,'(a,f13.3)')'syiel0 [MPa]:',syiel0/mpa

c$$$
c$$$  Determine if actively yielding
c$$$

      if (smises.gt.(one+toler)*syiel0) then
         write(imsg,'(a)') '*********************'
         write(imsg,'(a)') '**     PLASTIC     **'
         write(imsg,'(a)') '*********************'
         if (idiaw) write(idia,'(x,a1,x)',advance='no') 'P'
c        obtain plastic flow direction
         call vm_devi_flow(stress,sdeviator,shydro,flow,ntens,ndi)
         write(imsg,'(a)') "After vm_devi_flow"

c        solve for equivalent von mises stress
c        and equivalent plastic strain increment using newton iteration

c        syiel0: the yield stress size using the predictor stress
         syield=syiel0
         deqpl=zero

         do 50 kewton=1, newton ! NR iteration
c$$$
c$$$  Newton Raphson Method to find deqpl
c$$$
c     \partial(rhs)/\partial(deqpl) = -eg3 - \partial(syield)\partial(deqpl)
c                                   = -eg3 - hard(1)
c     x_1 = x_0 - f/f`
c
            f     = smises-eg3*deqpl-syield
            fp    = -eg3-hard(1)
            deqpl = deqpl - f / fp

            write(imsg,*)'Step:',kewton
            write(imsg,*)'deqpl:',deqpl
            write(imsg,*)'f  [MPa]:',f/mpa
            write(imsg,*)'fp [MPa]:',fp/mpg
            write(imsg,*)'eqplas+deqpl:',eqplas+deqpl
            if (abs(f).lt.toler*syiel0) goto 100
            call uhard(syield,hard,eqplas+deqpl,eqplasrt,time,dtime,
     $           temp,dtemp,noel,npt,layer,kspt,kstep,kinc,cmname,
     $           nstatv,statev,numfieldv,predef,dpred,numprops,props)
            write(imsg,'(a,f7.1)')'syield [MPa]',syield/mpa

 50      continue               ! End of NR iteration

         write(imsg,'(a)') "*** Warning - Plasticity NR
     $did not converge"
         stop -1
 100     continue

         write(imsg,'(a,i3)') " NR converged at kinc", kinc

c$$$
c$$$  Update stress, elastic and plastic strains and
c$$$  equivalent plastic strain
c$$$
         do 105 i=1,ndi
            stress(i) = flow(i)  * syield+shydro
            eplas(i)  = eplas(i) + three/two * flow(i) * deqpl
            eelas(i)  = eelas(i) - three/two * flow(i) * deqpl
 105     continue

         do 110 i=ndi+1,ntens
            stress(i) = flow(i) * syield
            eplas(i)  = eplas(i) + three * flow(i) * deqpl
            eelas(i)  = eelas(i) - three * flow(i) * deqpl
 110     continue
         eqplas=eqplas+deqpl
c$$$
c$$$  calculate plastic dissipation
c$$$
         spd = deqpl*(syiel0+syield)/two

c$$$
c$$$  In case elasto-plasticity, the jacobian
c$$$  is different from that in elasticity
c$$$
c$$$  Formulate the jacobian (material tangent) in epl domain
c$$$  First calculate effective moduli
c$$$
         mus = G * syield / smises
         write(imsg,'(a,f7.1)') 'syield    [MPa]',syield/mpa
         write(imsg,'(a,f7.1)') 'predictor [MPa]',smises/mpa
         write(imsg,*)          'mus        ',mus
         labs = kappa - 2. / 3. * mus
         ddsdde(:,:) = 0.d0
         do 120 i=1,ndi
         do 120 j=1,ndi
            ddsdde(i,j) = labs
 120     continue

         do 140 i=1,ndi
            do 130 j=1,ndi
               ddsdde(i,j) = labs
 130        continue
            ddsdde(i,i) = labs + 2*mus
 140     continue
         do i=ndi+1,ntens
            ddsdde(i,i) = mus
         enddo
         hs = hard(1)/ (1.+hard(1)/3./G - 3. * mus)
         write(imsg,*) 'hs:',hs
         do 160 i=1,ntens
         do 160 j=1,ntens
            ddsdde(i,j) = ddsdde(i,j) + flow(i) * flow(j) * hs
 160     continue
         if (isnan_in_marr(ddsdde,ntens,ntens)) then
            write(imsg,'(a)') 'NAN found in DDSDDE'
            stop -1
         endif
         write(imsg,'(a)')'-----------'
         write(imsg,'(a)')'DDSDDE'

         call w_mdim(imsg,ddsdde,ntens)
         call w_empty_lines(imsg,5)
      else
         write(imsg,'(a)') '*********************'
         write(imsg,'(a)') '** Remains ELASTIC **'
         write(imsg,'(a)') '*********************'
         if (idiaw) write(idia,'(x,a1,x)',advance='no') 'E'
      endif                     ! end of plasticity case
c$$$
c$$$  store elastic and (equivalent) plastic strains
c$$$  in state variable array
c$$$
      do 180 i=1, ntens
         statev(i) = eelas(i)
         statev(i+ntens) = eplas(i)
 180  continue
      statev(1+2*ntens) = eqplas

      if (idiaw) then
         write(idia,'(2e13.4,3f11.4)',advance='no')
     $        eqplas, deqplas,
     $        stress(1)/mpa,
     $        stress(2)/mpa, stress(3)/mpa
         write(idia,'(x,a1,x)',advance='no') '|'
         write(idia,'(x,a5,x)',advance='no')'Flow:'
         call w_dim(idia,flow,ntens,1.,.true.)
         close(idia)
      endif


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
