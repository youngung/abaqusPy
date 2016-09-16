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

      parameter(one=1d0,two=2d0,three=3d0,six=6d0,enumax=0.4999d0,
     $     newton=10,toler=1d-6)

c$$$  local arrays
c               eps^el        eps^pl        flow directions
      dimension eelas(ntens), eplas(ntens), flow(ntens),
     $     sdeviator(ntens), hard(3)

c$$$  Assuming real*8
      real*8 dtime,temp,dtemp,celent,smises,shydro,eg3,deqpl,f,fp,
     $     hard,mpa,gpa
      integer imsg,k,i,j,kewton
      mpa=1.e6
      gpa=1.e9

      imsg=7
      call w_empty_lines(imsg,3)
      write(imsg,'(a)') "Beginning of the UMAT"
      write(imsg,'(a,i5)') 'NOEL:',NOEL
      write(imsg,'(a,i5)') 'NPT:',NPT
      write(imsg,'(3a13)')'time1','time2', 'dtime'
      write(imsg,'(3e13.3)') (time(i),i=1,2), dtime
      write(imsg,*)
      write(imsg,'(a10,i5)')"kinc :",kinc
      write(imsg,'(a10,4i5)')"jstep:",jstep
      write(imsg,'(a10,i5)')"kspt :",kstp

c$$$  Elastic properties
c$$$      props(1) = 200.d9
c$$$      props(2) = 0.3
      emod = props(1)
      enu  = props(2)
      ! call emod_iso(emod,enu,ddsdde) ! isotropic elastic modulus
      call emod_iso_shell(emod,enu,ddsdde)
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
         write(msg,'(a)')'----------------------------------------------------
     $--------------------'
      endif

c$$$
c$$$  Recover elastic and plastic strains and rotate forward
c$$$
      write(*,*) statev
      call rotsig(statev(      1),drot,eelas,2,ndi,nshr)
      call rotsig(statev(ntens+1),drot,eplas,2,ndi,nshr)

      write(imsg,'(a)') "After rotsig"
      write(imsg,'(a)') 'eelas'
      call w_dim(imsg,eelas,ntens)
      write(imsg,'(a)') 'eplas'
      call w_dim(imsg,eplas,ntens)
      write(imsg,'(a)') 'dstrans'
      call w_dim(imsg,dstrans,ntens)

      eqplas=statev(1+2*ntens)  ! equivalent plastic strain
      write(imsg,'(a,e13.3)') 'Eqv pl strain:',eqplas
c$$$
c$$$  Caclulate predictor stress and elastic strain
c$$$
      do i=1,ntens
         do j=1,ntens
c           !multiply by total strain increment
            stress(j)=stress(j) + ddsdde(j,i)*dstran(i)
         enddo
         eelas(i) = eelas(i) + dstran(i)
      enddo
      write(imsg,'(a)') 'stress'
      !call w_dim(imsg,stress,ntens)
c$$$
c$$$  Calculate equivalent Von Mises stress
c$$$
      smises = vm(stress,ndi,ntens)
      !write(imsg,'(a,e13.3)') 'vm:', smises

      nvalue=nprops/2-1
      call UHARD(SYIEL0,HARD,EQPLAS,EQPLASRT,TIME,DTIME,TEMP,
     1     DTEMP,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,CMNAME,NSTATV,
     2     STATEV,NUMFIELDV,PREDEF,DPRED,NUMPROPS,PROPS)
      write(imsg,'(a)') "After uhard"
      write(imsg,'(a,f13.3)')'syiel0 [MPa]:',syiel0/mpa

c$$$
c$$$  Determine if actively yielding
c$$$

      if (smises.gt.(one+toler)*syiel0) then
c        obtain plastic flow direction
         call vm_devi_flow(stress,sdeviator,shydro,flow,ntens,ndi)


         write(imsg,'(a)') "After vm_devi_flow"

c        solve for equivalent von mises stress
c        and equivalent plastic strain increment using newton iteration

c        syiel0: the yield stress size using the predictor stress
         syield=syiel0
         deqpl=zero

         do kewton=1, newton     ! NR iteration
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

            write(imsg,'(a,i3)')'NR Step:',kewton
            write(imsg,'(a,e13.3)')'deqpl:',deqpl
            write(imsg,'(a,e13.3)') 'f :',f
            write(imsg,'(a,e13.3)') 'fp:',fp
            write(imsg,'(a,e13.3)')'eqplas+deqpl:',eqplas+deqpl
            if (abs(f).lt.toler*syiel0) goto 10
            call uhard(syield,hard,eqplas+deqpl,eqplasrt,time,dtime,temp,
     1           dtemp,noel,npt,layer,kspt,kstep,kinc,cmname,nstatv,
     2           statev,numfieldv,predef,dpred,numprops,props)
            write(imsg,'(a,e13.3)')'syield',syield

         enddo                  ! End of NR iteration

         write(imsg,'(a)') "*** Warning - Plasticity NR
     $did not converge"
 10      continue
         write(imsg,'(a,i)') " NR converged at kinc", kinc

c$$$
c$$$  Update stress, elastic and plastic strains and
c$$$  equivalent plastic strain
c$$$
         do i=1,ndi
            stress(i) = flow(i) * syield+shydro
            eplas(i)  = eplas(i) + three/two * flow(i) * deqpl
            eelas(i)  = eelas(i) - three/two * flow(i) * deqpl
         enddo

         do i=ndi+1,ntens
            stress(i) = flow(i) + syield
            eplas(i)  = eplas(i) + three * flow(i) * deqpl
            eelas(i)  = eelas(i) - three * flow(i) * deqpl
         enddo
         eqplas=eqplas+deqpl

c$$$
c$$$  calculate plastic dissipation
c$$$
         spd = deqpl*(syiel0+syield)/two

c$$$
c$$$  Formulate the jacobian (material tangent)
c$$$  First calculate effective moduli
c$$$
      else
         write(imsg,'(a)') '*********************'
         write(imsg,'(a)') '** Remains ELASTIC **'
         write(imsg,'(a)') '*********************'
      endif                     ! end of plasticity case

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
