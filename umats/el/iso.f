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

c$$$ Assuming real*8
      real*8 dtime,temp,dtemp,celent,ddsdde

c     Beginning of user-specific part
      dimension eelas(6), eplas(6), flow(6)
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,six=6.d0,
     $     enumax=0.4999d0,newton=10,toler=1d-6)

      integer imsg, k1, k2
c$$$
c$$$c     PROPS(1) - E   Young's modulus
c$$$c     PROPS(2) - Nu  Poisson's ratio

      imsg=7                    !! msg file ID
      if (ndi.ne.3) then
         write(imsg,*)'This UMAT may only be used for elements with
     $three direct (normal) stress components'
         call xit               ! exit
      endif
      emod =props(1)
      enu = props(2)
c$$$      emod=300d9
c$$$      enu=0.3
c     initialization
      do i=1,ntens
         do j=1,ntens
            ddsdde(i,j) = zero
         enddo
      enddo
c
c     construct elastic tensor (6x6) with assuming that
c     \gamma_ij = 2\varepsilon_ij is the engineering shear strain
c

c     Multiplier
      x = emod/(one+enu)/(one-two*enu)

c     off-diagonal terms
      do i=1,3
         do j=1,3
            ddsdde(i,j) = enu * x
         enddo
      enddo

      do i=1,3
c        overwrite the diganogal term
         ddsdde(i,i) = (one-enu)*x
c        digonal terms
         ddsdde(i+3,i+3) = (one-two*enu)/two * x
      enddo

!     call emod_iso(emod,enu,ddsdde) ! isotropic elastic modulus
      if (kinc.eq.1 .and. noel.eq.1) then
         do i=1,ntens
            write(imsg,'(6e13.2)') (ddsdde(i,j),j=1,ntens)
         enddo
         write(msg,'(a)')'----------------------------------------------------
     $--------------------'
      endif

c     updates stress
      do k1=1,ntens
         do k2=1, ntens
            stress(k2)=stress(k2)+ddsdde(k2,k1)*dstran(k1)
         enddo
      enddo

      return
      end subroutine umat
c$$$      
c$$$
c$$$      subroutine emod_iso(e,nu,c)
c$$$c     intent(in): e,nu
c$$$c     intent(out): c
c$$$      parameter(ndi=3,ntens=6)
c$$$      real*8 e, nu, c(ntens,ntens), x
c$$$      integer i,j
c$$$
c$$$c     initialization
c$$$      do i=1,ntens
c$$$         do j=1,ntens
c$$$            c(i,j) = 0.d0
c$$$         enddo
c$$$      enddo
c$$$c
c$$$c     construct elastic tensor (6x6) with assuming that
c$$$c     \gamma_ij = 2\varepsilon_ij is the engineering shear strain
c$$$c
c$$$
c$$$c     Multiplier
c$$$      x = e/(1.+nu)/(1.-2.*nu)
c$$$
c$$$c     off-diagonal terms
c$$$      do i=1,3
c$$$         do j=1,3
c$$$            c(i,i) = nu * x
c$$$         enddo
c$$$      enddo
c$$$
c$$$      do i=1,3
c$$$         c(i,i) = (1.-nu)*x     !! overwrite the diganogal term
c$$$         c(i+3,i+3) = (1.-2.*nu)/2. * x
c$$$      enddo
c$$$
c$$$      return
c$$$      end subroutine emod_iso
c$$$
c$$$c      include "/home/younguj/repo/abaqusPy/umats/lib/elast.f"
