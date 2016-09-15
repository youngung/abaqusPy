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
      real*8 dtime,temp,dtemp,celent

c     Beginning of user-specific part
      dimension eelas(6), eplas(6), flow(6)
      parameter(one=1d0,two=2d0,three=3d0,six=6d0,enumax=0.4999d0,
     $     newton=10,toler=1d-6)

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
      emod=props(1)
      enu = props(2)

      call emod_iso(emod,enu,ddsdde) !! isotropic elastic modulus

      if (kinc.eq.1 .and. noel.eq.1) then
         write(imsg,*) 'Elastic modulus'
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

      subroutine emod_iso(e,nu,c)
c     intent(in): e,nu
c     intent(out): c
      parameter(ndi=3,ntens=6)
      real*8 e, nu, c(ntens,ntens), x
      integer i,j

c     initialization
      do i=1,ntens
         do j=1,ntens
            c(i,j) = 0.d0
         enddo
      enddo
c
c     construct elastic tensor (6x6) with assuming that
c     \gamma_ij = 2\varepsilon_ij is the engineering shear strain
c

c     Multiplier
      x = e/(1.+nu)/(1.-2.*nu)

c     off-diagonal terms
      do i=1,3
         do j=1,3
            c(i,i) = nu * x
         enddo
      enddo

      do i=1,3
         c(i,i) = (1.-nu)*x     !! overwrite the diganogal term
         c(i+3,i+3) = (1.-2.*nu)/2. * x
      enddo

      return
      end subroutine emod_iso


c      include "/home/younguj/repo/abaqusPy/umats/lib/elast.f"
