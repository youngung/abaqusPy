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

c     updates stress
      do k1=1,ntens
         do k2=1, ntens
            stress(k2)=stress(k2)+ddsdde(k2,k1)*dstran(k1)
         enddo
      enddo

      return
      end subroutine umat
      include "../lib/elast.f"
