c$$$  ABAQUS 
      program main
c      INCLUDE 'ABA_PARAM.INC'
      implicit none  !! to test the namespace in UMAT
C
      CHARACTER*80 CMNAME
      integer NDI,NSHR,NTENS,NSTATV,NPROPS,NOEL,NPT,LAYER,KSPT,JSTEP,
     $     KINC,kstep      
      parameter(ntens=6,ndi=3,nshr=3,nstatv=1,nprops=1,noel=1,npt=1,
     $     layer=1,kspt=1,kinc=1,kstep=1)

      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1     DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2     STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3     PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4     JSTEP(4)
c$$$  implicit none - variables passed through umat
c-----------------------------------------------------------------------

      real*8 STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,DRPLDE,DRPLDT,
     $     STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,PROPS,COORDS,
     $     DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,syield

c$$$  local arrays
      dimension eelas(ntens),eplas(ntens),flow(ntens),sdeviator(ntens),
     $     hard(3),dphi(ntens),d2phi(ntens,ntens)
      character*255 fndia,cwd
      real*8 zero,one,two,three,six,enumax,toler
      real*8 eqplas,syiel0,hard,sdeviator,empa,gpa,emod,enu,f,fp,eg3,
     $     deqpl,smises,shydro,flow,G,ekappa,emus,hs,elabs,eelas,
     $     eplas,eqplasrt,dphi,d2phi,de,emx
      integer imsg,idia,k,i,j,NUMFIELDV,numprops,kewton,newton
      logical isnan,idiaw,isnan_in_marr
      parameter(zero=0.d0,one=1d0,two=2d0,three=3d0,six=6d0,
     $     enumax=0.4999d0,newton=10,toler=1d-6)
      empa=1.e6
      gpa=1.e9

      de  = 0.01
      emx = 0.1
      do while(eqplas.le.emx)
         call uhard(SYIELD,HARD,EQPLAS,EQPLASRT,TIME,DTIME,TEMP,
     1        DTEMP,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,CMNAME,NSTATV,
     2        STATEV,NUMFIELDV,PREDEF,DPRED,NUMPROPS,PROPS)
         write(*,'(2(f7.3,2x))') eqplas, syield/empa
         eqplas = eqplas + de
      enddo
         
      end program main

c     UHARD using Voce
      include "/home/younguj/repo/abaqusPy/umats/lib/uhard.f"      
