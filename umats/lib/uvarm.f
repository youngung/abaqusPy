c-----------------------------------------------------------------------
      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)
c      INCLUDE 'ABA_PARAM.INC'
      implicit none
C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

C     The dimensions of the variables FLGRAY, ARRAY and JARRAY
C     must be set equal to or greater than 15.

c     given arrays
      real*8 uvar,direct,t,time,array,coord
      integer laccfla,jarray,jmac,jmatyp
c     Arguments
c     direct : material orientation in terms of global coordinates
c     t      : material orientation w.r.t. element basis
c     time   : time stamps at steps n and n+1
c     dtime  : time increment
c     cmname : UMAT name
c     coord  : coordinates at this material point
c     jmac   : variable that must be passed into the GETVRM utility
c             routine to access an output variable.
c     jmatyp : Variable that must be passed into the GETVRM utility
c             routine to access an output variable.
c
c     local variables
      dimension aux33(3,3),bux33(3,3),aux6(6),directt(3,3)
      integer jrcd,aux33,aux6,directt,bux33

c     user coding to define UVAR
c-----------------------------------------------------------------------
      call tens33_trans(direct,directt)

c-----------------------------------------------------------------------
c     Getting stress tensor transformed to global coordinates
c-----------------------------------------------------------------------
      call getvrm('S',array,jarray,flgray,jrcd,jmac,jmatyp,
     $     matlayo,laccfla)
      jerror = jerror + jrcd
      UVAR(1:6) = ARRAY(1:6)
      call voigt2(UVAR(1:6),aux33)
      call rot_tensor(aux33,directt,bux33)
      call voigt1(bux33,UVAR(1:6))

c-----------------------------------------------------------------------
c     Getting strain tensor transformed to global coordinates
c-----------------------------------------------------------------------
      call getvrm('E',array,jarray,flgray,jrcd,jmac,jmatyp,
     $     matlayo,laccfla)
      jerror = jerror + jrcd
      UVAR(7:12) = ARRAY(1:6)
      call voigt4(UVAR(7:12),aux33)
      call rot_tensor(aux33,directt,bux33)
      call voigt2(bux33,UVAR(7:12))

      if(jerror.ne.0)then
        write(6,*) 'request error in uvarm for element number ',
     1      noel,'integration point number ',npt
      endif

      return
      end subroutine uvarm

c$$$      include '/home/younguj/repo/abaqusPy/lib/cnv.f' ! voigt
c$$$      include '/home/younguj/repo/abaqusPy/lib/algb.f'! rot_tensor
