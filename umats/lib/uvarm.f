c-----------------------------------------------------------------------
      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)
c      INCLUDE 'ABA_PARAM.INC'
      implicit none
C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      integer nuvarm
      DIMENSION UVAR(nuvarm),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

C     The dimensions of the variables FLGRAY, ARRAY and JARRAY
C     must be set equal to or greater than 15.

c     given arrays
      real*8 uvar,direct,t,time,array,coord,dtime
      integer matlayo,laccfla,jarray,jmac,jmatyp,jerror
      integer noel,npt,layer,kspt,kstep,kinc,ndi,nshr
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
c-----------------------------------------------------------------------
c     local variables
      dimension aux33(3,3),bux33(3,3),aux6(6),directt(3,3),aux4(4),
     $     aux3(3),bux3(3)
      real*8 aux33,aux6,directt,bux33,aux4,empa,aux3,bux3
      integer jrcd,idia
      logical idiaw

      idiaw=.false.
c      idiaw=.true.
c$$$      if (kinc.eq.30 .and. npt.eq.1 .and. noel.eq.1 .and.dtime.ne.0) then
c$$$         idiaw=.true.
c$$$      endif
      idia=0
      empa = 1d6

c     user coding to define UVAR

c-----------------------------------------------------------------------
c     direct converts global directions to that of material
c     directt = transpose(direct) so that it converts material
c     direction to global.
c      call tens33_trans(direct,directt)
      directt(:,:)=direct(:,:)
c-----------------------------------------------------------------------
      if (idiaw) then
         call fill_line(idia,'*',52)
         call w_chr(idia,'direct')
         call w_mdim(idia,direct,3,1d0)
         call w_chr(idia,'directt')
         call w_mdim(idia,directt,3,1d0)
      endif
c-----------------------------------------------------------------------
c     Getting solution dependent state variables
c-----------------------------------------------------------------------
      call getvrm('SDV',array,jarray,flgray,jrcd,jmac,jmatyp,
     $     matlayo,laccfla)
      idia=0
      jerror = jerror + jrcd
      if (idiaw) then
         call w_chr(idia,'array(1:3); elastic strain')
         call w_dim(idia,array(1:3),3,1d0,.true.)
         call w_chr(idia,'array(4:6): plastic strain')
         call w_dim(idia,array(1:3),3,1d0,.true.)
      endif
      aux3(:) = array(1:3)
      call rot_tensor_shell2(aux3,directt,bux3,1) ! elastic strain
      uvar(1:3)=bux3(:)
      aux3(:) = array(4:6)
      call rot_tensor_shell2(aux3,directt,bux3,1) ! plastic strain
      uvar(4:6)=bux3(:)
      if (idiaw) then
         call w_chr(idia,'elastic strain in global')
         call w_dim(idia,uvar(1:3),3,1d0,.true.)
         call w_chr(idia,'plastic strain in global')
         call w_dim(idia,uvar(4:6),3,1d0,.true.)
      endif
c-----------------------------------------------------------------------
c     Getting stress tensor transformed to global coordinates
c-----------------------------------------------------------------------
      call getvrm('S',array,jarray,flgray,jrcd,jmac,jmatyp,
     $     matlayo,laccfla)
      if (idiaw) then
         call w_chr(idia,'Stress tensor in local coordinate [MPa]')
         call w_dim(idia,array(1:4),4,1d0/empa,.true.)
      endif
      jerror = jerror + jrcd
      call rot_tensor_shell(array(1:4),directt,aux4,0)
      uvar(7)=aux4(1)
      uvar(8)=aux4(2)
      uvar(9)=aux4(4)
c-----------------------------------------------------------------------
c     Getting strain tensor transformed to global coordinates
c-----------------------------------------------------------------------
      call getvrm('E',array,jarray,flgray,jrcd,jmac,jmatyp,
     $     matlayo,laccfla)
      jerror = jerror + jrcd
      call rot_tensor_shell(array(1:4),directt,aux4,1)
      uvar(10)=aux4(1)
      uvar(11)=aux4(2)
      uvar(12)=aux4(4)
c      if (idiaw) then
      if (idiaw) then
         call w_chr(idia,'Elastic strain in global coordinates')
         call w_dim(idia,uvar(1:3),3,1d0,.true.)
         call w_chr(idia,'Plastic strain in global coordinates')
         call w_dim(idia,uvar(4:6),3,1d0,.true.)
         call w_chr(idia,'stress in global [MPa]')
         call w_dim(idia,uvar(7:9),3,1d0/empa,.true.)
         call fill_line(idia,'*',52)
         call w_empty_lines(idia,2)
      endif
c-----------------------------------------------------------------------
      if(jerror.ne.0)then
        write(6,*) 'request error in uvarm for element number ',
     1      noel,'integration point number ',npt
      endif

      call stop_debug(0)
      
      return
      end subroutine uvarm
c-----------------------------------------------------------------------
      subroutine rot_tensor_shell(a4,rot,b4,iopt)
c     Convert 4 dimensional plane-stress tensor
c     - convention:
c      vec  ij component
c        1:      11
c        2:      22
c        3:      33
c        4:      12
c     Arguments
c     a4  : a4 tensor
c     rot : rotation matrix (3,3)
c     b4  : b4 (rotated tensor)
c     iopt: 0; stress, 1: strain
      implicit none
      dimension a4(4),b4(4),rot(3,3),aux6(6),aux33(3,3),bux33(3,3)
      real*8 a4,b4,rot,aux6,aux33,bux33
      integer iopt
c      call w_chr(0,'a4:')
c      call w_dim(0,a4,4,1d0,.true.)
      call reduce_4to6(a4,aux6)
c      call w_chr(0,'aux6:')
c      call w_dim(0,aux6,6,1d0,.true.)
      if (iopt.eq.0) then
         call voigt2(aux6,aux33)
      elseif (iopt.eq.1) then
         call voigt4(aux6,aux33)
      endif
c      call w_chr(0,'aux33:')
c      call w_mdim(0,aux33,3,1d0)
      bux33(:,:) = 0d0
      call rot_tensor(aux33,rot,bux33)
c      call w_chr(0,'bux33:')
c      call w_mdim(0,bux33,3,1d0)
      if (iopt.eq.0) then
         call voigt1(bux33,aux6)
      elseif(iopt.eq.1) then
         call voigt3(bux33,aux6)
      endif
c      call w_chr(0,'aux6:')
c      call w_dim(0,aux6,6,1d0,.true.)
      call reduce_6to4(aux6,b4)
c      call w_chr(0,'b4:')
c      call w_dim(0,b4,4,1d0,.true.)
      return
      end subroutine rot_tensor_shell
c-----------------------------------------------------------------------
      subroutine rot_tensor_shell2(a3,rot,b3,iopt)
c     Convert 3 dimensional plane-stress tensor
c     - convention:
c      vec  ij component
c        1:      11
c        2:      22
c        3:      12
c     Arguments
c     a3  : a3 tensor
c     rot : rotation matrix (3,3)
c     b3  : b3 (rotated tensor)
c     iopt: 0; stress, 1: strain
      implicit none
      dimension a3(3),b3(3),rot(3,3),aux6(6),aux33(3,3),bux33(3,3)
      real*8 a3,b3,rot,aux6,aux33,bux33
      integer iopt
c      call w_chr(0,'--- shell2 ---')
c      call w_chr(0,'a3')
c      call w_dim(0,a3,3,1d0,.true.)
      call reduce_3to6(a3,aux6)
c      call w_chr(0,'aux6')
c      call w_dim(0,aux6,6,1d0,.true.)
      if (iopt.eq.0) then
         call voigt2(aux6,aux33)
      elseif (iopt.eq.1) then
         call voigt4(aux6,aux33)
      endif
c      call w_chr(0,'aux33')
c      call w_mdim(0,aux33,3,1d0)
      bux33(:,:)=0d0
      call rot_tensor(aux33,rot,bux33)
c      call w_chr(0,'bux33')
c      call w_mdim(0,bux33,3,1d0)
c      call w_chr(0,'** bux33 **')
      if (iopt.eq.0) then
         call voigt1(bux33,aux6)
      elseif(iopt.eq.1) then
         call voigt3(bux33,aux6)
      endif
c      call w_chr(0,'aux6')
c      call w_dim(0,aux6,6,1d0,.true.)
      call reduce_6to3(aux6,b3)
c      call w_chr(0,'b3')
c      call w_dim(0,b3,6,1d0,.true.)
c      call w_chr(0,'--- end of shell 2 ---')
      return
      end subroutine rot_tensor_shell2
