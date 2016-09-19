c-----------------------------------------------------------------------
c     calculate transformation matrix
      subroutine inplane_rot_matrix(psi,r)
      implicit none
      real*8 r(3,3), psi, c, s
      c=dcos(psi)
      s=dsin(psi)
      r(:,:) = 0.d0
      r(1,1) =  c
      r(1,2) = -s
      r(2,1) =  s
      r(2,2) =  c
      r(3,3) =  1.
      return
      end subroutine inplane_rot_matrix
c-----------------------------------------------------------------------
      subroutine inplane_rot(psi,a33,b33)
      implicit none
      real*8 psi,a33(3,3),b33(3,3),rot(3,3)
      integer i,j,k,l
      call inplane_rot_matrix(psi,rot)
      b33(:,:)=0.d0
      do 10 i=1,3
      do 10 j=1,3
      do 10 k=1,3
      do 10 l=1,3
         b33(i,j) = b33(i,j) + rot(i,k) * a33(k,l) * rot(j,l)
 10   continue
      return
      end subroutine
