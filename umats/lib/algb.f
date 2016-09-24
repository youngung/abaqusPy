c-----------------------------------------------------------------------
c     calculate in-plane rotation matrix
      subroutine inplane_rot_matrix(psi,r)
c     intent(in) psi
c     intent(out) r
      implicit none
      real*8 r(3,3), psi, c, s
      c=dcos(psi)
      s=dsin(psi)
      r(:,:) = 0.d0
      r(1,1) =  c
      r(1,2) =  s*(-1.d0)
      r(2,1) =  s
      r(2,2) =  c
      r(3,3) =  1.d0
      return
      end subroutine inplane_rot_matrix
c-----------------------------------------------------------------------
c     Apply in-plane rotation on array b33
      subroutine inplane_rot(psi,a33,b33)
c     intent(in)  psi, a33
c     intent(out) b33
c     psi
c     a33
c     b33
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
c-----------------------------------------------------------------------
c     Apply tensor inner dot
c     ci = aij x bj
      subroutine mult_array(aij,bj,ntens,ci)
c     intent(int) aij,bj,ntens
c     intent(out) ci
      implicit none
      integer i,j,ntens
      real*8 aij(ntens,ntens),bj(ntens),ci(ntens)
      ci(:) = 0.d0
      do 10 i=1,ntens
      do 10 j=1,ntens
         ci(i) = ci(i) + aij(i,j) * bj(j)
 10   continue
      return
      end subroutine mult_array
c-----------------------------------------------------------------------
c     Apply tensor inner dot for 2nd x 2nd
c     cij = aik x bkj
      subroutine mult_array2(aik,bkj,ntens,cij)
c     intent(int) aik,bkj,ntens
c     intent(out) cij
      implicit none
      integer i,j,k,ntens
      real*8 aik(ntens,ntens),bkj(ntens,ntens),cij(ntens,ntens)
      cij(:,:) = 0.d0
      do 10 i=1,ntens
      do 10 j=1,ntens
      do 10 k=1,ntens
         cij(i,j) = cij(i,j) + aik(i,k) * bkj(k,j)
 10   continue
      return
      end subroutine mult_array2
c-----------------------------------------------------------------------

c     Apply incremental update on array
c     ai = ai + d_ai
      subroutine add_array(ai,d_ai,ntens)
c     intent(in) ai, d_ai, ntens
c     intent(out) ai
c     ai   : old array
c     d_ai : increments
c     ntens: len
      implicit none
      integer ntens, i
      real*8 ai(ntens),d_ai(ntens)
      do 10 i=1,ntens
         ai(i) = ai(i) + d_ai(i)
 10   continue
      return
      end subroutine add_array
c-----------------------------------------------------------------------
c     Apply incremental update on array
c     ci = ai + d_ai
      subroutine add_array2(ai,d_ai,ci,ntens)
c     intent(in) ai, d_ai, ntens
c     intent(out) ci
c     ai   : old array
c     d_ai : increments
c     ci   : new array
c     ntens: len
      implicit none
      integer ntens, i
      real*8 ai(ntens),d_ai(ntens),ci(ntens)
      do 10 i=1,ntens
         ci(i) = ai(i) + d_ai(i)
 10   continue
      return
      end subroutine add_array2
c-----------------------------------------------------------------------
      subroutine get_idx(n,a)
      dimension a(n,n)
      integer n,i
      real*8 a
      a(:,:) = 0d0
      do 10 i=1,n
         a(i,i) = 1d0
 10   continue
      end subroutine get_idx
