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
c     psi
c     a33
c     b33
      implicit none
      real*8 psi,a33(3,3),b33(3,3),rot(3,3)
cf2py intent(in)  psi, a33
cf2py intent(out) b33
      call inplane_rot_matrix(psi,rot)
      call rot_tensor(a33,rot,b33)
      return
      end subroutine
c-----------------------------------------------------------------------
      subroutine rot_tensor(a33,rot,b33)
      implicit none
      dimension a33(3,3),rot(3,3),b33(3,3)
      real*8 a33,rot,b33
      integer i,j,k,l
      b33(:,:) = 0d0
      do 10 i=1,3
      do 10 j=1,3
      do 10 k=1,3
      do 10 l=1,3
         b33(i,j) = b33(i,j) + rot(i,k) * a33(k,l) * rot(j,l)
 10   continue
      return
      end subroutine rot_tensor
c-----------------------------------------------------------------------
c     Apply tensor inner dot
c     ci = aij x bj
      subroutine mult_array(aij,bj,ntens,ci)
      implicit none
      integer, intent(in) :: ntens
      dimension aij(ntens,ntens),bj(ntens),ci(ntens)
      real*8, intent(in)  :: aij,bj
      real*8, intent(out) :: ci
      integer i,j
cf2py intent(in) aij, bj, ntens
cf2py intent(out) ci
cf2py depend(ntens) aij, bj
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
c------------------------------------------------------------------------
      subroutine inner_dot_voigt(ntens,ndi,nshr,a,b,c)
      implicit none
      integer i
      integer, intent(in)::ndi,nshr,ntens
      dimension a(ntens), b(ntens)
      real*8, intent(in) ::a, b
      real*8, intent(out)::c
      c=0d0
      do 10 i=1,ndi
         c=c+a(i)*b(i)
 10   continue
      do 20 i=ndi+1, nshr+ndi
         c=c+a(i)*b(i)*2d0
 20   continue
      return
      end subroutine inner_dot_voigt
c------------------------------------------------------------------------
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

c-----------------------------------------------------------------------
      SUBROUTINE LU_INVERSE (A,N)
C *** INVERTS A MATRIX USING LU DECOMPOSITION
      DIMENSION A(N,N),Y(N,N),INDX(N)      ! MAY CHOKE SOME COMPILERS
C     DIMENSION A(5,5),Y(5,5),INDX(5)

C **************************************************************
C *** BLOCK ADDED 03/DEC/05 TO AVOID NUMERICALLY SINGULAR MATRIX
      AMAX=0.D0
      DO I=1,N
      DO J=1,N
        DUM=ABS(A(I,J))
        IF(DUM .GT. AMAX) AMAX=DUM
      ENDDO
      ENDDO
      DO I=1,N
      DO J=1,N
        A(I,J)=A(I,J)/AMAX      ! normalize the matrix
      ENDDO
      ENDDO
C **************************************************************

      DO I=1,N
        DO J=1,N
          Y(I,J)=0.
        ENDDO
        Y(I,I)=1.
      ENDDO

      CALL LUDCMP(A,N,N,INDX,D,ISINGULAR)
      IF(ISINGULAR.EQ.1) THEN
        WRITE(*,*) ' *** SINGULAR MATRIX IN LU_INVERSE !!'
        STOP -1
      ENDIF

      DO J=1,N
        CALL LUBKSB(A,N,N,INDX,Y(1,J))
      ENDDO

      DO I=1,N
      DO J=1,N
        A(I,J)=Y(I,J) /AMAX      ! renormalize the inverse
      ENDDO
      ENDDO

      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE ludcmp(a,n,np,indx,d,isingular)
      INTEGER n,np,indx(n),NMAX
      REAL d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k,isingular
      REAL aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
c
c        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
c
        if(aamax.eq.0.) then
        isingular=1
        return
        endif
c
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.

        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax

        if(a(j,j).eq.0.) a(j,j)=TINY

c       if(a(j,j).eq.0.) then
c       isingular=1
c       return
c       endif

        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
c
      isingular=0
c
      return
      END subroutine ludcmp
c-----------------------------------------------------------------------
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
