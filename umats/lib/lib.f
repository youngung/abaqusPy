c-----------------------------------------------------------------------
      subroutine w33f(imsg,aux)
      implicit none
      integer imsg,i,j
      real*8 aux(3,3)
      do i=1,3
         write(imsg,'(3(e13.3,2x))') (aux(i,j),j=1,3)
      enddo
      return
      end subroutine w33f
c-----------------------------------------------------------------------
      subroutine w33i(imsg,aux)
      implicit none
      integer imsg,i,j
      integer aux(3,3)
      do i=1,3
         write(imsg,'(3(i13.3,2x))') (aux(i,j),j=1,3)
      enddo
      return
      end subroutine w33i
c-----------------------------------------------------------------------
      subroutine w_empty_lines(imsg,n)
      implicit none
      integer imsg, n, i
      do  i=1,n
         write(imsg,*)
      enddo
      return
      end subroutine w_empty_lines
c-----------------------------------------------------------------------
      subroutine w_mdim(imsg,array,ndi)
      implicit none
      integer ndi, nshr, imsg, i, j
      real*8 array(ndi,ndi)
      do i=1,ndi
         do j=1,ndi
            write(imsg, '(e13.3)',advance='no') array(i,j)
         enddo
         write(imsg,*)
      enddo
      return
      end subroutine w_mdim
c-----------------------------------------------------------------------
      subroutine w_dim(imsg,array,ndi,fact,ibr)
c     imsg: file ID
c     array: 1-Dimensional array
c     ndi: size of the array
c     fact:multiplicative factor to scale the elements in the array
c     ibr: flag to insert line-breaker
      implicit none
      character*80 fmt
      integer ndi,nshr,imsg,i,j
      real*8 array(ndi),fact,mxval,get_mx,mxv
      logical ibr
      mxv = get_mx(array,ndi)
      if (mxv.le.1.) then
         fmt="(e9.2,x)"
      elseif ((mxv.gt.1.).and.(mxv.lt.1e5)) then
         fmt="(f5.2)"
      elseif ((mxv.ge.1e5)) then
         fmt="(e9.2,x)"
      else
         stop
      endif
      do i=1,ndi
         write(imsg,fmt,advance='no') array(i) * fact
      enddo
      if (ibr) write(imsg,*)
      return
      end subroutine w_dim
c-----------------------------------------------------------------------
      real*8 function get_mx(array,ndi)
      implicit none
      integer i,ndi
      real*8 array(ndi)
      get_mx=array(1)
      do i=2,ndi
         if (array(i).gt.get_mx) then
            get_mx=array(i)
         endif
      enddo
      return
      end function
