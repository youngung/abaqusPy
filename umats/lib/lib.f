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
c     imsg: file unit - if imsg =0 use std (*)
      implicit none
      integer imsg, n, i
      if (imsg.eq.0) then
         do i=1,n
            write(*,*)
         enddo
      else
         do i=1,n
            write(imsg,*)
         enddo
      endif
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
      subroutine w_dim(imsg,array,ntens,fact,ibr)
c     imsg: file ID
c     array: 1-Dimensional array
c     ntens: size of the array
c     fact:multiplicative factor to scale the elements in the array
c     ibr: flag to insert line-breaker
      implicit none
      character*80 fmt,clen
      dimension array(ntens),brray(ntens)
      integer ntens,nshr,imsg,i,j
      real*8 array,fact,mxval,get_mx,mxv,brray
      logical ibr
      mxv = get_mx(array,ntens)
      if (mxv.le.1.) then
         fmt="e9.2,x"
      elseif ((mxv.gt.1.).and.(mxv.lt.1e5)) then
         fmt="f5.2"
      elseif ((mxv.ge.1e5)) then
         fmt="e9.2,x"
      else
         write(*,*) 'Err: Unexpected case of max value'
         stop
      endif

      do i=1,ntens
         brray(i) = array(i) * fact
      enddo

      if (imsg.ne.0) then
         fmt = '('//trim(fmt)//')'
         write(*,*)'check fmt for imsg.ne.0'
         write(*,*)'fmt:',fmt
         write(imsg,fmt,advance='no') (brray(i),i=1,ntens)
         if (ibr) write(imsg,*)
      else
         write(*,*)'check fmt for imsg.eq.0'
         write(clen,"(i2)") ntens
         fmt = '('//trim(clen)//trim(fmt)//')'
         write(*,*)'fmt:',fmt
         write(*   ,fmt) (brray(i),i=1,ntens)
      endif
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
c-----------------------------------------------------------------------
      subroutine pjoin(pardir,fn,fullfn)
      implicit none
      character(len=255) :: path,fn,fullfn,pardir
      character(len=1) :: sep
      call get_environment_variable('PATH',path)
      sep=path(1:1)
      fullfn = trim(pardir) // sep // trim(fn)
      end subroutine
