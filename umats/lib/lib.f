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
      end function get_mx
c-----------------------------------------------------------------------
      real*8 function get_mmx(array,ncol,nrow)
      implicit none
      integer i,j,ncol,nrow
      real*8 array(ncol,nrow)
      get_mmx=array(1,1)
      do i=1,ncol
         do j=1,nrow
            if (array(i,j).gt.get_mmx) then
               get_mmx=array(i,j)
            endif
         enddo
      enddo
      return
      end function get_mmx
c-----------------------------------------------------------------------
      subroutine pjoin(pardir,fn,fullfn)
      implicit none
      character(len=255) :: path,fn,fullfn,pardir
      character(len=1) :: sep
      call get_environment_variable('PATH',path)
      sep=path(1:1)
      fullfn = trim(pardir) // sep // trim(fn)
      end subroutine
c-----------------------------------------------------------------------
