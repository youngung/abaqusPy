c-----------------------------------------------------------------------
      real*8 function get_mx(array,ndi)
      implicit none
      dimension array(ndi)
      real*8, intent(in):: array
      integer, intent(in):: ndi
      integer i

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
      integer, intent(in) :: ncol,nrow
      dimension array(ncol,nrow)
      real*8, intent(in) ::  array
      integer i,j
      get_mmx=array(1,1)
      do 5 i=1,ncol
      do 5 j=1,nrow
         if (array(i,j).gt.get_mmx) then
            get_mmx=array(i,j)
         endif
 5    continue
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
