c-----------------------------------------------------------------------
      logical function isnan_in_arr(array,ndim)
      implicit none
      integer i,ndim
      real*8 array(ndim)
      logical iexit
      isnan_in_arr=.false.
      do 5 i=1,ndim
         if (isnan(array(i))) then
            isnan_in_arr=.true.
            goto 10
         endif
 5    continue
 10   return
      end function isnan_in_arr
c-----------------------------------------------------------------------
      logical function isnan_in_marr(array,ncol,nrow)
      implicit none
      integer i,j,ncol,nrow
      real*8 array(ncol,nrow)
      isnan_in_marr=.false.
      do 5 i=1,ncol
      do 5 j=1,nrow
         if (isnan(array(i,j))) then
            isnan_in_marr=.true.
            goto 10
         endif
 5    continue
 10   return
      end function isnan_in_marr
