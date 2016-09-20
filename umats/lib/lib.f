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
      character*80 function get_fmt(val)
      implicit none
      real*8 val,aval
      aval = dabs(val)
      if (                       (aval.lt.1e-3)) then
         get_fmt="e10.2,x"
      elseif ((aval.ge.1e-3).and.(aval.lt.1e0 )) then
         get_fmt="f9.5,x"
      elseif ((aval.ge.1e0 ).and.(aval.lt.1e3 )) then
         get_fmt="f9.3,x"
      elseif ((aval.ge.1e3 ).and.(aval.lt.1e6 )) then
         get_fmt="f9.0,x"
      elseif ((aval.ge.1e6 ))                   then
         get_fmt="e10.2,x"
      else
         write(*,*) 'Err: Unexpected case of max value'
         stop
      endif

c     TEST
c$$$      write(*,*)'val:',aval
c$$$      write(*,*)'fmt:'
c$$$      write(*,*) get_fmt

      return
      end function
c-----------------------------------------------------------------------
      subroutine w_val(iunit,str,v)
c     iunit: file ID
c     chr : chracter
      implicit none
      integer iunit,nc
      character*80 fmt,get_fmt,ncc
      character(len=*) str
      real*8 v
      fmt = get_fmt(v)
      nc = len(str)
      write(ncc,'(i)') nc
      ncc=adjustl(ncc)
      fmt = trim('(a')//trim(ncc)//',x,'//trim(fmt)//')'
      if (iunit.eq.0) then
         write(*,    trim(fmt)) trim(str),v
      else
         write(iunit,trim(fmt)) trim(str),v
      endif
      return
      end subroutine w_val
c-----------------------------------------------------------------------
      subroutine w_mdim(iunit,array,ndi,fact)
c     iunit: file ID
c     array: 1-Dimensional array
c     ndi: size of the (ndi x ndi) array
c     fact:multiplicative factor to scale the elements in the array
c     ibr: flag to insert line-breaker
      implicit none
      character*80 fmt,clen,get_fmt
      integer ndi,iunit, i, j
      real*8 array(ndi,ndi),brray(ndi,ndi),get_mmx,mxv,fact
      do 10 i=1,ndi
      do 10 j=1, ndi
         brray(i,j) = array(i,j) * fact
 10   continue
      mxv = get_mmx(brray,ndi,ndi)
      fmt = get_fmt(mxv) ! common filter

      write(clen,"(i2)") ndi
      fmt = '('//trim(clen)//trim(fmt)//')'

      if (iunit.ne.0) then
         do i=1,ndi
            write(iunit,fmt) (brray(i,j),j=1,ndi)
         enddo
      else
         do i=1,ndi
            write(*   ,fmt) (brray(i,j),j=1,ndi)
         enddo
      endif
      return
      end subroutine w_mdim
c-----------------------------------------------------------------------
      subroutine w_dim(iunit,array,ndi,fact,ibr)
c     iunit: file ID
c     array: 1-Dimensional array
c     ndi: size of the array
c     fact:multiplicative factor to scale the elements in the array
c     ibr: flag to insert line-breaker
      implicit none
      character*80 fmt,clen,get_fmt
      dimension array(ndi),brray(ndi)
      integer ndi,iunit,i,j
      real*8 array,fact,mxval,get_mx,mxv,brray
      logical ibr
      do 10 i=1,ndi
         brray(i) = array(i) * fact
 10   continue
      mxv = get_mx(brray,ndi)
      fmt = get_fmt(mxv)

      if (iunit.ne.0) then
         fmt = '('//trim(fmt)//')'
         write(iunit,fmt,advance='no') (brray(i),i=1,ndi)
         if (ibr) write(iunit,*)
      else
         write(clen,"(i2)") ndi
         fmt = '('//trim(clen)//trim(fmt)//')'
         write(*   ,fmt) (brray(i),i=1,ndi)
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
