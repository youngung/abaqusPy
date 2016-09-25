c-----------------------------------------------------------------------
c     write an empty line(s)
      subroutine w_empty_lines(imsg,n)
c     Arguments
c     imsg: file unit - if imsg =0 use std (*)
c     n   : the number of lines to be emptied
      implicit none
      integer imsg,n,i
      if (imsg.eq.0) then
         do 5 i=1,n
            write(*,*)
 5       continue
      else
         do 10 i=1,n
            write(imsg,*)
 10      continue
      endif
      return
      end subroutine w_empty_lines
c-----------------------------------------------------------------------
c     Fill a line(s) with the given <chr> by repeating it for <n> times
      subroutine fill_line(imsg,chr,n)
c     imsg: file unit - if imsg =0 use std (*)
c     chr : symbol that will fill the lines
c     n   : the number of repetition
      implicit none
      character*1 chr
      integer imsg, n, i
      if (imsg.eq.0) then
         do 5 i=1,n
            write(*,   '(a1)',advance='no') chr
 5       continue
         write(*,*)
      else
         do 10 i=1,n
            write(imsg,'(a1)',advance='no') chr
 10      continue
         write(imsg,*)
      endif

      return
      end subroutine fill_line
c-----------------------------------------------------------------------
c     Get the suitable format to print the given value <val>
      character*80 function get_fmt(val)
c     Argument
c     val: the value
      implicit none
      real*8 val,aval
      logical is_inf
      aval = dabs(val)
      if     (is_inf(aval))                      then
         get_fmt='f10.0,x'
      elseif (                   (aval.lt.1e-3)) then
         get_fmt="e12.4,x"
      elseif ((aval.ge.1e-3).and.(aval.lt.1e0 )) then
         get_fmt="f9.5,x"
      elseif ((aval.ge.1e0 ).and.(aval.lt.1e3 )) then
         get_fmt="f9.3,x"
      elseif ((aval.ge.1e3 ).and.(aval.lt.1e6 )) then
         get_fmt="f9.0,x"
      elseif ((aval.ge.1e6 ))                    then
         get_fmt="e12.4,x"
      else
         write(*,*) 'Err: Unexpected case of max value'
         stop -1
      endif
      return
      end function get_fmt
c-----------------------------------------------------------------------
c     Given the character value <str> find a suitible format to write
c     this value
      character*80 function get_fmt_str(str)
c     Argument
c     str: the subjected character with any arbitrary size
      implicit none
      character(len=*) str
      character*80 ncc
      integer nc
      logical is_inf
      nc = len(str)
      write(ncc,'(i80)') nc
      ncc = adjustl(ncc)
      get_fmt_str='(a'//trim(ncc)//',x)'
      return
      end function get_fmt_str
c-----------------------------------------------------------------------
c     Given the integer value <ival> find a suitible format to write
c     this value.
      character*80 function get_fmt_int(ival)
c     Argument
c     ival: integer value
      implicit none
      integer ival,nc
      character*80 ncc
      character*20 nc_no
      write(ncc,'(i10)')ival
      ncc = adjustl(ncc)
      nc=len(trim(ncc)) ! number of character for this integer value
      write(nc_no,'(i10)') nc
      get_fmt_int='(i'//trim(adjustl(nc_no))//',x)'
      return
      end function get_fmt_int
c-----------------------------------------------------------------------
      subroutine w_val(iunit,str,v)
c     Arguments
c     iunit: file ID
c     str : chracter that describe the passed value
c     v   : the value
      implicit none
      integer iunit,nc
      character*80 fmt,get_fmt,ncc
      character(len=*) str
      real*8 v
      fmt = get_fmt(v)
      nc = len(str)
      write(ncc,'(i80)') nc
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
      subroutine w_ival(iunit,str,iv)
c     Arguments
c     iunit: file ID
c     str : chracter that describe the passed value
c     iv  : the integer value
      implicit none
      integer iunit,nc
      character*80 fmt,get_fmt_int,ncc
      character(len=*) str
      integer iv
      fmt = get_fmt_int(iv)
      nc = len(str)
      write(ncc,'(i10)') nc
      ncc=adjustl(ncc)
      fmt = trim('(a')//trim(ncc)//',x,'//trim(fmt)//')'
      if (iunit.eq.0) then
         write(*,    trim(fmt)) trim(str),iv
      else
         write(iunit,trim(fmt)) trim(str),iv
      endif
      return
      end subroutine w_ival
c-----------------------------------------------------------------------
c     print out a 2nd order square matrix array
      subroutine w_mdim(iunit,array,ndi,fact)
c     iunit: file ID
c     array: 2-Dimensional array
c     ndi: size of the (ndi x ndi) array
c     fact: multiplicative factor to scale the elements in the array
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
         do 15 i=1,ndi
            write(iunit,fmt) (brray(i,j),j=1,ndi)
 15      continue
      else
         do 20 i=1,ndi
            write(*   ,fmt) (brray(i,j),j=1,ndi)
 20      continue
      endif
      return
      end subroutine w_mdim
c-----------------------------------------------------------------------
c     print out a 2nd order non-square matrix array
      subroutine w_mndim(iunit,array,ndi,nci,fact)
c     iunit: file ID
c     array: 2-Dimensional array
c     ndi: size of the (ndi x ndi) array
c     fact: multiplicative factor to scale the elements in the array
      implicit none
      character*80 fmt,clen,get_fmt
      integer ndi,nci,iunit,i,j
      dimension array(ndi,nci),brray(ndi,nci)
      real*8 array,brray,get_mmx,mxv,fact
      do 10 i=1,ndi
      do 10 j=1,nci
         brray(i,j) = array(i,j) * fact
 10   continue
      mxv = get_mmx(brray,ndi,nci)
      fmt = get_fmt(mxv) ! common filter

      write(clen,"(i2)") nci
      fmt = '('//trim(clen)//trim(fmt)//')'

      if (iunit.ne.0) then
         do 15 i=1,ndi
            write(iunit,fmt) (brray(i,j),j=1,nci)
 15      continue
      else
         do 20 i=1,ndi
            write(*   ,fmt) (brray(i,j),j=1,nci)
 20      continue
      endif
      return
      end subroutine w_mndim
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
      subroutine w_chr(iunit,str)
c     iunit: file ID
c     chr : chracter
      implicit none
      integer iunit,nc
      character*80 fmt,get_fmt_str,ncc
      character(len=*) str
      fmt = get_fmt_str(str)
      if (iunit.eq.0) then
         write(*,    trim(fmt)) trim(str)
      else
         write(iunit,trim(fmt)) trim(str)
      endif
      return
      end subroutine w_chr
c-----------------------------------------------------------------------
c     print header
      subroutine print_head(i)
c     i: file unit (use std if 0, use imsg elsewhere)
      integer i
      call w_empty_lines(i,2)
      if (i.eq.0) then ! USE std
         write(*,*) '*------------------*'
         write(*,*) '|       UMAT       |'
         write(*,*) '|                  |'
         write(*,*) '*------------------*'
         write(*,*) 'Youngung Jeong'
         write(*,*) 'youngung.jeong@gmail.com'
      else
         write(i,*) '*------------------*'
         write(i,*) '|       UMAT       |'
         write(i,*) '*------------------*'
         write(i,*) 'Youngung Jeong'
         write(i,*) 'youngung.jeong@gmail.com'
      endif
      call w_empty_lines(i,1)
      return
      end subroutine print_head
c-----------------------------------------------------------------------
c     print footer
      subroutine print_foot(i)
c     i: file unit (use std if 0, use imsg elsewhere)
      integer i
      call w_empty_lines(i,2)
      if (i.eq.0) then ! USE std
         write(*,*) '*------------------*'
         write(*,*) '|       END        |'
         write(*,*) '*------------------*'
      else
         write(i,*) '*------------------*'
         write(i,*) '|       END        |'
         write(i,*) '*------------------*'
      endif
      call w_empty_lines(i,1)
      return
      end subroutine print_foot
c-----------------------------------------------------------------------
      subroutine stop_debug(iunit)
      implicit none
      integer iunit
      call w_empty_lines(iunit,3)
      call fill_line(iunit,'*',52)
      call w_chr(iunit, '--  Reached safely at debugging the point  --')
      call fill_line(iunit,'*',52)
      call w_empty_lines(iunit,3)
      stop -1
      end subroutine
c-----------------------------------------------------------------------

c$$$      program test
c$$$      character*80 get_fmt_int,get_fmt
c$$$      integer i,n
c$$$      real*8 val
c$$$
c$$$      do i=-9,9
c$$$         val = 10d0**i
c$$$         !write(*,*) i, trim(get_fmt(val))
c$$$         call w_val(0,'dum',val)
c$$$      enddo
c$$$
c$$$c$$$      call w_chr(0,'dum')
c$$$c$$$      call w_chr(0,'dum2439 33')
c$$$
c$$$c$$$      write(*,*) get_fmt_int(1)
c$$$c$$$      write(*,*) get_fmt_int(10)
c$$$c$$$      write(*,*) get_fmt_int(100)
c$$$c$$$      write(*,*) get_fmt_int(10000)
c$$$
c$$$      end program test
c$$$      include '/home/younguj/repo/abaqusPy/umats/lib/lib.f'
c$$$      include '/home/younguj/repo/abaqusPy/umats/lib/is.f'
