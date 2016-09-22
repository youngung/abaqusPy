c-----------------------------------------------------------------------
c     subroutine that stores/restores state variables.
c     This subroutine makes transactions between statev and other
c     variables passed as arguments.
      subroutine restore_statev(statev,nstatv,eqpl,stran_el,stran_pl,
     $     ntens,yldp,nyldp,iopt,verbose,iunit,

     $     iw,kinc,noel,npt,time,stress)
c-----------------------------------------------------------------------
c     Arguments
c     statev  : state variable array
c     nstatv  : len of statev
c     eqpl    : equivalent plastic strain
c     stran_el: cumulative elastic strain
c     stran_pl: cumulative plastic strain
c     ntens   : Len of stran_el and stran_pl
c     yldp    : parameters for yield surface
c     nyldp   : Len of yldp ! mostly static as for isotropic hardening...
c     iopt    : option to define the behavior of restore_statev
c                  0: read from statev
c                  1: save to statev
c     verbose : flag to be verbose or not
c     iunit   : file unit (if not zero) or std to which
c               inquiry stream will be written
c     kinc    : increment number
c     noel    : element number
c     npt     : integration point number
c     time    : time array passed to umat
c     stress  : stress is not stored to statev
c               (this is just to write to fn_stav)
c-----------------------------------------------------------------------
      implicit none
      integer nstatv,ntens
      dimension statev(nstatv),stran_el(ntens),stran_pl(ntens),
     $     yldp(nyldp),stress(ntens)
      real*8 statev,eqpl,stran_el,stran_pl,yldp,time,stress
      character*80 fmt,fnstv,ncc
      integer iopt,i,nyldp,iunit,kinc,noel,npt
      logical verbose,iw
      if (iopt.eq.0) then
         ! read from statev
         eqpl = statev(1)
         do 5 i=1,ntens
            stran_el(i) = statev(i+1)
            stran_pl(i) = statev(i+1+ntens)
 5       continue
         do 10 i=1,nyldp
            yldp(i)     = statev(i+1+ntens*2)
 10      continue
      elseif (iopt.eq.1) then
         ! save to statev
         ! read from statev
         statev(1) = eqpl
         do 15 i=1,ntens
            statev(i+1)       = stran_el(i)
            statev(i+1+ntens) = stran_pl(i)
 15      continue
         do 20 i=1,nyldp
            statev(i+1+ntens*2) = yldp(i)
 20      continue
      else
         write(*,*) 'Unexpected iopt given'
         stop
      endif
      if (verbose) then
         call w_empty_lines(iunit,2)
         call fill_line(iunit,'*',72)
         call w_chr(iunit,'inquiry request on state variables')
         call w_val(iunit,'eqpl    :',eqpl)
         call w_chr(iunit,'stran_el:')
         call w_dim(iunit,stran_el,ntens,1.d0,.true.)
         call w_chr(iunit,'stran_pl:')
         call w_dim(iunit,stran_pl,ntens,1.d0,.true.)
         call w_chr(iunit,'yldp    :')
         call w_dim(iunit,yldp,nyldp,1.d0,.true.)
         call fill_line(iunit,'*',72)
         call w_empty_lines(iunit,2)
      endif

      if (iw) then
         fnstv='/home/younguj/repo/abaqusPy/examples/one/statev.txt'
         open(425,position='append',file=fnstv)
         write(ncc,'(i)')ntens+nstatv
         fmt='('//trim(adjustl(ncc))//'(e12.4,1x))'
         write(425,fmt,advance='no') time
         write(425,fmt) (stress(i),i=1,ntens),(statev(i),i=1,nstatv)
         ! write(425,*)           ! line-breaker
         close(425)
      endif

      return
      end subroutine restore_statev
c-----------------------------------------------------------------------
