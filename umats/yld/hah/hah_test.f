c     read_alpha is located in
      program test
c     Arguments
c     iyld_choice - isotropic yield surface kernel
c     stress      - cauchy stress tensor
c     phi         - HAH yield surface
c     dphi        - HAH yield surface 1st derivatives
c     d2phi       - HAH yield surface 2nd derivatives
c     yldc        - isotropic yield surface constants
c     yldp        - HAH yield surface contants/state variables
c     nyldc       - Len of yldc
c     nyldp       - Len of yldp
c     ntens       - Len of stress tensor
      implicit none
      integer ntens,nyldc,nyldp
      parameter(ntens=3,nyldc=9,nyldp=30)
      dimension yldc(nyldc),yldp(nyldp),stress(ntens)
      dimension dphi(ntens),d2phi(ntens,ntens)
      real*8 yldc,yldp,stress,phi,dphi,d2phi
      real*8 phi_chi
      integer iyld_choice

c     local - microstructure deviator
      dimension emic(ntens)
      real*8 emic
c     local - Bauschinger parameters
      dimension gk(4)
      dimension e_ks(5),aux_ten(ntens)
      dimension f_ks(2)
      real*8 gk,e_ks,f_ks,eeq,aux_ten
c     local - Latent hardening parameters
      real*8 gL,ekL,eL
c     local - cross hardening parameters
      real*8 gS,c_ks,ss
c     local - gen
      real*8 ref,hydro
c     local - controlling
      integer imsg
      logical idiaw

c      idiaw=.false.
      idiaw=.true.
      imsg=0


c     create a dummy microstructure deviator
      aux_ten(:)=0d0
      aux_ten(1)=1d0
      call deviat(ntens,aux_ten,emic,hydro)
c     isotropic conditions
      gL = 1d0
      gS = 1d0
      e_ks(:)=1d0               ! k1,k2,k3,k4,k5
      f_ks(:)=0d0               ! f1, f2 that are functions of (k1,k2,k3,k4,k5)
c      call read_alpha(
c     $     '/home/younguj/repo/abaqusPy/umats/yld/alfas.txt',yldc)

      yldc(:8) = 1d0
      yldc(9) = 8d0

      call hah_io(1,nyldp,ntens,yldp,emic,gk,e_ks,f_ks,eeq,ref,
     $     gL,ekL,eL,gS,c_ks,ss)

      iyld_choice=2             ! yld2000-2d

      if (idiaw) then
         call fill_line(imsg,'*',72)
      endif

      stress(:)=0d0
      stress(1)=1d0
      if (idiaw) then
         call w_chr(imsg,'cauchy stress')
         call w_dim(imsg,stress,ntens,1d0,.true.)
         call w_chr(imsg,'yldc')
         call w_dim(imsg,yldc,nyldc,1d0,.true.)
         call fill_line(imsg,'*',72)
         call w_chr(imsg,'just before entering hah')
         call w_chr(imsg,'cauchy stress')
         call w_dim(imsg,stress,ntens,1d0,.true.)
      endif
      call hah(iyld_choice,stress,phi,dphi,d2phi,yldc,yldp,nyldc,
     $     nyldp,ntens)
      if (idiaw) then
         call w_chr(imsg,'right after exit hah')
c$$$      call w_ival(imsg,'iyld_choice:',iyld_choice)
         call w_val( imsg,'phi_chi    :',phi_chi)
         call w_val( imsg,'phi        :',phi)
         call fill_line(imsg,'*',72)
      endif

c      call hah_uten(iyld_choice,yldc,nyldc,yldp,nyldp,ntens)
      call hah_locus(iyld_choice,yldc,nyldc,yldp,nyldp,ntens)

      end program test
c--------------------------------------------------------------------------------
c     Test one stress state
      subroutine one(iyld_choice,nyldc,nyldp,ntens,yldc,yldp,stress,
     $     d2phi,dphi,phi)
      implicit none
c     Arguments passed into
      integer iyld_choice,nyldc,nyldp,ntens
      dimension yldc(nyldc),yldp(nyldp),stress(ntens),
     $     d2phi(ntens,ntens),dphi(ntens)
      real*8 yldc,yldp,stress,d2phi,dphi,phi

      call hah(iyld_choice,stress,phi,dphi,d2phi,yldc,yldp,nyldc,nyldp,
     $     ntens)

      return
      end subroutine one
c--------------------------------------------------------------------------------
      subroutine hah_uten(iyld_choice,yldc,nyldc,yldp,nyldp,ntens)
      implicit none
c     Arguments passed into
      integer ntens,nyldc,nyldp
      dimension yldc(nyldc),yldp(nyldp)
      real*8 yldc,yldp
      integer iyld_choice

c     local variables.
      dimension dphi(ntens),d2phi(ntens,ntens),s33lab(3,3),s33mat(3,3),
     $     s6mat(6),dphi33m(3,3),s3mat(3),
     $     dphi33l(3,3),s6lab(6),dphi6(6)
      real*8 dphi,d2phi,pi,th,s33lab,s33mat,s6mat,time0,time1,
     $     dphi33m,s3mat,dphi33l,rv,dphi6,s6lab,phim
      integer nth,i,j
      parameter(nth=10)

      pi=4.d0*datan(1.d0)

      call cpu_time(time0)


      if (ntens.ne.3) then
         write(*,*)' *********************************************'
         write(*,*)' Warning: case that ntens not equal 3 was not '
         write(*,*)' throughly considered in hah_test.hah_uten    '
         write(*,*)' *********************************************'
         call exit(-1)
      endif

      s6lab(:)=0d0
      s6lab(1)=1d0
      call voigt2(s6lab,s33lab)

      write(*,'(a7,5(4a7,x,a1,x),2a7)')'th',
     $     's11_l','s22_l','s33_l','s12_l','|',
     $     's11_m','s22_m','s33_m','s12_m','|',
     $     's11_l','s22_l','s33_l','s12_l','|',
     $     'e11_m','e22_m','e33_m','e12_m','|',
     $     'e11_l','e22_l','e33_l','e12_l','|',
     $     'rv','phim'

      do 10 j=1,nth
         th = pi/2d0 - pi/2d0/(nth-1)*(j-1)
         write(*,'(f7.2)',advance='no') th*180.d0/pi
c$$$
c$$$  si_lab (sij_l)
c$$$
         write(*,'(4f7.2,x,a1,x)',advance='no')
     $        (s33lab(i,i),i=1,3),s33lab(1,2),'|'
c$$$
c$$$  si_mat (sij_m)
c$$$

         call inplane_rot(th,s33lab,s33mat)
         write(*,'(4f7.2,x,a1,x)',advance='no')
     $        (s33mat(i,i),i=1,3),s33mat(1,2),'|'
c$$$
c$$$  si_lab (sij_l)
c$$$
         call inplane_rot(th*(-1d0),s33mat,s33lab)
         write(*,'(4f7.2,x,a1,x)',advance='no')
     $        (s33lab(i,i),i=1,3),s33lab(1,2),'|'

c$$$
c$$$  ei_mat
c$$$

         call voigt1(s33mat,s6mat)
         call reduce_6to3(s6mat,s3mat)
c         write(*,*)'s3mat:',s3mat
c         call exit(-1)

         call hah(iyld_choice,s3mat,phim,dphi,d2phi,
     $        yldc,yldp,nyldc,nyldp,ntens)

c         write(*,*)'dphi:',dphi
         if (ntens.eq.3) then
            call reduce_3to6(dphi,dphi6)
            dphi6(3) = -dphi(1)-dphi(2)
            call voigt4(dphi6,dphi33m) !! this is wrong.
         elseif (ntens.eq.6) then
            call voigt4(dphi,dphi33m) !! this is wrong.
         else
            call exit(-1)
         endif

c         call w_mdim(0,dphi33m,3,1d0)
c         call exit(-1)
!        dphi in material axes
         write(*,'(4f7.2,x,a1,x)',advance='no')
     $        (dphi33m(i,i),i=1,3),dphi33m(1,2),'|'
c$$$
c$$$  ei_lab - strain in the lab space
c$$$
         call inplane_rot(th*(-1.d0),dphi33m,dphi33l)
         rv =-dphi33l(2,2)/(dphi33l(1,1)+dphi33l(2,2))
         write(*,'(4f7.2,x,a1,x,2f7.2)',advance='no')
     $        (dphi33l(i,i),i=1,3),dphi33l(1,2),'|',rv,phim
         write(*,*)
 10   continue

      call cpu_time(time1)
c      write(*,'(a,f7.1)') 'Elapsed time: [\mu s]',
c     $     (time1-time0)*1e6

      return
      end subroutine hah_uten
c-----------------------------------------------------------------------
c     Subroutine to draw yield locus
      subroutine hah_locus(iyld_choice,yldc,nyldc,yldp,nyldp,ntens)
c     Arguments
c     iyld_choice: yield function choice
c     yldc       : yield function constants
c     nyldc      : Len of yldc
c     yldp       : yield function parameters
c     nyldp      : Len of yldp
c     ntens      : Len of tensor
      implicit none
c     Arguments passed into
      integer, intent(in) :: iyld_choice
      integer, intent(in) :: ntens,nyldc,nyldp
      dimension yldc(nyldc),yldp(nyldp)
      real*8, intent(in) :: yldc
      real*8 yldp
c     Local variables.
      dimension d2phi(ntens),smat(ntens),dphi(ntens)
      real*8 dphi,d2phi,pi,th,time0,time1,
     $     phim,q,smat
      integer nth,j
      parameter(nth=6)

      call cpu_time(time0)
      open(1,file='ys.txt',status='unknown')

c     pi and yield surface exponent q stored in yldc
      pi=4.d0*datan(1.d0)
      q = yldc(9)

c$$$      write(*,'(a11,x,(4a9,x,a1,x),2(4a11,x,a1,x,a11,x))')'th',
c$$$     $     's1_m', 's2_m', 's3_m', 's6_m', '|',
c$$$     $     'e11_m','e22_m','e33_m','e12_m','|','phim',
c$$$     $     's1_m', 's2_m', 's3_m', 's6_m', '|','phim'

c$$$      s6mat_ref(:)=0d0 !! reference stress state (cauchy)
c$$$      s6mat_ref(1)=1d0
c$$$      call hah(iyld_choice,s6mat_ref,ref,dphi,d2phi,yldc,yldp,nyldc,
c$$$     $     nyldp,ntens)

      do 10 j=1,nth
         th = 2*pi/(nth-1)*(j-1)
         write(*,'(f11.2,a)',advance='no') th*180.d0/pi,'|'
         if (ntens.eq.3) then
            smat(1)   = dcos(th)
            smat(2)   = dsin(th)
            smat(3:ntens)   = 0.
         elseif (ntens.eq.6) then
            smat(1)   = dcos(th)
            smat(2)   = dsin(th)
            smat(3:ntens)   = 0.
         else
            write(*,*) 'unexpected dimension of ntens'
            call exit(-1)
         endif
         call w_dim(0,smat,ntens,1d0,.false.)
         call w_chrc(0,'|')
c$$$         write(*,'(4f11.3,x,a1,x)',advance='no')
c$$$     $        (smat(i),i=1,ntens),'|'
c$$$
c$$$  si_mat
c$$$
c         call deviat(ntens,smat,sdev,hydro)
c         write(*,*)'sdev:',sdev
c         call exit(-1)

c         call reduce_6to3(s6mat,s3mat)
c         call voigt2(s6mat,s33mat)
c         write(*,*)
c         write(*,*)'ntens---:',ntens
         call hah(iyld_choice,smat,phim,dphi,d2phi,
     $        yldc,yldp,nyldc,nyldp,ntens)
c          write(*,*)'ntens-:',ntens
c         call exit(-1)
c         call w_vals(0,phim)
c         call w_chrc(0,'|')
c         write(*,*) dphi(1)
c         write(*,*) dphi(2)
c         write(*,*) dphi(3)
c         write(*,*) smat
c         call w_dim(0,smat,ntens,1d0,.false.)
         call w_dim(0,dphi,ntens,1d0,.false.)
c         call exit(-1)

c$$$         write(*,'(4f11.3,x,a1,x,f11.3,x)',advance='no')
c$$$     $        (dphi(i),i=1,3),dphi(6),'|',phim


c$$$c     size adjustment
c$$$         dum=1d0/phim
c$$$         sdev(:) = sdev(:) * dum
c$$$         smat = sdev(:) + hydro
c$$$c         write(*,*)'dum:',dum
c$$$c         write(*,*)'sdev:',sdev
c$$$         call hah(iyld_choice,smat,phim,dphi,d2phi,
c$$$     $        yldc,yldp,nyldc,nyldp,ntens)
c$$$         write(*,'(4f11.3,x,a1,x,f11.3)',advance='no')
c$$$     $        (smat(i),i=1,3),smat(3),'|',phim
c$$$c         call exit(-1)
c$$$
c$$$c         write(*,*)
c$$$
c$$$c         s6mat(:) = (s6mat(:)/phim)**(yldp(9))
c$$$
c$$$c         call hah(iyld_choice,s6mat,phim,dphi,d2phi,
c$$$c     $        yldc,yldp,nyldc,nyldp,ntens)
c$$$
c$$$c         write(*,'(a,4f11.3,x,f11.7)',advance='no'), '|',
c$$$c     $        (s6mat(i),i=1,3),dphi(6),phim
c$$$
c$$$c         write(*,*)'phim:',phim
c$$$c          s6mat(:) = s6mat(:)    !!## /phim**8

         write(*,*)

c$$$         if (s6mat(1).eq.0 .and. s6mat(2).eq.0) then
c$$$            write(*,*)'Something went wrong in hah_test.hah_locus'
c$$$            call exit(-1)
c$$$         endif

c$$$c         call hah(iyld_choice,s6mat,phim,dphi,d2phi,
c$$$c     $        yldc,yldp,nyldc,nyldp,ntens)
c$$$
c$$$         call reduce_3to6(aux3,dphim)
c$$$         call voigt4(dphim,dphi33m)
c$$$!        dphi in material axes
c$$$
c$$$         call voigt2(s6mat,s33mat)
c$$$         write(*,'(4e11.3,x,a1,x)',advance='no')
c$$$     $        (s33mat(i,i),i=1,3),s33mat(1,2),'|'
c$$$
c$$$         write(*,'(4e11.3,x,a1,f11.7)',advance='no')
c$$$     $        (dphi33m(i,i),i=1,3),dphi33m(1,2),'|',phim
c$$$c         write(*,*)
c$$$
c$$$c         write(1,'(f11.2,2e13.5,f11.7)')
c$$$c     $        th*180d0/pi, s6mat(1), s6mat(2), phim

 10   continue

      call cpu_time(time1)
      write(*,'(a,f7.1)') 'Elapsed time: [\mu s]',
     $     (time1-time0)*1e6


      close(1)

      return
      end subroutine hah_locus
