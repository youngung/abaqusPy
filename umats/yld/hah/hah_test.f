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
      dimension dphi_chi(ntens),d2phi_chi(ntens,ntens)
      real*8 phi_chi,dphi_chi,d2phi_chi
      integer iyld_choice,i

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

      aux_ten(:)=0d0
      aux_ten(1)=1d0
      gL = 1d0
      gS = 1d0
      e_ks(:)=1d0               ! k1,k2,k3,k4,k5
      f_ks(:)=0d0               ! f1, f2 that are functions of (k1,k2,k3,k4,k5)

      call deviat(ntens,aux_ten,emic,hydro)
      call hah_io(1,nyldp,ntens,yldp,emic,gk,e_ks,f_ks,eeq,ref,
     $     gL,ekL,eL,gS,c_ks,ss)

      iyld_choice=2             ! yld2000-2d

      if (idiaw) then
         call fill_line(imsg,'*',72)
      endif

c      call read_alpha(
c     $     '/home/younguj/repo/abaqusPy/umats/yld/alfas.txt',yldc)

      yldc(:8) = 1d0
      yldc(9) = 8d0

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


      call hah_uten(iyld_choice,yldc,nyldc,yldp,nyldp,ntens)
c      call hah_locus(iyld_choice,yldc,nyldc,yldp,nyldp,ntens)

      end program test
c--------------------------------------------------------------------------------

      subroutine hah_uten(iyld_choice,yldc,nyldc,yldp,nyldp,ntens)
      implicit none
c     Arguments passed into
      integer ntens,nyldc,nyldp
      dimension yldc(nyldc),yldp(nyldp)
      real*8 yldc,yldp
      integer iyld_choice

c     local variables.
      dimension dphi(ntens),d2phi(ntens),s33lab(3,3),s33mat(3,3),
     $     s6mat(6),aux3(3),aux33(3,3),dphi33m(3,3),s3mat(3),
     $     dphi33l(3,3),s6lab(6)
      real*8 phi,dphi,d2phi,pi,th,s33lab,s33mat,s6mat,time0,time1,
     $     aux3,aux33,phim,dphim,dphi33m,s3mat,dphi33l,rv,
     $     s6lab
      integer nth,i,j
      parameter(nth=5)

      pi=4.d0*datan(1.d0)

      call cpu_time(time0)

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
c$$$  si_lab
c$$$
         write(*,'(4f7.2,x,a1,x)',advance='no')
     $        (s33lab(i,i),i=1,3),s33lab(1,2),'|'

c$$$
c$$$  si_mat
c$$$
         call inplane_rot(th,s33lab,s33mat)
         write(*,'(4f7.2,x,a1,x)',advance='no')
     $        (s33mat(i,i),i=1,3),s33mat(1,2),'|'
c$$$
c$$$  si_lab
c$$$
         call inplane_rot(th*(-1d0),s33mat,s33lab)
         write(*,'(4f7.2,x,a1,x)',advance='no')
     $        (s33lab(i,i),i=1,3),s33lab(1,2),'|'
c$$$
c$$$  ei_mat
c$$$
         call voigt1(s33mat,s6mat)

         call reduce_6to3(s6mat,s3mat)
         call hah(iyld_choice,s6mat,phim,dphi,d2phi,
     $        yldc,yldp,nyldc,nyldp,ntens)

         call reduce_3to6(aux3,dphim)
         call voigt4(dphim,dphi33m)
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
      write(*,'(a,f7.1)') 'Elapsed time: [\mu s]',
     $     (time1-time0)*1e6

      return
      end subroutine hah_uten
c-----------------------------------------------------------------------
      subroutine hah_locus(iyld_choice,yldc,nyldc,yldp,nyldp,ntens)
      implicit none
c     Arguments passed into
      integer ntens,nyldc,nyldp
      dimension yldc(nyldc),yldp(nyldp)
      real*8 yldc,yldp
      integer iyld_choice

c     local variables.
      dimension dphi(ntens),d2phi(ntens),s33lab(3,3),s33mat(3,3),
     $     s6mat(6),aux3(3),aux33(3,3),dphi33m(3,3),s3mat(3),
     $     dphi33l(3,3),s6lab(6),s6mat_ref(ntens)
      real*8 phi,dphi,d2phi,pi,th,s33lab,s33mat,s6mat,time0,time1,
     $     aux3,aux33,phim,dphim,dphi33m,s3mat,dphi33l,rv,
     $     s6lab,q,ref,s6mat_ref
      integer nth,i,j
      parameter(nth=9)

      call cpu_time(time0)

      open(1,file='ys.txt',status='unknown')

c     pi and yield surface exponent q stored in yldc
      pi=4.d0*datan(1.d0)
      q = yldc(9)


      s6lab(:)=0d0
      s6lab(1)=2d0
      call voigt2(s6lab,s33lab)

      write(*,'(a11,3(4a11,x,a1,x),a11)')'th',
     $     's1_m','s2_m','s3_m','s6_m','|',
     $     's11_m','s22_m','s33_m','s12_m','|',
     $     'e11_m','e22_m','e33_m','e12_m','|','phim'

      s6mat_ref(:)=0d0 !! reference stress state (cauchy)
      s6mat_ref(1)=1d0
      call hah(iyld_choice,s6mat_ref,ref,dphi,d2phi,
     $     yldc,yldp,nyldc,nyldp,ntens)

c      ref = ref ** q
c      write(*,'(a,f11.2)')'ref:',ref

      do 10 j=1,nth
         th = 2*pi/(nth-1)*(j-1)

         write(*,'(f11.2,a)',advance='no') th*180.d0/pi,'|'

         s6mat(1)   = dcos(th)
         s6mat(2)   = dsin(th)
         s6mat(3:6) = 0d0

c         do i=1,6
c            s6mat(i) = s6mat(i)/ref
c         enddo
         write(*,'(4e11.1,x,a1,x)',advance='no')
     $        (s6mat(i),i=1,3),s6mat(6),'|'

c$$$
c$$$  si_mat
c$$$

c$$$
c$$$  ei_mat
c$$$
         call voigt1(s33mat,s6mat)

         call reduce_6to3(s6mat,s3mat)
         call hah(iyld_choice,s6mat,phim,dphi,d2phi,
     $        yldc,yldp,nyldc,nyldp,ntens)

c         write(*,*)'phim:',phim

         s6mat(:) = s6mat(:)/phim**8


         if (s6mat(1).eq.0 .and. s6mat(2).eq.0) then
            write(*,*)'something went wrong'
            call exit(-1)
         endif

c         call hah(iyld_choice,s6mat,phim,dphi,d2phi,
c     $        yldc,yldp,nyldc,nyldp,ntens)

         call reduce_3to6(aux3,dphim)
         call voigt4(dphim,dphi33m)
!        dphi in material axes

         call voigt2(s6mat,s33mat)
         write(*,'(4e11.3,x,a1,x)',advance='no')
     $        (s33mat(i,i),i=1,3),s33mat(1,2),'|'

         write(*,'(4e11.3,x,a1,f11.7)',advance='no')
     $        (dphi33m(i,i),i=1,3),dphi33m(1,2),'|',phim
c         write(*,*)

c         write(1,'(f11.2,2e13.5,f11.7)')
c     $        th*180d0/pi, s6mat(1), s6mat(2), phim

 10   continue

      call cpu_time(time1)
      write(*,'(a,f7.1)') 'Elapsed time: [\mu s]',
     $     (time1-time0)*1e6


      close(1)

      return
      end subroutine hah_locus
