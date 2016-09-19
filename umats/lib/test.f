      program test_vm
      implicit none
      integer nth,i,j,ii,jj
      parameter(nth=15)
      real*8 th,s33(3,3),e33(3,3),s6(6),pi,s6lab(6),
     $     s6mat(6),s33lab(3,3),s33mat(3,3),phim,dphim(6),d2phim(6,6),
     $     phil,dphil(6),d2phil(6,6),dphi33m(3,3),dphi33l(3,3),
     $     dphi33ld(3,3),am
      pi=4.d0*datan(1.d0)

c     Uniaxial tensin stress state referred in the lab axes
      s6lab(:)=0.
      s6lab(1)=1.
      call voigt2(s6lab,s33lab)

      write(*,'(a7,5(4a7,x,a1,x))')'th',
     $     's11_l','s22_l','s33_l','s12_l','|',
     $     's11_m','s22_m','s33_m','s12_m','|',
     $     's11_l','s22_l','s33_l','s12_l','|',
     $     'e11_m','e22_m','e33_m','e12_m','|',
     $     'e11_l','e22_l','e33_l','e12_l','|'

      do 10 j=1,nth
         th = pi/2.d0 - pi/2.d0/(nth-1) * (j-1)
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
         call inplane_rot(th*(-1.),s33mat,s33lab)
         write(*,'(4f7.2,x,a1,x)',advance='no')
     $        (s33lab(i,i),i=1,3),s33lab(1,2),'|'
c$$$
c$$$  ei_mat
c$$$
         call voigt1(s33mat,s6mat)
c        yield stress is written in the material axes
         call vm_gen(s6mat,phim,dphim,d2phim)
!        shear strains: =1/2.shear
         call voigt4(dphim,dphi33m)
         !! dphi in material axes
         write(*,'(4f7.2,x,a1,x)',advance='no')
     $        (dphi33m(i,i),i=1,3),dphi33m(1,2),'|'
c$$$
c$$$  ei_lab - strain in the lab space
c$$$
         call inplane_rot(th*(-1.d0),dphi33m,dphi33l)
         write(*,'(4f7.2,x,a1,x)',advance='no')
     $        (dphi33l(i,i),i=1,3),dphi33l(1,2),'|'
         write(*,*)

 10   continue
      return
      end program
c-----------------------------------------------------------------------
c$$$!     pal
      include "/home/younguj/repo/abaqusPy/umats/lib/vm.f"
!     mac
c      include "/Users/yj/repo/abaqusPy/umats/lib/vm.f"
