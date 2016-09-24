c-----------------------------------------------------------------------
c     to test various yield functions
      program test_yld
      implicit none
      integer nth,i,j,ii,jj
      parameter(nth=10)
      real*8 th,s33(3,3),e33(3,3),s6(6),pi,s6lab(6),s6mat(6),s3mat(3),
     $     s33lab(3,3),s33mat(3,3),phim,dphim(6),d2phim(6,6),phil,
     $     dphil(6),d2phil(6,6),dphi33m(3,3),dphi33l(3,3),
     $     dphi33ld(3,3),yldp_hill(6),yldp_yld2000(9),aux33(3,3),
     $     aux3(3),bux3(3)
      pi=4.d0*datan(1.d0)

c     Uniaxial tensin stress state referred in the lab axes
      s6lab(:)=0.
      s6lab(1)=1.d0
      call voigt2(s6lab,s33lab)

c-----------------------------------------------------------------------
      yldp_hill(1)=0.5d0
      yldp_hill(2)=0.5d0
      yldp_hill(3)=0.5d0
      yldp_hill(4)=1.5d0
      yldp_hill(5)=1.5d0
      yldp_hill(6)=1.5d0

      yldp_yld2000(1)=1d0
      yldp_yld2000(2)=1d0
      yldp_yld2000(3)=1d0
      yldp_yld2000(4)=1d0
      yldp_yld2000(5)=1d0
      yldp_yld2000(6)=1d0
      yldp_yld2000(7)=1d0
      yldp_yld2000(8)=1d0
      yldp_yld2000(9)=8d0

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

c         call vm_gen(    s6mat,phim,dphim,d2phim)
c         call hill48_gen(s6mat,phim,dphim,d2phim,yldp_hill)

         call reduce_6to3(s6mat,s3mat)
         call yld2000_2d(s3mat,phim,aux3,aux33,yldp_yld2000)
         call reduce_3to6(aux3,dphim)

!        shear strains: =1/2.shear
         call voigt4(dphim,dphi33m)
!        dphi in material axes
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
      subroutine reduce_6to3(a6,a3)
      dimension a6(6),a3(3)
      real*8 a6,a3
      a3(1) = a6(1)
      a3(2) = a6(2)
      a3(3) = a6(6)
      return
      end subroutine reduce_6to3
c-----------------------------------------------------------------------
      subroutine reduce_3to6(a3,a6)
      dimension a6(6),a3(3)
      real*8 a6,a3
      a6(1) = a3(1)
      a6(2) = a3(2)
      a6(6) = a3(3)
      return
      end subroutine reduce_3to6
c-----------------------------------------------------------------------
c$$$!     pal
      include "/home/younguj/repo/abaqusPy/umats/yld/vm.f"
      include "/home/younguj/repo/abaqusPy/umats/yld/hill48.f"
      include "/home/younguj/repo/abaqusPy/umats/yld/yld2000_2d.f"
      include "/home/younguj/repo/abaqusPy/umats/lib/lib_write.f"
      include "/home/younguj/repo/abaqusPy/umats/lib/is.f"
      include "/home/younguj/repo/abaqusPy/umats/lib/lib.f"
c$$$!     mac
c$$$      include "/Users/yj/repo/abaqusPy/umats/yld/vm.f"
c$$$      include "/Users/yj/repo/abaqusPy/umats/yld/hill48.f"
c$$$      include "/Users/yj/repo/abaqusPy/umats/yld/yld2000_2d.f"
