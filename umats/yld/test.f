c-----------------------------------------------------------------------
c     to test various yield functions
      program test_yld
      implicit none
      integer nth,nyldc
      parameter(nth=15)
      dimension yldc(9),rs(nth),ys(nth),locus(nth,2)
      real*8 yldc,rs,ys,locus
      real*8 toler
      parameter(toler=1e-12)
      integer i
      logical verbose
      verbose=.false.
!     von Mises
      write(*,*) 'testing von Mises'
      nyldc=1
      call uten(0,yldc(:nyldc),nyldc,nth,rs,ys,verbose)
      do 5 i=1,nth
         if (abs(ys(i)-1d0)>toler .or. abs(rs(i)-1d0)>toler) then
            write(*,*) ys(i)
            write(*,*)'Something went wrong in VM'
         endif
 5    continue
!     Hill 48
      write(*,*) 'testing hill48'
      yldc(1)=0.5d0
      yldc(2)=0.5d0
      yldc(3)=0.5d0
      yldc(4)=1.5d0
      yldc(5)=1.5d0
      yldc(6)=1.5d0
      nyldc=6
      call uten(1,yldc(:nyldc),nyldc,nth,rs,ys,verbose) !! Hill48
      do 10 i=1,nth
         if (abs(ys(i)-1d0)>toler .or. abs(rs(i)-1d0)>toler) then
            write(*,*) ys(i)
            write(*,*)'Something went wrong in Hill48'
         endif
 10   continue
!     yld2000-2d isotropic
      write(*,*) 'testing yld2000-2d isotropic'
      yldc(1)=1.d0
      yldc(2)=1.d0
      yldc(3)=1.d0
      yldc(4)=1.d0
      yldc(5)=1.d0
      yldc(6)=1.d0
      yldc(7)=1.d0
      yldc(8)=1.d0
      yldc(9)=8.0000
      nyldc=9
      call uten(2,yldc(:nyldc),nyldc,nth,rs,ys,verbose)
      do 15 i=1,nth
         if (abs(ys(i)-1d0)>toler .or. abs(rs(i)-1d0)>toler) then
            write(*,*) ys(i),rs(i)
            write(*,*)'Something went wrong in yld2000-2d'
         endif
 15   continue
!     yld2000-2d anisotropic
      write(*,*) 'testing yld2000-2d anisotropic'
      yldc(1)=0.4865
      yldc(2)=1.3783
      yldc(3)=0.7536
      yldc(4)=1.0246
      yldc(5)=1.0363
      yldc(6)=0.9036
      yldc(7)=1.2321
      yldc(8)=1.4858
      yldc(9)=8.0000
      nyldc=9

!     yld2000-2d anisotropic uten
      verbose=.true.
      call uten(2,yldc(:nyldc),nyldc,nth,rs,ys,verbose) !! Yld2000-2d
      verbose=.false.
      call inplane_locus(0,yldc(:nyldc),nyldc,nth,locus,verbose)
      call inplane_locus(1,yldc(:nyldc),nyldc,nth,locus,verbose)
      call inplane_locus(2,yldc(:nyldc),nyldc,nth,locus,verbose)

      end program test_yld
c-----------------------------------------------------------------------
      subroutine uten(iyld,yldc,nyldc,nth,rs,ys,verbose)
      implicit none
      integer nth,i,j,iyld,nyldc
      dimension ys(nth),rs(nth)
      real*8 th,pi,s6lab(6),s6mat(6),s3mat(3),
     $     s33lab(3,3),s33mat(3,3),phim,dphim(6),d2phim(6,6),
     $     dphi33m(3,3),dphi33l(3,3),
     $     yldc(nyldc),aux33(3,3),
     $     aux3(3),rs,ys,rv
      logical verbose
      pi=4d0*datan(1d0)
c     Uniaxial tensin stress state referred in the lab axes
      s6lab(:)=0d0
      s6lab(1)=1d0
      call voigt2(s6lab,s33lab)
c-----------------------------------------------------------------------
      if (verbose)
     $     write(*,'(a5,5(4a6,x,a1,x),2a6)')'th',
     $     's11_l','s22_l','s33_l','s12_l','|',
     $     's11_m','s22_m','s33_m','s12_m','|',
     $     's11_l','s22_l','s33_l','s12_l','|',
     $     'e11_m','e22_m','e33_m','e12_m','|',
     $     'e11_l','e22_l','e33_l','e12_l','|',
     $     'R','Phi'
      do 10 j=1,nth
         th = pi/2.d0 - pi/2.d0/(nth-1) * (j-1)
         if (verbose) write(*,'(f5.1)',advance='no') th*180.d0/pi
c$$$
c$$$  si_lab
c$$$
         if (verbose) write(*,'(4f6.2,x,a1,x)',advance='no')
     $        (s33lab(i,i),i=1,3),s33lab(1,2),'|'
c$$$
c$$$  si_mat
c$$$
         call inplane_rot(th,s33lab,s33mat)
         if (verbose) write(*,'(4f6.2,x,a1,x)',advance='no')
     $        (s33mat(i,i),i=1,3),s33mat(1,2),'|'
c$$$
c$$$  si_lab
c$$$
         call inplane_rot(th*(-1.),s33mat,s33lab)
         if (verbose) write(*,'(4f6.2,x,a1,x)',advance='no')
     $        (s33lab(i,i),i=1,3),s33lab(1,2),'|'
c$$$
c$$$  ei_mat
c$$$
         call voigt1(s33mat,s6mat)
c        yield stress is written in the material axes
         if (iyld.eq.0) then
            call vm_gen(    s6mat,phim,dphim,d2phim)
         elseif (iyld.eq.1) then
            call hill48_gen(s6mat,phim,dphim,d2phim,yldc)
         elseif (iyld.eq.2) then
            call reduce_6to3(s6mat,s3mat)
            call yld2000_2d(s3mat,phim,aux3,aux33,yldc)
            call reduce_3to6(aux3,dphim)
         endif
!        shear strains: =1/2.shear
         call voigt4(dphim,dphi33m)
!        dphi in material axes
         if (verbose) write(*,'(4f6.2,x,a1,x)',advance='no')
     $        (dphi33m(i,i),i=1,3),dphi33m(1,2),'|'
c$$$
c$$$  ei_lab - strain in the lab space
c$$$
         call inplane_rot(th*(-1.d0),dphi33m,dphi33l)
         rv =-dphi33l(2,2)/(dphi33l(1,1)+dphi33l(2,2))
         rs(j)=rv
         if (verbose) write(*,'(4f6.2,x,a1,x,2f6.2)',advance='no')
     $        (dphi33l(i,i),i=1,3),dphi33l(1,2),'|',rv,phim
         ys(j)=phim
         if (verbose) write(*,*)
 10   continue
      return
      end subroutine
c-----------------------------------------------------------------------
      subroutine inplane_locus(iyld,yldc,nyldc,nth,locus,verbose)
      implicit none
      integer nth,j,iyld,nyldc
      dimension locus(nth,2)
      real*8 th,pi,s6mat(6),s3mat(3),
     $     phim,dphim(6),d2phim(6,6),
     $     yldc(nyldc),aux33(3,3),aux3(3),locus,q
      logical verbose

      pi=4d0*datan(1d0)

      if (iyld.eq.0) then
         open(unit=1,file='locus_vm.txt')
      elseif (iyld.eq.1) then
         open(unit=1,file='locus_hill48.txt')
      elseif (iyld.eq.2) then
         open(unit=1,file='locus_yld2000_2d.txt')
      endif

      do 10 j=1,nth
         th = 2*pi/(nth-1)*(j-1)
         s6mat(1)   = dcos(th)
         s6mat(2)   = dsin(th)
         s6mat(3:6) = 0d0

         if (iyld.eq.0) then
            q=2d0
            call vm_gen(    s6mat,phim,dphim,d2phim)
         elseif (iyld.eq.1) then
            q=2d0
            call hill48_gen(s6mat,phim,dphim,d2phim,yldc)
         elseif (iyld.eq.2) then
            q=yldc(9)
            call reduce_6to3(s6mat,s3mat)
            call yld2000_2d(s3mat,phim,aux3,aux33,yldc)
            call reduce_3to6(aux3,dphim)
         endif

         if (verbose) then
            write(*,'(f7.1 ,x,a2,x)',advance='no')
     $           th*180d0/pi,'|'
            write(*,'(3f7.2,x,a2,x)',advance='no')
     $           s6mat(1),s6mat(2),phim,'|'
         endif

         s6mat(:) = (s6mat(:)/phim)

         if (iyld.eq.0) then
            q=2d0
            call vm_gen(    s6mat,phim,dphim,d2phim)
         elseif (iyld.eq.1) then
            q=2d0
            call hill48_gen(s6mat,phim,dphim,d2phim,yldc)
         elseif (iyld.eq.2) then
            q=yldc(9)
            call reduce_6to3(s6mat,s3mat)
            call yld2000_2d(s3mat,phim,aux3,aux33,yldc)
            call reduce_3to6(aux3,dphim)
         endif

         if (verbose) then
            write(*,'(3f7.2)',advance='no') s6mat(1),s6mat(2),phim
            write(*,*)
         endif
         write(1,'(2f9.4)') s6mat(1),s6mat(2)
         locus(j,1:2) = s6mat(1:2)

 10   continue
      if (verbose) write(*,*)'-----------'

      return
      end subroutine inplane_locus
c-----------------------------------------------------------------------

c$$$!     pal
c      include "/home/younguj/repo/abaqusPy/umats/yld/vm.f"
c      include "/home/younguj/repo/abaqusPy/umats/yld/hill48.f"
c      include "/home/younguj/repo/abaqusPy/umats/yld/yld2000_2d.f"
c      include "/home/younguj/repo/abaqusPy/umats/lib/lib_write.f"
c      include "/home/younguj/repo/abaqusPy/umats/lib/is.f"
c      include "/home/younguj/repo/abaqusPy/umats/lib/lib.f"
c      include "/home/younguj/repo/abaqusPy/umats/lib/cnv.f"
c$$$!     mac
c$$$      include "/Users/yj/repo/abaqusPy/umats/yld/vm.f"
c$$$      include "/Users/yj/repo/abaqusPy/umats/yld/hill48.f"
c$$$      include "/Users/yj/repo/abaqusPy/umats/yld/yld2000_2d.f"
