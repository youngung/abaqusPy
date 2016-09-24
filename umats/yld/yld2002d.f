c-----------------------------------------------------------------------
      subroutine yld2000_2d(cauchy,phi,dphi,d2phi)
      implicit none
      integer ntens
      parameter(ntens=3)
      dimension cauchy(ntens),dphi(ntens),d2phi(ntens,ntens),
     $     c(2,ntens,ntens),x1(ntens),x2(ntens),l1(ntens,ntens),
     $     l2(ntens,ntens),dh(2,ntens),hs(2)
      real*8 cauchy,dphi,d2phi,phi,c,x1,x2,pp,ppp,h,l1,l2,a,dh,
     $     time0,time1,hershey1,hershey2,hs

      call fill_line(0,'*',52)
      call get_idx(3,c(1,:,:))
      call get_idx(3,c(2,:,:))
      call c2x(cauchy,c(1,:,:),x1) ! x1: linearly transformed stress
      call c2x(cauchy,c(2,:,:),x2) ! x2: linearly transformed stress

      a = 1d0                   ! yield surface exponent
      call hershey(x1,x2,a,hs(1),hs(2))
      phi = (0.5d0*(hs(1)+hs(2)))**(1d0/a)

      write(*,*) 'phi:', phi


!     calculate l1/l2 matrix using c1,c2
      call calc_l(c(1,:,:),c(2,:,:),l1,l2)
      call calc_dh1(x1,x2,l1,l2,a,dh(1,:)) ! calculate dphi`/dsig in Section A1.1
      call calc_dh2(x1,x2,l1,l2,a,dh(2,:)) ! calculate dphi`/dsig in Section A1.1


      call w_chr(0,'x1')
      call w_dim(0,x1,3,1d0,.true.)
      call w_chr(0,'x2')
      call w_dim(0,x2,3,1d0,.true.)


      write(*,*)'h1:',hs(1)
      write(*,*)'h2:',hs(2)


      call w_chr(0,'dh1')
      call w_dim(0,dh(1,:),3,1d0,.true.)
      call w_chr(0,'dh2')
      call w_dim(0,dh(2,:),3,1d0,.true.)

c-----------------------------------------------------------------------
      dphi(:) = 1d0/2d0/a * (phi**(1d0-a)) * (dh(1,:)+dh(2,:))
c-----------------------------------------------------------------------
      return
      end subroutine yld2000_2d
c-----------------------------------------------------------------------
      subroutine calc_l(c1,c2,l1,l2)
c     Arguments
c     c1: c1 matrix
c     c2: c2 matrix
c     l1: l1 matrix
c     l2: l2 matrix
c     intent(in) c1,c2
c     intent(out) l1,l2
      implicit none
      dimension c1(3,3),c2(3,3),l1(3,3),l2(3,3),t(3,3)
      real*8 c1,c2,l1,l2,t

      t(:,:)= 0d0
      t(1,1)= 2d0/3d0
      t(1,2)=-1d0/3d0
      t(2,2)= 2d0/3d0
      t(2,1)=-1d0/3d0
      t(3,3)= 1d0

      call mult_array(c1,t,3,l1)
      call mult_array(c2,t,3,l2)

      return
      end subroutine calc_l
c-----------------------------------------------------------------------
c     calculate dphi`/d(sigma) with phi` being the first
c     homogeneous function in eq (18)
      subroutine calc_dh2(x1,x2,l1,l2,a,dh)
c     Arugments
c     x1: linearly transformed stress using c1
c     c2: linearly transformed stress using c2
      implicit none
      dimension x1(3),x2(3),l1(3,3),l2(3,3),dh1(3),dh2(3),dh(3),dhl(2,3)
      real*8 x1,x2,l1,l2,dh1,dh2,a,dhl,dh

      call calc_dhdxab_2(x1,a,dh1)
      call calc_dhdxab_2(x2,a,dh2)

      call mult_array(dh1,l1,3,dhl(1,:))
      call mult_array(dh2,l2,3,dhl(2,:))
      dh(:) = dhl(1,:) + dhl(2,:)
      return
      end subroutine calc_dh2
c-----------------------------------------------------------------------
c     calculate dphi`/d(sigma) with phi` being the first
c     homogeneous function in eq (18)
      subroutine calc_dh1(x1,x2,l1,l2,a,dh)
c     Arguments
c     x1: linearly transformed stress using c1 (in)
c     x2: linearly transformed stress using c2 (in)
c     l1: l1 matrix (in)
c     l2: l2 matrix (in)
c     a : yield function exponent
c     dh: dphi`/dsigma with sigma being cauchy stress
      implicit none
      dimension x1(3),x2(3),l1(3,3),l2(3,3),dh1(3),dh2(3),dh(3),dhl(2,3)
      real*8 x1,x2,l1,l2,dh1,dh2,a,dhl,dh

      call calc_dhdxab(x1,a,dh1) ! Calculate 1st term on RHS in eq (A1.3)
      call calc_dhdxab(x2,a,dh2) ! Calculate 2nd term on RHS in eq (A1.3) 
      call mult_array(dh1,l1,3,dhl(1,:))
      call mult_array(dh2,l2,3,dhl(2,:))
      dh(:) = dhl(1,:) + dhl(2,:)
      return
      end subroutine calc_dh1
c-----------------------------------------------------------------------
c     eq (A1.3)
      subroutine calc_dhdxab(x,a,dhdxab)
c     Arguments
c     x: linearly transformed stress using a particular c matrix
c     a: yield surface exponent
c     dhdxab: (intent:out)
      implicit none
      dimension x(3),dhdx(2),dhdxab(3),dxdxij(2,3),xp(2)
      real*8 x,a,dhdx,dhdxab,dxdxij,delta,xp
      integer i
      call calc_delta(x,delta)  ! depending on delta
      if     (delta.ne.0) then
c        calc dprin/dx`_i in (A1.4) and principal stress
         call calc_dp1_dxij(x,dxdxij,delta,xp)
         call calc_dhdxi(xp,a,dhdx) ! calc dphi`/dxp (2)
         do 10 i=1,3
            dhdxab(i) = dhdx(1) * dxdxij(1,i) + dhdx(2) * dxdxij(2,i)
 10      continue
      else ! singular case.
         dhdxab(:) = 0d0
      endif
      return
      end subroutine calc_dhdxab
c-----------------------------------------------------------------------
c     eq (A1.7)
      subroutine calc_dhdxab_2(x,a,dhdxab)
c     x: 
      implicit none
      dimension x(3),dhdx(2),dhdxab(3),dxdxij(2,3),xp(2)
      real*8 x,a,dhdx,dhdxab,dxdxij,delta,xp
      integer i
      call calc_delta(x,delta) ! equivalent
      if     (delta.ne.0) then
         call calc_dp1_dxij_2(x,dxdxij,delta,xp)
         call calc_dhdxi_2(xp,a,dhdx)
         do 10 i=1,3
            dhdxab(i) = dhdx(1) * dxdxij(1,i) + dhdx(2) * dxdxij(2,i)
 10      continue
      else ! singular case. (A1.10)
         call calc_dhdxi_2(x,a,dhdx)
         dhdxab(1) = dhdx(1)
         dhdxab(2) = dhdx(1)
         dhdxab(3) = 0d0
      endif
      return
      end subroutine calc_dhdxab_2
c-----------------------------------------------------------------------
c     eq (A1.5)
c     round(princi)/round(transformed stress)
      subroutine calc_dp1_dxij(xi,dpdxij,delta,xp)
c     xi: transformed stress tensor
c     dpdxij
c     xp: principal stress (2)
c     intent(in) xi,delta
c     intent(out) dpdxij,xp
      
      implicit none
      dimension xi(3),dpdxij(2,3),xp(2)
      real*8 xi,delta,dpdxij,x1,x2,x3,xp
      call princ(xi(1),xi(2),xi(3),xp)

      x1=xi(1)
      x2=xi(2)
      x3=xi(3)

      dpdxij(1,1) = 0.5d0 * (1d0 + (x1 - x2)/dsqrt(delta))
      dpdxij(1,2) = 0.5d0 * (1d0 - (x1 - x2)/dsqrt(delta))
      dpdxij(1,3) = 2d0   *               x3/dsqrt(delta)
c     **
      dpdxij(2,1) = 0.5d0 * (1d0 - (x1 - x2)/dsqrt(delta))
      dpdxij(2,2) = 0.5d0 * (1d0 + (x1 - x2)/dsqrt(delta))
      dpdxij(2,3) =-2d0   *               x3/dsqrt(delta)
      return
      end subroutine calc_dp1_dxij
c-----------------------------------------------------------------------
c     eq (A1.9)=(A1.5)
      subroutine calc_dp1_dxij_2(x,dpdxij,delta,xp)
      implicit none
      dimension x(3),dpdxij(2,3),xp(2)
      real*8 x,delta,dpdxij,xp
      call calc_dp1_dxij(x,dpdxij,delta,xp)
      return
      end subroutine calc_dp1_dxij_2
c-----------------------------------------------------------------------
c     eq (A1.4)
      subroutine calc_dhdxi(xp,a,dhdx)
c     xp: principal stress (2)
c     a : exponent a
      implicit none
      dimension xp(2),dhdx(2)
      real*8 xp,a,x1,x2,dhdx
      x1=xp(1)
      x2=xp(2)
      dhdx(1) = a*(x1-x2)**(a-1d0)
      dhdx(2) =-a*(x1-x2)**(a-1d0)
      return
      end subroutine calc_dhdxi
c-----------------------------------------------------------------------
c     eq (A1.8)
      subroutine calc_dhdxi_2(xp,a,dhdx)
c     Arguments
c     xp: principal stress components
c     a : yield surface exponent
      implicit none
      dimension xp(2),dhdx(2)
      real*8 xp,a,x1,x2,dhdx
      x1=xp(1)
      x2=xp(2)
      dhdx(1) =    a*dabs(2d0*x2+x1)**(a-1d0)*sign(1d0,2d0*x2+x1)+
     $         2d0*a*dabs(2d0*x1+x2)**(a-1d0)*sign(1d0,2d0*x1+x2)
      dhdx(2) =2d0*a*dabs(2d0*x2+x1)**(a-1d0)*sign(1d0,2d0*x2+x1)+
     $             a*dabs(2d0*x1+x2)**(a-1d0)*sign(1d0,2d0*x1+x2)
      return
      end subroutine calc_dhdxi_2
c-----------------------------------------------------------------------
c     eq (A1.2) delta used in principal stress calculation
      subroutine calc_delta(x,delta)
c     Arguments
c     x: linearly transformed stress using a particular c matrix (intent: in)
c     delta: to be calculated (intent: out)
      implicit none
      dimension x(3)
      real*8 x,delta
      integer i
      delta = (x(1)-x(2))**2+4d0*x(3)**2d0
      return
      end subroutine calc_delta
c-----------------------------------------------------------------------
      subroutine hershey(x1,x2,a,h1,h2)
c     Arguments
c     xp1, xp2: principal stres components of the two linearly
c              transformed stresses, respectively.
c     a: exponent
      implicit none
      dimension x1(3),x2(3),xp1(2),xp2(2)
      real*8 x1,x2,xp1,xp2,a,h1,h2

      call princ(x1(1),x1(2),x1(3),xp1)
      call princ(x2(1),x2(2),x2(3),xp2)

      h1 = dabs(xp1(1) - xp1(2))**a
      h2 = dabs(2d0*xp2(2)+xp2(1))**a +dabs(2d0*xp2(1)+xp2(2))**a
      end subroutine hershey
c-----------------------------------------------------------------------
      subroutine princ(xx,yy,xy,x)
      implicit none
      dimension x(2)
      real*8 xx,yy,xy,f,x
      f  = dsqrt((xx-yy)**2 + 4d0 * xy**2)
      x(1) = 0.5d0*(xx+yy+f)
      x(2) = 0.5d0*(xx+yy-f)
      return
      end subroutine princ
c-----------------------------------------------------------------------
c     Convert cauchy stress to X (eq 14 )
      subroutine c2x(cauchy,c,x)
c     Arguments
c     cauchy : cauchy stress
c     c      : c matrix
c     x      : x tensor
      implicit none
      dimension cauchy(3),c(3,3),x(3),sdev(3)
      real*8 cauchy,c,x,sdev
      call cauchy2sdev(cauchy,sdev)
      call cmat(c,sdev,x)
      return
      end subroutine c2x
c-----------------------------------------------------------------------
      subroutine cauchy2sdev(c,s)
      implicit none
c     Arguments
c     c : c matrix          (in)
c     s : deviatoric stress (out)
      dimension t(3,3),c(3),s(3)
      real*8 t,c,s
      t(:,:)= 0d0
      t(1,1)= 2d0/3d0
      t(1,2)=-1d0/3d0
      t(2,2)= 2d0/3d0
      t(2,1)=-1d0/3d0
      t(3,3)= 1d0
      call mult_array(t,c,3,s)
      return
      end subroutine cauchy2sdev
c-----------------------------------------------------------------------
c     Convert deviatoric stress to x tensor
      subroutine cmat(c,sdev,x)
      implicit none
      dimension c(3,3),sdev(3),x(3)
      real*8 c, sdev,x
      call mult_array(c,sdev,3,x)
      return
      end subroutine cmat
c-----------------------------------------------------------------------
      program test
      implicit none
      dimension cauchy(3),dphi(3),d2phi(3,3),c(2,3,3),
     $     x1(3),x2(3)
      real*8 cauchy,dphi,d2phi,phi,c,x1,x2,pp,ppp,
     $     hershey1,hershey2,a,time0,time1

      cauchy(:)=0d0
      cauchy(1)=1d0
      dphi(:)=0d0
      d2phi(:,:)=0d0
      call w_chr(0,'cauchy stress')
      call w_dim(0,cauchy,3,1d0,.true.)

      call cpu_time(time0)
      call yld2000_2d(cauchy,phi,dphi,d2phi)
      call cpu_time(time1)
      write(*,'(a,f7.2,a)')'Elapsed time:',
     $     (time1-time0)/1d-3,'mu seconds'

      call w_val(0,'phi:',phi)
      call w_chr(0,'dphi:')
      call w_dim(0,dphi,3,1d0,.true.)
      call w_chr(0,'d2phi:')
      call w_mdim(0,d2phi,3,1d0)


      end program
c-----------------------------------------------------------------------
c     Palmetto
!     include '/home/younguj/repo/abaqusPy/umats/lib/algb.f'
c     include '/home/younguj/repo/abaqusPy/umats/lib/lib_write.f'
c     include '/home/younguj/repo/abaqusPy/umats/lib/lib.f'
c     include '/home/younguj/repo/abaqusPy/umats/lib/is.f'

c     Mac
      include '/Users/yj/repo/abaqusPy/umats/lib/algb.f'
      include '/Users/yj/repo/abaqusPy/umats/lib/lib_write.f'
      include '/Users/yj/repo/abaqusPy/umats/lib/lib.f'
      include '/Users/yj/repo/abaqusPy/umats/lib/is.f'
