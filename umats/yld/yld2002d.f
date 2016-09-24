c-----------------------------------------------------------------------
      subroutine yld2000_2d(cauchy,psi,dpsi,d2psi)
      implicit none
      integer ntens
      parameter(ntens=3)
      dimension cauchy(ntens),dpsi(ntens),d2psi(ntens,ntens),
     $     c(2,ntens,ntens),x1(ntens),x2(ntens),l1(ntens,ntens),
     $     l2(ntens,ntens),dh(2,ntens),phis(2),dphis(2,ntens)
      real*8 cauchy,dpsi,d2psi,psi,c,x1,x2,pp,ppp,h,l1,l2,a,dh,
     $     time0,time1,hershey1,hershey2,phis,dphis
      integer i

      dphis(:,:)=0d0

      call fill_line(0,'*',52)
      call get_idx(3,c(1,:,:))
      call get_idx(3,c(2,:,:))
      call w_chr(0,'c1')
      call w_mdim(0,c(1,:,:),3,1d0)
      call w_chr(0,'c2')
      call w_mdim(0,c(2,:,:),3,1d0)
      call calc_x(cauchy,c(1,:,:),x1) ! x1: linearly transformed stress (prime)
      call calc_x(cauchy,c(2,:,:),x2) ! x2: linearly transformed stress (double prime)

      a = 2d0                   ! yield surface exponent
      call hershey(x1,x2,a,phis(1),phis(2))
      psi = (0.5d0*(phis(1)+phis(2)))**(1d0/a)
      call w_chr(0,'phis:')
      call w_dim(0,phis,2,1d0,.true.)
      call w_val(0,'psi :',psi)

      call fill_line(0,'*',52)
      call calc_dphi_dcauchy(cauchy,c(1,:,:),a,dphis(1,:),0)
      call fill_line(0,'-',52)
      call calc_dphi_dcauchy(cauchy,c(2,:,:),a,dphis(2,:),1)
      call fill_line(0,'*',52)

      dpsi(:) = 0d0
      do i=1,ntens
         dpsi(i) = 1d0/2d0/a * ((phis(1)+phis(2))/2d0)**(1d0/a-1d0) *
     $        (dphis(1,i) + dphis(2,i))
      enddo
      call w_dim(0,dpsi,3,1d0,.true.)
      stop

      return
      end subroutine yld2000_2d
c-----------------------------------------------------------------------
      subroutine calc_dphi_dcauchy(cauchy,c,a,dphi,iopt)
      implicit none
      dimension cauchy(3),c(3,3),dphi(3),x(3),xp(2),dphi_dp(2),
     $     dp_dx(2,3),dphi_dx(3),l(3,3)
      real*8 cauchy,c,delta,x,dphi,xp,dp_dx,dphi_dp,a,dphi_dx,l
      integer iopt,i,j

      call calc_x(cauchy,c,x)
      call calc_delta(x,delta)
      call calc_princ(x,xp)
      dphi_dx(:) = 0d0
      if (delta.ne.0) then
         if (iopt.eq.0) call calc_a14(xp,a,   dphi_dp) ! eq A1.4
         if (iopt.eq.1) call calc_a18(xp,a,dphi_dp)    ! eq A1.8
         call calc_a15(x,delta,dp_dx)                  ! eq A1.9 = eq A1.5

         call w_chr(0,'dphi_dp')
         call w_dim(0,  dphi_dp,2,  1d0,.true.)
         call w_chr(0, 'dp_dx')
         call w_mndim(0,dphi_dp,2,3,1d0)
!        eq A1.3
            do 5 j=1,3
            do 5 i=1,2
               dphi_dx(j) = dphi_dx(j) + dphi_dp(i) * dp_dx(i,j)
 5          continue
      else
         write(*,*) 'delta is zero'
         stop
         if     (iopt.eq.0) then
            write(*,*)
!        eq A1.6
         elseif (iopt.eq.1) then
            write(*,*)
!        eq A1.10
         endif
      endif

      ! dphi_dx
      dphi(:)=0d0
      call calc_l(c,l)
      call w_chr(0,'L:')
      call w_mdim(0,l,3,1d0)
      do 10 j=1,3
      do 10 i=1,3
         dphi(j) = dphi(j) + dphi_dx(i) * l(i,j)
 10   continue

      call w_chr(0,'dphi')
      call w_dim(0,dphi,3,1d0,.true.)

      return
      end subroutine calc_dphi_dcauchy
c-----------------------------------------------------------------------
      subroutine calc_a15(x,delta,dp_dx)
      implicit none
      dimension x(3),dp_dx(2,3)
      real*8 x,delta,dp_dx,deno
      deno = dsqrt(delta)
      dp_dx(1,1) = 0.5d0 * ( 1d0 + (x(1)-x(2))/deno)
      dp_dx(1,2) = 0.5d0 * ( 1d0 - (x(1)-x(2))/deno)
      dp_dx(1,3) = 2.0d0 * x(3) / deno

      dp_dx(2,1) = 0.5d0 * ( 1d0 - (x(1)-x(2))/deno)
      dp_dx(2,2) = 0.5d0 * ( 1d0 + (x(1)-x(2))/deno)
      dp_dx(2,3) =-2.0d0 * x(3) / deno
      return
      end subroutine
c-----------------------------------------------------------------------
      subroutine calc_a14(xp,a,dphi_dp)
      implicit none
      dimension xp(2),dphi_dp(2)
      real*8 xp,dphi_dp,a
      dphi_dp(1) = a*(xp(1)-xp(2))**(a-1d0)
      dphi_dp(2) =-a*(xp(1)-xp(2))**(a-1d0)
      end subroutine calc_a14
c-----------------------------------------------------------------------
      subroutine calc_a18(xp,a,dphi_dp)
      implicit none
      dimension xp(2),dphi_dp(2)
      real*8 xp,dphi_dp,sgn1,sgn2,a
      sgn1=sign(1d0, 2d0*xp(2)+xp(1))
      sgn2=sign(1d0, 2d0*xp(1)+xp(2))

      dphi_dp(1) =       a * dabs(2d0*xp(2) + xp(1))**(a-1d0)*sgn1+
     $             2d0 * a * dabs(2d0*xp(1) + xp(2))**(a-1d0)*sgn2
      dphi_dp(2) = 2d0 * a * dabs(2d0*xp(2) + xp(1))**(a-1d0)*sgn1+
     $                   a * dabs(2d0*xp(1) + xp(2))**(a-1d0)*sgn2
      end subroutine calc_a18
c-----------------------------------------------------------------------
      subroutine calc_l(c,l)
c     Arguments
c     c: c matrix
c     l: l matrix
c     intent(in) c
c     intent(out) l
      implicit none
      dimension c(3,3),l(3,3),t(3,3)
      real*8 c,l,t

      t(:,:)= 0d0
      t(1,1)= 2d0/3d0
      t(1,2)=-1d0/3d0
      t(2,2)= 2d0/3d0
      t(2,1)=-1d0/3d0
      t(3,3)= 1d0
      call mult_array2(c,t,3,l)
      return
      end subroutine calc_l
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

      call calc_princ(x1,xp1)
      call calc_princ(x2,xp2)

      h1 = dabs(xp1(1) - xp1(2))**a
      h2 = dabs(2d0*xp2(2)+xp2(1))**a +dabs(2d0*xp2(1)+xp2(2))**a
      end subroutine hershey
c-----------------------------------------------------------------------
      subroutine calc_princ(s,xp)
      implicit none
      dimension xp(2),s(3)
      real*8 xx,yy,xy,f,x,s,xp
      xx=s(1)
      yy=s(2)
      xy=s(3)
      f  = dsqrt((xx-yy)**2 + 4d0 * xy**2)
      xp(1) = 0.5d0*(xx+yy+f)
      xp(2) = 0.5d0*(xx+yy-f)
      return
      end subroutine calc_princ
c-----------------------------------------------------------------------
c     Convert cauchy stress to linearly transformed X (eq 14 )
      subroutine calc_x(cauchy,c,x)
c     Arguments
c     cauchy : cauchy stress
c     c      : c matrix
c     x      : x tensor
      implicit none
      dimension cauchy(3),c(3,3),x(3),sdev(3)
      real*8 cauchy,c,x,sdev
      call cauchy2sdev(cauchy,sdev)
      call mult_array(c,sdev,3,x)
      return
      end subroutine calc_x
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

      call w_val(0,'phi:',phi)
      call w_chr(0,'dphi:')
      call w_dim(0,dphi,3,1d0,.true.)
      call w_chr(0,'d2phi:')
      call w_mdim(0,d2phi,3,1d0)

      call w_empty_lines(0,3)
      write(*,'(a,f7.2,a)')'Elapsed time:',
     $     (time1-time0)/1d-3,'mu seconds'

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






c-----------------------------------------------------------------------

c$$$c-----------------------------------------------------------------------
c$$$c     calculate dphi`/d(sigma) with phi` being the first
c$$$c     homogeneous function in eq (18)
c$$$      subroutine calc_dh2(x1,x2,l1,l2,a,dh)
c$$$c     Arugments
c$$$c     x1: linearly transformed stress using c1
c$$$c     c2: linearly transformed stress using c2
c$$$      implicit none
c$$$      dimension x1(3),x2(3),l1(3,3),l2(3,3),dh1(3),dh2(3),dh(3),dhl(2,3)
c$$$      real*8 x1,x2,l1,l2,dh1,dh2,a,dhl,dh
c$$$
c$$$      call calc_dhdxab_2(x1,a,dh1)
c$$$      call calc_dhdxab_2(x2,a,dh2)
c$$$
c$$$      call mult_array(dh1,l1,3,dhl(1,:))
c$$$      call mult_array(dh2,l2,3,dhl(2,:))
c$$$      dh(:) = dhl(1,:) + dhl(2,:)
c$$$      return
c$$$      end subroutine calc_dh2
c$$$c-----------------------------------------------------------------------
c$$$c     calculate dphi`/d(sigma) with phi` being the first
c$$$c     homogeneous function in eq (18)
c$$$      subroutine calc_dh1(x1,x2,l1,l2,a,dh)
c$$$c     Arguments
c$$$c     x1: linearly transformed stress using c1 (in)
c$$$c     x2: linearly transformed stress using c2 (in)
c$$$c     l1: l1 matrix (in)
c$$$c     l2: l2 matrix (in)
c$$$c     a : yield function exponent
c$$$c     dh: dphi`/dsigma with sigma being cauchy stress
c$$$      implicit none
c$$$      dimension x1(3),x2(3),l1(3,3),l2(3,3),dh1(3),dh2(3),dh(3),dhl(2,3)
c$$$      real*8 x1,x2,l1,l2,dh1,dh2,a,dhl,dh
c$$$
c$$$      call calc_dhdxab(x1,a,dh1) ! Calculate 1st term on RHS in eq (A1.3)
c$$$      call calc_dhdxab(x2,a,dh2) ! Calculate 2nd term on RHS in eq (A1.3)
c$$$      call mult_array(dh1,l1,3,dhl(1,:))
c$$$      call mult_array(dh2,l2,3,dhl(2,:))
c$$$      dh(:) = dhl(1,:) + dhl(2,:)
c$$$      return
c$$$      end subroutine calc_dh1
c$$$c-----------------------------------------------------------------------
c$$$c     eq (A1.3)
c$$$      subroutine calc_dhdxab(x,a,dhdxab)
c$$$c     Arguments
c$$$c     x: linearly transformed stress using a particular c matrix
c$$$c     a: yield surface exponent
c$$$c     dhdxab: (intent:out)
c$$$      implicit none
c$$$      dimension x(3),dhdx(2),dhdxab(3),dxdxij(2,3),xp(2)
c$$$      real*8 x,a,dhdx,dhdxab,dxdxij,delta,xp
c$$$      integer i
c$$$      call calc_delta(x,delta)  ! depending on delta
c$$$      if     (delta.ne.0) then
c$$$c        calc dprin/dx`_i in (A1.4) and principal stress
c$$$         call calc_dp1_dxij(x,dxdxij,delta,xp)
c$$$         call calc_dhdxi(xp,a,dhdx) ! calc dphi`/dxp (2)
c$$$         do 10 i=1,3
c$$$            dhdxab(i) = dhdx(1) * dxdxij(1,i) + dhdx(2) * dxdxij(2,i)
c$$$ 10      continue
c$$$      else ! singular case.
c$$$         dhdxab(:) = 0d0
c$$$      endif
c$$$      return
c$$$      end subroutine calc_dhdxab
c$$$c-----------------------------------------------------------------------
c$$$c     eq (A1.7)
c$$$      subroutine calc_dhdxab_2(x,a,dhdxab)
c$$$c     x:
c$$$      implicit none
c$$$      dimension x(3),dhdx(2),dhdxab(3),dxdxij(2,3),xp(2)
c$$$      real*8 x,a,dhdx,dhdxab,dxdxij,delta,xp
c$$$      integer i
c$$$      call calc_delta(x,delta) ! equivalent
c$$$      if     (delta.ne.0) then
c$$$         call calc_dp1_dxij_2(x,dxdxij,delta,xp)
c$$$         call calc_dhdxi_2(xp,a,dhdx)
c$$$         do 10 i=1,3
c$$$            dhdxab(i) = dhdx(1) * dxdxij(1,i) + dhdx(2) * dxdxij(2,i)
c$$$ 10      continue
c$$$      else ! singular case. (A1.10)
c$$$         call calc_dhdxi_2(x,a,dhdx)
c$$$         dhdxab(1) = dhdx(1)
c$$$         dhdxab(2) = dhdx(1)
c$$$         dhdxab(3) = 0d0
c$$$      endif
c$$$      return
c$$$      end subroutine calc_dhdxab_2
c$$$c-----------------------------------------------------------------------
c$$$c     eq (A1.5)
c$$$c     round(princi)/round(transformed stress)
c$$$      subroutine calc_dp1_dxij(xi,dpdxij,delta,xp)
c$$$c     xi: transformed stress tensor
c$$$c     dpdxij
c$$$c     xp: principal stress (2)
c$$$c     intent(in) xi,delta
c$$$c     intent(out) dpdxij,xp
c$$$
c$$$      implicit none
c$$$      dimension xi(3),dpdxij(2,3),xp(2)
c$$$      real*8 xi,delta,dpdxij,x1,x2,x3,xp
c$$$      call princ(xi(1),xi(2),xi(3),xp)
c$$$
c$$$      x1=xi(1)
c$$$      x2=xi(2)
c$$$      x3=xi(3)
c$$$
c$$$      dpdxij(1,1) = 0.5d0 * (1d0 + (x1 - x2)/dsqrt(delta))
c$$$      dpdxij(1,2) = 0.5d0 * (1d0 - (x1 - x2)/dsqrt(delta))
c$$$      dpdxij(1,3) = 2d0   *               x3/dsqrt(delta)
c$$$c     **
c$$$      dpdxij(2,1) = 0.5d0 * (1d0 - (x1 - x2)/dsqrt(delta))
c$$$      dpdxij(2,2) = 0.5d0 * (1d0 + (x1 - x2)/dsqrt(delta))
c$$$      dpdxij(2,3) =-2d0   *               x3/dsqrt(delta)
c$$$      return
c$$$      end subroutine calc_dp1_dxij
c$$$c-----------------------------------------------------------------------
c$$$c     eq (A1.9)=(A1.5)
c$$$      subroutine calc_dp1_dxij_2(x,dpdxij,delta,xp)
c$$$      implicit none
c$$$      dimension x(3),dpdxij(2,3),xp(2)
c$$$      real*8 x,delta,dpdxij,xp
c$$$      call calc_dp1_dxij(x,dpdxij,delta,xp)
c$$$      return
c$$$      end subroutine calc_dp1_dxij_2
c$$$c-----------------------------------------------------------------------
c$$$c     eq (A1.4)
c$$$      subroutine calc_dhdxi(xp,a,dhdx)
c$$$c     xp: principal stress (2)
c$$$c     a : exponent a
c$$$      implicit none
c$$$      dimension xp(2),dhdx(2)
c$$$      real*8 xp,a,x1,x2,dhdx
c$$$      x1=xp(1)
c$$$      x2=xp(2)
c$$$      dhdx(1) = a*(x1-x2)**(a-1d0)
c$$$      dhdx(2) =-a*(x1-x2)**(a-1d0)
c$$$      return
c$$$      end subroutine calc_dhdxi
c$$$c-----------------------------------------------------------------------
c$$$c     eq (A1.8)
c$$$      subroutine calc_dhdxi_2(xp,a,dhdx)
c$$$c     Arguments
c$$$c     xp: principal stress components
c$$$c     a : yield surface exponent
c$$$      implicit none
c$$$      dimension xp(2),dhdx(2)
c$$$      real*8 xp,a,x1,x2,dhdx
c$$$      x1=xp(1)
c$$$      x2=xp(2)
c$$$      dhdx(1) =    a*dabs(2d0*x2+x1)**(a-1d0)*sign(1d0,2d0*x2+x1)+
c$$$     $         2d0*a*dabs(2d0*x1+x2)**(a-1d0)*sign(1d0,2d0*x1+x2)
c$$$      dhdx(2) =2d0*a*dabs(2d0*x2+x1)**(a-1d0)*sign(1d0,2d0*x2+x1)+
c$$$     $             a*dabs(2d0*x1+x2)**(a-1d0)*sign(1d0,2d0*x1+x2)
c$$$      return
c$$$      end subroutine calc_dhdxi_2
c$$$c-----------------------------------------------------------------------

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
c$$$c     Convert deviatoric stress to x tensor
      subroutine cmat(c,sdev,x)
      implicit none
      dimension c(3,3),sdev(3),x(3)
      real*8 c, sdev,x
      call mult_array(c,sdev,3,x)
      return
      end subroutine cmat
