c-----------------------------------------------------------------------
c     Derivatives associated with Homogeneous Anisotropic Hardening
c     model

c     General references
c     [1] Lee et al., IJP 29, (2012) p 13-41
c     [2] Lee et al., Comput. Methods. Appl. Mech. Engrg. 247 (2012) p73-92
c     [3] Lee et al., Comput. Methods. Appl. Mech. Engrg. 286 (2015) p63-86
c     [4] Manual of abaqusPy

c-----------------------------------------------------------------------
      subroutine hah_deriv(ntens,nyldp,s,yldp,
     $     phi_h,  dphi_h,
     $     psi_sdp,dpsi_dsdp,
     $     psi_sp ,dpsi_dsp,
     $     phi_iso, phi_hah)
c     Arguments
c     ntens     : Len of stress tensor
c     nyldp     : Len of yldp
c     s         : deviatoric stress
c     yldp      : yield surface parameters (including state variables)
c     phi_h     : phi_h
c     dphi_h    : derivative of phi_h
c     psi_sdp   : psi_sdp
c     dpsi_dsdp :dpsi_dsdp
c     psi_sp    : psi_sp
c     dpsi_dsp   :dpsi_dsp
c     phi_iso   : isotropic yield surface
c     phi_hah   : Homogeneous aniostropic yield surface
      implicit none
      integer, intent(in) :: ntens,ndi,nshr,nyldp
      dimension s(ntens),yldp(nyldp),dphi_h(ntens),dpsi_dsdp(ntens),
     $     dpsi_dsp(ntens)
      real*8, intent(in) :: s,yldp,phi_h,dphi_h,psi_sdp,dpsi_dsdp,
     $     psi_sp,dpsi_dsp,phi_iso,phi_hah
      dimension sdev(6)
      real*8 sdev
      dimension dh_ds(6,6),cel(6,6),deps_dsig(6),aux6(6),bux6(6),
     $     dhs_ds(6),idx(6,6),dsc_ds(6,6),dfk_ds(2,6)
      real*8 dh_ds,deps_dsig,aux6,bux6,dhs_ds,dsc_ds,sdh,shhat,dfk_ds
      integer idx

      dimension dsig_ds(6,6)
      integer i,j, ind

      dimension emic(ntens),demic(ntens),dgR,gk(4),e_ks(5),f_ks(2),
     $     ekrs(5),target(ntens),dso_ds(6,6),dsdp_ds(6,6),dgb_de(4),
     $     dfks_dgk(2),dfk_ds(2),dphih_dsig(6)
      real*8 dgR,gk,e_ks,f_ks,eeq,ref,gL,ekL,eL,gS,c_ks,ss,ekrs,target,
     $     dso_ds,dsdp_ds,dfks_dgk,dgb_de,dfk_ds,sgn,dphih_dsig,dgs_deps,
     $     dgs_ds(6)

      dimension dphib_ds(6), dphib_dsig(6)
      real*8 dphib_ds,dphib_dsig

      if (ntens.eq.3) then
         sdev(1:2) = s(1:2)
         sdev(3)   = -s(1)-s(2)
         sdev(4:5) = 0d0
         sdev(6)   = s(3)
      elseif (ntens.eq.6) then
         sdev(:)=s(:)
      endif
      call emod_iso(Cel,0.3,200e9,3,3)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Calculate demic/deps and dgr/deps
      call micro_dev_deriv(6,3,3,sdev,nyldp,yldp) ! in <microd.f>
      call hah_io(0,nyldp,ntens,yldp,emic,demic,dgR,gk,e_ks,f_ks,eeq,
     $     ref,gL,ekL,eL,gS,c_ks,ss,ekrs,target)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     dsig/ds
      call calc_dsig_ds(6,dsig_ds)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     deps/dsig: see Eq 27 (its inverse actually)
      call calc_deps_dsig(ntens,dphi_red,Cel,deps_dsig)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     dh/ds: Eq. 28 in Ref [4]
      do 5 j=1,3
      do 5 i=1,3
         aux6(j)   = aux6(j)  +deps_dsig(i)  *dsig_ds(i,j)
         aux6(j+3) = aux6(j+3)+deps_dsig(i+3)*dsig_ds(i+3,j+3)*2d0
 5    continue
c     Outer product to obtain dh_ds
      do 10 i=1,6
      do 10 j=1,6
         dh_ds(i,j) = dh_de(i) * aux6(j)
 10   continue
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Eq 25 in Ref [4]
      idx(:,:)=0
      do 15 i=1,6
         idx(i,i)=1
 15   continue
      do 20 i=1,3
         dhs_ds(i)  =dhs_ds(i)+ dh_ds(i,  j)  *s(j)  +
     $        h(j)  *idx(j,i)
         dhs_ds(i+3)=dhs_ds(i)+ dh_ds(i+3,j+3)*s(j+3)+
     $        h(j+3)*idx(j+3,i+3)
         dhs_ds(i+3)=dhs_ds(i+3) * 2d0
 20   continue
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Eq 33 in Ref [4]
      sdh = 0d0
      do 25 i=1,3
         sdh = sdh + dhs_ds(i)  *emic(i)
         sdh = sdh + dhs_ds(i+3)*emic(i+3)*2d0
 25   continue
c     sshat
      sshat=0d0
      do 30 i=1,3
         sshat = sshat + s(i)   * emic(i)
         sshat = sshat + s(i+3) * emic(i+3) * 2d0
 30   continue
      do 35 i=1,6
      do 35 j=1,6
         dsc_ds(i,j) = 8d0/3d0 * (sdh+sshat*dh_ds(i,j))
 35   continue
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Eq 35 in Ref [4]
      do 40 i=1,6
      do 40 j=1,6
         dso_ds(i,j) = idx(i,j) - dsc_ds(i,j)
 40   continue
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Eq 37 in Ref [4]
      do 45 i=1,6
      do 45 j=1,6
         dsdp_ds(i,j) = dsc_ds(i,j) + dso_ds(i,j) / gL
 45   continue
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Eq 1 in Ref [4]
      call calc_bau(6,3,3,emic,target6,gs,e_ks,yldc(9),phi_iso,phi_hah,
     $     debar,fks,gs_new,dfks_dgk,dgb_de)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      do 50 j=1,3
      do 50 i=1,3
         dfk_ds(1,j) = dfks_dgk(1)   * dgb_de(1) * deps_dsig(i) *
     $        dsig_ds(i,j)
         dfk_ds(2,j) = dfks_dgk(2)   * dgb_de(2) * deps_dsig(i) *
     $        dsig_ds(i,j)
         dfk_ds(1,j+3) = dfks_dgk(1) * dgb_de(1) * deps_dsig(i) *
     $        dsig_ds(i+3,j+3) * 2d0
         dfk_ds(2,j+3) = dfks_dgk(2) * dgb_de(2) * deps_dsig(i) *
     $        dsig_ds(i+3,j+3) * 2d0
 50   continue
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Eq 21 (or 24) in Ref [4]
      sgn = sign(1d0,sdh)
      if (sgn.ge.0) ind=2
      if (sng.lt.0) ind=1
      do 55 i=1,6
         dphib_ds(i)  = 2d0*sgn*dfk_ds(ind,i)*sdh+fks(ind)*dhs_ds(i)
 55   continue
c     Eq 17 in Ref [4]
      do 60 i=1,6
      do 60 j=1,3
         dphib_dsig(i)=dphib_dsig(i)+dphib_ds(j)  *ds_dsig(j,  i)
         dphib_dsig(i)=dphib_dsig(i)+dphib_ds(j+3)*ds_dsig(j+3,i)*2d0
 60   continue
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Eq 16 in Ref [4]
      do 62 i=1,6
      do 62 j=1,6
      do 62 k=1,6
         dphih_dsig(i) = 1d0 / phi_h * (
     $        psi_sdp* dpsi_dsdp(j)*dsdp_ds(j,k) +
     $        psi_sp * dpsi_dsp( j)*dsp_ds(j,ik)) * ds_dsig(k,i)
 62   continue
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Eq 41 in Ref [4]
      call corss_hardening(6,3,3,emic,tensor_ref,ks,ss,gS,dgS_deps)
      dgs_ds(:)=0d0
      do 65 i=1,3
         dgs_ds(i)  =dgs_ds(i)  +dgs_deps*deps_dsig(j)*dsig_ds(j,i)
         dgs_ds(i+3)=dgs_ds(i+3)+dgs_deps*deps_dsig(j)*dsig_ds(j,i+3)*2
 65   continue
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Eq 40 in Ref [4]
      do 70 i=1,6
      do 70 j=1,6
         dsp_ds(i,j) = - 4d0 * dgs_ds(i) * so(j) +
     $        4d0* (1d0-gS) * dso_ds(i,j)
 70   continue
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c$$$c     Eq 16 in Ref [4]
c$$$      aux6(:)=0d0
c$$$      do 75 i=1,6
c$$$      do 75 j=1,6
c$$$         aux6(i) = 2d0 * psi_sdp * dpsi_dsdp(j) * dsdp_ds(j,i) +
c$$$     $             2d0 * psi_sp  * dpsi_dsp(j)  * dsp_ds(j,i)
c$$$ 75   continue
c$$$      do 80 i=1,6
c$$$      do 80 j=1,6
c$$$         dphih_dsig
c$$$ 80   continue
c$$$c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end subroutine hah_deriv
c-----------------------------------------------------------------------
c     calculate \frac{\partial\varepsilon){\partial\sigma}
      subroutine calc_deps_dsig(ntens,dphi_red,Cel,deps_dsig)
      integer, intent(in) :: ntens
      dimension Cel(6,6),dphi_red(ntens),dsig_deps(6),deps_dsig(6),
     $     aux33(3,3),aux33_inv(3,3)
      real*8, intent(in) :: Cel
      real*8, intent(out):: deps_dsig
      real*8 dphi_red,dsig_deps
      integer i,j

      dsig_deps(:) = 0d0
      do 10 i=1,3
      do 10 j=1,3
         dsig_deps(i)  =dsig_deps(i)  -Cel(i,  j)  *dphi(j)
         dsig_deps(i+3)=dsig_deps(i+3)-Cel(i+3,j+3)*dphi(j+3)*2d0
 10   continue
      call voigt2(dsig_deps,aux33)
      aux33_inv(:,:)=aux33(:,:)
      call lu_inverse(aux33_inv,3)
      call voigt1(aux33_inv,deps_dsig)
      return
      end subroutine calc_deps_dsig
c-----------------------------------------------------------------------
      subroutine calc_dsig_ds(ntens,d)
      implicit none
      integer, intent(in) :: ntens
      dimension d(ntens,ntens)
      real*8, intent(out) :: d
      real*8 onethird
      parameter (onethird=1d0/3d0)
      integer i,j,ndi,nshr
cf2py intent(in) ntens
cf2py intent(out) d
cf2py depends(ntens) d
      d(:,:)=0d0
      if (ntens.eq.3) then
         ndi=2
         nshr=1
      elseif (ntens.eq.6) then
         ndi=3
         nshr=3
      else
         write(*,*)'ntens should be either 3 or 6 in',
     $        ' microd/calc_dsig_ds'
         call exit(-1)
      endif
      do 1 i=1,ntens
         d(i,i) = 1d0
 1    continue
      do 5 j=1,ndi
      do 5 i=1,ndi
         d(i,j) = d(i,j) -  onethird
 5    continue
      return
      end subroutine calc_dsig_ds
c$$$      subroutine main_deriv(nyldp,ntens,yldp,sdev)
c$$$      implicit none
c$$$      integer, intent(in) nyldp,ntens
c$$$      dimension yldp(nyldp),sdev(ntens)
c$$$      real*8, intent(in) yldp,sdev
c$$$
c$$$      return
c$$$      end subroutine main_deriv
c$$$c-----------------------------------------------------------------------
c$$$c     Calculate derivatives between various stress components
c$$$c     dsc_ds, dso_ds,dsig_ds,dsp_dcauchy
c$$$      subroutine deriv_stress_components(iopt,emic,gS,f1,f2,
c$$$     $     dsc_ds,dso_ds,dsig_ds,dsp_dcauchy)
c$$$      implicit none
c$$$      integer, intent(in) :: iopt
c$$$      dimension emic(3,3),dsc_ds(3,3,3,3),dso_ds(3,3,3,3),
c$$$     $     dsig_ds(3,3,3,3),kro(3,3)
c$$$      real*8, intent(in) :: emic,gS,f1,f2
c$$$      real*8, intent(out) :: dsc_ds,dso_ds,dsig_ds,dsp_dcauchy
c$$$      real*8 temp
c$$$      integer kro,i,j,krok,onethird
c$$$      parameter(onethird=1d0/3d0)
c$$$      real H
c$$$      H=8d0/3d0
c$$$      kro(:,:)=0
c$$$c     Kronecker
c$$$      kro(1,1) = 1
c$$$      kro(2,2) = 1
c$$$      kro(3,3) = 1
c$$$
c$$$
c$$$c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c$$$c     Eq 48-a in Ref [3]
c$$$c     dsc_ds: \partial{s_c} / \partial{sdev}
c$$$      temp = H*0.5d0
c$$$      do 10 i=1,3
c$$$      do 10 j=1,3
c$$$      do 10 k=1,3
c$$$      do 10 l=1,3
c$$$         dsc_ds(i,j,k,l) = 0d0
c$$$      do 10 m=1,3
c$$$      do 10 n=1,3
c$$$         krok = (kro(m,k) * kro(n,l) + kro(m,l) * kro(n,k))
c$$$         dsc_ds(i,j,k,l) = dsc_ds(i,j,k,l) + temp * (
c$$$     $        krok * emic(m,n) * emic(i,j))
c$$$ 10   continue
c$$$c     Eq 48-b
c$$$      temp = 0.5d0
c$$$      do 20 i=1,3
c$$$      do 20 j=1,3
c$$$         dso_ds(i,j,:,:)=0d0
c$$$      do 20 k=1,3
c$$$      do 20 l=1,3
c$$$         dso_ds(i,j,k,l,) = dso_ds(i,j,k,l) +
c$$$     $        (kro(i,k)*kro(j,l)+kro(i,l)*kro(j,k))*temp-
c$$$     $        dsc_ds(i,j,k,l)
c$$$ 20   continue
c$$$c     Eq 49-a
c$$$      do 30 i=1,3
c$$$      do 30 j=1,3
c$$$         dsig_ds(i,j,:,:)=0d0
c$$$      do 30 k=1,3
c$$$      do 30 l=1,3
c$$$         dsig_ds(i,j,k,l)=dsig_ds(i,j,k,l) +
c$$$     $        kro(i,k)*kro(j,l)-onethird*kro(i,j)*kro(k,l)
c$$$ 30   continue
c$$$c     Eq 49-b
c$$$      temp = 4d0*(1d0-gS)
c$$$      do 40 i=1,3
c$$$      do 40 j=1,3
c$$$      do 40 k=1,3
c$$$      do 40 l=1,3
c$$$         dsp_dcauchy(i,j,k,l) = 0d0
c$$$      do 40 m=1,3
c$$$      do 40 n=1,3
c$$$         dsp_dcauchy(i,j,k,l) = dsp_dcauchy(i,j,k,l) +
c$$$     $        temp * dso_ds(i,j,m,n) * dsig_ds(m,n,k,l)
c$$$ 40   continue
c$$$c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c$$$c     Eq 47-a in Ref [3]
c$$$
c$$$      return
c$$$      end subroutine
