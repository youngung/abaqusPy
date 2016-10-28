c-----------------------------------------------------------------------
c     Derivatives associated with Homogeneous Anisotropic Hardening
c     model

c     General references
c     [1] Lee et al., IJP 29, (2012) p 13-41
c     [2] Lee et al., Comput. Methods. Appl. Mech. Engrg. 247 (2012) p73-92
c     [3] Lee et al., Comput. Methods. Appl. Mech. Engrg. 286 (2015) p63-86
c     [4] Manual of abaqusPy

c-----------------------------------------------------------------------
      subroutine hah_deriv(nyldp,sdev,yldp,
     $     phi_h,  dphi_h,
     $     psi_sdp,dpsi_dsdp,
     $     psi_sp ,dpsi_dsp,
     $     phi_iso, phi_hah,dphi_hah)
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
c     dphi_hah  : derivative of the Homogeneous aniostropic yield surface
      implicit none
      integer, intent(in) :: ntens,ndi,nshr,nyldp
      dimension yldp(nyldp),dphi_h(6),dpsi_dsdp(6),dpsi_dsp(6),
     $     dphi_hah(6),sdev(6)
      real*8, intent(in) :: s,yldp,phi_h,dphi_h,psi_sdp,dpsi_dsdp,
     $     psi_sp,dpsi_dsp,phi_iso,phi_hah,dphi_hah
      real*8,intent(in):: sdev
      dimension dh_ds(6,6),cel(6,6),deps_dsig(6),dsig_deps(6),aux6(6),
     $     bux6(6),cux6(6),dhs_ds(6),idx(6,6),dsc_ds(6,6),dfk_ds(2,6),
     $     dsig_ds(6,6),emic(6),dh_de(6),dgR,gk(4),e_ks(5),f_ks(2),
     $     ekrs(5),target(6),dso_ds(6,6),dsdp_ds(6,6),dgb_de(4),
     $     dfks_dgk(2),dfk_ds(2),dphih_dsig(6),dgs_ds(6),aux33(3,3),
     $     bux33(3,3),dsp_ds(6,6),dphib_ds(6),dphib_dsig(6),dfks_ds(6)
      real*8 dh_ds,deps_dsig,aux6,bux6,cux6,dhs_ds,dsc_ds,sdh,shhat,
     $     dfk_ds,dgR,gk,e_ks,f_ks,eeq,ref0,ref1,gL,ekL,eL,gS,c_ks,ss,
     $     ekrs,target,dso_ds,dsdp_ds,dfks_dgk,dgb_de,dfk_ds,sgn,
     $     dphih_dsig,dgs_deps,dgs_ds,aux33,bux33,dsp_ds,dsig_ds,
     $     dphib_ds,dphib_dsig,dfks_ds
      integer idx,i,j, ind

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call emod_iso(Cel,0.3,200e9,3,3) !! to be improved later.
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Step 1
c     Calculate dh/deps <dh_de> (and dgr/deps) and saved them to yldp
      call micro_dev_deriv(6,3,3,nyldp,sdev,yldp) ! in <microd.f>
      call hah_io(0,nyldp,ntens,yldp,emic,dh_de,dgR,gk,e_ks,f_ks,eeq,
     $     ref0,ref1,gL,ekL,eL,gS,c_ks,ss,ekrs,target)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Step 2 - calculate dh_ds
      call calc_dh_ds(Cel,dphi_hah,dh_de,dh_ds,deps_dsig)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Step 3 - calculate d(h:s)/ds
      call calc_dhs_ds(dphi_hah,sdev,emic,dh_de,dh_ds,dhs_ds)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Steps 4,5,6,7
      call calc_stress_deriv(dhs_ds,dh_ds,deps_dsig,emic,sdev,eL,gL,ekL,
     $     ref0,ref1,c_kS,sS,gS,dsdp_ds,dsp_ds)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Step 8 - see Eq 30 in Ref [4]
      aux6(:)=0d0
      bux6(:)=0d0
      do 10 i=1,6
      do 10 j=1,3
         aux6(i) = aux6(i) + dpsi_dsdp(j)  *dsdp_ds(j,i)
         aux6(i) = aux6(i) + dpsi_dsdp(j+3)*dsdp_ds(j+3,i)*2d0
         bux6(i) = bux6(i) + dpsi_dsp(j)   *dps_ds(j,i)
         bux6(i) = bux6(i) + dpsi_dsp(j+3) *dps_ds(j+3,i)*2d0
 10   continue
      cux6(:)=psi_sdp * aux6(:) + psi_sp* bux6(:)
      dphi_h(:)=0d0
      do 15 i=1,6
      do 15 j=1,3
         dphi_h(i) = dphi_h(i) + cux6(j)  *ds_dsig(j,i)
         dphi_h(i) = dphi_h(i) + cux6(j+3)*ds_dsig(j+3,i)*2d0
 15   continue
      dphi_h(:)=dphi_h(:) / phi_h
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Step 9
      call calc_dfk_ds(emic,sdev,deps_dsig,dsig_ds,gs,e_ks,yldc(9),
     $     ref0,ref1,dfks_ds)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Step 10
      call calc_dphib_dsig(emic,sdev,dhs_ds,dfks_ds,f_ks,dphib_dsig)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end subroutine hah_deriv


c     Small operations below

c-----------------------------------------------------------------------
c     Equation 39 in Ref [4]
      subroutine calc_dhs_ds(dphi_dsig,sdev,emic,dh_de,dh_ds,dhs_ds)
c     Arguments
c     dphi_dsig : Derivative of over all yield surface with respect to cauchy stress
c     sdev      : Deviatoric stress tensor
c     emic      : Microstructure deviator
c     dh_de     : dh/de
c     dh_ds     : dh/ds
c     dhs_ds    : d(h:s)/ds (to be obtained)
c     dsig/ds
      implicit none
      dimension dphi_dsig(6),sdev(6),emic(6),dh_de(6),dh_ds(6,6),
     $     dhs_ds(6)
      real*8, intent(in):: dphi_dsig,sdev,emic,dh_de,dh_ds
      real*8, intent(out):: dhs_ds
c     local
      dimension dsig_ds(6,6),aux6(6)
      real*8 dsig_ds,aux6
      integer i,j
c     dsig/ds
      call calc_dsig_ds(6,dsig_ds)
      do 10 i=1,6
      do 10 j=1,3
         aux6(i) = aux6(i) + dh_ds(i,j)  *sdev(j)       + emic(i)
         aux6(i) = aux6(i) + dh_ds(i,j+3)*sdev(j+3)*2d0 + emic(i)
 10   continue
      end subroutine calc_dhs_ds
c-----------------------------------------------------------------------
c     Eq 41. in Ref [4]
      subroutine calc_dh_ds(Cel,dphi_dsig,dh_de,dh_ds,deps_dsig)
c     Arguments
c     Cel       : elastic constants
c     dphi_dsig : Derivative of over all yield surface with respect to cauchy stress
c     dh_de     : dh_de
c     dh_ds     : dh_ds (to be obtained)
      implicit none
      dimension cel(6,6),dphi_dsig(6),dh_ds(6,6),deps_dsig(6)
      real*8 intent(in):: cel,dphi_dsig
      real*8 intent(out):: dh_ds,deps_dsig
c     locals
      dimension aux(6),de_ds(6),dsig_ds(6,6)
      real*8 aux,de_ds,dsig_ds
      integer i,j,k,l

c     dsig/ds
      call calc_dsig_ds(6,dsig_ds)
c     debar/dsig
      call calc_deps_dsig(Cel,dphi_dsig,deps_dsig)
c     dh_de/de_dsig/dsig_ds

c     calc de_ds
      de_ds(:)=0d0
      do 10 i=1,6
      do 10 j=1,3
         de_ds(i) = de_ds(j)   + deps_dsig(j)  * dsig_ds(j,  i)
         de_ds(i) = de_ds(j+3) + deps_dsig(j+3)* dsig_ds(j+3,i) * 2d0
 10   continue
c     tensor product
      do 20 j=1,6
      do 20 i=1,6
         dh_ds(i,j) = dh_de(i)*de_ds(j)
 20   continue
      return
      end subroutine calc_dh_ds
c-----------------------------------------------------------------------
c     deps/dsig: see Eq 40 in Ref [4] (its inverse actually)
      subroutine calc_deps_dsig(Cel,dphi_dsig,deps_dsig)
c     Arguments
c     Cel      : Elastic constants (6x6)
c     dphi_dsig: Derivative of over all yield surface
c     deps_dsig: to be calculated in this subroutine
      implicit none
      dimension Cel(6,6), dphi_dsig(6),deps_dsig(6)
      real*8, intent(in) :: Cel, dphi_dsig
      real*8, intent(out) ::  deps_dsig
c     locals
      dimension aux33(3,3)
      real*8 aux33
      integer i,j
      deps_dsig(:)=0d0
      dsig_deps(:)=0d0
      do 1 i=1,6
      do 1 j=1,3
!        normal components
         dsig_deps(i)=dsig_deps(i)-Cel(i,j)  *dphi_dsig(j)
!        shear components
         dsig_deps(i)=dsig_deps(i)-Cel(i,j+3)*dphi_dsig(j+3)*2d0
 1    continue
      call voigt2(dsig_deps,aux33)
      aux33_inv(:,:)=aux33(:,:)
      call lu_inverse(aux33_inv,3)
      call voigt1(aux33_inv,deps_dsig)
      return
      end subroutine calc_deps_dsig
c-----------------------------------------------------------------------
      subroutine calc_dsig_ds(ntens,dsig_ds)
      implicit none
      integer, intent(in) :: ntens
      dimension dsig_ds(ntens,ntens)
      real*8, intent(out) :: dsig_ds
      real*8 onethird
      parameter (onethird=1d0/3d0)
      integer i,j,ndi,nshr
cf2py intent(in) ntens
cf2py intent(out) dsig_ds
cf2py depends(ntens) dsig_ds
      dsig_ds(:,:)=0d0
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
         dsig_ds(i,i) = 1d0
 1    continue
      do 5 j=1,ndi
      do 5 i=1,ndi
         dsig_ds(i,j) = dsig_ds(i,j) +  onethird
 5    continue
      return
      end subroutine calc_dsig_ds
c-----------------------------------------------------------------------
      subroutine calc_ds_dsig(ds_dsig)
      implicit none
      dimension ds_dsig(6,6)
      real*8, intent(out) :: ds_dsig

      end subroutine calc_ds_dsig
c-----------------------------------------------------------------------
c     Calculate ds``/ds, dsp/ds
      subroutine calc_stress_deriv(dhs_ds,dh_ds,deps_dsig,
     $     emic,sdev,eL,gL,ekL,
     $     ref0,ref1,c_kS,sS,gS,dsdp_ds,dsp_ds)
      dimension dhs_ds(6),dh_ds(6,6),emic(6),sdev(6),deps_dsig(6)
      dimension dsdp_ds(6,6), dsp_ds(6,6),dsc_ds(6,6),aux6(6),sc(6),
     $     so(6),dsig_ds(6,6),dgs_ds(6)
      real*8, intent(in):: dhs_ds,dh_ds,emic,sdev,deps_dsig,eL,gL,ekL,
     $     ref0,ref1,c_kS,sS,gS
      real*8, intent(out)::dsdp_ds,dsp_ds
c     locals
      real*8 dsc_ds,sh,tempval,sc,so,dsig_ds,dgL_deps,aux6,eta,eta2,
     $     dgs_deps,dgs_ds

      call calc_dsig_ds(6,dsig_ds)

      eta=1d0/gL
      eta2=eta/gL

c     calc h:s
      sh=0d0
      do 5 i=1,3
         sh = sh + emic(i)  *sdev(i)
         sh = sh + emic(i+3)*sdev(i+3)*2d0
 5    continue
c     calc dsc_ds, equation 5 in Ref [4]
      dsc_ds(:)=0d0
      do 10 i=1,6
      do 10 j=1,6
         dsc_ds(i,j) = dsc_ds(i,j) + dhs_ds(i) * emic(j) + sh*dh_ds(i,j)
 10   continue
      dsc_ds(:,:)=dsc_ds(:,:)*8d0/3d0
c     calc dso/ds
      dso_ds(:,:) = dsc_ds(:,:)
      do 15 i=1,6
         dso_ds(i,i)=1d0-dso_ds(i,i)
 15   continue
c     calculate so/sc
      call hah_decompose(6,3,3,sdev,emic,sc,so)
c     Using latent_update obtain dgL_deps
      call latent_update(6,emic,sdev,eL,Gl,ekL,ref0,ref1,dgL_deps)
c     calculate dsig/ds : so
      aux6(:)=0d0
      do 20 i=1,6
      do 20 j=1,3
         aux6(i) = aux6(i) + dsig_ds(i,j) * so(j)
         aux6(i) = aux6(i) + dsig_ds(i,j+3) * so(j+3)*2d0
 20   continue
c     calculate ds``/ds  - see equation 9 in Ref [4]
      dsdp_ds(:,:)=0d0
      do 25 i=1,6
      do 25 j=1,6
         dsdp_ds(i,j) = dsc_ds(i,j) + eta*dso_ds(i,j) -eta2*dgL_deps*
     $        deps_dsig(i) * aux6(j)
 25   continue
c     Calculate dgs_deps
      call cross_hardening(6,emic,sdev,c_ks,ss,gs,dgs_deps)
c     Calculate dgs_ds
      dgs_ds(:)=0d0
      do 30 i=1,6
      do 30 j=1,3
         dgs_ds(i)=dgs_ds(i)+deps_dsig(j)*dsig_ds(j,i)
         dgs_ds(i)=dgs_ds(i)+deps_dsig(j+3)*dsig_ds(j+3,i)*2d0
 30   continue
      dgs_ds(:) = dgs_ds(:)*dgs_deps
c     calculate dsp_ds - see Eq 12 in Ref [4]
      do 35 i=1,6
      do 35 j=1,6
         dps_ds(i,j) = -4d0 * dgs_ds(i) * so(j) + 4d0*(1d0-gs)*
     $        dso_ds(i,j)
 35   continue
      return
      end subroutine calc_stress_deriv
c-----------------------------------------------------------------------
c     Calculate \partial f_k /\partial s
c     Refer to eq 23 in Ref [4]
      subroutine calc_dfk_ds(emic,sdev,deps_dsig,dsig_ds,
     $     gs,e_ks,q,ys_iso,ys_hah,dfk_ds)
      implicit none
      dimension emic(6),sdev(6),deps_dsig(6,6),dsig_ds(6,6),
     $     gs(4),e_ks(5),dfk_ds(6)
      real*8,intent(in):: emic,e_ks,q,ys_iso,ys_hah,sdev,deps_dsig,
     $     dsig_ds
      real**, intent(out)::dfk_ds
c
c     locals
      integer i,j,ind
      dimension dfks_dgs(2), dgs_deps(2), dfk_ds(6), gs_new(4),fks(2)
      real*8 dfks_dgs, dgs_deps, dot_prod, sh,gs_new,fks

c     calculate s:h
      sh = dot_prod(6,3,3,emic,sdev)

c     obtain dfk/dgk and dgks/deps
c     gs_new and fks are dummy in this subroutine
      call calc_bau(6,3,3,emic,sdev,gs,e_ks,q,ys_iso,
     $     ys_hah,fks,gs_new,dfks_dgs,dgs_deps)

c     obtain dfk_ds
      if (sign(1d0,sh).ge.0) then
         ind=1
      else
         ind=2
      endif
      do 10 i=1,6
      do 10 j=1,3
         dfk_ds(i) = dfk_ds(i) + dfk_dgs(ind) * dgs_deps(ind) *
     $        deps_dsig(j)*dsig_ds(j,i)
         dfk_ds(i) = dfk_ds(i) + dfk_dgs(ind) * dgs_deps(ind) *
     $        deps_dsig(j+3)*dsig_ds(j+3,i) * 2d0
 10   continue
      end subroutine calc_dfk_ds
c-----------------------------------------------------------------------
      subroutine calc_dphib_dsig(emic,sdev,dhs_ds,dfk_ds,f_ks,
     $     dphib_dsig)
c     Arguments
c     emic
c     sdev
c     dhs_ds
c     dfk_ds
c     f_ks
c     dphib_ds
      dimension emic(6), sdev(6), dhs_ds(6),dfk_ds(6),dphib_ds(6),
     $     dphib_dsig(6),f_ks(2)
      real*8, intent(in)::emic,sdev,dhs_ds,dfk_ds,f_ks
      real*8, intent(out)::dphib_ds,dphib_dsig
c     locals
      real*8 sh, dot_prod
      integer i, ind

      call calc_dsig_ds(6,dsig_ds)

c     calculate s:h
      sh = dot_prod(6,3,3,emic,sdev)
c     obtain dfk_ds
      if (sign(1d0,sh).ge.0) then
         ind=1
      else
         ind=2
      endif
c     Eq 38.
      dphib_dsig(:) = 2.d0 *  sign(1d0,sh) * (dfk_ds(:)*sh +f_ks(ind)*
     $     dhs_ds(:))
c     Eq. 31.
      do 10 i=1,6
      do 10 j=1,3

 10   continue
      return
      end subroutine calc_dphib_dsig
