c-----------------------------------------------------------------------
c     Derivatives associated with Homogeneous Anisotropic Hardening
c     model

c     General references
c     [1] Lee et al., IJP 29, (2012) p 13-41
c     [2] Lee et al., Comput. Methods. Appl. Mech. Engrg. 247 (2012) p73-92
c     [3] Lee et al., Comput. Methods. Appl. Mech. Engrg. 286 (2015) p63-86
c     [4] Manual of abaqusPy

c-----------------------------------------------------------------------
      subroutine hah_deriv(nyldp,sdev,yldp,q,
     $     phi_h,  dphi_h,
     $     psi_sdp,dpsi_dsdp,
     $     psi_sp ,dpsi_dsp,
     $     phi_hah,dphi_hah)
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
c     phi_hah   : Homogeneous aniostropic yield surface
c     dphi_hah  : derivative of the Homogeneous aniostropic yield surface
      implicit none
      integer, intent(in) :: nyldp
      integer ntens
      parameter (ntens=6)
      dimension yldp(nyldp),dphi_h(6),dpsi_dsdp(6),
     $     dpsi_dsp(6),dphi_hah(6),sdev(6),dphi_hah_fin(6)
      real*8, intent(in) :: yldp,phi_h,psi_sdp,dpsi_dsdp,
     $     psi_sp,dpsi_dsp,phi_hah,dphi_hah,q
      real*8,intent(in):: sdev
      dimension dh_dsig(6,6),cel(6,6),deps_dsig(6),aux6(6),
     $     bux6(6),cux6(6),
     $     emic(6),dh_de(6),gk(4),e_ks(5),f_ks(2),
     $     ekrs(5),target(6),dsdp_ds(6,6),
     $     dphi_b(6),
     $     dfks_dsig(6),dhs_dsig(6),dps_ds(6,6),ds_dsig(6,6),gS(4)
      real*8 dh_dsig,deps_dsig,cel,aux6,bux6,cux6,
     $     dgR,gk,e_ks,f_ks,eeq,ref0,ref1,gL,ekL,eL,gS,c_ks,ss,
     $     ekrs,target,dsdp_ds,
     $     dphi_b,phi_b,dfks_dsig,dh_de,dhs_dsig,
     $     dphi_hah_fin,emic,dps_ds,dphi_h,ds_dsig
      integer i,j

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call emod_iso(Cel,0.3,200e9,3,3) !! to be improved later.
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Step 1
c     Calculate dh/deps <dh_de> (and dgr/deps) and saved them to yldp
      call micro_dev_deriv(6,nyldp,sdev,yldp) ! in <microd.f>
      call hah_io(0,nyldp,ntens,yldp,emic,dh_de,dgR,gk,e_ks,f_ks,eeq,
     $     ref0,ref1,gL,ekL,eL,gS,c_ks,ss,ekrs,target)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Step 2 - calculate dh_dsig and deps_dsig
      call calc_dh_dsig(Cel,dphi_hah,dh_de,dh_dsig,deps_dsig)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Step 3 - calculate d(h:s)/dsig
      call calc_dhs_dsig(sdev,emic,dh_dsig,dhs_dsig)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Steps 4,5,6
      call calc_stress_deriv(emic,dsdp_ds)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Step 7 - see Eq 26 in Ref [4]
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
c     Step 8
      call calc_dfk_dsig(emic,sdev,deps_dsig,gs,e_ks,q,
     $     ref0,ref1,dfks_dsig,phi_b)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Step 9
      call calc_dphib_dsig(emic,sdev,dhs_dsig,dfks_dsig,f_ks,dphi_b)
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     Step 10 in Eq 23
      dphi_hah_fin(:)=(phi_h**(q-1d0)*dphi_h(:)+
     $     phi_b**(q-1d0)*dphi_b(:))*phi_hah**(1d0-q)
      return
      end subroutine hah_deriv
c     Small operations below
c-----------------------------------------------------------------------
c     Equation 38 in Ref [4]
c     step 3
      subroutine calc_dhs_dsig(sdev,emic,dh_dsig,
     $     dhs_dsig)
c     Arguments
c     sdev      : Deviatoric stress tensor
c     emic      : Microstructure deviator
c     dh_dsig   : dh/dsig
c     dhs_dsig  : d(h:s)/dsig (to be obtained)
      implicit none
      dimension sdev(6),emic(6),dh_dsig(6,6),
     $     dhs_dsig(6)
      real*8, intent(in):: sdev,emic,dh_dsig
      real*8, intent(out):: dhs_dsig
c     local
      dimension ds_dsig(6,6),aux6(6)
      real*8 ds_dsig,aux6
      integer i,j
c     calc dh_dsig: s
      aux6(:)=0d0
      do 5 i=1,6
      do 5 j=1,3
         aux6(i) = aux6(i) + dh_dsig(i,j)  *sdev(j)
         aux6(i) = aux6(i) + dh_dsig(i,j+3)*sdev(j+3) *2d0
 5    continue
      call calc_ds_dsig(6,ds_dsig)
      dhs_dsig(:)=0d0
      do 10 i=1,6
      do 10 j=1,3
         dhs_dsig(i)=dhs_dsig(i)+aux6(i)+emic(j)  *ds_dsig(j,i)
         dhs_dsig(i)=dhs_dsig(i)+aux6(i)+emic(j+3)*ds_dsig(j+3,i)*2d0
 10   continue
      return
      end subroutine calc_dhs_dsig
c-----------------------------------------------------------------------
c     Eq 41. in Ref [4]
c     Step 2 - calculate dh_dsig
      subroutine calc_dh_dsig(Cel,dphi_dsig,dh_de,dh_dsig,deps_dsig)
c     Arguments
c     Cel       : elastic constants
c     dphi_dsig : Derivative of over all yield surface with respect to cauchy stress
c     dh_de     : dh_de
c     dh_dsig   : dh_dsig (to be obtained)
c     deps_dsig : deps_dsig (to be obtained)
      implicit none
      dimension cel(6,6),dphi_dsig(6),dh_dsig(6,6),deps_dsig(6),dh_de(6)
      real*8, intent(in):: cel,dphi_dsig,dh_de
      real*8, intent(out):: dh_dsig,deps_dsig
c     locals
      integer i,j
c     debar/dsig
      call calc_deps_dsig(Cel,dphi_dsig,deps_dsig)
c     dh_dsig = dh_de/de_dsig
c     tensor product
      do 20 j=1,6
      do 20 i=1,6
         dh_dsig(i,j) = dh_de(i)*deps_dsig(j)
 20   continue
      return
      end subroutine calc_dh_dsig
c-----------------------------------------------------------------------
c     deps/dsig: see Eq 40 in Ref [4] (its inverse actually)
      subroutine calc_deps_dsig(Cel,dphi_dsig,deps_dsig)
c     Arguments
c     Cel      : Elastic constants (6x6)
c     dphi_dsig: Derivative of over all yield surface
c     deps_dsig: to be calculated in this subroutine
      implicit none
      dimension Cel(6,6), dphi_dsig(6),deps_dsig(6)
c     locals
      dimension aux33(3,3),dsig_deps(6),aux33_inv(3,3)
      real*8, intent(in) :: Cel, dphi_dsig
      real*8, intent(out) ::  deps_dsig
c     locals
      real*8 aux33,dsig_deps,aux33_inv
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
c     Calculate {\partial s}/{\partial/sigma}
      subroutine calc_ds_dsig(ntens,ds_dsig)
      implicit none
      integer, intent(in) :: ntens
      dimension ds_dsig(ntens,ntens)
      real*8, intent(out) :: ds_dsig
      real*8 onethird
      parameter (onethird=1d0/3d0)
      integer i,j,ndi,nshr
cf2py intent(in) ntens
cf2py intent(out) ds_dsig
cf2py depends(ntens) ds_dsig
      if (ntens.eq.3) then
         ndi=2
         nshr=1
      elseif (ntens.eq.6) then
         ndi=3
         nshr=3
      else
         write(*,*)'ntens should be either 3 or 6 in',
     $        ' microd/calc_ds_dsig'
         call exit(-1)
      endif
      ds_dsig(:,:)=0d0
      do 1 i=1,ntens
         ds_dsig(i,i) = 1d0
 1    continue
      do 5 j=1,ndi
      do 5 i=1,ndi
         ds_dsig(i,j) = ds_dsig(i,j) -  onethird
 5    continue
      return
      end subroutine calc_ds_dsig
c-----------------------------------------------------------------------
c     Calculate ds``/ds, dsp/ds
      subroutine calc_stress_deriv(emic,dsdp_ds)
      dimension emic(6)
      dimension dsdp_ds(6,6),dsp_ds(6,6),dsc_ds(6,6),dso_ds(6,6)
      real*8, intent(in):: emic
      real*8, intent(out)::dsdp_ds
c     locals
      real*8 dsc_ds,eta,dso_ds,dsp_ds

      eta=1d0/gL
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     calc dsc_ds, equation 5 in Ref [4]
      dsc_ds(:,:)=0d0
      do 10 i=1,6
      do 10 j=1,6
         dsc_ds(i,j) = emic(i)*emic(j)
 10   continue
      dsc_ds(:,:)=dsc_ds(:,:)*8d0/3d0
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     calc dso/ds
      dso_ds(:,:) = dsc_ds(:,:)
      do 15 i=1,6
         dso_ds(i,i)=1d0-dsc_ds(i,i)
 15   continue
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     calculate ds``/ds  - see equation 9 in Ref [4]
      dsdp_ds(:,:)=0d0
      do 25 i=1,6
      do 25 j=1,6
         dsdp_ds(i,j) = dsc_ds(i,j) + eta*dso_ds(i,j)
 25   continue
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     calculate dsp_ds - see Eq 11 in Ref [4]
      do 35 i=1,6
      do 35 j=1,6
         dsp_ds(i,j) =  4d0*(1d0-gs)* dso_ds(i,j)
 35   continue
      return
      end subroutine calc_stress_deriv
c-----------------------------------------------------------------------
c     Calculate \partial f_k /\partial s and phi_b
c     Refer to Eqs 23 and 29 in Ref [4]
      subroutine calc_dfk_dsig(emic,sdev,deps_dsig,gs,e_ks,q,
     $     ys_iso,ys_hah,dfk_dsig,phi_b)
      implicit none
      dimension emic(6),sdev(6),deps_dsig(6),gs(4),e_ks(5),dfk_dsig(6)
      real*8, intent(in):: emic,e_ks,q,ys_iso,ys_hah,sdev,deps_dsig
      real*8, intent(out)::dfk_dsig,phi_b
c     locals
      integer ind
      dimension dfks_dgs(2),dgs_deps(2),gs_new(4),
     $     fks(2)
      real*8 dfks_dgs,dgs_deps,dot_prod,sh,gs_new,fks,gs

c     calculate s:h
      sh = dot_prod(6,3,3,emic,sdev)

c     obtain dfk/dgk and dgks/deps
c     gs_new and fks are dummy in this subroutine
      call calc_bau(6,3,3,emic,sdev,gs,e_ks,q,ys_iso,ys_hah,fks,gs_new,
     $     dfks_dgs,dgs_deps)

c     obtain dfk_dsig
      if (sign(1d0,sh).ge.0) then
         ind=1
      else
         ind=2
      endif
      dfk_dsig(:)=0d0
      dfk_dsig(:)=dfks_dgs(ind)*dgs_deps(ind)*deps_dsig(:)
      phi_b = fks(ind) * sh * 2d0
      return
      end subroutine calc_dfk_dsig
c-----------------------------------------------------------------------
c     Eq. 33 in Ref [4]
      subroutine calc_dphib_dsig(emic,sdev,dhs_dsig,dfk_dsig,f_ks,
     $     dphib_dsig)
c     Arguments
c     emic
c     sdev
c     dhs_dsig
c     dfk_dsig
c     f_ks
c     dphib_dsig
      dimension emic(6), sdev(6), dhs_dsig(6),dfk_dsig(6),
     $     dphib_dsig(6),f_ks(2)
      real*8, intent(in)::emic,sdev,dhs_dsig,dfk_dsig,f_ks
      real*8, intent(out)::dphib_dsig
c     locals
      real*8 sh, dot_prod
      integer ind
c     calculate s:h
      sh = dot_prod(6,3,3,emic,sdev)
c     obtain dfk_dsig
      if (sign(1d0,sh).ge.0) then
         ind=1
      else
         ind=2
      endif
c     Eq. 33
      dphib_dsig(:)=2.d0*sign(1d0,sh)*(dfk_dsig(:)*sh+f_ks(ind)*
     $     dhs_dsig(:))
      return
      end subroutine calc_dphib_dsig
c-----------------------------------------------------------------------
