c-----------------------------------------------------------------------
c     Derivatives associated with Homogeneous Anisotropic Hardening
c     model

c     General references
c     [1] Manual of abaqusPy

c-----------------------------------------------------------------------
      subroutine hah_deriv(nyldp,cauchy,yldp,q,phi_h,dphi_h,psi_spp,
     $     dpsi_dspp,psi_sp,dpsi_dsp,phi_hah,dphi_hah,dphi_hah_fin)
c     Arguments
c     nyldp     : Len of yldp
c     cauchy    : cauchy stress
c     yldp      : yield surface parameters (including state variables)
c     q         : yield surface exponent
c     phi_h     : phi_h
c     dphi_h    : derivative of phi_h, i.e., dphi_h/dsig
c     psi_spp   :  psi_spp,  i.e., psi(s``)
c     dpsi_dspp : dpsi_dspp, i.e., dpsi(s``)/ds``
c     psi_sp    : psi_sp, i.e., psi(sp)
c     dpsi_dsp  : dpsi_dsp, i.e., dpsi(sp)/dsp
c     phi_hah   : Homogeneous aniostropic yield surface
c     dphi_hah  : derivative of the Homogeneous aniostropic yield surface
      implicit none
      integer ntens,ndi,nshr
      parameter (ntens=6,ndi=3,nshr=3)
      integer, intent(in) :: nyldp
      dimension cauchy(ntens),yldp(nyldp),dphi_h(ntens),
     $     dpsi_dspp(ntens),dpsi_dsp(ntens),dphi_hah(ntens),
     $     dphi_hah_fin(ntens)
      real*8, intent(in) :: cauchy,yldp,q,phi_h,dphi_h,psi_spp,
     $     dpsi_dspp,psi_sp,dpsi_dsp,phi_hah,dphi_hah
      real*8, intent(out) :: dphi_hah_fin

c     locals
      dimension cel(ntens,ntens),dphi_dsig(ntens),dsp_dsig(ntens),
     $     dspp_dsig(ntens),sdev(ntens)
      real*8  cel,dphi_dsig,dsp_dsig,dspp_dsig,sdev,hydro

c     ! calculate deviatoric stress
      call deviat(ntens,cauchy,sdev,hydro)
      call emod_iso(200e9,0.3,Cel,ndi,nshr) !! to be improved later.
      call calc_sdev_deriv(ntens,ndi,nshr,nyldp,yldp,sdev,
     $     dphi_dsig,dsp_dsig,dspp_dsig)
      return
      end subroutine hah_deriv

c$$$c     Calculate the stress derivatives Eq 14 in Ref [1]
c$$$      subroutine calc_sdev_deriv(ntens,nyldp,yldp,sdev)
c$$$      end subroutine
c-----------------------------------------------------------------------
c     Calculate d(sp)/d(Cauchy)
c     Calculate d(spp)/d(Cauchy)

c     Eq 14 in Ref [1]
      subroutine calc_sdev_deriv(ntens,ndi,nshr,nyldp,yldp,sdev,
     $     dphi_dsig,dsp_dsig,dspp_dsig)
c     Arguments
c     ntens
c     ndi
c     nshr
c     nyldp
c     yldp
c     sdev
c     dphi_dsig
      integer, intent(in) :: ntens,ndi,nshr,nyldp
      dimension yldp(nyldp),sdev(ntens),dphi_dsig(ntens)
      real*8,  intent(in) :: yldp,sdev,dphi_dsig
      dimension dsh_ds(ntens), ds_dsig(ntens,ntens)
      real*8 dsh_ds
      integer i
      dimension hhat(ntens),dhhat_de(ntens),gk(4),e_ks(5),f_ks(2),
     $     ekrs(5),target(ntens),hds_dsig(ntens),deps_dsig(ntens)
      real*8 hhat,dhhat_de,gk,e_ks,f_ks,ekrs,target,hds_dsig,eeq,dgR,
     $     ref0,ref1,gL,ekL,eL,gS,c_ks,ss,deps_dsig
      dimension dsc_dsig(ntens,ntens),dso_dsig(ntens,ntens),so(ntens),
     $     sc(ntens),dsp_dsig(ntens,ntens),dspp_dsig(ntens,ntens),
     $     Cel(ntens,ntens)
      real*8 ds_dsig,dsc_dsig,dso_dsig,so,sc,dgL_deps,dsp_dsig,dspp_dsig
     $     ,Cel

c     Restore from yldp
      call hah_io(0,nyldp,ntens,yldp,hhat,dhhat_de,dgR,gk,e_ks,
     $     f_ks,eeq,ref0,ref1,gL,ekL,eL,gS,c_ks,ss,ekrs,target)
c     Calculate ds/dsig
      call calc_ds_dsig(ntens,ds_dsig)
c     Calculate s:h
      call inner_dot_voigt(ntens,ndi,nshr,sdev,hhat,sh)
c     Calculate dsh_ds
      do 10 i=1,ndi
         dsh_ds(i)=hhat(i)
 10   continue
      do 20 i=ndi+1,ntens
         dsh_ds(i)=hhat(i) * 2d0
 20   continue
c     Calculate dh/de and dgr/de and save them to yldp
      call micro_dev_deriv(ntens,nyldp,sdev,yldp)
c     Restore dhhat_de from yldp
      call hah_io(0,nyldp,ntens,yldp,hhat,dhhat_de,dgR,gk,e_ks,
     $     f_ks,eeq,ref0,ref1,gL,ekL,eL,gS,c_ks,ss,ekrs,target)
c     Calculate deps_dsig
c     ---
      call calc_deps_dsig(Cel,dphi_dsig,deps_dsig)
c     calculate hhat_k ds_k/dsig_j -> hds_dsig

      hds_dsig(:)=0d0
      do 40 i=1,ntens
      do 30 j=1,ndir
         hds_dsig(i) = hds_dsig(i) + hhat(j) * ds_dsig(j,i)
 30   continue
      do 40 j=ndir+1,ntens
         hds_dsig(i) = hds_dsig(i) + hhat(j) * ds_dsig(j,i) * 2d0
 40   continue

c     Calculate dsc_dsig
      dsc_dsig(:,:)=0d0
      do 50 i=1,ntens
      do 50 j=1,ntens
         dsc_dsig(i,j)=dsh_ds(i)*hds_dsig(j)+sh*dhhat_de(i)*deps_dsig(j)
         dso_dsig(i,j)=ds_dsig(i,j)-dsc_dsig(i,j)
 50   continue

      call cross_hardening(ntens,hhat,target,c_ks,ss,gS,dgS_de)
      call hah_decompose(ntens,ndi,nshr,sdev,hhat,sc,so)

      do 60 i=1,ntens
      do 60 j=1,ntens
         dsp_dsig(i,j) = -4d0 * so(i) * dgS_de * deps_dsig(j) +
     $        4d0 * (1d0 - gS) * dso_dsig(i,j)
 60   continue

      call latent_update(ntens,nyldp,yldp,dgL_deps)

c     Calculate ds``/dsig
      do 70 i=1,ntens
      do 70 j=1,ntens
         dspp_dsig(i,j) = dsc_dsig(i,j)-
     $        (gL**(-2d0))*dgL_deps*deps_dsig(i)*so(j)+
     $        (gL**(-1d0))*dso_dsig(i,j)
 70   continue

      return
      end subroutine calc_sdev_deriv
c-----------------------------------------------------------------------
      subroutine calc_dsh_ds(ntens,ndi,nshr,sdev,hhat,dsh_ds)
      integer, intent(in)::ntens,ndi,nshr
      dimension sdev(ntens), hhat(ntens)
      real*8, intent(in):: sdev,hhat
      real*8, intent(out):: dsh_ds
      real*8 sh
      integer i,j

      if (ntens.ne.6) then
         write(*,*) 'Subr calc_dsh_ds needs ntens.eq.6'
         call exit(-1)
      endif

      return
      end subroutine calc_dsh_ds


c$$$      subroutine hah_deriv(nyldp,sdev,yldp,q,
c$$$     $     phi_h,  dphi_h,
c$$$     $     psi_sdp,dpsi_dsdp,
c$$$     $     psi_sp ,dpsi_dsp,
c$$$     $     phi_hah,dphi_hah,dphi_hah_fin)
c$$$c     Arguments
c$$$c     ntens     : Len of stress tensor
c$$$c     nyldp     : Len of yldp
c$$$c     s         : deviatoric stress
c$$$c     yldp      : yield surface parameters (including state variables)
c$$$c     phi_h     : phi_h
c$$$c     dphi_h    : derivative of phi_h
c$$$c     psi_sdp   : psi_sdp
c$$$c     dpsi_dsdp :dpsi_dsdp
c$$$c     psi_sp    : psi_sp
c$$$c     dpsi_dsp   :dpsi_dsp
c$$$c     phi_hah   : Homogeneous aniostropic yield surface
c$$$c     dphi_hah  : derivative of the Homogeneous aniostropic yield surface
c$$$      implicit none
c$$$      integer, intent(in) :: nyldp
c$$$      integer ntens
c$$$      parameter (ntens=6)
c$$$      dimension yldp(nyldp),dphi_h(6),dpsi_dsdp(6),
c$$$     $     dpsi_dsp(6),dphi_hah(6),sdev(6),dphi_hah_fin(6)
c$$$      real*8, intent(in) :: yldp,phi_h,psi_sdp,dpsi_dsdp,
c$$$     $     psi_sp,dpsi_dsp,phi_hah,dphi_hah,q
c$$$      real*8,intent(in):: sdev
c$$$      dimension dh_dsig(6,6),cel(6,6),deps_dsig(6),aux6(6),
c$$$     $     bux6(6),cux6(6),
c$$$     $     emic(6),dh_de(6),gk(4),e_ks(5),f_ks(2),
c$$$     $     ekrs(5),target(6),dsdp_ds(6,6),
c$$$     $     dphi_b(6),
c$$$     $     dfks_dsig(6),dhs_dsig(6),dps_ds(6,6),ds_dsig(6,6),gS(4)
c$$$      real*8 dh_dsig,deps_dsig,cel,aux6,bux6,cux6,
c$$$     $     dgR,gk,e_ks,f_ks,eeq,ref0,ref1,gL,ekL,eL,gS,c_ks,ss,
c$$$     $     ekrs,target,dsdp_ds,
c$$$     $     dphi_b,phi_b,dfks_dsig,dh_de,dhs_dsig,
c$$$     $     dphi_hah_fin,emic,dps_ds,dphi_h,ds_dsig
c$$$      integer i,j
c$$$
c$$$c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c$$$      call emod_iso(Cel,0.3,200e9,3,3) !! to be improved later.
c$$$c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c$$$c     Step 1
c$$$c     Calculate dh/deps <dh_de> (and dgr/deps) and saved them to yldp
c$$$      call micro_dev_deriv(6,nyldp,sdev,yldp) ! in <microd.f>
c$$$      call hah_io(0,nyldp,ntens,yldp,emic,dh_de,dgR,gk,e_ks,f_ks,eeq,
c$$$     $     ref0,ref1,gL,ekL,eL,gS,c_ks,ss,ekrs,target)
c$$$c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c$$$c     Step 2 - calculate dh_dsig and deps_dsig
c$$$      call calc_dh_dsig(Cel,dphi_hah,dh_de,dh_dsig,deps_dsig)
c$$$c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c$$$c     Step 3 - calculate d(h:s)/dsig
c$$$      call calc_dhs_dsig(sdev,emic,dh_dsig,dhs_dsig)
c$$$c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c$$$c     Steps 4,5,6
c$$$      call calc_stress_deriv(emic,dsdp_ds)
c$$$c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c$$$c     Step 7 - see Eq 26 in Ref [1]
c$$$      aux6(:)=0d0
c$$$      bux6(:)=0d0
c$$$      do 10 i=1,6
c$$$      do 10 j=1,3
c$$$         aux6(i) = aux6(i) + dpsi_dsdp(j)  *dsdp_ds(j,i)
c$$$         aux6(i) = aux6(i) + dpsi_dsdp(j+3)*dsdp_ds(j+3,i)*2d0
c$$$         bux6(i) = bux6(i) + dpsi_dsp(j)   *dps_ds(j,i)
c$$$         bux6(i) = bux6(i) + dpsi_dsp(j+3) *dps_ds(j+3,i)*2d0
c$$$ 10   continue
c$$$      cux6(:)=psi_sdp * aux6(:) + psi_sp* bux6(:)
c$$$      dphi_h(:)=0d0
c$$$      do 15 i=1,6
c$$$      do 15 j=1,3
c$$$         dphi_h(i) = dphi_h(i) + cux6(j)  *ds_dsig(j,i)
c$$$         dphi_h(i) = dphi_h(i) + cux6(j+3)*ds_dsig(j+3,i)*2d0
c$$$ 15   continue
c$$$      dphi_h(:)=dphi_h(:) / phi_h
c$$$c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c$$$c     Step 8
c$$$      call calc_dfk_dsig(emic,sdev,deps_dsig,gs,e_ks,q,
c$$$     $     ref0,ref1,dfks_dsig,phi_b)
c$$$c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c$$$c     Step 9
c$$$      call calc_dphib_dsig(emic,sdev,dhs_dsig,dfks_dsig,f_ks,dphi_b)
c$$$c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c$$$c     Step 10 in Eq 23
c$$$      dphi_hah_fin(:)=(phi_h**(q-1d0)*dphi_h(:)+
c$$$     $     phi_b**(q-1d0)*dphi_b(:))*phi_hah**(1d0-q)
c$$$      return
c$$$      end subroutine hah_deriv
c$$$c     Small operations below
c$$$c-----------------------------------------------------------------------
c$$$c     Equation 38 in Ref [1]
c$$$c     step 3
c$$$      subroutine calc_dhs_dsig(sdev,emic,dh_dsig,
c$$$     $     dhs_dsig)
c$$$c     Arguments
c$$$c     sdev      : Deviatoric stress tensor
c$$$c     emic      : Microstructure deviator
c$$$c     dh_dsig   : dh/dsig
c$$$c     dhs_dsig  : d(h:s)/dsig (to be obtained)
c$$$      implicit none
c$$$      dimension sdev(6),emic(6),dh_dsig(6,6),
c$$$     $     dhs_dsig(6)
c$$$      real*8, intent(in):: sdev,emic,dh_dsig
c$$$      real*8, intent(out):: dhs_dsig
c$$$c     local
c$$$      dimension ds_dsig(6,6),aux6(6)
c$$$      real*8 ds_dsig,aux6
c$$$      integer i,j
c$$$c     calc dh_dsig: s
c$$$      aux6(:)=0d0
c$$$      do 5 i=1,6
c$$$      do 5 j=1,3
c$$$         aux6(i) = aux6(i) + dh_dsig(i,j)  *sdev(j)
c$$$         aux6(i) = aux6(i) + dh_dsig(i,j+3)*sdev(j+3) *2d0
c$$$ 5    continue
c$$$      call calc_ds_dsig(6,ds_dsig)
c$$$      dhs_dsig(:)=0d0
c$$$      do 10 i=1,6
c$$$      do 10 j=1,3
c$$$         dhs_dsig(i)=dhs_dsig(i)+aux6(i)+emic(j)  *ds_dsig(j,i)
c$$$         dhs_dsig(i)=dhs_dsig(i)+aux6(i)+emic(j+3)*ds_dsig(j+3,i)*2d0
c$$$ 10   continue
c$$$      return
c$$$      end subroutine calc_dhs_dsig
c$$$c-----------------------------------------------------------------------
c$$$c     Eq 41. in Ref [1]
c$$$c     Step 2 - calculate dh_dsig
c$$$      subroutine calc_dh_dsig(Cel,dphi_dsig,dh_de,dh_dsig,deps_dsig)
c$$$c     Arguments
c$$$c     Cel       : elastic constants
c$$$c     dphi_dsig : Derivative of over all yield surface with respect to cauchy stress
c$$$c     dh_de     : dh_de
c$$$c     dh_dsig   : dh_dsig (to be obtained)
c$$$c     deps_dsig : deps_dsig (to be obtained)
c$$$      implicit none
c$$$      dimension cel(6,6),dphi_dsig(6),dh_dsig(6,6),deps_dsig(6),dh_de(6)
c$$$      real*8, intent(in):: cel,dphi_dsig,dh_de
c$$$      real*8, intent(out):: dh_dsig,deps_dsig
c$$$c     locals
c$$$      integer i,j
c$$$c     debar/dsig
c$$$      call calc_deps_dsig(Cel,dphi_dsig,deps_dsig)
c$$$c     dh_dsig = dh_de/de_dsig
c$$$c     tensor product
c$$$      do 20 j=1,6
c$$$      do 20 i=1,6
c$$$         dh_dsig(i,j) = dh_de(i)*deps_dsig(j)
c$$$ 20   continue
c$$$      return
c$$$      end subroutine calc_dh_dsig
c$$$c-----------------------------------------------------------------------
c$$$c     deps/dsig: see Eq 40 in Ref [1] (its inverse actually)
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
c$$$c-----------------------------------------------------------------------
c$$$c     Calculate {\partial s}/{\partial/sigma}
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
      elseif (ntens.eq.4) then
         ndi=3
         nshr=1
      elseif (ntens.eq.6) then
         ndi=3
         nshr=3
      else
         write(*,*)'ntens should be either 3 or 4 or 6 in',
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
c$$$c-----------------------------------------------------------------------
c$$$c     Calculate ds``/ds, dsp/ds
c$$$      subroutine calc_stress_deriv(emic,dsdp_ds)
c$$$      dimension emic(6)
c$$$      dimension dsdp_ds(6,6),dsp_ds(6,6),dsc_ds(6,6),dso_ds(6,6)
c$$$      real*8, intent(in):: emic
c$$$      real*8, intent(out)::dsdp_ds
c$$$c     locals
c$$$      real*8 dsc_ds,eta,dso_ds,dsp_ds
c$$$
c$$$      eta=1d0/gL
c$$$c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c$$$c     calc dsc_ds, equation 5 in Ref [1]
c$$$      dsc_ds(:,:)=0d0
c$$$      do 10 i=1,6
c$$$      do 10 j=1,6
c$$$         dsc_ds(i,j) = emic(i)*emic(j)
c$$$ 10   continue
c$$$      dsc_ds(:,:)=dsc_ds(:,:)*8d0/3d0
c$$$c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c$$$c     calc dso/ds
c$$$      dso_ds(:,:) = dsc_ds(:,:)
c$$$      do 15 i=1,6
c$$$         dso_ds(i,i)=1d0-dsc_ds(i,i)
c$$$ 15   continue
c$$$c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c$$$c     calculate ds``/ds  - see equation 9 in Ref [1]
c$$$      dsdp_ds(:,:)=0d0
c$$$      do 25 i=1,6
c$$$      do 25 j=1,6
c$$$         dsdp_ds(i,j) = dsc_ds(i,j) + eta*dso_ds(i,j)
c$$$ 25   continue
c$$$c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c$$$c     calculate dsp_ds - see Eq 11 in Ref [1]
c$$$      do 35 i=1,6
c$$$      do 35 j=1,6
c$$$         dsp_ds(i,j) =  4d0*(1d0-gs)* dso_ds(i,j)
c$$$ 35   continue
c$$$      return
c$$$      end subroutine calc_stress_deriv
c$$$c-----------------------------------------------------------------------
c$$$c     Calculate \partial f_k /\partial s and phi_b
c$$$c     Refer to Eqs 23 and 29 in Ref [1]
c$$$      subroutine calc_dfk_dsig(emic,sdev,deps_dsig,gs,e_ks,q,
c$$$     $     ys_iso,ys_hah,dfk_dsig,phi_b)
c$$$      implicit none
c$$$      dimension emic(6),sdev(6),deps_dsig(6),gs(4),e_ks(5),dfk_dsig(6)
c$$$      real*8, intent(in):: emic,e_ks,q,ys_iso,ys_hah,sdev,deps_dsig
c$$$      real*8, intent(out)::dfk_dsig,phi_b
c$$$c     locals
c$$$      integer ind
c$$$      dimension dfks_dgs(2),dgs_deps(2),gs_new(4),
c$$$     $     fks(2)
c$$$      real*8 dfks_dgs,dgs_deps,dot_prod,sh,gs_new,fks,gs
c$$$
c$$$c     calculate s:h
c$$$      sh = dot_prod(6,3,3,emic,sdev)
c$$$
c$$$c     obtain dfk/dgk and dgks/deps
c$$$c     gs_new and fks are dummy in this subroutine
c$$$      call calc_bau(6,3,3,emic,sdev,gs,e_ks,q,ys_iso,ys_hah,fks,gs_new,
c$$$     $     dfks_dgs,dgs_deps)
c$$$
c$$$c     obtain dfk_dsig
c$$$      if (sign(1d0,sh).ge.0) then
c$$$         ind=1
c$$$      else
c$$$         ind=2
c$$$      endif
c$$$      dfk_dsig(:)=0d0
c$$$      dfk_dsig(:)=dfks_dgs(ind)*dgs_deps(ind)*deps_dsig(:)
c$$$      phi_b = fks(ind) * sh * 2d0
c$$$      return
c$$$      end subroutine calc_dfk_dsig
c$$$c-----------------------------------------------------------------------
c$$$c     Eq. 33 in Ref [1]
c$$$      subroutine calc_dphib_dsig(emic,sdev,dhs_dsig,dfk_dsig,f_ks,
c$$$     $     dphib_dsig)
c$$$c     Arguments
c$$$c     emic
c$$$c     sdev
c$$$c     dhs_dsig
c$$$c     dfk_dsig
c$$$c     f_ks
c$$$c     dphib_dsig
c$$$      dimension emic(6), sdev(6), dhs_dsig(6),dfk_dsig(6),
c$$$     $     dphib_dsig(6),f_ks(2)
c$$$      real*8, intent(in)::emic,sdev,dhs_dsig,dfk_dsig,f_ks
c$$$      real*8, intent(out)::dphib_dsig
c$$$c     locals
c$$$      real*8 sh, dot_prod
c$$$      integer ind
c$$$c     calculate s:h
c$$$      sh = dot_prod(6,3,3,emic,sdev)
c$$$c     obtain dfk_dsig
c$$$      if (sign(1d0,sh).ge.0) then
c$$$         ind=1
c$$$      else
c$$$         ind=2
c$$$      endif
c$$$c     Eq. 33
c$$$      dphib_dsig(:)=2.d0*sign(1d0,sh)*(dfk_dsig(:)*sh+f_ks(ind)*
c$$$     $     dhs_dsig(:))
c$$$      return
c$$$      end subroutine calc_dphib_dsig
c$$$c-----------------------------------------------------------------------
