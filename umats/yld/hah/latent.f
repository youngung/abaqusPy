c------------------------------------------------------------------------
c     Latent hardening subroutines to be wrapped by f2py.
c     This module should be also used in Abaqus UMAT (under development...)
c     Find more information of Abaqus UMAT in
c       www.github.com/youngung/abaqusPy
c
c     General references
c     [1] Barlat et al. IJP 58, 2014 p201-218
c     [2] Jeong et al., IJP, 2016 (in press)

c     Youngung Jeong
c     youngung.jeong@gmail.com
c------------------------------------------------------------------------
c     calculate incremental gL (dgL_deps)
      subroutine latent_update(ntens,emic,target,eL,gL,ekL,
     $     ys_0,ys_hah,dgL_deps)
c     Arguments
c     ntens  : Len of <emic> and <target>
c     target : the target direction, with which the microstructure
c               deviator aims to realign.
c     eL     : L parameter
c     gL     : gL parameter at step n
c     ekL    : kL constant
c     ys_0   : \bar{\sigma}(0) - the initial yield surface (without HAH effect)
c     ys_hah : \bar{\sigma}(\bar{\varepsilon}) - current yield surface
c     emic   : microstructure deviator
c     dgL_deps
      implicit none
c     Arguments of subroutine
      integer, intent(in):: ntens
      dimension emic(ntens),target(ntens)
      real*8, intent(in) :: emic,target,eL,gL,ekL,ys_0,ys_hah
      real*8, intent(out):: dgL_deps
c     Local variables
      real*8 coschi,term
cf2py intent(in) eL,gL,ekL,ys_0,ys_hah,emic,target,ntens
cf2py intent(out) dgL_deps
cf2py depend(ntens) emic,target

      call calc_coschi(ntens,target,emic,coschi)
c     Eq 16 in Ref [1]
      term = dsqrt(eL * (1d0-coschi*coschi) + coschi*coschi)-1d0
      dgL_deps = ekL *( (ys_hah-ys_0) / ys_hah * term  + 1d0 - gL )

      return
      end subroutine latent_update
c------------------------------------------------------------------------
c     Latent hardening effect accounted and saved to phi_h
c     returns:  sqrt(phi(sp)**2 + phi(sdp)**2)
      subroutine latent(iyld_law,ntens,ndi,nshr,nyldp,nyldc,cauchy,yldp,
     $     yldc,dpsi_sdp,dpsi_sp,psi_sdp,psi_sp,phi_h)
c     Arguments
c-----------------------------------------------------------------------
c     iyld_law : type of yield function
c     ntens    : Len of sdev, emic
c     ndi      : Number of normal components
c     nshr     : Number of shear components
c     nyldp    : Len of yldp
c     nyldc    : Len of yldc
c     cauchy   : cauchy stress
c     yldp     : yldp
c     yldc     : yldc
c     dpsi_spd
c     dpsi_sp
c     psi_spd
c     psi_sp
c     phi_h    : sqrt(psi(sp)**2 + psi(sdp)**2)
      implicit none
c     Arguments passed
      integer,intent(in)    :: iyld_law,ntens,ndi,nshr,nyldp,nyldc
      dimension cauchy(ntens), yldp(nyldp),yldc(nyldc)
      real*8, intent(inout) :: yldp,yldc
      real*8, intent(in)    :: cauchy
      real*8, intent(out)   :: phi_h
c     locals
      dimension sdev(ntens),so(ntens),sc(ntens),sdp(ntens),emic(ntens),
     $     demic(ntens),krs(4),target(ntens),dpsi_sdp(ntens),
     $     dpsi_sp(ntens),d2psi(ntens)
      real*8 sdev,so,sc,sdp,emic,demic,dgr,krs,target,d2psi,dpsi_sdp,
     $     dpsi_sp
c     variables to be stored from yldp
      dimension gk(4),e_ks(5),f_ks(2),sp(ntens)
      real*8 gk,e_ks,f_ks,eeq,ref0,ref1,gL,ekL,eL,gS,c_ks,ss,sp,psi_sdp,
     $     psi_sp,hydro
      logical idiaw
      integer imsg
      imsg=0
c      idiaw=.true.
      idiaw=.false.

c**   restore parameters from yldp
      call hah_io(0,nyldp,ntens,yldp,emic,demic,dgr,gk,e_ks,f_ks,eeq,
     $     ref0,ref1,gL,ekL,eL,gS,c_ks,ss,krs,target)

      if (idiaw) then
         call w_val(imsg,'gL:',gL)
         call w_val(imsg,'gS:',gS)
      endif

c**   Apply linear transformation to stress deviator
c     1. deviator
      call deviat(ntens,cauchy,sdev,hydro)
c     2. Obtain orthogonal / collinear components
      call hah_decompose(ntens,ndi,nshr,sdev,emic,sc,so)
c     3. Transform to allow extension along so (double prime s)
      sdp(:) = 1d0*sc(:) + so(:)/gL !! orthogonal distortion
c     4. sp = 4(1-g_s) s_o
      sp(:) = 4d0*(1d0-gS) * so(:)
      if (idiaw) then
         call w_chr(imsg,'cauchy')
         call w_dim(imsg,cauchy,ntens,1d0,.false.)
         call w_chr(imsg,'sdev')
         call w_dim(imsg,sdev,ntens,1d0,.false.)
         call w_chr(imsg,'emic')
         call w_dim(imsg,emic,ntens,1d0,.false.)
         call w_chr(imsg,'sc')
         call w_dim(imsg,sc,ntens,1d0,.false.)
         call w_chr(imsg,'so')
         call w_dim(imsg,so,ntens,1d0,.false.)
         call w_chr(imsg,'sdp')
         call w_dim(imsg,sdp,ntens,1d0,.false.)
         call w_chr(imsg,'sp')
         call w_dim(imsg,sp, ntens,1d0,.false.)
      endif
      if (iyld_law.eq.2) then
         !! psi_sdp, dpsi_sdp
         call yld2000_2d_dev(sdp,psi_sdp,dpsi_sdp,d2psi,yldc)
         !! psi_sp,  dpsi_sp
         call yld2000_2d_dev(sp, psi_sp, dpsi_sp, d2psi,yldc)
         phi_h = dsqrt(psi_sdp*psi_sdp + psi_sp*psi_sp)
      else
         call w_chr(imsg,'** iyld_law not expected in latent')
         call exit(-1)
      endif
c**   Calculate derivative (i.e., \frac{\partial\phi_h}{\partial\cauchy})
      return
      end subroutine latent
