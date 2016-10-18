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
c     calculate incremental gL (dgL)
      subroutine latent_update(ntens,ndi,nshr,emic,target,eL,gL,ekL,
     $     ys_iso,ys_hah,debar,dgL)
c     Arguments
c     ntens  : Len of <emic> and <target>
c     ndi    : Number of normal components
c     nshr   : Number of shear components
c     target : the target direction, with which the microstructure
c               deviator aims to realign.
c     eL     : L parameter
c     gL     : gL parameter at step n
c     ekL    : kL constant
c     ys_iso : \bar{\sigma}(0) - the initial yield surface (without HAH effect)
c     ys_hah : \bar{\sigma}(\bar{\varepsilon}) - current yield surface
c     emic   : microstructure deviator
c     debar  : incremental equivalent strain
c               (used as multiplier when updating the state variables...)
c     dgL    : incremental gL
      implicit none
c     Arguments of subroutine
      integer, intent(in):: ntens,ndi,nshr
      dimension emic(ntens),target(ntens)
      real*8, intent(in) :: emic,target,eL,gL,ekL,ys_iso,ys_hah,debar
      real*8, intent(out):: dgL
c     Local variables
      real*8 cos2chi,term
cf2py intent(in) eL,gL,ekL,ys_iso,ys_hah,emic,target,ntens,debar
cf2py intent(out) dgL
cf2py depend(ntens) emic,target

      call calc_cos2chi(ntens,ndi,nshr,target,emic,cos2chi)
c     Eq 16 --
      term = dsqrt(eL * (1d0-cos2chi) + cos2chi)-1d0
      dgL = ekL *( (ys_hah-ys_iso) / ys_hah * term  + 1d0 - gL )
      dgL = dgL * debar

      return
      end subroutine latent_update

c------------------------------------------------------------------------
c     Latent hardening effect accounted and saved to phi
c     returns:  (sqrt(phi(sp)**2 + phi(sdp)**2)) ** q
      subroutine latent(iyld_law,ntens,ndi,nshr,nyldp,nyldc,cauchy,yldp,
     $     yldc,phi)
c     Arguments
c     iyld_law : type of yield function
c     ntens    : Len of sdev, emic
c     ndi      : Number of normal components
c     nshr     : Number of shear components
c     nyldp    : Len of yldp
c     nyldc    : Len of yldc
c     cauchy   : cauchy stress
c     yldp     : yldp
c     yldc     : yldc
c     phi      : (phi1**2+phi2**2)**(1/2)
      implicit none
c     Arguments passed
      integer,intent(in) :: iyld_law,ntens,ndi,nshr,nyldp,nyldc
      dimension cauchy(ntens),yldp(nyldp),yldc(nyldc)
      real*8, intent(inout) ::  yldp,yldc
      real*8, intent(in)  :: cauchy
      real*8, intent(out) :: phi
c     locals
      dimension sdev(ntens),so(ntens), sc(ntens), sdp(ntens),
     $     emic(ntens),dphi(ntens),d2phi(ntens)
      real*8 sdev,so,sc,sdp,emic,dphi,d2phi
c     variables to be stored from yldp
      dimension gk(4),e_ks(5),f_ks(2),sp(ntens),phis(2)
      real*8 gk,e_ks,f_ks,eeq,ref,gL,ekL,eL,gS,c_ks,ss,sp,phis,hydro
      logical idiaw
      integer imsg
      imsg=0
c      idiaw=.true.
      idiaw=.false.

c**   restore parameters from yldp
      call hah_io(0,nyldp,ntens,yldp,emic,gk,e_ks,f_ks,eeq,ref,gL,ekL,
     $     eL,gS,c_ks,ss)

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
      sdp = 1d0*sc(:) + so(:)/gL !! orthogonal distortion
c     4. sp = 4(1-g_s) s_o
      sp(:) = 4d0*(1d0-gS) * so(:)/gL
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
         call yld2000_2d_dev(sdp,phis(1),dphi,d2phi,yldc)
         call yld2000_2d_dev(sp, phis(2),dphi,d2phi,yldc)
         phi = dsqrt(phis(1)**2 + phis(2)**2)
      else
         call w_chr(imsg,'** iyld_law not expected in latent')
         call exit(-1)
      endif

      return
      end subroutine latent
