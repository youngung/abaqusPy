c-----------------------------------------------------------------------
c     Return mapping subroutine to find stress at n-1 step
c     and integrate all state variables during the given incremental
c     step defined by dstran (rate-independent ... yet?)
      subroutine return_mapping(Cel,spr,phi_n,eeq_n,dphi_n,dstran,
     $     stran_el,stran_pl,ntens,idiaw,hrdp,nhrdp,hrdc,nhrdc,
     $     ihard_law,iyld_law,yldc,nyldc,yldp_ns,nyldp,

c          variables to be updated.
     $     spd,snew,statev,nstatv,ddsdde,failnr
     $     )
c-----------------------------------------------------------------------
c***  Arguments
c     Cel    : elastic moduli
c     spr    : predictor stress at k=0
c     phi_n  : yield surface at the given step n
c     eeq_n  : accumulative plastic equivalent strain at step n
c     dphi_n : dphi at step n
c     dstran : total incremental strain given between steps n and n+1
c     stran_el: total elastic strain at step n
c     stran_pl: total plastic strain at step n
c     ntens   : len of stress/strain tensor (also applied to Cel)
c     idiaw   : whether or not turn on the diagnostic streams
c     hrdp    : hardening state variable arrays associated with hrdp
c               (e.g., equivalent plastic strain for isotropic hardening)
c     hrdc    : hardening constants (invariable)
c     nhrdp   : Len of hrdp
c     nhrdc   : Len of hrdc
c     ihard_law: hardening law (refer to uhard.f for more details)
c-----------------------------------------------------------------------
c***  Intents of Arguments
c     intent(in) Cel, spr, phi_n, eeq_n, dphi_n, dstran, stran_el,
c                stran_pl,ntens,idiaw, hrdp,nhrdp,hrdc,nhrdc,ihard_law
c     intent(out) -- states at n+1 step
c-----------------------------------------------------------------------
      implicit none
      character*255 fndia
      character*20 chr
      integer ntens,mxnr,nhrdc,nhrdp,ihard_law,iyld_law,nyldc,nyldp,
     $     nstatv
      parameter(mxnr=10)
c-----------------------------------------------------------------------
      dimension spr(ntens),dphi_n(ntens),snew(ntens),
     $     spr_ks(mxnr,ntens),statev(nstatv),
     $     dstran(ntens),

     $     stran_el(ntens),
     $     dstran_el(ntens),dstran_el_ks(mxnr,ntens),
     $     stran_el_ks(mxnr,ntens),
     $     stran_pl(ntens),
     $     dstran_pl(ntens),dstran_pl_ks(mxnr,ntens),
     $     stran_pl_ks(mxnr,ntens),

     $     aux_n(ntens),em_k(ntens),Cel(ntens,ntens),eeq_ks(mxnr),
     $     enorm_ks(mxnr,ntens),fo_ks(mxnr),fp_ks(mxnr),dlamb_ks(mxnr),
     $     dphi_ks(mxnr,ntens),d2phi_ks(mxnr,ntens,ntens),phi_ks(mxnr),
     $     dh_ks(mxnr),h_flow_ks(mxnr),

     $     hrdc(nhrdc),hrdp(nhrdp),
     $     yldc(nyldc),yldp_ns(0:1,nyldp),
     $     ddsdde(ntens,ntens)
c-----------------------------------------------------------------------
      real*8 Cel,spr,dphi_n,dstran,stran_el,dstran_el,dstran_el_ks
     $     ,stran_el_k,stran_el_ks,stran_pl,dstran_pl,dstran_pl_ks,
     $     stran_pl_k,stran_pl_ks,yldc,yldp_ns,statev,snew,ddsdde,
     $     spd
      real*8 seq_k,spr_ks       ! eq stress at nr-step k, stress predic at nr-step k
      real*8 enorm_ks           ! m_(n+alpha)
      real*8 fo_ks,fp_ks        ! Fobjective, Jacobian for NR
      real*8 dlamb_k,dlamb_ks,phi_n
      real*8 dphi_ks,d2phi_ks
      real*8 delta_eeq,eeq_n,aux_n,eeq_k,eeq_ks,empa,gpa
      real*8 h_flow_ks,dh_ks,phi_ks,em_k,tolerance,tol_val
      real*8 hrdc,hrdp
      integer k,idia,imsg
      parameter(tolerance=1d-4)
      logical idiaw,ibreak,failnr
      failnr=.false.             ! in case NR fails

      if (ntens.ne.3) then
         call fill_line(0,'*',72)
         write(*,*)'ntens:',ntens
         write(*,*)'Err: unexpected dimension of tensor given',ntens
         call fill_line(0,'*',72)
         stop -1
      endif

      empa=1d6
      gpa =1d9

c***  non-zero initial might be very risky.
c     delta_eeq = (dstran(1)**2+dstran(2)**2+dstran(2)**2)/3.d0 !! initial guess
      delta_eeq = 0d0           ! initial guess on equivalent strain rate contribution to dstran
      dlamb_ks(1) = delta_eeq
      spr_ks(1,:) = spr(:)       !! stress predictor
      dphi_ks(1,:) = dphi_n(:)
      phi_ks(1) = phi_n

      dstran_el_ks(1,:) = dstran(:)
      dstran_pl_ks(1,:) = 0d0   !! or using delta_eeq ...

c------------------------------------------------------------------------
c     iv. return mapping (loop over k)
      k=1

c      idia=315 ! write to diagnostic file
      idia=0   ! write to stdo
c      idia=7   ! write to Abaqus msg file

      if (idiaw) then
         call w_empty_lines(idia,3)
         call w_chr(idia,'Enter NR--')
      endif

      ibreak=.false.
      do while (k.le.mxnr)

c         s_k(:) = spr_ks(k,:)    ! predictor stress at current k
         em_k(:) = dphi_ks(k,:) ! yield normal at current k
         eeq_ks(k) = eeq_n + dlamb_ks(k) ! assumed plastic strain at current k

c***  Hardening state variable updates according to NR step k
         hrdp(1) = eeq_ks(k)
c***  --------------------------------
         call uhard(ihard_law,hrdp,nhrdp,hrdc,nhrdc,h_flow_ks(k),
     $        dh_ks(k),empa)

         if (k.eq.1) tol_val=h_flow_ks(1)*tolerance

         stran_el_ks(k,:) = dstran_el_ks(k,:) + stran_el(:)
         stran_pl_ks(k,:) = dstran_pl_ks(k,:) + stran_pl(:)

         if (idiaw) then
            call w_empty_lines(idia,2)
            call fill_line(idia,'*',72)
            call w_val(idia,'I-NR: ', float(k)*1d0)
            call w_chr(idia,'Spr_k')
            call w_dim(idia,spr_ks(k,:),ntens,1d0/empa,.true.)
            call w_chr(idia,'m_k')
            call w_dim(idia,em_k,ntens,1d0,.true.)
            call w_val(idia,'dlamb_ks(k) :',dlamb_ks(k))
            call w_val(idia,'eeq_ks(k)   :',eeq_ks(k))
            call w_val(idia,'phi_k [MPa] :',phi_ks(k)/empa)
            call w_val(idia,'h_flow_ks [MPa] :',h_flow_ks(k)/empa)
            call w_val(idia,'hf(k)/hf(1): [%]',h_flow_ks(k)/h_flow_ks(1)*1d2)
            call w_val(idia,'dh_ks [MPa] :',dh_ks(k)/empa)
         endif
c-----------------------------------------------------------------------
c        f   = yield - hardening             (objective function)
         fo_ks(k) = phi_ks(k) - h_flow_ks(k)
         if (abs(fo_ks(k)).le.tol_val)then
            goto 100
         else
c           Find Fp
c           ** Use values pertaining to n+1 step (assuming that current eeq_ks(k) is correct)

            call calc_fp(dphi_ks(k,:),Cel,dh_ks(k),ntens,fp_ks(k))
         endif

         if (idiaw) then
            call w_val(idia,'h_flow_k [MPa]:',h_flow_ks(k)/empa)
            call w_val(idia,'dh(k)[MPa]:',  dh_ks(k)/empa)
            call w_val(idia,'fo_ks(k)[MPa]:',fo_ks(k)/empa)
            call w_val(idia,'fp_ks(k)[GPa]:',fp_ks(k)/gpa)
            call w_val(idia,'fo/tol_val %   :',fo_ks(k)/tol_val*1d2)
         endif
c------------------------------------------------------------------------
c         2.  Update the multiplier^(k+1)  (dlamb)
c             dlamb^(k+1) = dlamb^k - fo_ks(k)/fp_ks(k)
         dlamb_ks(k+1) = dlamb_ks(k) + fo_ks(k)/fp_ks(k)

c     new plastic strain increment
         dstran_pl_ks(k+1,:) = dlamb_ks(k+1) * dphi_ks(k,:) ! backward
c     new elastic strain increment?
         dstran_el_ks(k+1,:) = dstran(:)  - dstran_pl_ks(k+1,:)
c     new plastic acc strain
         stran_pl_ks(k+1,:) = stran_pl(:) + dstran_pl_ks(k+1,:)
c        new elastic acc strain
c     Update dE^(el)^(k+1) and update the predictor stress.
         stran_el_ks(k+1,:) = stran_el(:) + dstran_el_ks(k+1,:)

c     find the new predictor stress for next NR step
c     Using  dE = dE^(el)^(k+1) + dlamb^(k+1),
         call mult_array(cel,stran_el_ks(k+1,:),ntens,spr_ks(k+1,:))
         if (idiaw) then
            call w_val(idia,'dlamb_ks(k+1)',dlamb_ks(k+1))
            call w_empty_lines(idia,2)
            call w_chr(idiaw,'dstran')
            call w_dim(idia,dstran(:),ntens,1d0,.false.)
            call w_empty_lines(idia,2)
            call w_chr(idia,'dstran_pl_k')
            call w_dim(idia,dstran_pl_ks(k,:),ntens,1d0,.false.)
            call w_chr(idia,'dstran_pl_k+1')
            call w_dim(idia,dstran_pl_ks(k+1,:),ntens,1d0,.false.)
            call w_empty_lines(idia,2)
            call w_chr(idia,'dstran_el_k')
            call w_dim(idia,dstran_el_ks(k,:),ntens,1d0,.false.)
            call w_chr(idia,'dstran_el_k+1')
            call w_dim(idia,dstran_el_ks(k+1,:),ntens,1d0,.false.)
            call w_empty_lines(idia,2)
            call w_chr(idia,' stran_pl_k')
            call w_dim(idia, stran_pl_ks(k,:),ntens,1d0,.false.)
            call w_chr(idia,' stran_pl_k+1')
            call w_dim(idia, stran_pl_ks(k+1,:),ntens,1d0,.false.)
            call w_empty_lines(idia,2)
            call w_chr(idia,' stran_el_k')
            call w_dim(idia, stran_el_ks(k,:),ntens,1d0,.false.)
            call w_chr(idia,' stran_el_k+1')
            call w_dim(idia, stran_el_ks(k+1,:),ntens,1d0,.false.)
            call w_chr(idia,' old predictor stress [MPa]')
            call w_dim(idia, spr_ks(k,:),ntens,1/empa,.false.)
            call w_chr(idia,' new predictor stress [MPa]')
            call w_dim(idia, spr_ks(k+1,:),ntens,1/empa,.false.)
         endif

c------------------------------------------------------------------------
c        3. Find normal of the updated predictor stress (s_(n+1)^(k+1))
         call update_yldp(iyld_law,yldp_ns,nyldp,dlamb_ks(k+1))
         call yld(iyld_law,yldp_ns(1,:),yldc,nyldp,nyldc,
     $        spr_ks(k+1,:),phi_ks(k+1),dphi_ks(k+1,:),
     $        d2phi_ks(k+1,:,:),ntens)
         k=k+1
      enddo ! end of do while loop for NR procedure

c     case when k exceeds mxnr
      call w_chr(idia,'Warning: NR procedure failed to converge')
      if (idiaw) then
         call fill_line(idia,'===',72)
      endif
      failnr=.true.
      return ! to get out before updating state variables.

c-----------------------------------------------------------------------
 100  continue ! successful NR run
      call w_chr(idia,'NR procedure converged')
      if (idiaw) call fill_line(idia,'===',72)
c***  update state variables
      call restore_statev(statev,nstatv,eeq_n+dlamb_ks(k),stran_el_ks(k+1,:),
     $     stran_pl_ks(k,:),ntens,yldp_ns(1,:),nyldp,1,.false.,idia)
c***  new stress
      snew(:)=spr_ks(k,:)
c$$$c***  plastic dissipation
c$$$      spd = spd +
c***  new jacobian
      call calc_epl_jacob(Cel,dphi_ks(k,:),dh_ks(k),ntens,ddsdde)

c     if (idiaw) close(idia)
      return
      end subroutine return_mapping
c-----------------------------------------------------------------------
c     Calculate fp using the below formula
c     fp  = r(s^eq_(n+1)^k)/r(s_(n+1)^k) : -C^el : r(s^eq_(n+1)^k / r(s_(n+1)^k) + H`)
c     fp = dphi_i C_ij dphi_j + H
      subroutine calc_fp(dphi,Cel,dh,ntens,fp)
c     intent(in) dphi,Cel,dh,ntens
c     intent(out) fp
c     dphi: round(s^eq)/round(s)
c     Cel : elastic modulus
c     dh  : round(s^flow)/round(e^eq)
c     ntens: len of dphi (also applied to Cel)
c     fp   : slope
      implicit none
      integer ntens
      dimension s(ntens),Cel(ntens,ntens),dphi(ntens)
      real*8 s,seq,Cel,dphi,fp,dh
      integer i,j
      fp=0.d0
      do 10 i=1,ntens
      do 10 j=1,ntens
         fp=fp+dphi(i) * Cel(i,j) * dphi(j) +dh
 10   continue
      return
      end subroutine calc_fp
c-----------------------------------------------------------------------
c     calculate elasto-plastic consistent tangent modulus
c     Following J. W. Yoon et. al., 1999, IJP 15, p35-67
      subroutine calc_epl_jacob(Cel,dphi,dh,ntens,jacob)
      implicit none
      integer ntens
      dimension Cel(ntens,ntens),dphi(ntens),jacob(ntens,ntens)
      real*8 Cel,dphi,jacob,dh

c     local varaiables
      dimension a(ntens),b(ntens)
      real*8 deno,a,b
      integer i,j,k

      deno=0d0
      a(:)=0d0
      b(:)=0d0
      jacob(:,:)=0d0

      do 5 i=1,ntens
      do 5 j=1,ntens
         deno = deno + dphi(i)*cel(i,j)*dphi(j)
         a(i) = a(i) + cel(i,j) * dphi(j)
         b(i) = b(i) + cel(i,j) * dphi(j)
 5    continue
      deno = deno + dh

      do 10 i=1,ntens
      do 10 j=1,ntens
         jacob(i,j) = cel(i,j) - a(i) * b(j) / deno
 10   continue


      write(*,*)'------------------------------------------'
      write(*,*)'dphi:',dphi
      write(*,*)'------------------------------------------'
      write(*,*)'dh:',dh
      write(*,*)'deno:',deno
      write(*,*)'a(i)'
      do 20 i=1,ntens
         write(*,'(3e13.3)') a(i)
 20   continue
      write(*,*)'b(i)'
      do 30 i=1,ntens
         write(*,'(3e13.3)') b(i)
 30   continue
      write(*,*)'------------------------------------------'

      write(*,*)'jacob - Cel'
      do 40 i=1,ntens
      do 40 j=1,ntens
         write(*,'(3e13.3)') (jacob(i,j) -cel(i,j))/1d9
 40   continue
      write(*,*)'------------------------------------------'


      return
      end subroutine calc_epl_jacob
