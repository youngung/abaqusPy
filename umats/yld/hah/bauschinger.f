c------------------------------------------------------------------------
c     Bauschinger subroutines to be wrapped by f2py.
c     This module should be also used in Abaqus UMAT
c     Find more information of Abaqus UMAT in
c       www.github.com/youngung/abaqusPy
c
c     General references
c     [1] Barlat et al. IJP 58, 2014 p201-218
c     [2] Jeong et al., IJP, 2016 (in press)

c     Youngung Jeong
c     youngung.jeong@gmail.com
c------------------------------------------------------------------------
c     subroutine that calculates fk parameter for HAH
c     Refer to Eq 7 in Ref. [1]
c     The equation in Ref. [1] has a typo: missing exponent q
c     to the value of g_k
      subroutine calc_fk(gk,H,q,fk,dfk)
c     Arguments
c     gk
c     H
c     q  : the exponent
c     fk :f value
c     dfk: \partial(fk)/\partial(gk)
      implicit none
      real*8, intent(in) ::gk,H,q
      real*8, intent(out)::fk,dfk
      real*8, c
cf2py intent(in) gk,H,q
cf2py intent(out) fk,dfk
      c = dsqrt(6d0*H)/4d0
      fk  = (c * (1d0/(gk**q)-1d0)) ** (1d0/q)
c     d(fk)/d(gk) - to be used for derivative of yield surface
      dfk = c**(1d0/q) * (- gk**(-q-1d0) * fk **(1d0-q))
      return
      end subroutine
c------------------------------------------------------------------------
c$$$c     subroutine that calculates gk parameter for HAH
c$$$c         --   Eq 8 in Ref. [1]
c$$$      subroutine calc_gk(dk,e_ik,gk)
c$$$c     Arguments
c$$$c     dk
c$$$c     e_ik
c$$$c     gk
c$$$      implicit none
c$$$      real*8 dk,e_ik,gk
c$$$cf2py intent(in) dk,e_ik
c$$$cf2py intent(out) gk
c$$$      gk = dk / e_ik
c$$$      return
c$$$      end subroutine calc_gk
c------------------------------------------------------------------------
c     subroutine that calculates gk parameter for HAH
c         --   Eqs 7/8 in Ref. [1]
      subroutine calc_bau(e_mic,target,gs,e_ks,q,ys_iso,ys_hah,debar,
     $     fks,gs_new)
c     Arguments
c     e_mic  : microstructure deviator
c     target : the target direction towards which the microstructure
c              deviator aims to realign.
c     gs     : state variables that quantifies the <flattening>
c     e_ks   : k1,k2,k3,k4,k5 parameters that controls the evolunary
c              behavior of gs
c     q      : exponent
c     ys_iso : \bar{\sigma}(0) - the initial yield surface (without HAH effect)
c     ys_hah : \bar{\sigma}(\bar{\varepsilon}) - current yield surface
c     debar  : incremental equivalent strain
c              (used as multiplier when updating the state variables...)
c     fks    : fks state variables for n+1 step
c     gs_new : gs state variables for n+1 step
c     Note
c     Both arguments should be 'hat' properties.
c     Use <bauschinger_lib.f / hat> to convert a tensor
c     to its hatted property.
      implicit none
c     Arguments passed into
      integer ntens,ndi,nshr
      parameter(ntens=6,ndi=3,nshr=3)
      dimension e_mic(ntens)
      dimension target(ntens)
      dimension gs(4)
      dimension e_ks(5)
      dimension fks(2),dfks(2) ! dfks : dfk/dgk
      dimension gs_new(4)
      real*8 e_mic,target,gs,e_ks,fks,dfks,gs_new
      real*8 q,ys_iso,ys_hah,debar
c     locals
      dimension dgs(4)
      real*8 dot_prod,dd,dgs
      integer k
      integer i
cf2py intent(in) e_mic,target,gs,e_ks,ys_iso,ys_hah
cf2py intent(out) fks,gs_new
      dd = dot_prod(ntens,ndi,nshr,e_mic,target)
c***  index k depends on the sign of (mic:target)
      if (sign(1d0,dd).ge.0) then
         k = 1
      else
         k = 2
      endif
c***  Eqs 14&15 in Ref. [1]
      dgs(:)=0d0
      if (k.eq.1) then
         dgs(1) = e_ks(2) * (e_ks(3) * ys_iso/ys_hah - gs(1))
         dgs(2) = e_ks(1) * (  gs(3) - gs(2)) / gs(2)
         dgs(4) = e_ks(5) * (e_ks(4) - gs(4))
      elseif(k.eq.2) then
         dgs(1) = e_ks(1) * ( gs(4) - gs(1)) / gs(1)
         dgs(2) = e_ks(2) * (e_ks(3) * ys_iso/ys_hah - gs(2))
         dgs(3) = e_ks(5) * (e_ks(4) - gs(3))
      else
         call exit(-1)
      endif
      do 10 i=1,4
         gs_new(i) = gs(i) + dgs(i) * debar
 10   continue
      do 20 i=1,2
         call calc_fk(gs_new(i),8d0/3d0,q,fks(i),dfks(i))
 20   continue
      return
      end subroutine calc_bau
c------------------------------------------------------------------------
c     phib1 = f1**q  *  |dotp - |dotp||**q
c     phib1 = f2**q  *  |dotp + |dotp||**q
c     phib2
      subroutine bauschinger(ntens,ndi,nshr,emic,sdev,f_ks,q,phib1,
     $     phib2)
c     Arguments
c     ntens: Len of sdev, emic
c     ndi   : Number of normal components
c     nshr  : Number of shear components
c     emic : microstructure deviator
c     sdev : deviatoric stress
c     f_ks : f1, f2 parameters
c     q    : yield function exponent
      implicit none
      integer, intent(in) :: ntens,ndi,nshr
      dimension f_ks(2),emic(ntens),sdev(ntens)
      real*8, intent(in)  :: emic,sdev,f_ks,q
      real*8, intent(out) :: phib1,phib2
c     locals
      real*8 dotp,dot_prod
cf2py intent(in) f_ks,q,emic,sdev,ntens
cf2py intent(out) f1,f2
      dotp = dot_prod(ntens,ndi,nshr,emic,sdev)
      phib1 = (f_ks(1)**q) * (dabs(dotp-dabs(dotp))**q)
      phib2 = (f_ks(2)**q) * (dabs(dotp+dabs(dotp))**q)
      return
      end subroutine bauschinger
c------------------------------------------------------------------------
