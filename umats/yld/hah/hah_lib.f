c-----------------------------------------------------------------------
c     Library subroutines for Homogeneous Anisotropic Hardening model
c     General references

c     [1] Barlat et al., IJP 58, 2014 p201-218
c     [2] Jeong et al., IJP, 2016 (in press)

c     Youngung Jeong
c     youngung.jeong@gmail.com
c-----------------------------------------------------------------------
      subroutine hah_decompose(tensor,ntens,emic,
     $     tensor_colin,tensor_ortho)
c     Arguments
c     tensor       : subjected tensor
c     ntens        : Len of <tensor>
c     emic         : microstructure deviator
c     tensor_colin : colinear component of <tensor>
c     tensor_ortho : orthogonal component of <tensor>
      implicit none
      integer, intent(in):: ntens
      dimension tensor(ntens)
      dimension tensor_colin(ntens)
      dimension tensor_ortho(ntens)
      dimension emic(ntens)
      real*8,intent(in) ::tensor,emic
      real*8,intent(out)::tensor_colin,tensor_ortho
      real*8 dot_prod,H,dd
      integer imsg
      logical idiaw
cf2py intent(in) tensor, ntens, emic
cf2py intent(out) tensor_colin, tensor_ortho
cf2py depend(ntens) tensor,emic,tensor_colin,tensor_ortho

      idiaw=.false.
      imsg = 0
      if (idiaw) then
         call w_empty_lines(imsg,2)
         call fill_line(imsg,'#',52)
         call w_chr(imsg,'Enter HAH_DEDCOMPOSE')
         call w_chr(imsg,'Given <tensor>')
         call w_dim(imsg,tensor,ntens,1d0,.false.)
         call w_chr(imsg,'Given microstructure deviator <emic>')
         call w_dim(imsg,emic,ntens,1d0,.false.)
      endif

      H = 8d0/3d0
      dd = dot_prod(tensor,emic,ntens)
      if (idiaw) call w_val(imsg,'dd',dd)
c     Eqs 4&5 in Ref [1]
      tensor_colin(:) = H * dd * emic(:)
      tensor_ortho(:) = tensor(:) - tensor_colin(:)

      if (idiaw) then
         call w_chr(imsg,'Exit HAH_DEDCOMPOSE')
         call fill_line(imsg,'#',52)
         call w_empty_lines(imsg,1)
      endif

      return
      end subroutine hah_decompose
c-----------------------------------------------------------------------
c     Return a 'hat' property of the given tensor <tensor6>
c     Refer to eq 1 in Ref [1]
      subroutine hat(H,tensor6,tensor6_hat)
      implicit none
      dimension tensor6(6),tensor6d(6)
      dimension tensor6_hat(6)
      real*8, intent(in) :: tensor6
      real*8, intent(out):: tensor6_hat
      real*8 H,tensor6d,hydro,dotproducts,denom
      integer i
cf2py intent(in)  H,tensor6
cf2py intent(out) tensor6_hat
c***  it should be a deviator
      call deviat6(tensor6,tensor6d,hydro)
c***  calculate dotproducts
      dotproducts=0d0
      do 5 i=1,6
         dotproducts=dotproducts + tensor6d(i)*tensor6d(i)
 5    continue
c***  denominator in eq 1.
      denom = dsqrt(H * dotproducts)
c***  normalization.
      tensor6_hat(:) = tensor6d(:) / denom
      return
      end subroutine hat
c------------------------------------------------------------------------
c     Calculate cos2chi value using the two given hat tensors, <a> and <b>
      subroutine calc_cos2chi(a,b,ntens,val)
c     Arguments
c     a   : tensor in 6d (it should be a hat property)
c     b   : tensor in 6d (it should be a hat property)
c     val : value to be returned
      implicit none
      integer,intent(in):: ntens
      dimension a(ntens),b(ntens)
      real*8, intent(in)::a,b
      real*8, intent(out)::val
      real*8  H,dot_prod
cf2py intent(in) a,b,ntens
cf2py intent(out) val
      H = 8d0/3d0
      val = dot_prod(a,b,ntens) * H
      if (abs(val).gt.1d0) then
         write(*,*)'Something went wrong in calc_cos2chi'
      endif
      return
      end subroutine calc_cos2chi
c------------------------------------------------------------------------
      real*8 function dot_prod(a,b,n)
c     Arguments
c     a   : tensor in n-dimension
c     b   : tensor in n-dimension
c     n   : Len of a and b
      implicit none
      dimension a(n)
      dimension b(n)
      real*8, intent(in):: a,b
      integer, intent(in):: n
      integer i
      dot_prod=0d0
      do 10 i=1,n
         dot_prod = dot_prod + a(i) * b(i)
 10   continue
      end function dot_prod
c------------------------------------------------------------------------
      subroutine pi_proj(sdev,s1,s2)
c     Arguments
c     sdev : Deviatoric stress
c     s1   : x coordinate of stress projected on pi-plane
c     s2   : y coordinate of stress projected on pi-plane
      implicit none
      real*8, intent(in)::  sdev(6)
      real*8, intent(out):: s1,s2
cf2py intent(in) sdev
cf2py intent(out) s1,s2
      s1 = 2*sdev(1)/sqrt(6d0)-sdev(2)/sqrt(6d0)-sdev(3)/sqrt(6d0)
      s2 =                     sdev(2)/sqrt(2d0)-sdev(3)/sqrt(2d0)
      return
      end subroutine pi_proj
c-----------------------------------------------------------------------
c     subroutine that converts back and forth
c     Used to convert/restore yldp of HAH case.
c     When yldp is used for HAH model, this subroutine must be used.
      subroutine hah_io(iopt,nyldp,ntens,yldp,emic,gk,e_ks,f_ks,eeq,ref,
     $     gL,ekL,eL,gS,c_ks,ss)
c     Arguments
c     iopt  :  if 0, state variables <- yldp
c              if 1, state variables -> yldp
c     nyldp : Len of yldp
c     ntens : Len of stress tensor and <emic>
c     yldp  : yield surface state variables
c     emic  : microstructure deviator
c     eeq   : equivalent plastic strain
c     ref   : yield surface reference sizen
c             used to renormalized yield surface
c     gk    : gk parameters
c     e_ks  : ks parameters for Bauschinger effect
c     f_ks  : f1 and f2 parameters for Bauschinger effect
c     gL    : gL latent hardening parameter
c     ekL   : kL latent hardening parameter
c     eL    : L  latent hardening parameter
c     gS    : gS parameter for cross hardening
c     c_ks  : ks parameter for cross hardening
c     ss    : S  parameter for cross hardening
      implicit none
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c**   Arguments passed
      integer, intent(in)::iopt,nyldp,ntens
      dimension yldp(nyldp),emic(ntens),gk(4),e_ks(5),f_ks(2)
      real*8, intent(inout):: yldp,emic,gk,e_ks,f_ks,eeq,ref,gL,ekL,eL,
     $     gS,c_ks,ss
c     local
      integer i
      if (iopt.eq.0) then       ! state variables <- yldp
c     HAH yield surface/state variables (and a few constants if any)
c***  equivalent plastic strain (cumulative)
         eeq     = yldp(1)
c***  Reference size
         ref     = yldp(2)
c***  Microstructure deviator
         do 5 i=1,ntens
            emic(i) = yldp(i+2)
 5       continue
c***  Bauschinger effect
         gk(:)   = yldp(ntens+2:ntens+5)   ! state variables
         e_ks(:) = yldp(ntens+5:ntens+9)  ! k1,k2,k3,k4,k5 constants
         do 10 i=1,2
            f_ks(i) = yldp(ntens+9+i) ! f_k state parameters
 10      continue
c***  Latent hardening
         gL      = yldp(ntens+12)
         ekL     = yldp(ntens+13)
         eL      = yldp(ntens+14)
c***  cross hardening
         gS      = yldp(ntens+15)
         c_ks    = yldp(ntens+16)
         ss      = yldp(ntens+17)
c     diagnose
         if (dabs(gL).lt.1e-3) then
            write(*,*)'gL is too small'
            call exit(-1)
         endif
      elseif (iopt.eq.1) then   ! state variables -> yldp
c     HAH yield surface/state variables (and a few constants if any)
c***  equivalent plastic strain (cumulative)
         yldp(1)   = eeq
c***  Reference size
         yldp(2)   = ref
c***  Microstructure deviator
         do 15 i=1,ntens
            yldp(i+2) = emic(i)
 15      continue
c***  Bauschinger effect
         yldp(ntens+2:ntens+5) = gk(:)      ! state variables
         yldp(ntens+5:ntens+9) = e_ks(:)  ! k1,k2,k3,k4,k5 constants
         do 20 i=1,2
            yldp(ntens+9+i) = f_ks(i) ! f_k state parameters
 20      continue
c***  Latent hardening
         yldp(ntens+12) = gL
         yldp(ntens+13) = ekL
         yldp(ntens+14) = eL
c***  cross hardening
         yldp(ntens+15) = gS
         yldp(ntens+16) = c_ks
         yldp(ntens+17) = ss
      endif
c     yield surface constants
      return
      end subroutine hah_io
c-----------------------------------------------------------------------
c     Calculate referece size of yield surface and save it to yldp
      subroutine hah_calc_ref(ntens,nyldp,nyldc,yldp,yldc,iyld_choice)
      implicit none
c     Arguments
c     ntens       : Len of stress
c     nyldp       : Len of yldp
c     nyldc       : Len of yldc
c     yldp        : Yield function parameters
c     yldc        : Yield function constants
c     iyld_choice : yield surface choice
      integer, intent(in) :: ntens,nyldp,nyldc,iyld_choice
      dimension yldp(nyldp),yldc(nyldc)
      real*8 yldp,yldc
c     hah_io
      dimension emic(ntens),gk(4),e_ks(5),f_ks(2)
      real*8 emic,gk,e_ks,f_ks,eeq,ref,gL,ekL,eL,gS,c_ks,ss
c     locals
      dimension cauchy_ref(ntens)
      real*8 cauchy_ref,phi

      call hah_io(0,nyldp,ntens,yldp,emic,gk,e_ks,f_ks,eeq,ref,
     $     gL,ekL,eL,gS,c_ks,ss)
      cauchy_ref(:)=0d0
      cauchy_ref(1)=1d0
c     cauchy_ref(2)=1d0
c     returns:  (sqrt(phi(sp)**2 + phi(sdp)**2)) ** q
      call latent(iyld_choice,ntens,nyldp,nyldc,cauchy_ref,yldp,yldc,
     $     phi)
      ref = phi**(1d0/yldp(9))
c      call w_val(imsg,'ref',ref)
c     save ref to yldp
      call hah_io(1,nyldp,ntens,yldp,emic,gk,e_ks,f_ks,eeq,ref,
     $     gL,ekL,eL,gS,c_ks,ss)
      return
      end subroutine hah_calc_ref
c-----------------------------------------------------------------------
