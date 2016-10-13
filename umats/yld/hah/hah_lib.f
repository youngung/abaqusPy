c-----------------------------------------------------------------------
c     General references
c     [1] Barlat et al. IJP 58, 2014 p201-218
c     [2] Jeong et al., IJP, 2016 (in press)
      subroutine hah_decompose(tensor,ntens,emic,
     $     tensor_colin,tensor_ortho)
c     Arguments
c     tensor       : subjected tensor
c     ntens        : Len of <tensor>
c     emic         : microstructure deviator
c     tensor_colin : colinear component of <tensor>
c     tensor_ortho : orthogonal component of <tensor>
      implicit none
      integer ntens
      dimension tensor(ntens)
      dimension tensor_colin(ntens)
      dimension tensor_ortho(ntens)
      dimension emic(ntens)
      real*8 tensor,tensor_colin,tensor_ortho,dot_prod,H,dd,emic
      integer i,imsg
      logical idiaw
cf2py intent(in) tensor, ntens, emic
cf2py intent(out) tensor_colin, tensor_ortho

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
      real*8 H, tensor6,tensor6_hat,tensor6d,hpress,dotproducts,
     $     denom
      integer i
cf2py intent(in) H,tensor6
cf2py intent(out) tensor6_hat

c***  it should be a deviator
      call deviat(tensor6,tensor6d,hpress)

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
      integer ntens
      dimension a(ntens),b(ntens)
      real*8 a,b,val,H,dot_prod
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
      real*8 a,b
      integer n,i
      dot_prod=0d0
      do 10 i=1,n
         dot_prod = dot_prod + a(i) * b(i)
 10   continue
      end function dot_prod
c------------------------------------------------------------------------
      subroutine pi_proj(sd,s1,s2)
c     Arguments
c     sd : Deviatoric stress
c     s1 : x coordinate of stress projected on pi-plane
c     s2 : y coordinate of stress projected on pi-plane
      implicit none
      real*8 sd(6),s1,s2
cf2py intent(in) sd
cf2py intent(out) s1,s2
      s1 = 2*sd(1)/sqrt(6d0)-sd(2)/sqrt(6d0)-sd(3)/sqrt(6d0)
      s2 =                   sd(2)/sqrt(2d0)-sd(3)/sqrt(2d0)
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
      integer iopt,nyldp,ntens
      dimension yldp(nyldp),emic(ntens),gk(4),e_ks(5),f_ks(2)
      real*8 yldp,emic,gk,e_ks,f_ks,eeq,ref,gL,ekL,eL,gS,c_ks,ss
c     local
      integer i

      if (iopt.eq.0) then       ! state variables <- yldp

c     HAH yield surface/state variables (and a few constants if any)
c***  equivalent plastic strain (cumulative)
         eeq     = yldp(1)
c***  reference size
         ref     = yldp(2)

c***  microstructure deviator
         do i=1,ntens
            emic(i) = yldp(i+2)
         enddo
c***  Bauschinger effect
         gk(:)   = yldp(ntens+2:ntens+5)   ! state variables
         e_ks(:) = yldp(ntens+5:ntens+9)  ! k1,k2,k3,k4,k5 constants
         do i=1,2
            f_ks(i) = yldp(ntens+9+i) ! f_k state parameters
         enddo
c***  Latent hardening
         gL      = yldp(ntens+12)
         ekL     = yldp(ntens+13)
         eL      = yldp(ntens+14)
c***  cross hardening
         gS      = yldp(ntens+15)
         c_ks    = yldp(ntens+16)
         ss      = yldp(ntens+17)
      elseif (iopt.eq.1) then   ! state variables -> yldp

c     HAH yield surface/state variables (and a few constants if any)
c***  equivalent plastic strain (cumulative)
         yldp(1)   = eeq
         yldp(2)   = ref
c***  microstructure deviator
         do i=1,ntens
            yldp(i+2) = emic(i)
         enddo
c***  Bauschinger effect
         yldp(ntens+2:ntens+5) = gk(:)      ! state variables
         yldp(ntens+5:ntens+9) = e_ks(:)  ! k1,k2,k3,k4,k5 constants
         do i=1,2
            yldp(ntens+9+i) = f_ks(i) ! f_k state parameters
         enddo
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
