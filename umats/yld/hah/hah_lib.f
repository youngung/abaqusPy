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
cf2py intent(in) tensor, ntens, emic
cf2py intent(out) tensor_colin, tensor_ortho

      imsg = 0
      call w_empty_lines(imsg,2)
      call fill_line(imsg,'#',52)
      call w_chr(imsg,'Enter HAH_DEDCOMPOSE')
      call w_chr(imsg,'Given <tensor>')
      call w_dim(imsg,tensor,ntens,1d0,.false.)
      call w_chr(imsg,'Given microstructure deviator <emic>')
      call w_dim(imsg,emic,ntens,1d0,.false.)

      H = 8d0/3d0

      dd = dot_prod(tensor,emic,ntens)
      
      call w_val(imsg,'dd',dd)

c     Eqs 4&5 in Ref [1]
      tensor_colin(:) = H * dd * emic(:)
      tensor_ortho(:) = tensor(:) - tensor_colin(:)

      call fill_line(imsg,'#',52)
      call w_chr(imsg,'Exit HAH_DEDCOMPOSE')
      call w_empty_lines(imsg,1)

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
      subroutine calc_cos2chi(a,b,ntens,val)
c     Arguments
c     a   : tensor in 6d
c     b   : tensor in 6d
c     val : value to be returned
      implicit none
      integer ntens
      dimension a(ntens),b(ntens)
      real*8 a,b,val,H,dot_prod
cf2py intent(in) a,b,ntens
cf2py intent(out) val
      H = 8d0/3d0
      val = dot_prod(a,b,ntens) * H
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
c$$$      subroutine deviat(s,sd,sm)
c$$$c     Arguments
c$$$c     s : tensor
c$$$c     sd: deviator of <s>
c$$$c     sm: hydrostatic pressure
c$$$      implicit none
c$$$      real*8 s(6),sd(6),sm
c$$$cf2py intent(in) s
c$$$cf2py intent(out) sd,sm
c$$$      sm = (s(1)+s(2)+s(3))/3.
c$$$      sd(1) = s(1)-sm
c$$$      sd(2) = s(2)-sm
c$$$      sd(3) = s(3)-sm
c$$$      sd(4) = s(4)
c$$$      sd(5) = s(5)
c$$$      sd(6) = s(6)
c$$$      return
c$$$      end subroutine deviat
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
      subroutine hah_io(yldp,nyldp,
     $     eeq,

     $     ntens,

     $     emic,

     $     gk,e_ks,f_ks,

     $     gL,ekL,eL,

     $     gS,c_ks,ss,

     $     iopt)
c     Arguments
c     yldp  : yield surface state variables
c     nyldp : Len of yldp
c     eeq   : equivalent plastic strain
c     emic   : microstructure deviator
c     gk    : gk parameters
c     gL    : gL
c     ekL   : kL
c     eL    : L parameters
c     gS    : gS parameter
c     c_ks  : ks parameter
c     ss    : S parameter
c     iopt  :  if 0, state variables <- yldp
c              if 1, state variables -> yldp
      implicit none
      integer iopt,nyldp
      dimension yldp(nyldp)
      real*8 yldp

c     local - microstructure deviator
      integer ntens
      dimension emic(ntens)
      real*8 emic
c     local - Bauschinger parameters
      dimension gk(4)
      dimension e_ks(5)
      dimension f_ks(2)
      real*8 gk,e_ks,f_ks,eeq
c     local - Latent hardening parameters
      real*8 gL,ekL,eL
c     local - cross hardening parameters
      real*8 gS,c_ks,ss

c     local
      integer i

      if (iopt.eq.0) then       ! state variables <- yldp

c     HAH yield surface/state variables (and a few constants if any)
c***  equivalent plastic strain (cumulative)
         eeq     = yldp(1)
c***  microstructure deviator
         do i=1,ntens
            emic(i) = yldp(i+1)
         enddo
c***  Bauschinger effect
         gk(:)   = yldp(ntens+1:ntens+4)   ! state variables
         e_ks(:) = yldp(ntens+4:ntens+8)  ! k1,k2,k3,k4,k5 constants
         do i=1,2
            f_ks(i) = yldp(ntens+8+i) ! f_k state parameters
         enddo
c***  Latent hardening
         gL      = yldp(ntens+11)
         ekL     = yldp(ntens+12)
         eL      = yldp(ntens+13)
c***  cross hardening
         gS      = yldp(ntens+14)
         c_ks    = yldp(ntens+15)
         ss      = yldp(ntens+16)
      elseif (iopt.eq.1) then   ! state variables -> yldp

c     HAH yield surface/state variables (and a few constants if any)
c***  equivalent plastic strain (cumulative)
         yldp(1)   = eeq
c***  microstructure deviator
         do i=1,ntens
            yldp(i+1) = emic(i)
         enddo
c***  Bauschinger effect
         yldp(ntens+1:ntens+4) = gk(:)      ! state variables
         yldp(ntens+4:ntens+8) = e_ks(:)  ! k1,k2,k3,k4,k5 constants
         do i=1,2
            yldp(ntens+8+i) = f_ks(i) ! f_k state parameters
         enddo
c***  Latent hardening
         yldp(ntens+11) = gL
         yldp(ntens+12) = ekL
         yldp(ntens+13) = eL
c***  cross hardening
         yldp(ntens+14) = gS
         yldp(ntens+15) = c_ks
         yldp(ntens+16) = ss
      endif

c     yield surface constants

      return
      end subroutine hah_io
c-----------------------------------------------------------------------
