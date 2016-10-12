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
      integer ntens
      dimension tensor(ntens)
      dimension tensor_colin(ntens)
      dimension tensor_ortho(ntens)
      dimension emic(ntens)
      real*8 tensor,tensor_colin,tensor_ortho,dot_prod,H,dd
      integer i
cf2py intent(in) tensor, ntens, emic
cf2py intent(out) tensor_colin, tensor_ortho

      H = 8d0/3d0

      dd = dot_prod(tensor,emic)

c     Eqs 4&5 in Ref [1]
      tensor_colin(:) = H * dd * emic(:)
      tensor_ortho(:) = tensor(:) - tensor_colin(:)

      return
      end subroutine hah_decompose
