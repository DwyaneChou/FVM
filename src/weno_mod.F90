    module weno_coef_mod
      use constants_mod
      implicit none
      
    contains
    subroutine calc_weno_coef(coef,stencil_width,nPointsOnEdge,quad_pos)
      real   (r_kind), dimension(:,:), intent(out) :: coef
      integer(i_kind)                , intent(in ) :: stencil_width
      integer(i_kind)                , intent(in ) :: nPointsOnEdge
      real   (r_kind), dimension(:,:), intent(in ) :: quad_pos
      
    end subroutine calc_weno_coef
    
    end module weno_coef_mod