    module weno_coef_mod
      use constants_mod
      implicit none
      
    contains
    subroutine calc_weno_coef(coef,stencil_width,nPointsOnEdge,quad_pos)
      real(r_kind),dimension(:,:) :: coef
    end subroutine calc_weno_coef
    
    end module weno_coef_mod