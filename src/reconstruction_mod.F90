    module reconstruction_mod
      use constants_mod, only: i_kind,r_kind
      use parameters_mod, only: ids,ide,jds,jde,ifs,ife,stencil_width
      implicit none
      
      integer(i_kind) :: maxRecCells
      integer(i_kind) :: maxRecTerms
      
      integer(i_kind),dimension(:,:,:), allocatable :: nGstRecCells ! number of cells for ghost point reconstruction
      
    contains
      subroutine init_reconstruction
        
        maxRecCells = stencil_width**2
        maxRecTerms = maxRecCells
        
        allocate(nGstRecCells(ids:ide,jds:jde,ifs:ife))
        
      end subroutine init_reconstruction
      
    end module reconstruction_mod
    