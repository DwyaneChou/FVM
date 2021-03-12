module spatial_operators_mod
  use constants_mod
  use parameters_mod
  implicit none
      
  integer(i_kind) :: maxRecCells
  
  integer(i_kind), dimension(:,:,:,:), allocatable :: iGstCell ! x index of ghost reconstruction cells
  integer(i_kind), dimension(:,:,:,:), allocatable :: jGstCell ! y index of ghost reconstruction cells
  
  logical, dimension(:,:,:), allocatable :: inDomain
  
  integer(i_kind),dimension(:,:,:), allocatable :: nGstRecCells ! number of cells for ghost point reconstruction
  
    contains
    subroutine init_spatial_operator
      integer(i_kind) :: i,j,iPatch
      integer(i_kind) :: iCOS ! indices of Cells On Stencils
      integer(i_kind) :: iRec,jRec
        
      maxRecCells = stencil_width**2
        
      allocate(nGstRecCells(ids:ide,jds:jde,ifs:ife))
      
      allocate(iGstCell  (maxRecCells,ids:ide,jds:jde,ifs:ife))
      allocate(jGstCell  (maxRecCells,ids:ide,jds:jde,ifs:ife))
    
      allocate(inDomain  (ims:ime,jms:jme,ifs:ife))
        
      inDomain(ims:ime,jms:jme,ifs:ife) = .false. 
      inDomain(ids:ide,jds:jde,ifs:ife) = .true.
      
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            ! Ghost interpolation cells
            iCOS = 0
            do jRec = -recBdy,recBdy
              do iRec = -recBdy,recBdy
                if(inDomain(i+iRec,j+jRec,iPatch))then
                  iCOS = iCOS + 1
                  iGstCell(iCOS,i,j,iPatch) = i + iRec
                  jGstCell(iCOS,i,j,iPatch) = j + jRec
                endif
              enddo
            enddo
            nGstRecCells(i,j,iPatch) = iCOS
          enddo
        enddo
      enddo
      
      !$OMP PARALLEL
      !$OMP END PARALLEL
      
    end subroutine init_spatial_operator
    
    subroutine fill_halo
      
      real(r_kind), dimension(nVar,ims:ime,jms:jme,ifs:ife) :: qC
      !real(r_kind), dimension(:,:,:,:), allocatable :: qC
      
      integer(i_kind) :: iVar,i,j,iPatch,iPOC
      integer(i_kind) :: iRec,jRec
      integer(i_kind) :: m,n
      
      real(r_kind), dimension(maxRecCells) :: u
      
      !allocate(qC(nVar,ims:ime,jms:jme,ifs:ife))
      
      qC = 321321
      
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            m = nGstRecCells(i,j,iPatch)
            do iVar = 1,nVar
              do iPOC = 1,m
                iRec = iGstCell(iPOC,i,j,iPatch)
                jRec = jGstCell(iPOC,i,j,iPatch)
                u(iPOC) = qC(iVar,iRec,jRec,iPatch)
              enddo
              print*,''
              print*,u(1:m)
              print*,''
              stop
            enddo
          enddo
        enddo
      enddo
    end subroutine fill_halo
    
end module spatial_operators_mod

