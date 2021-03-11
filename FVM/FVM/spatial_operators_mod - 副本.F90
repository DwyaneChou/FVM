module spatial_operators_mod
  use constants_mod
  use mesh_mod
  use parameters_mod
  use stat_mod
  use tend_mod
  !use projection_mod
  !use reconstruction_mod
  !use math_mod
  implicit none
  
  private
  
  public init_spatial_operator,spatial_operator
  
  logical, dimension(:,:,:), allocatable :: inDomain
  
  integer(i_kind), dimension(:,:,:,:), allocatable :: iGstCell ! x index of ghost reconstruction cells
  integer(i_kind), dimension(:,:,:,:), allocatable :: jGstCell ! y index of ghost reconstruction cells
  
    contains
    subroutine init_spatial_operator
      integer(i_kind) :: i,j,iPatch
      integer(i_kind) :: iCOS ! indices of Cells On Stencils
      integer(i_kind) :: iRec,jRec
    
      allocate(inDomain  (ims:ime,jms:jme,ifs:ife))
      
      allocate(iGstCell  (maxRecCells,ims:ime,jms:jme,ifs:ife))
      allocate(jGstCell  (maxRecCells,ims:ime,jms:jme,ifs:ife))
      
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
      
      !$OMP PARALLEL DO PRIVATE(i,j)
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
              
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine init_spatial_operator
    
    subroutine spatial_operator(stat,tend)
      type(stat_field), intent(inout) :: stat
      type(tend_field), intent(inout) :: tend
      
      call fill_halo(stat%q)
      
    end subroutine spatial_operator
    
    subroutine fill_halo(q)
      real(r_kind), dimension(nVar,ims:ime,jms:jme,ifs:ife), intent(inout) :: q
      
      integer(i_kind) :: iVar,i,j,iPatch,iPOC
      integer(i_kind) :: iRec,jRec
      integer(i_kind) :: m,n
      
      real(r_kind), dimension(maxRecCells) :: u
      
      real(r_kind) :: uc,vc,us,vs
      
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            m = nGstRecCells(i,j,iPatch)
            n = nGhostPointsOnCell(i,j,iPatch)
            do iVar = 1,nVar
              do iPOC = 1,m
                iRec = iGstCell(iPOC,i,j,iPatch)
                jRec = jGstCell(iPOC,i,j,iPatch)
                u(iPOC) = q(iVar,iRec,jRec,iPatch)
              enddo
              print*,''
              print*,u(1:m)
              print*,''
              
              do iPOC = 1,m
                iRec = iGstCell(iPOC,i,j,iPatch)
                jRec = jGstCell(iPOC,i,j,iPatch)
                !print*,q(iVar,iRec,jRec,iPatch)
                u(iPOC) = q(iVar,iRec,jRec,iPatch)
                print*,u(iPOC)
              enddo
              stop
            enddo
          enddo
        enddo
      enddo
    end subroutine fill_halo
    
end module spatial_operators_mod

