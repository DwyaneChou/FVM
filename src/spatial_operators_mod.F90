module spatial_operators_mod
  use constants_mod
  use mesh_mod
  use parameters_mod
  use stat_mod
  use tend_mod
  use projection_mod
  
  implicit none
  
  private
  
  public init_spatial_operator,spatial_operator
  
  logical, dimension(:,:,:), allocatable :: inDomain
  
  integer(i_kind), dimension(:,:,:), allocatable :: iCenCell ! center cell index on reconstruction stencil
  
  integer(i_kind), dimension(:,:,:,:), allocatable :: iRecCell ! x index of reconstruction cells
  integer(i_kind), dimension(:,:,:,:), allocatable :: jRecCell ! y index of reconstruction cells
  
  integer(i_kind), dimension(:,:,:,:), allocatable :: iGstCell ! x index of ghost reconstruction cells
  integer(i_kind), dimension(:,:,:,:), allocatable :: jGstCell ! y index of ghost reconstruction cells
  
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: xRel ! relative x coordinate of reconstruction cells
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: yRel ! relative y coordinate of reconstruction cells
  
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: xGst ! relative x coordinate of Ghost reconstruction cells
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: yGst ! relative y coordinate of Ghost reconstruction cells
  
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: recMatrixL
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: recMatrixR
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: recMatrixB
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: recMatrixT
  
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: A
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: invA
  
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: gstMatrixL
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: gstMatrixR
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: gstMatrixB
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: gstMatrixT
  
  real   (r_kind) :: dV
  
  real   (r_kind), dimension(:,:,:,:  ), allocatable :: qC
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: qL
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: qR
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: qB
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: qT
  
  real(r_kind), dimension(:,:,:,:,:), allocatable :: FL    ! Reconstructed F_(i-1/2,j)
  real(r_kind), dimension(:,:,:,:,:), allocatable :: FR    ! Reconstructed F_(i+1/2,j)
  real(r_kind), dimension(:,:,:,:,:), allocatable :: FB    ! Reconstructed F_(i,j-1/2)
  real(r_kind), dimension(:,:,:,:,:), allocatable :: FT    ! Reconstructed F_(i,j+1/2)
      
  real(r_kind), dimension(:,:,:,:,:), allocatable :: GL    ! Reconstructed G_(i-1/2,j)
  real(r_kind), dimension(:,:,:,:,:), allocatable :: GR    ! Reconstructed G_(i+1/2,j)
  real(r_kind), dimension(:,:,:,:,:), allocatable :: GB    ! Reconstructed G_(i,j-1/2)
  real(r_kind), dimension(:,:,:,:,:), allocatable :: GT    ! Reconstructed G_(i,j+1/2)
  
  real(r_kind), dimension(:,:,:,:), allocatable :: Fe    ! F on edges of each cell
  real(r_kind), dimension(:,:,:,:), allocatable :: Ge    ! H on edges of each cell
  
  real(r_kind), dimension(:,:,:,:,:), allocatable :: FeP   ! F on points on edges of each cell
  real(r_kind), dimension(:,:,:,:,:), allocatable :: GeP   ! H on points on edges of each cell
  
  real(r_kind), dimension(:), allocatable :: xL
  real(r_kind), dimension(:), allocatable :: xR
  real(r_kind), dimension(:), allocatable :: xB
  real(r_kind), dimension(:), allocatable :: xT
  
  real(r_kind), dimension(:), allocatable :: yL
  real(r_kind), dimension(:), allocatable :: yR
  real(r_kind), dimension(:), allocatable :: yB
  real(r_kind), dimension(:), allocatable :: yT
  
    contains
    subroutine init_spatial_operator
      integer(i_kind) :: i,j,iPatch
      integer(i_kind) :: iCOS ! indices of Cells On Stencils
      integer(i_kind) :: iRec,jRec
      integer(i_kind) :: iQP,jQP
      integer(i_kind) :: iR,jR
      
      integer(i_kind) :: nxp,nyp
      integer(i_kind) :: invstat
      integer(i_kind) :: nRC
    
      real(r_kind) :: recCoef
      real(r_kind) :: recdx  
      real(r_kind) :: recdy  
      real(r_kind) :: recdV  
    
      allocate(inDomain  (ims:ime,jms:jme,ifs:ife))
      
      allocate(iCenCell  (ims:ime,jms:jme,ifs:ife))
      
      allocate(iRecCell  (maxRecCells,ims:ime,jms:jme,ifs:ife))
      allocate(jRecCell  (maxRecCells,ims:ime,jms:jme,ifs:ife))
      
      allocate(iGstCell  (maxRecCells,ims:ime,jms:jme,ifs:ife))
      allocate(jGstCell  (maxRecCells,ims:ime,jms:jme,ifs:ife))
      
      allocate(xRel(4,maxRecCells,ids:ide,jds:jde,ifs:ife))
      allocate(yRel(4,maxRecCells,ids:ide,jds:jde,ifs:ife))
      
      allocate(xGst(4,maxRecCells,ids:ide,jds:jde,ifs:ife))
      allocate(yGst(4,maxRecCells,ids:ide,jds:jde,ifs:ife))
      
      allocate(recMatrixL(nPointsOnEdge,maxRecTerms,ids:ide,jds:jde,ifs:ife))
      allocate(recMatrixR(nPointsOnEdge,maxRecTerms,ids:ide,jds:jde,ifs:ife))
      allocate(recMatrixB(nPointsOnEdge,maxRecTerms,ids:ide,jds:jde,ifs:ife))
      allocate(recMatrixT(nPointsOnEdge,maxRecTerms,ids:ide,jds:jde,ifs:ife))
      
      allocate(A   (maxRecCells,maxRecCells,ids:ide,jds:jde,ifs:ife))
      allocate(invA(maxRecCells,maxRecCells,ids:ide,jds:jde,ifs:ife))
      
      allocate(gstMatrixL(nPointsOnEdge,maxRecCells,ids:ide,jds:jde,ifs:ife))
      allocate(gstMatrixR(nPointsOnEdge,maxRecCells,ids:ide,jds:jde,ifs:ife))
      allocate(gstMatrixB(nPointsOnEdge,maxRecCells,ids:ide,jds:jde,ifs:ife))
      allocate(gstMatrixT(nPointsOnEdge,maxRecCells,ids:ide,jds:jde,ifs:ife))
      
      allocate(qC(nVar,              ims:ime,jms:jme,ifs:ife))
      allocate(qL(nVar,nPointsOnEdge,ims:ime,jms:jme,ifs:ife))
      allocate(qR(nVar,nPointsOnEdge,ims:ime,jms:jme,ifs:ife))
      allocate(qB(nVar,nPointsOnEdge,ims:ime,jms:jme,ifs:ife))
      allocate(qT(nVar,nPointsOnEdge,ims:ime,jms:jme,ifs:ife))
      
      allocate(FL(nVar,nPointsOnEdge,ims:ime,jms:jme,ifs:ife))
      allocate(FR(nVar,nPointsOnEdge,ims:ime,jms:jme,ifs:ife))
      allocate(GB(nVar,nPointsOnEdge,ims:ime,jms:jme,ifs:ife))
      allocate(GT(nVar,nPointsOnEdge,ims:ime,jms:jme,ifs:ife))
      
      allocate(Fe(nVar,ids:ide+1,jds:jde  ,ifs:ife))
      allocate(Ge(nVar,ids:ide  ,jds:jde+1,ifs:ife))
      
      allocate(FeP(nVar,nPointsOnEdge,ids:ide+1,jds:jde  ,ifs:ife))
      allocate(GeP(nVar,nPointsOnEdge,ids:ide  ,jds:jde+1,ifs:ife))
      
      allocate(xL(nPointsOnEdge))
      allocate(xR(nPointsOnEdge))
      allocate(xB(nPointsOnEdge))
      allocate(xT(nPointsOnEdge))
      
      allocate(yL(nPointsOnEdge))
      allocate(yR(nPointsOnEdge))
      allocate(yB(nPointsOnEdge))
      allocate(yT(nPointsOnEdge))
      
      dV = dx * dy
    
      ! set reconstruction cells on each stencil
      print*,'Set reconstruction cells on each stencil'
      print*,''
      inDomain(ims:ime,jms:jme,ifs:ife) = .false. 
      inDomain(ids:ide,jds:jde,ifs:ife) = .true.
      
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            ! Reconstruction cells
            iCOS = 0
            do jRec = -recBdy,recBdy
              do iRec = -recBdy,recBdy
                iCOS = iCOS + 1
                iRecCell(iCOS,i,j,iPatch) = i + iRec
                jRecCell(iCOS,i,j,iPatch) = j + jRec
              enddo
            enddo
            nRecCells(i,j,iPatch) = iCOS
            
            locPolyDegree(i,j,iPatch) = min( maxval(iRecCell(1:iCOS,i,j,iPatch)) - minval(iRecCell(1:iCOS,i,j,iPatch)), maxval(jRecCell(1:iCOS,i,j,iPatch)) - minval(jRecCell(1:iCOS,i,j,iPatch)) )
            
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
            nGstCells(i,iCOS,iPatch) = iCOS
            
          enddo
        enddo
      enddo
      
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            nRecTerms(i,j,iPatch) = ( locPolyDegree(i,j,iPatch) + 1 ) * ( locPolyDegree(i,j,iPatch) + 2 ) / 2
          enddo
        enddo
      enddo
      
      ! initilize polyCoordCoef
      print*,'Initilize polyCoordCoef'
      print*,''
      
      xL = - dx / 2.
      xR =   dx / 2.
      do iQP = 1,nPointsOnEdge
        xB(iQP) = xL(1) + dx * quad_pos_1d(iQP)
        xT(iQP) = xB(iQP)
      enddo
      
      yB = - dy / 2.
      yT =   dy / 2.
      do jQP = 1,nPointsOnEdge
        yL(jQP) = yB(1) + dy * quad_pos_1d(jQP)
        yR(jQP) = yL(jQP)
      enddo
      
      recCoef = 1.
      recdx   = 1. / ( dx * recCoef )
      recdy   = 1. / ( dy * recCoef )
      recdV   = 1. / ( recCoef**2 )
      
      !$OMP PARALLEL DO PRIVATE(i,j,iCOS,iRec,jRec)
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            do iCOS = 1,nRecCells(i,j,iPatch)
              iRec = iRecCell(iCOS,i,j,iPatch)
              jRec = jRecCell(iCOS,i,j,iPatch)
              
              if(iRec==i.and.jRec==j)iCenCell(i,j,iPatch) = iCOS
              
              xRel(:,iCOS,i,j,iPatch) = ( x(ccs:cce,iRec,jRec,iPatch) - x(cc,i,j,iPatch) ) * recdx
              yRel(:,iCOS,i,j,iPatch) = ( y(ccs:cce,iRec,jRec,iPatch) - y(cc,i,j,iPatch) ) * recdy
              
              call calc_polynomial_square_integration(locPolyDegree(i,j,iPatch),xRel(1,iCOS,i,j,iPatch),xRel(2,iCOS,i,j,iPatch),&
                                                                                yRel(1,iCOS,i,j,iPatch),yRel(4,iCOS,i,j,iPatch),polyCoordCoef(iCOS,1:nRecTerms(i,j,iPatch),i,j,iPatch))
            enddo
            
            ! Calculate reconstruction matrix on edge
            call calc_polynomial_matrix(locPolyDegree(i,j,iPatch),nPointsOnEdge,nRecTerms(i,j,iPatch),xL*recdx,yL*recdy,recMatrixL(:,1:nRecTerms(i,j,iPatch),i,j,iPatch))
            call calc_polynomial_matrix(locPolyDegree(i,j,iPatch),nPointsOnEdge,nRecTerms(i,j,iPatch),xR*recdx,yR*recdy,recMatrixR(:,1:nRecTerms(i,j,iPatch),i,j,iPatch))
            call calc_polynomial_matrix(locPolyDegree(i,j,iPatch),nPointsOnEdge,nRecTerms(i,j,iPatch),xB*recdx,yB*recdy,recMatrixB(:,1:nRecTerms(i,j,iPatch),i,j,iPatch))
            call calc_polynomial_matrix(locPolyDegree(i,j,iPatch),nPointsOnEdge,nRecTerms(i,j,iPatch),xT*recdx,yT*recdy,recMatrixT(:,1:nRecTerms(i,j,iPatch),i,j,iPatch))
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      !$OMP PARALLEL DO PRIVATE(i,j,nRC,nxp,nyp,iCOS,jR,iR,invstat)
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            nRC = nGstCells(i,j,iPatch)
            !nxp = min( maxval(iGstCell(1:nRC,i,j,iPatch)) - minval(iGstCell(1:nRC,i,j,iPatch)) + 1, stencil_width )
            !nyp = min( maxval(jGstCell(1:nRC,i,j,iPatch)) - minval(jGstCell(1:nRC,i,j,iPatch)) + 1, stencil_width )
            nxp = maxval(iGstCell(1:nRC,i,j,iPatch)) - minval(iGstCell(1:nRC,i,j,iPatch)) + 1
            nyp = maxval(jGstCell(1:nRC,i,j,iPatch)) - minval(jGstCell(1:nRC,i,j,iPatch)) + 1
            
            do iCOS = 1,nGstCells(i,j,iPatch)
              iRec = iGstCell(iCOS,i,j,iPatch)
              jRec = jGstCell(iCOS,i,j,iPatch)
              
              xGst(:,iCOS,i,j,iPatch) = ( x(ccs:cce,iRec,jRec,iPatch) - x(cc,i,j,iPatch) ) * recdx
              yGst(:,iCOS,i,j,iPatch) = ( y(ccs:cce,iRec,jRec,iPatch) - y(cc,i,j,iPatch) ) * recdy
            enddo
            
            iCOS = 0
            do jR = 1,nyp
              do iR = 1,nxp
                iCOS = iCOS + 1
                call calc_rectangle_poly_integration(nxp,nyp,xGst(1,iCOS,i,j,iPatch),xGst(2,iCOS,i,j,iPatch),&
                                                             yGst(1,iCOS,i,j,iPatch),yGst(4,iCOS,i,j,iPatch),A(iCOS,1:nRC,i,j,iPatch))
              enddo
            enddo
            
            call BRINV(nRC,A(1:nRC,1:nRC,i,j,iPatch),invA(1:nRC,1:nRC,i,j,iPatch),invstat)
            if(invstat==0)then
              print*,'Inverse A dost not exist'
              print*,'i,j,iPatch,nRC,nxp,nyp are'
              print*,i,j,iPatch,nRC,nxp,nyp
              stop 'Check BRINV for Special treamtment on boundary cells'
            endif
            
            ! Calculate reconstruction matrix on edge
            call calc_rectangle_poly_matrix(nxp,nyp,nPointsOnEdge,xL*recdx,yL*recdy,gstMatrixL(:,1:nRC,i,j,iPatch)) ! need to be modified to ghost position
            call calc_rectangle_poly_matrix(nxp,nyp,nPointsOnEdge,xR*recdx,yR*recdy,gstMatrixR(:,1:nRC,i,j,iPatch)) ! need to be modified to ghost position
            call calc_rectangle_poly_matrix(nxp,nyp,nPointsOnEdge,xB*recdx,yB*recdy,gstMatrixB(:,1:nRC,i,j,iPatch)) ! need to be modified to ghost position
            call calc_rectangle_poly_matrix(nxp,nyp,nPointsOnEdge,xT*recdx,yT*recdy,gstMatrixT(:,1:nRC,i,j,iPatch)) ! need to be modified to ghost position
            
            gstMatrixL(:,1:nRC,i,j,iPatch) = matmul( gstMatrixL(:,1:nRC,i,j,iPatch), invA(1:nRC,1:nRC,i,j,iPatch) )
            gstMatrixR(:,1:nRC,i,j,iPatch) = matmul( gstMatrixR(:,1:nRC,i,j,iPatch), invA(1:nRC,1:nRC,i,j,iPatch) )
            gstMatrixB(:,1:nRC,i,j,iPatch) = matmul( gstMatrixB(:,1:nRC,i,j,iPatch), invA(1:nRC,1:nRC,i,j,iPatch) )
            gstMatrixT(:,1:nRC,i,j,iPatch) = matmul( gstMatrixT(:,1:nRC,i,j,iPatch), invA(1:nRC,1:nRC,i,j,iPatch) )
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
    end subroutine init_spatial_operator
    
    subroutine spatial_operator(stat,tend)
      type(stat_field), target, intent(in   ) :: stat
      type(tend_field), target, intent(inout) :: tend
      
      
    end subroutine spatial_operator
    
end module spatial_operators_mod

