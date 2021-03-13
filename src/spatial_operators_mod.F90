module spatial_operators_mod
  use constants_mod
  use mesh_mod
  use parameters_mod
  use stat_mod
  use tend_mod
  use projection_mod
  use reconstruction_mod
  use math_mod
  implicit none
  
  private
  
  public init_spatial_operator,spatial_operator
  
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
  
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: recMatrixQ ! reconstruction for gaussian quadrature points
  
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: recMatrixDx
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: recMatrixDy
  
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: polyMatrixL
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: polyMatrixR
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: polyMatrixB
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: polyMatrixT
  
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: Apoly
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: invApoly
  
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: A
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: invA
  
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: gstMatrix
  
  !real(r_kind) :: recCoef
  !real(r_kind) :: recdx
  !real(r_kind) :: recdy
  !real(r_kind) :: recdV
  real(r_kind) :: dV
  
  real   (r_kind), dimension(:,:,:,:  ), allocatable :: qC
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: qL
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: qR
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: qB
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: qT
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: qQ ! Just for phit quadrature on cell
  
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
  
  real(r_kind), dimension(:,:,:,:), allocatable :: src   ! source term
  
  real(r_kind), dimension(:,:,:,:), allocatable :: ghs    ! sqrtG * zs * gravity
  real(r_kind), dimension(:,:,:  ), allocatable :: ghsC   ! sqrtG * zs * gravity on Cell
  
  real(r_kind), dimension(:,:,:,:), allocatable :: phit    ! phi + phis
  real(r_kind), dimension(:,:,:  ), allocatable :: phitC   ! phi + phis on cell
  
  real(r_kind), dimension(:,:,:,:), allocatable :: dphitdx    ! d phi_t / dx
  real(r_kind), dimension(:,:,:,:), allocatable :: dphitdy    ! d phi_t / dy
  
    contains
    subroutine init_spatial_operator
      integer(i_kind) :: i,j,iPatch
      integer(i_kind) :: iCOS ! indices of Cells On Stencils
      integer(i_kind) :: iPOC ! indices of points on cell
      integer(i_kind) :: iRec,jRec
      integer(i_kind) :: iQP,jQP
      integer(i_kind) :: iR,jR
      
      integer(i_kind) :: nxp,nyp
      integer(i_kind) :: invstat
      integer(i_kind) :: nRC,nRT
      
      integer :: pg
  
      real(r_kind), dimension(:), allocatable :: xL
      real(r_kind), dimension(:), allocatable :: xR
      real(r_kind), dimension(:), allocatable :: xB
      real(r_kind), dimension(:), allocatable :: xT
      
      real(r_kind), dimension(:), allocatable :: yL
      real(r_kind), dimension(:), allocatable :: yR
      real(r_kind), dimension(:), allocatable :: yB
      real(r_kind), dimension(:), allocatable :: yT
      
      real(r_kind), dimension(:), allocatable :: xg
      real(r_kind), dimension(:), allocatable :: yg
      
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
      
      allocate(recMatrixQ(nQuadPointsOnCell,maxRecTerms,ids:ide,jds:jde,ifs:ife))

      allocate(recMatrixDx(nQuadPointsOnCell,maxRecTerms,ids:ide,jds:jde,ifs:ife))
      allocate(recMatrixDy(nQuadPointsOnCell,maxRecTerms,ids:ide,jds:jde,ifs:ife))
      
      allocate(polyMatrixL(nPointsOnEdge,maxRecTerms,ids:ide,jds:jde,ifs:ife))
      allocate(polyMatrixR(nPointsOnEdge,maxRecTerms,ids:ide,jds:jde,ifs:ife))
      allocate(polyMatrixB(nPointsOnEdge,maxRecTerms,ids:ide,jds:jde,ifs:ife))
      allocate(polyMatrixT(nPointsOnEdge,maxRecTerms,ids:ide,jds:jde,ifs:ife))
      
      allocate(Apoly   (maxRecCells,maxRecCells,ids:ide,jds:jde,ifs:ife))
      allocate(invApoly(maxRecCells,maxRecCells,ids:ide,jds:jde,ifs:ife))
      
      allocate(A   (maxRecCells,maxRecCells,ids:ide,jds:jde,ifs:ife))
      allocate(invA(maxRecCells,maxRecCells,ids:ide,jds:jde,ifs:ife))
      
      allocate(gstMatrix(maxGhostPointsOnCell,maxRecCells,ids:ide,jds:jde,ifs:ife))
      
      allocate(qC(nVar,                  ims:ime,jms:jme,ifs:ife))
      allocate(qL(nVar,nPointsOnEdge    ,ims:ime,jms:jme,ifs:ife))
      allocate(qR(nVar,nPointsOnEdge    ,ims:ime,jms:jme,ifs:ife))
      allocate(qB(nVar,nPointsOnEdge    ,ims:ime,jms:jme,ifs:ife))
      allocate(qT(nVar,nPointsOnEdge    ,ims:ime,jms:jme,ifs:ife))
      allocate(qQ(nVar,nQuadPointsOnCell,ims:ime,jms:jme,ifs:ife))
      
      allocate(FL(nVar,nPointsOnEdge,ims:ime,jms:jme,ifs:ife))
      allocate(FR(nVar,nPointsOnEdge,ims:ime,jms:jme,ifs:ife))
      allocate(GB(nVar,nPointsOnEdge,ims:ime,jms:jme,ifs:ife))
      allocate(GT(nVar,nPointsOnEdge,ims:ime,jms:jme,ifs:ife))
      
      allocate(Fe(nVar,ids:ide+1,jds:jde  ,ifs:ife))
      allocate(Ge(nVar,ids:ide  ,jds:jde+1,ifs:ife))
      
      allocate(FeP(nVar,nPointsOnEdge,ids:ide+1,jds:jde  ,ifs:ife))
      allocate(GeP(nVar,nPointsOnEdge,ids:ide  ,jds:jde+1,ifs:ife))
      
      allocate(src(nVar,ids:ide,jds:jde,ifs:ife))
      
      allocate(ghs (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
      allocate(ghsC(              ims:ime,jms:jme,ifs:ife))
      
      allocate(phit (nQuadPointsOnCell,ims:ime,jms:jme,ifs:ife))
      allocate(phitC(                  ims:ime,jms:jme,ifs:ife))
      
      allocate(dphitdx(nQuadPointsOnCell,ims:ime,jms:jme,ifs:ife))
      allocate(dphitdy(nQuadPointsOnCell,ims:ime,jms:jme,ifs:ife))
      
      allocate(xL(nPointsOnEdge))
      allocate(xR(nPointsOnEdge))
      allocate(xB(nPointsOnEdge))
      allocate(xT(nPointsOnEdge))
      
      allocate(yL(nPointsOnEdge))
      allocate(yR(nPointsOnEdge))
      allocate(yB(nPointsOnEdge))
      allocate(yT(nPointsOnEdge))
      
      allocate(xg(nPointsOnCell))
      allocate(yg(nPointsOnCell))
      
      dV = dx * dy
    
      ! set reconstruction cells on each stencil
      print*,'Set reconstruction cells on each stencil'
      print*,''
      
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
            nGstRecCells(i,j,iPatch) = iCOS
            
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
      
      !$OMP PARALLEL DO PRIVATE(i,j,iCOS,iRec,jRec,nRC,nRT) COLLAPSE(3)
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            nRC = nRecCells(i,j,iPatch)
            nRT = nRecTerms(i,j,iPatch)
            do iCOS = 1,nRC
              iRec = iRecCell(iCOS,i,j,iPatch)
              jRec = jRecCell(iCOS,i,j,iPatch)
              
              if(iRec==i.and.jRec==j)iCenCell(i,j,iPatch) = iCOS
              
              xRel(:,iCOS,i,j,iPatch) = ( x(ccs:cce,iRec,jRec,iPatch) - x(cc,i,j,iPatch) ) * recdx
              yRel(:,iCOS,i,j,iPatch) = ( y(ccs:cce,iRec,jRec,iPatch) - y(cc,i,j,iPatch) ) * recdy
              
              call calc_polynomial_square_integration(locPolyDegree(i,j,iPatch),xRel(1,iCOS,i,j,iPatch),xRel(2,iCOS,i,j,iPatch),&
                                                                                yRel(1,iCOS,i,j,iPatch),yRel(4,iCOS,i,j,iPatch),polyCoordCoef(iCOS,1:nRT,i,j,iPatch))
            enddo
            
            ! Calculate reconstruction matrix on edge
            call calc_polynomial_matrix(locPolyDegree(i,j,iPatch),nPointsOnEdge,nRT,xL*recdx,yL*recdy,recMatrixL(:,1:nRT,i,j,iPatch))
            call calc_polynomial_matrix(locPolyDegree(i,j,iPatch),nPointsOnEdge,nRT,xR*recdx,yR*recdy,recMatrixR(:,1:nRT,i,j,iPatch))
            call calc_polynomial_matrix(locPolyDegree(i,j,iPatch),nPointsOnEdge,nRT,xB*recdx,yB*recdy,recMatrixB(:,1:nRT,i,j,iPatch))
            call calc_polynomial_matrix(locPolyDegree(i,j,iPatch),nPointsOnEdge,nRT,xT*recdx,yT*recdy,recMatrixT(:,1:nRT,i,j,iPatch))
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      !$OMP PARALLEL DO PRIVATE(j,i,iPOC,xg,yg,nRC,nRT) COLLAPSE(3)
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            nRC = nRecCells(i,j,iPatch)
            nRT = nRecTerms(i,j,iPatch)
            do iPOC = 1,nQuadPointsOncell
              !xg(iPOC) = ( x(cqs+iPOC-1,i,j,iPatch) - x(cc,i,j,iPatch) ) * recdx
              !yg(iPOC) = ( y(cqs+iPOC-1,i,j,iPatch) - y(cc,i,j,iPatch) ) * recdy
              xg(iPOC) = x(cqs+iPOC-1,i,j,iPatch) - x(cc,i,j,iPatch)
              yg(iPOC) = y(cqs+iPOC-1,i,j,iPatch) - y(cc,i,j,iPatch)
            enddo
            
            call calc_polynomial_deriv_matrix(locPolyDegree(i,j,iPatch),nQuadPointsOnCell,nRT,&
                                              xg(1:nQuadPointsOnCell),yg(1:nQuadPointsOnCell),&
                                              recMatrixDx(:,1:nRT,i,j,iPatch)                ,&
                                              recMatrixDy(:,1:nRT,i,j,iPatch))
            
            call calc_polynomial_matrix(locPolyDegree(i,j,iPatch),nQuadPointsOnCell,nRT,xg(1:nQuadPointsOnCell)*recdx,yg(1:nQuadPointsOnCell)*recdy,recMatrixQ(:,1:nRT,i,j,iPatch))
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
            
      !$OMP PARALLEL DO PRIVATE(i,j,nRC,nxp,nyp,iCOS,jR,iR,iRec,jRec,invstat) COLLAPSE(3)
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            nRC = nRecCells(i,j,iPatch)
            nxp = maxval(iRecCell(1:nRC,i,j,iPatch)) - minval(iRecCell(1:nRC,i,j,iPatch)) + 1
            nyp = maxval(jRecCell(1:nRC,i,j,iPatch)) - minval(jRecCell(1:nRC,i,j,iPatch)) + 1
            
            iCOS = 0
            do jR = 1,nyp
              do iR = 1,nxp
                iCOS = iCOS + 1
                call calc_rectangle_poly_integration(nxp,nyp,xRel(1,iCOS,i,j,iPatch),xRel(2,iCOS,i,j,iPatch),&
                                                             yRel(1,iCOS,i,j,iPatch),yRel(4,iCOS,i,j,iPatch),Apoly(iCOS,1:nRC,i,j,iPatch))
              enddo
            enddo
            
            call BRINV(nRC,Apoly(1:nRC,1:nRC,i,j,iPatch),invApoly(1:nRC,1:nRC,i,j,iPatch),invstat)
            if(invstat==0)then
              print*,'Inverse Apoly dost not exist'
              print*,'i,j,iPatch,nRC,nxp,nyp are'
              print*,i,j,iPatch,nRC,nxp,nyp
              stop 'Check BRINV for Special treamtment on boundary cells'
            endif
            
            ! Calculate reconstruction matrix on edge
            call calc_rectangle_poly_matrix(nxp,nyp,nPointsOnEdge,xL*recdx,yL*recdy,polyMatrixL(:,1:nRC,i,j,iPatch))
            call calc_rectangle_poly_matrix(nxp,nyp,nPointsOnEdge,xR*recdx,yR*recdy,polyMatrixR(:,1:nRC,i,j,iPatch))
            call calc_rectangle_poly_matrix(nxp,nyp,nPointsOnEdge,xB*recdx,yB*recdy,polyMatrixB(:,1:nRC,i,j,iPatch))
            call calc_rectangle_poly_matrix(nxp,nyp,nPointsOnEdge,xT*recdx,yT*recdy,polyMatrixT(:,1:nRC,i,j,iPatch))
            
            polyMatrixL(:,1:nRC,i,j,iPatch) = matmul( polyMatrixL(:,1:nRC,i,j,iPatch), invApoly(1:nRC,1:nRC,i,j,iPatch) )
            polyMatrixR(:,1:nRC,i,j,iPatch) = matmul( polyMatrixR(:,1:nRC,i,j,iPatch), invApoly(1:nRC,1:nRC,i,j,iPatch) )
            polyMatrixB(:,1:nRC,i,j,iPatch) = matmul( polyMatrixB(:,1:nRC,i,j,iPatch), invApoly(1:nRC,1:nRC,i,j,iPatch) )
            polyMatrixT(:,1:nRC,i,j,iPatch) = matmul( polyMatrixT(:,1:nRC,i,j,iPatch), invApoly(1:nRC,1:nRC,i,j,iPatch) )
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      !$OMP PARALLEL DO PRIVATE(i,j,nRC,nxp,nyp,iCOS,jR,iR,iRec,jRec,invstat,xg,yg,pg) COLLAPSE(3)
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            nRC = nGstRecCells(i,j,iPatch)
            pg  = nGhostPointsOnCell(i,j,iPatch)
            if(pg/=0)then
              !nxp = min( maxval(iGstCell(1:nRC,i,j,iPatch)) - minval(iGstCell(1:nRC,i,j,iPatch)) + 1, stencil_width )
              !nyp = min( maxval(jGstCell(1:nRC,i,j,iPatch)) - minval(jGstCell(1:nRC,i,j,iPatch)) + 1, stencil_width )
              nxp = maxval(iGstCell(1:nRC,i,j,iPatch)) - minval(iGstCell(1:nRC,i,j,iPatch)) + 1
              nyp = maxval(jGstCell(1:nRC,i,j,iPatch)) - minval(jGstCell(1:nRC,i,j,iPatch)) + 1
              
              do iCOS = 1,nRC
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
              xg(1:pg) = ( x(cgs:cgs+pg-1,i,j,iPatch) - x(cc,i,j,iPatch) ) * recdx
              yg(1:pg) = ( y(cgs:cgs+pg-1,i,j,iPatch) - y(cc,i,j,iPatch) ) * recdy
              call calc_rectangle_poly_matrix(nxp,nyp,pg,xg,yg,gstMatrix(1:pg,1:nRC,i,j,iPatch))
              
              gstMatrix(1:pg,1:nRC,i,j,iPatch) = matmul( gstMatrix(1:pg,1:nRC,i,j,iPatch), invA(1:nRC,1:nRC,i,j,iPatch) )
            endif
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      ghs = zs * gravity
      
      do iPatch = ifs,ife
        do j = jms,jme
          do i = ims,ime
            ghsC(i,j,iPatch) = cell_quadrature(ghs(cqs:cqe,i,j,iPatch))
          enddo
        enddo
      enddo
      
      ! Calculate metric tensor for cube boundary condition
      sqrtGL_adj    = sqrtGL
      sqrtGR_adj    = sqrtGR
      sqrtGB_adj    = sqrtGB
      sqrtGT_adj    = sqrtGT
      
      matrixAL_adj  = matrixAL
      matrixAR_adj  = matrixAR
      matrixAB_adj  = matrixAB
      matrixAT_adj  = matrixAT
      
      matrixIAL_adj = matrixIAL
      matrixIAR_adj = matrixIAR
      matrixIAB_adj = matrixIAB
      matrixIAT_adj = matrixIAT
      
      ! Panel 1
      call restore_bdy_field(sqrtGR_adj(:,ids-1,jds:jde,1),sqrtGR(:,ide,jds:jde,4), 1) ! Left
      call restore_bdy_field(sqrtGL_adj(:,ide+1,jds:jde,1),sqrtGL(:,ids,jds:jde,2), 1) ! Right
      call restore_bdy_field(sqrtGT_adj(:,ids:ide,jds-1,1),sqrtGT(:,ids:ide,jde,6), 1) ! below
      call restore_bdy_field(sqrtGB_adj(:,ids:ide,jde+1,1),sqrtGB(:,ids:ide,jds,5), 1) ! over
      ! Panel 2
      call restore_bdy_field(sqrtGR_adj(:,ids-1,jds:jde,2),sqrtGR(:,ide,jds:jde,1), 1) ! Left
      call restore_bdy_field(sqrtGL_adj(:,ide+1,jds:jde,2),sqrtGL(:,ids,jds:jde,3), 1) ! Right
      call restore_bdy_field(sqrtGT_adj(:,ids:ide,jds-1,2),sqrtGR(:,ide,jds:jde,6),-1) ! below
      call restore_bdy_field(sqrtGB_adj(:,ids:ide,jde+1,2),sqrtGR(:,ide,jds:jde,5), 1) ! over
      ! Panel 3
      call restore_bdy_field(sqrtGR_adj(:,ids-1,jds:jde,3),sqrtGR(:,ide,jds:jde,2), 1) ! Left
      call restore_bdy_field(sqrtGL_adj(:,ide+1,jds:jde,3),sqrtGL(:,ids,jds:jde,4), 1) ! Right
      call restore_bdy_field(sqrtGT_adj(:,ids:ide,jds-1,3),sqrtGB(:,ids:ide,jds,6),-1) ! below
      call restore_bdy_field(sqrtGB_adj(:,ids:ide,jde+1,3),sqrtGT(:,ids:ide,jde,5),-1) ! over
      ! Panel 4
      call restore_bdy_field(sqrtGR_adj(:,ids-1,jds:jde,4),sqrtGR(:,ide,jds:jde,3), 1) ! Left
      call restore_bdy_field(sqrtGL_adj(:,ide+1,jds:jde,4),sqrtGL(:,ids,jds:jde,1), 1) ! Right
      call restore_bdy_field(sqrtGT_adj(:,ids:ide,jds-1,4),sqrtGL(:,ids,jds:jde,6), 1) ! below
      call restore_bdy_field(sqrtGB_adj(:,ids:ide,jde+1,4),sqrtGL(:,ids,jds:jde,5),-1) ! over
      ! Panel 5
      call restore_bdy_field(sqrtGR_adj(:,ids-1,jds:jde,5),sqrtGT(:,ids:ide,jde,4),-1) ! Left
      call restore_bdy_field(sqrtGL_adj(:,ide+1,jds:jde,5),sqrtGT(:,ids:ide,jde,2), 1) ! Right
      call restore_bdy_field(sqrtGT_adj(:,ids:ide,jds-1,5),sqrtGT(:,ids:ide,jde,1), 1) ! below
      call restore_bdy_field(sqrtGB_adj(:,ids:ide,jde+1,5),sqrtGT(:,ids:ide,jde,3),-1) ! over
      ! Panel 6
      call restore_bdy_field(sqrtGR_adj(:,ids-1,jds:jde,6),sqrtGB(:,ids:ide,jds,4), 1) ! Left
      call restore_bdy_field(sqrtGL_adj(:,ide+1,jds:jde,6),sqrtGB(:,ids:ide,jds,2),-1) ! Right
      call restore_bdy_field(sqrtGT_adj(:,ids:ide,jds-1,6),sqrtGB(:,ids:ide,jds,3),-1) ! below
      call restore_bdy_field(sqrtGB_adj(:,ids:ide,jde+1,6),sqrtGB(:,ids:ide,jds,1), 1) ! over
      
      do jQP = 1,2
        do iQP = 1,2
          ! MatrixA
          ! Panel 1
          call restore_bdy_field(matrixAR_adj(iQP,jQP,:,ids-1,jds:jde,1),matrixAR(iQP,jQP,:,ide,jds:jde,4), 1) ! Left
          call restore_bdy_field(matrixAL_adj(iQP,jQP,:,ide+1,jds:jde,1),matrixAL(iQP,jQP,:,ids,jds:jde,2), 1) ! Right
          call restore_bdy_field(matrixAT_adj(iQP,jQP,:,ids:ide,jds-1,1),matrixAT(iQP,jQP,:,ids:ide,jde,6), 1) ! below
          call restore_bdy_field(matrixAB_adj(iQP,jQP,:,ids:ide,jde+1,1),matrixAB(iQP,jQP,:,ids:ide,jds,5), 1) ! over
          ! Panel 2
          call restore_bdy_field(matrixAR_adj(iQP,jQP,:,ids-1,jds:jde,2),matrixAR(iQP,jQP,:,ide,jds:jde,1), 1) ! Left
          call restore_bdy_field(matrixAL_adj(iQP,jQP,:,ide+1,jds:jde,2),matrixAL(iQP,jQP,:,ids,jds:jde,3), 1) ! Right
          call restore_bdy_field(matrixAT_adj(iQP,jQP,:,ids:ide,jds-1,2),matrixAR(iQP,jQP,:,ide,jds:jde,6),-1) ! below
          call restore_bdy_field(matrixAB_adj(iQP,jQP,:,ids:ide,jde+1,2),matrixAR(iQP,jQP,:,ide,jds:jde,5), 1) ! over
          ! Panel 3
          call restore_bdy_field(matrixAR_adj(iQP,jQP,:,ids-1,jds:jde,3),matrixAR(iQP,jQP,:,ide,jds:jde,2), 1) ! Left
          call restore_bdy_field(matrixAL_adj(iQP,jQP,:,ide+1,jds:jde,3),matrixAL(iQP,jQP,:,ids,jds:jde,4), 1) ! Right
          call restore_bdy_field(matrixAT_adj(iQP,jQP,:,ids:ide,jds-1,3),matrixAB(iQP,jQP,:,ids:ide,jds,6),-1) ! below
          call restore_bdy_field(matrixAB_adj(iQP,jQP,:,ids:ide,jde+1,3),matrixAT(iQP,jQP,:,ids:ide,jde,5),-1) ! over
          ! Panel 4
          call restore_bdy_field(matrixAR_adj(iQP,jQP,:,ids-1,jds:jde,4),matrixAR(iQP,jQP,:,ide,jds:jde,3), 1) ! Left
          call restore_bdy_field(matrixAL_adj(iQP,jQP,:,ide+1,jds:jde,4),matrixAL(iQP,jQP,:,ids,jds:jde,1), 1) ! Right
          call restore_bdy_field(matrixAT_adj(iQP,jQP,:,ids:ide,jds-1,4),matrixAL(iQP,jQP,:,ids,jds:jde,6), 1) ! below
          call restore_bdy_field(matrixAB_adj(iQP,jQP,:,ids:ide,jde+1,4),matrixAL(iQP,jQP,:,ids,jds:jde,5),-1) ! over
          ! Panel 5
          call restore_bdy_field(matrixAR_adj(iQP,jQP,:,ids-1,jds:jde,5),matrixAT(iQP,jQP,:,ids:ide,jde,4),-1) ! Left
          call restore_bdy_field(matrixAL_adj(iQP,jQP,:,ide+1,jds:jde,5),matrixAT(iQP,jQP,:,ids:ide,jde,2), 1) ! Right
          call restore_bdy_field(matrixAT_adj(iQP,jQP,:,ids:ide,jds-1,5),matrixAT(iQP,jQP,:,ids:ide,jde,1), 1) ! below
          call restore_bdy_field(matrixAB_adj(iQP,jQP,:,ids:ide,jde+1,5),matrixAT(iQP,jQP,:,ids:ide,jde,3),-1) ! over
          ! Panel 6
          call restore_bdy_field(matrixAR_adj(iQP,jQP,:,ids-1,jds:jde,6),matrixAB(iQP,jQP,:,ids:ide,jds,4), 1) ! Left
          call restore_bdy_field(matrixAL_adj(iQP,jQP,:,ide+1,jds:jde,6),matrixAB(iQP,jQP,:,ids:ide,jds,2),-1) ! Right
          call restore_bdy_field(matrixAT_adj(iQP,jQP,:,ids:ide,jds-1,6),matrixAB(iQP,jQP,:,ids:ide,jds,3),-1) ! below
          call restore_bdy_field(matrixAB_adj(iQP,jQP,:,ids:ide,jde+1,6),matrixAB(iQP,jQP,:,ids:ide,jds,1), 1) ! over
          
          ! MatrixIA
          ! Panel 1
          call restore_bdy_field(matrixIAR_adj(iQP,jQP,:,ids-1,jds:jde,1),matrixIAR(iQP,jQP,:,ide,jds:jde,4), 1) ! Left
          call restore_bdy_field(matrixIAL_adj(iQP,jQP,:,ide+1,jds:jde,1),matrixIAL(iQP,jQP,:,ids,jds:jde,2), 1) ! Right
          call restore_bdy_field(matrixIAT_adj(iQP,jQP,:,ids:ide,jds-1,1),matrixIAT(iQP,jQP,:,ids:ide,jde,6), 1) ! below
          call restore_bdy_field(matrixIAB_adj(iQP,jQP,:,ids:ide,jde+1,1),matrixIAB(iQP,jQP,:,ids:ide,jds,5), 1) ! over
          ! Panel 2
          call restore_bdy_field(matrixIAR_adj(iQP,jQP,:,ids-1,jds:jde,2),matrixIAR(iQP,jQP,:,ide,jds:jde,1), 1) ! Left
          call restore_bdy_field(matrixIAL_adj(iQP,jQP,:,ide+1,jds:jde,2),matrixIAL(iQP,jQP,:,ids,jds:jde,3), 1) ! Right
          call restore_bdy_field(matrixIAT_adj(iQP,jQP,:,ids:ide,jds-1,2),matrixIAR(iQP,jQP,:,ide,jds:jde,6),-1) ! below
          call restore_bdy_field(matrixIAB_adj(iQP,jQP,:,ids:ide,jde+1,2),matrixIAR(iQP,jQP,:,ide,jds:jde,5), 1) ! over
          ! Panel 3
          call restore_bdy_field(matrixIAR_adj(iQP,jQP,:,ids-1,jds:jde,3),matrixIAR(iQP,jQP,:,ide,jds:jde,2), 1) ! Left
          call restore_bdy_field(matrixIAL_adj(iQP,jQP,:,ide+1,jds:jde,3),matrixIAL(iQP,jQP,:,ids,jds:jde,4), 1) ! Right
          call restore_bdy_field(matrixIAT_adj(iQP,jQP,:,ids:ide,jds-1,3),matrixIAB(iQP,jQP,:,ids:ide,jds,6),-1) ! below
          call restore_bdy_field(matrixIAB_adj(iQP,jQP,:,ids:ide,jde+1,3),matrixIAT(iQP,jQP,:,ids:ide,jde,5),-1) ! over
          ! Panel 4
          call restore_bdy_field(matrixIAR_adj(iQP,jQP,:,ids-1,jds:jde,4),matrixIAR(iQP,jQP,:,ide,jds:jde,3), 1) ! Left
          call restore_bdy_field(matrixIAL_adj(iQP,jQP,:,ide+1,jds:jde,4),matrixIAL(iQP,jQP,:,ids,jds:jde,1), 1) ! Right
          call restore_bdy_field(matrixIAT_adj(iQP,jQP,:,ids:ide,jds-1,4),matrixIAL(iQP,jQP,:,ids,jds:jde,6), 1) ! below
          call restore_bdy_field(matrixIAB_adj(iQP,jQP,:,ids:ide,jde+1,4),matrixIAL(iQP,jQP,:,ids,jds:jde,5),-1) ! over
          ! Panel 5
          call restore_bdy_field(matrixIAR_adj(iQP,jQP,:,ids-1,jds:jde,5),matrixIAT(iQP,jQP,:,ids:ide,jde,4),-1) ! Left
          call restore_bdy_field(matrixIAL_adj(iQP,jQP,:,ide+1,jds:jde,5),matrixIAT(iQP,jQP,:,ids:ide,jde,2), 1) ! Right
          call restore_bdy_field(matrixIAT_adj(iQP,jQP,:,ids:ide,jds-1,5),matrixIAT(iQP,jQP,:,ids:ide,jde,1), 1) ! below
          call restore_bdy_field(matrixIAB_adj(iQP,jQP,:,ids:ide,jde+1,5),matrixIAT(iQP,jQP,:,ids:ide,jde,3),-1) ! over
          ! Panel 6
          call restore_bdy_field(matrixIAR_adj(iQP,jQP,:,ids-1,jds:jde,6),matrixIAB(iQP,jQP,:,ids:ide,jds,4), 1) ! Left
          call restore_bdy_field(matrixIAL_adj(iQP,jQP,:,ide+1,jds:jde,6),matrixIAB(iQP,jQP,:,ids:ide,jds,2),-1) ! Right
          call restore_bdy_field(matrixIAT_adj(iQP,jQP,:,ids:ide,jds-1,6),matrixIAB(iQP,jQP,:,ids:ide,jds,3),-1) ! below
          call restore_bdy_field(matrixIAB_adj(iQP,jQP,:,ids:ide,jde+1,6),matrixIAB(iQP,jQP,:,ids:ide,jds,1), 1) ! over
        enddo
      enddo
      
    end subroutine init_spatial_operator
    
    subroutine spatial_operator(stat,tend)
      type(stat_field), target, intent(inout) :: stat
      type(tend_field), target, intent(inout) :: tend
      
      integer(i_kind) :: iVar,i,j,iPatch,iPOC,iEOC
      
      real(r_kind) :: avg
  
      real(r_kind), dimension(nPointsOnEdge) :: eigL,eigR
      
      call fill_halo(stat%q,qQ(1,:,:,:,:))
      
      qC = stat%q
      
      do iVar = 1,nVar
        call reconstruction(qC(iVar,:,:,:  ),&
                            qL(iVar,:,:,:,:),&
                            qR(iVar,:,:,:,:),&
                            qB(iVar,:,:,:,:),&
                            qT(iVar,:,:,:,:),&
                            qQ(iVar,:,:,:,:))
        
        !call reconstruction(qC(iVar,:,:,:  ),&
        !                    qL(iVar,:,:,:,:),&
        !                    qR(iVar,:,:,:,:),&
        !                    qB(iVar,:,:,:,:),&
        !                    qT(iVar,:,:,:,:))
        
        !print*,maxval(qC(iVar,:,:,:  )),minval(qC(iVar,:,:,:  ))
        !print*,maxval(qL(iVar,:,:,:,:)),minval(qL(iVar,:,:,:,:))
        !print*,maxval(qR(iVar,:,:,:,:)),minval(qR(iVar,:,:,:,:))
        !print*,maxval(qB(iVar,:,:,:,:)),minval(qB(iVar,:,:,:,:))
        !print*,maxval(qT(iVar,:,:,:,:)),minval(qT(iVar,:,:,:,:))
        !print*,maxval(qQ(iVar,:,:,:,:)),minval(qQ(iVar,:,:,:,:))
        !print*,''
      enddo
      
      do iPatch = ifs,ife
        do j = jms,jme
          do i = ims,ime
            phit(:,i,j,iPatch) = qQ(1,:,i,j,iPatch) / sqrtG(cqs:cqe,i,j,iPatch) + ghs(cqs:cqe,i,j,iPatch)
            phitC(i,j,iPatch) = cell_quadrature(phit(:,i,j,iPatch))
          enddo
        enddo
      enddo
      
      call reconstruction(phitC       ,&
                          dqdx=dphitdx,&
                          dqdy=dphitdy)
      
      !do iPatch = ifs,ife
      !  do j = jds,jde
      !    do i = ids,ide
      !      print*,qC(3,  i,j,iPatch)
      !      print*,Gaussian_quadrature_1d(qL(3,:,i,j,iPatch))
      !      print*,Gaussian_quadrature_1d(qR(3,:,i,j,iPatch))
      !      print*,Gaussian_quadrature_1d(qB(3,:,i,j,iPatch))
      !      print*,Gaussian_quadrature_1d(qT(3,:,i,j,iPatch))
      !      print*,''
      !    enddo
      !  enddo
      !enddo
      
      call fill_bdy_flux(qL,qR,qB,qT)
      
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids-1,ide+1
            do iPOC = 1,nPointsOnEdge
              FL(:,iPOC,i,j,iPatch) = calc_F(sqrtGL(iPOC,i,j,iPatch),matrixIGL(:,:,iPOC,i,j,iPatch),qL(:,iPOC,i,j,iPatch))
              FR(:,iPOC,i,j,iPatch) = calc_F(sqrtGR(iPOC,i,j,iPatch),matrixIGR(:,:,iPOC,i,j,iPatch),qR(:,iPOC,i,j,iPatch))
            enddo
          enddo
        enddo
        do j = jds-1,jde+1
          do i = ids,ide
            do iPOC = 1,nPointsOnEdge
              GB(:,iPOC,i,j,iPatch) = calc_G(sqrtGB(iPOC,i,j,iPatch),matrixIGB(:,:,iPOC,i,j,iPatch),qB(:,iPOC,i,j,iPatch))
              GT(:,iPOC,i,j,iPatch) = calc_G(sqrtGT(iPOC,i,j,iPatch),matrixIGT(:,:,iPOC,i,j,iPatch),qT(:,iPOC,i,j,iPatch))
            enddo
          enddo
        enddo
      enddo
      
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide+1
            do iPOC = 1,nPointsOnEdge
              eigL(iPOC) = calc_eigenvalue_x(sqrtGR(iPOC,i-1,j,iPatch),matrixIGR(:,:,iPOC,i-1,j,iPatch),qR(:,iPOC,i-1,j,iPatch))
              eigR(iPOC) = calc_eigenvalue_x(sqrtGL(iPOC,i  ,j,iPatch),matrixIGL(:,:,iPOC,i  ,j,iPatch),qL(:,iPOC,i  ,j,iPatch))
              
              FeP(:,iPOC,i,j,iPatch) = 0.5 * ( FL(:,iPOC,i,j,iPatch) + FR(:,iPOC,i-1,j,iPatch) - max(eigL(iPOC),eigR(iPOC)) * ( qL(:,iPOC,i,j,iPatch) - qR(:,iPOC,i-1,j,iPatch) ) )
            enddo
            do iVar = 1,nVar
              Fe(iVar,i,j,iPatch) = Gaussian_quadrature_1d(FeP(iVar,:,i,j,iPatch))
            enddo
          enddo
        enddo
      enddo
        
      do iPatch = ifs,ife
        do j = jds,jde+1
          do i = ids,ide
            do iPOC = 1,nPointsOnEdge
              eigL(iPOC) = calc_eigenvalue_y(sqrtGT(iPOC,i,j-1,iPatch),matrixIGT(:,:,iPOC,i,j-1,iPatch),qT(:,iPOC,i,j-1,iPatch))
              eigR(iPOC) = calc_eigenvalue_y(sqrtGB(iPOC,i,j  ,iPatch),matrixIGB(:,:,iPOC,i,j  ,iPatch),qB(:,iPOC,i,j  ,iPatch))
              
              GeP(:,iPOC,i,j,iPatch) = 0.5 * ( GB(:,iPOC,i,j,iPatch) + GT(:,iPOC,i,j-1,iPatch) - max(eigL(iPOC),eigR(iPOC)) * ( qB(:,iPOC,i,j,iPatch) - qT(:,iPOC,i,j-1,iPatch) ) )
            enddo
            do iVar = 1,nVar
              Ge(iVar,i,j,iPatch) = Gaussian_quadrature_1d(GeP(iVar,:,i,j,iPatch))
            enddo
          enddo
        enddo
      enddo
      
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            src(:,i,j,iPatch) = calc_src(sqrtG(cqs:cqe,i,j,iPatch),matrixG(:,:,cqs:cqe,i,j,iPatch),matrixIG(:,:,cqs:cqe,i,j,iPatch),&
                                         qQ(:,:,i,j,iPatch),dphitdx(:,i,j,iPatch),dphitdy(:,i,j,iPatch)                            ,&
                                         tanx(cqs:cqe,i,j,iPatch),tany(cqs:cqe,i,j,iPatch)                                         ,&
                                         Coriolis(cqs:cqe,i,j,iPatch),delta(cqs:cqe,i,j,iPatch),iPatch)
          enddo
        enddo
      enddo
      
      !print*,maxval(Fe),minval(Fe)
      !print*,maxval(Ge),minval(Ge)
      
      !print*,maxval(FL),minval(FL)
      !print*,maxval(FR),minval(FR)
      !print*,maxval(GB),minval(GB)
      !print*,maxval(GT),minval(GT)
      
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            do iVar = 1,nVar
              tend%q(iVar,i,j,iPatch) = - ( Fe(iVar,i+1,j,iPatch) - Fe(iVar,i,j,iPatch) ) / dx - ( Ge(iVar,i,j+1,iPatch) - Ge(iVar,i,j,iPatch) ) / dy + src(iVar,i,j,iPatch)
              !tend%q(iVar,i,j,iPatch) = - ( Fe(iVar,i+1,j,iPatch) - Fe(iVar,i,j,iPatch) ) / dx
              !tend%q(iVar,i,j,iPatch) = - ( Ge(iVar,i,j+1,iPatch) - Ge(iVar,i,j,iPatch) ) / dy
              !tend%q(iVar,i,j,iPatch) = src(iVar,i,j,iPatch)
              
              !print*,iVar,i,j,iPatch
              !print*,stat%q(iVar,i,j,iPatch),tend%q(iVar,i,j,iPatch), tend%q(iVar,i,j,iPatch) / stat%q(iVar,i,j,iPatch)
              !print*,-( Fe(iVar,i+1,j,iPatch) - Fe(iVar,i,j,iPatch) ) / dx
              !print*,-( Ge(iVar,i,j+1,iPatch) - Ge(iVar,i,j,iPatch) ) / dy
              !print*,src(iVar,i,j,iPatch)
              !print*,''
              
            enddo
          enddo
        enddo
      enddo
      
      !print*,maxval(tend%q/stat%q),minval(tend%q/stat%q)
      !
      !call check_halo(stat%q)
      !call check_tend(tend%q)
      !
      !stop 'spatial_operator'
      
    end subroutine spatial_operator
    
    subroutine fill_halo(q,qQ)
      real(r_kind), dimension(             nVar,ims:ime,jms:jme,ifs:ife), intent(inout) :: q
      real(r_kind), dimension(nQuadPointsOnCell,ims:ime,jms:jme,ifs:ife), intent(inout) :: qQ
      
      real(r_kind), dimension(nVar,maxGhostPointsOnCell,ims:ime,jms:jme,ifs:ife) :: qg
      real(r_kind), dimension(nVar,nQuadPointsOnCell   ,ims:ime,jms:jme,ifs:ife) :: tgq ! value on triangle gaussian quadrature points
      
      integer(i_kind) :: iVar,i,j,iPatch,iPOC
      integer(i_kind) :: iRec,jRec
      integer(i_kind) :: ig,jg,pg,ng
      integer(i_kind) :: m,n
      integer(i_kind) :: igp,iqp
      
      real(r_kind), dimension(maxRecCells) :: u
      
      real(r_kind) :: uc,vc,us,vs
      
      !$OMP PARALLEL
      !$OMP DO PRIVATE(i,j,m,n,iVar,iPOC,iRec,jRec,u) COLLAPSE(3)
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            m = nGstRecCells(i,j,iPatch)
            n = nGhostPointsOnCell(i,j,iPatch)
            if(n>0)then
              do iVar = 1,nVar
                do iPOC = 1,m
                  iRec = iGstCell(iPOC,i,j,iPatch)
                  jRec = jGstCell(iPOC,i,j,iPatch)
                  u(iPOC) = q(iVar,iRec,jRec,iPatch)
                enddo
                qg(iVar,1:n,i,j,iPatch) = matmul(gstMatrix(1:n,1:m,i,j,iPatch),u(1:m))
              enddo
            endif
          enddo
        enddo
      enddo
      !$OMP END DO
      
      !$OMP DO PRIVATE(i,j,iPOC,ig,jg,pg,ng,iVar,igp,iqp,uc,vc,us,vs) COLLAPSE(3)
      do iPatch = ifs,ife
        do j = jms,jme
          do i = ims,ime
            !if( sum(ghost_n(:,i,j,iPatch)) /= 0 )then
            !if( ghost_p(1,i,j,iPatch) /= 0 )then
            if( .not.inDomain(i,j,iPatch) )then
              do iPOC = 1,nQuadPointsOnCell
                ig = ghost_i(iPOC,i,j,iPatch)
                jg = ghost_j(iPOC,i,j,iPatch)
                pg = ghost_p(iPOC,i,j,iPatch)
                ng = ghost_n(iPOC,i,j,iPatch)
                do iVar = 1,nVar
                  tgq(iVar,iPOC,i,j,iPatch) = qg(iVar,ng,ig,jg,pg)
                enddo
                
                igp = cgs + ng   - 1
                iqp = cqs + iPOC - 1
                
                uc = tgq(2,iPOC,i,j,iPatch) / tgq(1,iPOC,i,j,iPatch)
                vc = tgq(3,iPOC,i,j,iPatch) / tgq(1,iPOC,i,j,iPatch)
                call contravProjPlane2Sphere(us, vs, uc, vc, matrixA (:,:,igp,ig,jg,pg    ))
                call contravProjSphere2Plane(uc, vc, us, vs, matrixIA(:,:,iqp,i ,j ,iPatch))
                
                tgq(1,iPOC,i,j,iPatch) = tgq(1,iPOC,i,j,iPatch) / sqrtG(igp,ig,jg,pg) * sqrtG(iqp,i,j,iPatch)
                tgq(2,iPOC,i,j,iPatch) = tgq(1,iPOC,i,j,iPatch) * uc
                tgq(3,iPOC,i,j,iPatch) = tgq(1,iPOC,i,j,iPatch) * vc
              enddo
                
              do iVar = 1,nVar
                q(iVar,i,j,iPatch) = cell_quadrature( tgq(iVar,:,i,j,iPatch) )
              enddo
              qQ(:,i,j,iPatch) = tgq(1,:,i,j,iPatch)
            endif
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
    end subroutine fill_halo
    
    subroutine reconstruction(q,qL,qR,qB,qT,qQ,dqdx,dqdy)
      real(r_kind), dimension(                  ims:ime,jms:jme,ifs:ife), intent(in   )          :: q
      real(r_kind), dimension(nPointsOnEdge    ,ims:ime,jms:jme,ifs:ife), intent(  out),optional :: qL
      real(r_kind), dimension(nPointsOnEdge    ,ims:ime,jms:jme,ifs:ife), intent(  out),optional :: qR
      real(r_kind), dimension(nPointsOnEdge    ,ims:ime,jms:jme,ifs:ife), intent(  out),optional :: qB
      real(r_kind), dimension(nPointsOnEdge    ,ims:ime,jms:jme,ifs:ife), intent(  out),optional :: qT
      real(r_kind), dimension(nQuadPointsOnCell,ims:ime,jms:jme,ifs:ife), intent(inout),optional :: qQ
      real(r_kind), dimension(nQuadPointsOnCell,ims:ime,jms:jme,ifs:ife), intent(  out),optional :: dqdx ! x derivative on quadrature points
      real(r_kind), dimension(nQuadPointsOnCell,ims:ime,jms:jme,ifs:ife), intent(  out),optional :: dqdy ! y derivative on quadrature points
      
      real(r_kind), dimension(maxRecCells            ) :: u
      real(r_kind), dimension(maxRecCells,maxRecTerms) :: coordMtx
      real(r_kind), dimension(            maxRecTerms) :: polyCoef
      
      integer(i_kind) :: iVar,i,j,iPatch,iCOS
      integer(i_kind) :: iRec,jRec
      integer(i_kind) :: ic
      integer(i_kind) :: m,n
      
      !$OMP PARALLEL DO PRIVATE(j,i,m,n,iCOS,iRec,jRec,u,coordMtx,ic,polyCoef) COLLAPSE(3)
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            m = nRecCells(i,j,iPatch)
            n = nRecTerms(i,j,iPatch)
            do iCOS = 1,m
              iRec = iRecCell(iCOS,i,j,iPatch)
              jRec = jRecCell(iCOS,i,j,iPatch)
              
              u       (iCOS    ) = q(iRec,jRec,iPatch)
              coordMtx(iCOS,1:n) = polyCoordCoef(iCOS,1:n,i,j,iPatch)
            enddo
            ic = iCenCell(i,j,iPatch)
            
            polyCoef(1:n) = WLS_ENO(coordMtx(1:m,1:n),u(1:m),dx,m,n,ic)
            
            if(present(qL  )) qL  (:,i,j,iPatch) = matmul(recMatrixL (:,1:n,i,j,iPatch),polyCoef(1:n))
            if(present(qR  )) qR  (:,i,j,iPatch) = matmul(recMatrixR (:,1:n,i,j,iPatch),polyCoef(1:n))
            if(present(qB  )) qB  (:,i,j,iPatch) = matmul(recMatrixB (:,1:n,i,j,iPatch),polyCoef(1:n))
            if(present(qT  )) qT  (:,i,j,iPatch) = matmul(recMatrixT (:,1:n,i,j,iPatch),polyCoef(1:n))
                                                                     
            if(present(qQ  )) qQ  (:,i,j,iPatch) = matmul(recMatrixQ (:,1:n,i,j,iPatch),polyCoef(1:n))
            
            if(present(dqdx)) dqdx(:,i,j,iPatch) = matmul(recMatrixDx(:,1:n,i,j,iPatch),polyCoef(1:n))
            if(present(dqdy)) dqdy(:,i,j,iPatch) = matmul(recMatrixDy(:,1:n,i,j,iPatch),polyCoef(1:n))
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
    end subroutine reconstruction
    
    function calc_F(sqrtG,matrixIG,q)
      real(r_kind), dimension(nVar) :: calc_F
      real(r_kind)                 , intent(in) :: sqrtG
      real(r_kind), dimension(2,2) , intent(in) :: matrixIG
      real(r_kind), dimension(nVar), intent(in) :: q
      
      real(r_kind) :: phi,u,v
      real(r_kind) :: IG11,IG21
      
      IG11 = matrixIG(1,1)
      IG21 = matrixIG(2,1)
      
      phi = q(1) / sqrtG
      u   = q(2) / q(1)
      v   = q(3) / q(1)
      
      calc_F(1) = q(2)
      calc_F(2) = q(2) * u! + 0.5 * IG11 * sqrtG * phi**2
      calc_F(3) = q(2) * v! + 0.5 * IG21 * sqrtG * phi**2
    end function calc_F
    
    function calc_G(sqrtG,matrixIG,q)
      real(r_kind), dimension(nVar) :: calc_G
      real(r_kind)                 , intent(in) :: sqrtG
      real(r_kind), dimension(2,2) , intent(in) :: matrixIG
      real(r_kind), dimension(nVar), intent(in) :: q
      
      real(r_kind) :: phi,u,v
      real(r_kind) :: IG12,IG22
      
      IG12 = matrixIG(1,2)
      IG22 = matrixIG(2,2)
      
      phi = q(1) / sqrtG
      u   = q(2) / q(1)
      v   = q(3) / q(1)
      
      calc_G(1) = q(3)
      calc_G(2) = q(3) * u! + 0.5 * IG12 * sqrtG * phi**2
      calc_G(3) = q(3) * v! + 0.5 * IG22 * sqrtG * phi**2
    end function calc_G
    
    function calc_eigenvalue_x(sqrtG,matrixIG,q)
      real(r_kind) :: calc_eigenvalue_x
      real(r_kind), dimension(nVar), intent(in) :: q
      real(r_kind)                 , intent(in) :: sqrtG
      real(r_kind), dimension(2,2) , intent(in) :: matrixIG
      
      real(r_kind) :: IG11
      
      real(r_kind) :: phi
      real(r_kind) :: u
      
      phi = q(1) / sqrtG
      u   = q(2) / q(1)
      
      IG11 = matrixIG(1,1)
      
      calc_eigenvalue_x = max( abs( u + sqrt(IG11*phi) ), abs( u - sqrt(IG11*phi) ) )
      
    end function calc_eigenvalue_x
    
    function calc_eigenvalue_y(sqrtG,matrixIG,q)
      real(r_kind) :: calc_eigenvalue_y
      real(r_kind), dimension(nVar), intent(in) :: q
      real(r_kind)                 , intent(in) :: sqrtG
      real(r_kind), dimension(2,2) , intent(in) :: matrixIG
      
      real(r_kind) :: IG22
      
      real(r_kind) :: phi
      real(r_kind) :: v
      
      phi = q(1) / sqrtG
      v   = q(3) / q(1)
      
      IG22 = matrixIG(2,2)
      
      calc_eigenvalue_y = max( abs( v + sqrt(IG22*phi) ), abs( v - sqrt(IG22*phi) ) )
      
    end function calc_eigenvalue_y
    
    function calc_src(sqrtG,matrixG,matrixIG,q,dphitdx,dphitdy,x,y,Coriolis,delta,iPatch)
      real(r_kind), dimension(nVar) :: calc_src
      real(r_kind), dimension(     nQuadPointsOnCell), intent(in) :: sqrtG
      real(r_kind), dimension(2, 2,nQuadPointsOnCell), intent(in) :: matrixG
      real(r_kind), dimension(2, 2,nQuadPointsOnCell), intent(in) :: matrixIG
      real(r_kind), dimension(nVar,nQuadPointsOnCell), intent(in) :: q
      real(r_kind), dimension(     nQuadPointsOnCell), intent(in) :: dphitdx
      real(r_kind), dimension(     nQuadPointsOnCell), intent(in) :: dphitdy
      real(r_kind), dimension(     nQuadPointsOnCell), intent(in) :: x ! tan(x) actually
      real(r_kind), dimension(     nQuadPointsOnCell), intent(in) :: y ! tan(y) actually
      real(r_kind), dimension(     nQuadPointsOnCell), intent(in) :: Coriolis
      real(r_kind), dimension(     nQuadPointsOnCell), intent(in) :: delta
      integer(i_kind) :: iPatch
      
      real(r_kind), dimension(nVar,nQuadPointsOnCell) :: psi_M
      real(r_kind), dimension(nVar,nQuadPointsOnCell) :: psi_C
      real(r_kind), dimension(nVar,nQuadPointsOnCell) :: psi_B
      
      real(r_kind), dimension(nQuadPointsOnCell) :: phi
      real(r_kind), dimension(nQuadPointsOnCell) :: u
      real(r_kind), dimension(nQuadPointsOnCell) :: v
      real(r_kind), dimension(nQuadPointsOnCell) :: phiu
      real(r_kind), dimension(nQuadPointsOnCell) :: phiv
      real(r_kind), dimension(nQuadPointsOnCell) :: G11
      real(r_kind), dimension(nQuadPointsOnCell) :: G12
      real(r_kind), dimension(nQuadPointsOnCell) :: G21
      real(r_kind), dimension(nQuadPointsOnCell) :: G22
      real(r_kind), dimension(nQuadPointsOnCell) :: IG11
      real(r_kind), dimension(nQuadPointsOnCell) :: IG12
      real(r_kind), dimension(nQuadPointsOnCell) :: IG21
      real(r_kind), dimension(nQuadPointsOnCell) :: IG22
      
      integer :: iVar
      
      phi  = q(1,:) / sqrtG
      u    = q(2,:) / q(1,:)
      v    = q(3,:) / q(1,:)
      phiu = phi * u
      phiv = phi * v
      
      G11  = matrixG(1,1,:)
      G12  = matrixG(1,2,:)
      G21  = matrixG(2,1,:)
      G22  = matrixG(2,2,:)
           
      IG11 = matrixIG(1,1,:)
      IG12 = matrixIG(1,2,:)
      IG21 = matrixIG(2,1,:)
      IG22 = matrixIG(2,2,:)
      
      psi_M(1,:) = 0
      psi_M(2,:) = 2. * sqrtG / (delta**2) * (-x * y**2 * phiu * u + y * ( 1. + y**2 ) * phiu * v )
      psi_M(3,:) = 2. * sqrtG / (delta**2) * ( x * ( 1. + x**2 ) * phiu * v - x**2 * y * phiv * v )
      
      psi_C(1,:) = 0
      !psi_C(2,:) = Coriolis * (  G12 * phiu + G22 * phiv )
      !psi_C(3,:) = Coriolis * ( -G11 * phiu - G12 * phiv )
      !psi_C(2,:) = sqrtG**2 * Coriolis * ( -IG12 * phiu + IG11 * phiv )
      !psi_C(3,:) = sqrtG**2 * Coriolis * ( -IG22 * phiu + IG12 * phiv )
      
      if(iPatch==6)then
        psi_C(2,:) = -sqrtG * 2. * Omega / (delta**2) * ( -x*y*phiu + (1.+y*y)*phiv )
        psi_C(3,:) = -sqrtG * 2. * Omega / (delta**2) * ( -(1.+x*x)*phiu + x*y*phiv )
      elseif(iPatch==5)then
        psi_C(2,:) = sqrtG * 2. * Omega / (delta**2) * ( -x*y*phiu + (1.+y*y)*phiv )
        psi_C(3,:) = sqrtG * 2. * Omega / (delta**2) * ( -(1.+x*x)*phiu + x*y*phiv )
      else
        psi_C(2,:) = sqrtG * 2. * Omega / (delta**2) * y * ( -x*y*phiu + (1.+y*y)*phiv )
        psi_C(3,:) = sqrtG * 2. * Omega / (delta**2) * y * ( -(1.+x*x)*phiu + x*y*phiv )
      endif
      
      !print*,iPatch
      !print*,psi_C(2,:)
      !print*,Coriolis * (  G12 * phiu + G22 * phiv )
      !print*,sqrtG**2 * Coriolis * ( -IG12 * phiu + IG11 * phiv )
      !!print*,psi_C(3,:)
      !!print*,sqrtG**2 * Coriolis * ( -IG22 * phiu + IG12 * phiv )
      !!print*,Coriolis * ( -G11 * phiu - G12 * phiv )
      !print*,''
      
      psi_B(1,:) = 0
      psi_B(2,:) = - sqrtG * phi * ( IG11 * dphitdx + IG12 * dphitdy )
      psi_B(3,:) = - sqrtG * phi * ( IG21 * dphitdx + IG22 * dphitdy )
      
      !iVar = 2
      !print*,cell_quadrature( psi_M(iVar,:) )
      !print*,cell_quadrature( psi_C(iVar,:) )
      !print*,cell_quadrature( psi_B(iVar,:) )
      !print*,''
      
      do iVar = 1,nVar
        calc_src(iVar) = cell_quadrature( psi_M(iVar,:) + psi_C(iVar,:) + psi_B(iVar,:) )
        !calc_src(iVar) = cell_quadrature( psi_M(iVar,:) )
        !calc_src(iVar) = cell_quadrature( psi_C(iVar,:) )
        !calc_src(iVar) = cell_quadrature( psi_B(iVar,:) )
      enddo
    end function calc_src
    
    subroutine fill_bdy_flux(qL,qR,qB,qT)
      real(r_kind), dimension(nVar,nPointsOnEdge,ims:ime,jms:jme,ifs:ife),intent(inout) :: qL
      real(r_kind), dimension(nVar,nPointsOnEdge,ims:ime,jms:jme,ifs:ife),intent(inout) :: qR
      real(r_kind), dimension(nVar,nPointsOnEdge,ims:ime,jms:jme,ifs:ife),intent(inout) :: qB
      real(r_kind), dimension(nVar,nPointsOnEdge,ims:ime,jms:jme,ifs:ife),intent(inout) :: qT
      
      integer(i_kind) :: iVar,iPOE,i,j,iPatch
      integer(i_kind) :: dir
      
      do iVar = 1,nVar
        ! Panel 1
        call restore_bdy_field(qR(iVar,:,ids-1,jds:jde,1),qR(iVar,:,ide,jds:jde,4), 1) ! Left
        call restore_bdy_field(qL(iVar,:,ide+1,jds:jde,1),qL(iVar,:,ids,jds:jde,2), 1) ! Right
        call restore_bdy_field(qT(iVar,:,ids:ide,jds-1,1),qT(iVar,:,ids:ide,jde,6), 1) ! below
        call restore_bdy_field(qB(iVar,:,ids:ide,jde+1,1),qB(iVar,:,ids:ide,jds,5), 1) ! over
        ! Panel 2
        call restore_bdy_field(qR(iVar,:,ids-1,jds:jde,2),qR(iVar,:,ide,jds:jde,1), 1) ! Left
        call restore_bdy_field(qL(iVar,:,ide+1,jds:jde,2),qL(iVar,:,ids,jds:jde,3), 1) ! Right
        call restore_bdy_field(qT(iVar,:,ids:ide,jds-1,2),qR(iVar,:,ide,jds:jde,6),-1) ! below
        call restore_bdy_field(qB(iVar,:,ids:ide,jde+1,2),qR(iVar,:,ide,jds:jde,5), 1) ! over
        ! Panel 3
        call restore_bdy_field(qR(iVar,:,ids-1,jds:jde,3),qR(iVar,:,ide,jds:jde,2), 1) ! Left
        call restore_bdy_field(qL(iVar,:,ide+1,jds:jde,3),qL(iVar,:,ids,jds:jde,4), 1) ! Right
        call restore_bdy_field(qT(iVar,:,ids:ide,jds-1,3),qB(iVar,:,ids:ide,jds,6),-1) ! below
        call restore_bdy_field(qB(iVar,:,ids:ide,jde+1,3),qT(iVar,:,ids:ide,jde,5),-1) ! over
        ! Panel 4
        call restore_bdy_field(qR(iVar,:,ids-1,jds:jde,4),qR(iVar,:,ide,jds:jde,3), 1) ! Left
        call restore_bdy_field(qL(iVar,:,ide+1,jds:jde,4),qL(iVar,:,ids,jds:jde,1), 1) ! Right
        call restore_bdy_field(qT(iVar,:,ids:ide,jds-1,4),qL(iVar,:,ids,jds:jde,6), 1) ! below
        call restore_bdy_field(qB(iVar,:,ids:ide,jde+1,4),qL(iVar,:,ids,jds:jde,5),-1) ! over
        ! Panel 5
        call restore_bdy_field(qR(iVar,:,ids-1,jds:jde,5),qT(iVar,:,ids:ide,jde,4),-1) ! Left
        call restore_bdy_field(qL(iVar,:,ide+1,jds:jde,5),qT(iVar,:,ids:ide,jde,2), 1) ! Right
        call restore_bdy_field(qT(iVar,:,ids:ide,jds-1,5),qT(iVar,:,ids:ide,jde,1), 1) ! below
        call restore_bdy_field(qB(iVar,:,ids:ide,jde+1,5),qT(iVar,:,ids:ide,jde,3),-1) ! over
        ! Panel 6
        call restore_bdy_field(qR(iVar,:,ids-1,jds:jde,6),qB(iVar,:,ids:ide,jds,4), 1) ! Left
        call restore_bdy_field(qL(iVar,:,ide+1,jds:jde,6),qB(iVar,:,ids:ide,jds,2),-1) ! Right
        call restore_bdy_field(qT(iVar,:,ids:ide,jds-1,6),qB(iVar,:,ids:ide,jds,3),-1) ! below
        call restore_bdy_field(qB(iVar,:,ids:ide,jde+1,6),qB(iVar,:,ids:ide,jds,1), 1) ! over
      enddo
      
      do iPatch = ifs,ife
        call convert_bdy_flux(sqrtGR(:,ids-1,jds:jde,iPatch),sqrtGR_adj(:,ids-1,jds:jde,iPatch),matrixAR(:,:,:,ids-1,jds:jde,iPatch),matrixAR_adj(:,:,:,ids-1,jds:jde,iPatch),matrixIAR(:,:,:,ids-1,jds:jde,iPatch),matrixIAR_adj(:,:,:,ids-1,jds:jde,iPatch),qR(:,:,ids-1,jds:jde,iPatch))
        call convert_bdy_flux(sqrtGL(:,ide+1,jds:jde,iPatch),sqrtGL_adj(:,ide+1,jds:jde,iPatch),matrixAL(:,:,:,ide+1,jds:jde,iPatch),matrixAL_adj(:,:,:,ide+1,jds:jde,iPatch),matrixIAL(:,:,:,ide+1,jds:jde,iPatch),matrixIAL_adj(:,:,:,ide+1,jds:jde,iPatch),qL(:,:,ide+1,jds:jde,iPatch))
        call convert_bdy_flux(sqrtGT(:,ids:ide,jds-1,iPatch),sqrtGT_adj(:,ids:ide,jds-1,iPatch),matrixAT(:,:,:,ids:ide,jds-1,iPatch),matrixAT_adj(:,:,:,ids:ide,jds-1,iPatch),matrixIAT(:,:,:,ids:ide,jds-1,iPatch),matrixIAT_adj(:,:,:,ids:ide,jds-1,iPatch),qT(:,:,ids:ide,jds-1,iPatch))
        call convert_bdy_flux(sqrtGB(:,ids:ide,jde+1,iPatch),sqrtGB_adj(:,ids:ide,jde+1,iPatch),matrixAB(:,:,:,ids:ide,jde+1,iPatch),matrixAB_adj(:,:,:,ids:ide,jde+1,iPatch),matrixIAB(:,:,:,ids:ide,jde+1,iPatch),matrixIAB_adj(:,:,:,ids:ide,jde+1,iPatch),qB(:,:,ids:ide,jde+1,iPatch))
      enddo
      
    end subroutine fill_bdy_flux
    
    subroutine convert_bdy_flux(sqrtG,sqrtG_adj,matrixA,matrixA_adj,matrixIA,matrixIA_adj,q)
      real(r_kind), dimension(     nPointsOnEdge,nx), intent(in   ) :: sqrtG
      real(r_kind), dimension(     nPointsOnEdge,nx), intent(in   ) :: sqrtG_adj
      real(r_kind), dimension(2,2, nPointsOnEdge,nx), intent(in   ) :: matrixA
      real(r_kind), dimension(2,2, nPointsOnEdge,nx), intent(in   ) :: matrixA_adj
      real(r_kind), dimension(2,2, nPointsOnEdge,nx), intent(in   ) :: matrixIA
      real(r_kind), dimension(2,2, nPointsOnEdge,nx), intent(in   ) :: matrixIA_adj
      real(r_kind), dimension(nVar,nPointsOnEdge,nx), intent(inout) :: q
      
      real(r_kind) :: phi
      real(r_kind) :: uc,vc,us,vs
      
      integer(i_kind) :: i,j,iPatch,iPOE
      
      do i = 1,nx
        do iPOE = 1,nPointsOnEdge
          phi = q(1,iPOE,i) / sqrtG_adj(iPOE,i)
          uc  = q(2,iPOE,i) / q(1,iPOE,i)
          vc  = q(3,iPOE,i) / q(1,iPOE,i)
          call contravProjPlane2Sphere(us, vs, uc, vc, matrixA_adj(:,:,iPOE,i))
          call contravProjSphere2Plane(uc, vc, us, vs, matrixIA   (:,:,iPOE,i))
          q(1,iPOE,i) = sqrtG(iPOE,i) * phi
          q(2,iPOE,i) = q(1,iPOE,i) * uc
          q(3,iPOE,i) = q(1,iPOE,i) * vc
        enddo
      enddo
      
    end subroutine convert_bdy_flux
    
    subroutine restore_bdy_field(q1,q2,dir)
      real   (r_kind), dimension(nPointsOnEdge,nx), intent(inout) :: q1
      real   (r_kind), dimension(nPointsOnEdge,nx), intent(in   ) :: q2
      integer(i_kind)                             , intent(in   ) :: dir
      
      real   (r_kind), dimension(nPointsOnEdge*nx) :: qp
      
      integer(i_kind) :: iPOE,i,j
      
      if( dir == 1 )then
        j = 0
        do i = 1,nx
          do iPOE = 1,nPointsOnEdge
            j = j + 1
            qp(j) = q2(iPOE,i)
          enddo
        enddo
      elseif( dir == -1 )then
        j = 0
        do i = nx,1,-1
          do iPOE = nPointsOnEdge,1,-1
            j = j + 1
            qp(j) = q2(iPOE,i)
          enddo
        enddo
      endif
      
      j = 0
      do i = 1,nx
        do iPOE = 1,nPointsOnEdge
          j = j + 1
          q1(iPOE,i) = qp(j)
        enddo
      enddo
      
    end subroutine restore_bdy_field
    
    subroutine check_halo(q)
      real(r_kind), dimension(nVar,ims:ime,jms:jme,ifs:ife),intent(in) :: q
      
      real(r_kind), dimension(ims:ime,jms:jme,ifs:ife) :: phi
      real(r_kind), dimension(ims:ime,jms:jme,ifs:ife) :: uc,vc
      real(r_kind), dimension(ims:ime,jms:jme,ifs:ife) :: us,vs
      
      integer(i_kind) :: iVar,i,j,iPatch
      
      open(123,file='check_halo.txt')
      do iPatch = ifs,ife
        do j = jms,jme
          do i = ims,ime
            phi(i,j,iPatch) = q(1,i,j,iPatch) / sqrtGC(i,j,iPatch)
            uc (i,j,iPatch) = q(2,i,j,iPatch) / q(1,i,j,iPatch)
            vc (i,j,iPatch) = q(3,i,j,iPatch) / q(1,i,j,iPatch)
            call contravProjPlane2Sphere(us(i,j,iPatch), vs(i,j,iPatch), uc(i,j,iPatch), vc(i,j,iPatch), matrixA(:,:,cc,i,j,iPatch))
          enddo
        enddo
      enddo

      do iPatch = ifs,ife
        do j = jms,jme
          write(123,'(38e35.16)')(phi(i,j,iPatch),i=ims,ime)
        enddo
      enddo

      do iPatch = ifs,ife
        do j = jms,jme
          write(123,'(38e35.16)')(us(i,j,iPatch),i=ims,ime)
        enddo
      enddo

      do iPatch = ifs,ife
        do j = jms,jme
          write(123,'(38e35.16)')(vs(i,j,iPatch),i=ims,ime)
        enddo
      enddo
      
      close(123)
    
    end subroutine check_halo
    
    subroutine check_tend(q)
      real(r_kind), dimension(nVar,ims:ime,jms:jme,ifs:ife),intent(in) :: q
      
      integer(i_kind) :: iVar,i,j,iPatch
      
      open(123,file='check_tend.txt')
      
      do iVar = 1,nVar
        do iPatch = ifs,ife
          do j = jds,jde
            write(123,'(30e35.16)')(q(iVar,i,j,iPatch),i=ids,ide)
          enddo
        enddo
      enddo
      
      close(123)
      
    end subroutine check_tend
    
end module spatial_operators_mod

