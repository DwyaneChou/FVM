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
  
  public qC, fill_halo, reconstruction ! just for diag
      
  integer(i_kind),dimension(:,:,:  ), allocatable :: nRecCells ! number of cells for reconstruction
  integer(i_kind),dimension(:,:,:  ), allocatable :: nGstRecCells ! number of cells for ghost point reconstruction
  integer(i_kind),dimension(:,:,:,:), allocatable :: nWENOCells ! number of cells for WENO reconstruction
  integer(i_kind),dimension(:,:,:  ), allocatable :: nRecTerms
  
  integer(i_kind),dimension(:,:,:), allocatable :: locPolyDegree ! degree of local reconstruction polynomial
  
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: polyCoordCoef
  
  integer(i_kind), dimension(:,:,:), allocatable :: iCenCell ! center cell index on reconstruction stencil
  
  integer(i_kind), dimension(:,:,:,:), allocatable :: iRecCell ! x index of reconstruction cells
  integer(i_kind), dimension(:,:,:,:), allocatable :: jRecCell ! y index of reconstruction cells
  
  integer(i_kind), dimension(:,:,:,:), allocatable :: iGstCell ! x index of ghost reconstruction cells
  integer(i_kind), dimension(:,:,:,:), allocatable :: jGstCell ! y index of ghost reconstruction cells
  
  integer(i_kind), dimension(:,:,:,:,:), allocatable :: iWENOCell ! x index of weno reconstruction cells
  integer(i_kind), dimension(:,:,:,:,:), allocatable :: jWENOCell ! y index of weno reconstruction cells
  
  real   (r_kind), dimension(:,:,:,:,:,:), allocatable :: WENOPoly ! polynomial coefficients for WENO
  real   (r_kind), dimension(:,:,:,:,:,:), allocatable :: invWENOPoly ! polynomial coefficients for WENO
  
  integer(i_kind), dimension(:,:,:,:), allocatable :: iCenWENO ! center cell index on reconstruction stencil for WENO2D
  
  real   (r_kind), dimension(:,:), allocatable :: r ! optimal coefficients for WENO 2D
  
  real   (r_kind), dimension(:,:), allocatable :: existWENOTerm
  
  integer(i_kind), dimension(:,:,:,:), allocatable :: rematch_idx_3_to_3
  integer(i_kind), dimension(:,:,:,:), allocatable :: rematch_idx_5_to_5
  integer(i_kind), dimension(:,:,:,:), allocatable :: rematch_idx_7_to_7
  integer(i_kind), dimension(:      ), allocatable :: rematch_idx_2_to_3
  integer(i_kind), dimension(:      ), allocatable :: rematch_idx_3_to_5
  integer(i_kind), dimension(:      ), allocatable :: rematch_idx_5_to_7
  
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: recMatrixL
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: recMatrixR
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: recMatrixB
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: recMatrixT
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: recMatrixQ ! reconstruction for gaussian quadrature points
  
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: polyMatrixL
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: polyMatrixR
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: polyMatrixB
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: polyMatrixT
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: polyMatrixQ
  
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: recMatrixDx
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: recMatrixDy
  
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: gstMatrix
  
  real   (r_kind), dimension(:,:,:,:), allocatable :: dh
  
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
  
  real(r_kind), dimension(:,:,:,:), allocatable :: Fe    ! F on edges of each cell
  real(r_kind), dimension(:,:,:,:), allocatable :: Ge    ! H on edges of each cell
  
  real(r_kind), dimension(:,:,:,:,:), allocatable :: FeP   ! F on points on edges of each cell
  real(r_kind), dimension(:,:,:,:,:), allocatable :: GeP   ! H on points on edges of each cell
  
  real(r_kind), dimension(:,:,:,:), allocatable :: src   ! source term
  
  real(r_kind), dimension(:,:,:,:), allocatable :: phit    ! phi + phis
  real(r_kind), dimension(:,:,:  ), allocatable :: phitC   ! phi + phis on cell
  
  real(r_kind), dimension(:,:,:,:), allocatable :: dphitdx    ! d phi_t / dx
  real(r_kind), dimension(:,:,:,:), allocatable :: dphitdy    ! d phi_t / dy
  
    contains
    subroutine init_spatial_operator
  
      real(r_kind), dimension(:,:,:,:,:), allocatable :: A
      real(r_kind), dimension(:,:,:,:,:), allocatable :: invA
  
      real(r_kind), dimension(:,:,:,:,:), allocatable :: Apoly
      real(r_kind), dimension(:,:,:,:,:), allocatable :: invApoly
  
      real(r_kind), dimension(:,:,:,:,:), allocatable :: xRel ! relative x coordinate of reconstruction cells
      real(r_kind), dimension(:,:,:,:,:), allocatable :: yRel ! relative y coordinate of reconstruction cells
      
      real(r_kind), dimension(:,:,:,:,:,:), allocatable :: xRelWENO ! relative x coordinate of reconstruction cells for WENO2D
      real(r_kind), dimension(:,:,:,:,:,:), allocatable :: yRelWENO ! relative y coordinate of reconstruction cells for WENO2D
      
      real(r_kind), dimension(:,:,:,:,:), allocatable :: xGst ! relative x coordinate of Ghost reconstruction cells
      real(r_kind), dimension(:,:,:,:,:), allocatable :: yGst ! relative y coordinate of Ghost reconstruction cells
  
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
      
      real(r_kind), dimension(:), allocatable :: xq
      real(r_kind), dimension(:), allocatable :: yq
      
      real(r_kind), dimension(:), allocatable :: existPolyTerm
      
      integer(i_kind) :: i,j,k,iPatch
      integer(i_kind) :: iCOS ! indices of Cells On Stencils
      integer(i_kind) :: iPOC ! indices of points on cell
      integer(i_kind) :: iStencil,jStencil
      integer(i_kind) :: iRec,jRec
      integer(i_kind) :: iQP,jQP
      integer(i_kind) :: iR,jR
      integer(i_kind) :: irs,ire,jrs,jre
      
      integer(i_kind) :: nxp,nyp
      integer(i_kind) :: invstat
      integer(i_kind) :: nRC,nRT
      
      integer(i_kind) :: pg
      
      integer(i_kind) :: iidx(2)
      integer(i_kind) :: jidx(2)
      integer(i_kind) :: xdir
      integer(i_kind) :: ydir
        
      allocate(nRecCells   (         ids:ide,jds:jde,ifs:ife))
      allocate(nGstRecCells(         ids:ide,jds:jde,ifs:ife))
      allocate(nRecTerms   (         ids:ide,jds:jde,ifs:ife))
      
      allocate(iRecCell  (maxRecCells,ims:ime,jms:jme,ifs:ife))
      allocate(jRecCell  (maxRecCells,ims:ime,jms:jme,ifs:ife))
      
      allocate(iGstCell  (maxRecCells,ims:ime,jms:jme,ifs:ife))
      allocate(jGstCell  (maxRecCells,ims:ime,jms:jme,ifs:ife))
      
      allocate(xRel(4,maxRecCells,ids:ide,jds:jde,ifs:ife))
      allocate(yRel(4,maxRecCells,ids:ide,jds:jde,ifs:ife))
      
      allocate(xRelWENO(nStencil_all,4,maxRecCells,ids:ide,jds:jde,ifs:ife))
      allocate(yRelWENO(nStencil_all,4,maxRecCells,ids:ide,jds:jde,ifs:ife))
      
      allocate(xGst(4,maxRecCells,ids:ide,jds:jde,ifs:ife))
      allocate(yGst(4,maxRecCells,ids:ide,jds:jde,ifs:ife))
      
      if(trim(reconstruct_scheme)=='WLS-ENO')then
        allocate(iCenCell  (                              ims:ime,jms:jme,ifs:ife))
        allocate(dh        (                  maxRecCells,ids:ide,jds:jde,ifs:ife))
        allocate(recMatrixL(nPointsOnEdge    ,maxRecTerms,ids:ide,jds:jde,ifs:ife))
        allocate(recMatrixR(nPointsOnEdge    ,maxRecTerms,ids:ide,jds:jde,ifs:ife))
        allocate(recMatrixB(nPointsOnEdge    ,maxRecTerms,ids:ide,jds:jde,ifs:ife))
        allocate(recMatrixT(nPointsOnEdge    ,maxRecTerms,ids:ide,jds:jde,ifs:ife))
        allocate(recMatrixQ(nQuadPointsOnCell,maxRecTerms,ids:ide,jds:jde,ifs:ife))
        
        allocate(locPolyDegree(                        ids:ide,jds:jde,ifs:ife))
        allocate(polyCoordCoef(maxRecCells,maxRecTerms,ids:ide,jds:jde,ifs:ife))
      elseif(trim(reconstruct_scheme)=='Polynomial')then
        allocate(polyMatrixL(nPointsOnEdge    ,maxRecTerms,ids:ide,jds:jde,ifs:ife))
        allocate(polyMatrixR(nPointsOnEdge    ,maxRecTerms,ids:ide,jds:jde,ifs:ife))
        allocate(polyMatrixB(nPointsOnEdge    ,maxRecTerms,ids:ide,jds:jde,ifs:ife))
        allocate(polyMatrixT(nPointsOnEdge    ,maxRecTerms,ids:ide,jds:jde,ifs:ife))
        allocate(polyMatrixQ(nQuadPointsOnCell,maxRecTerms,ids:ide,jds:jde,ifs:ife))
      
        allocate(Apoly   (maxRecCells,maxRecCells,ids:ide,jds:jde,ifs:ife))
        allocate(invApoly(maxRecCells,maxRecCells,ids:ide,jds:jde,ifs:ife))
      
        allocate(existPolyTerm(maxRecCells))
      elseif(trim(reconstruct_scheme)=='WENO2D')then
        allocate(iCenWENO    (nStencil_all,ids:ide,jds:jde,ifs:ife))
        allocate(nWENOCells  (nStencil_all,ids:ide,jds:jde,ifs:ife))
      
        allocate(iWENOCell (nStencil_all,maxRecCells,ims:ime,jms:jme,ifs:ife))
        allocate(jWENOCell (nStencil_all,maxRecCells,ims:ime,jms:jme,ifs:ife))
        
        allocate(WENOPoly   (nStencil_all,maxRecCells,maxRecCells,ims:ime,jms:jme,ifs:ife))
        allocate(invWENOPoly(nStencil_all,maxRecCells,maxRecCells,ims:ime,jms:jme,ifs:ife))
      
        allocate(r(nStencil,nStencil))
        
        allocate(existWENOTerm(nStencil_all,maxRecCells))
        
        allocate(rematch_idx_3_to_3(9 ,ids:ide,jds:jde,ifs:ife))
        allocate(rematch_idx_5_to_5(25,ids:ide,jds:jde,ifs:ife))
        allocate(rematch_idx_7_to_7(49,ids:ide,jds:jde,ifs:ife))
        allocate(rematch_idx_2_to_3(3 ))
        allocate(rematch_idx_3_to_5(9 ))
        allocate(rematch_idx_5_to_7(25))
        
        allocate(polyMatrixL(nPointsOnEdge    ,maxRecTerms,ids:ide,jds:jde,ifs:ife))
        allocate(polyMatrixR(nPointsOnEdge    ,maxRecTerms,ids:ide,jds:jde,ifs:ife))
        allocate(polyMatrixB(nPointsOnEdge    ,maxRecTerms,ids:ide,jds:jde,ifs:ife))
        allocate(polyMatrixT(nPointsOnEdge    ,maxRecTerms,ids:ide,jds:jde,ifs:ife))
        allocate(polyMatrixQ(nQuadPointsOnCell,maxRecTerms,ids:ide,jds:jde,ifs:ife))
      endif
      
      allocate(recMatrixDx(nQuadPointsOnCell,maxRecTerms,ids:ide,jds:jde,ifs:ife))
      allocate(recMatrixDy(nQuadPointsOnCell,maxRecTerms,ids:ide,jds:jde,ifs:ife))
      
      allocate(A   (maxRecCells,maxRecCells,ids:ide,jds:jde,ifs:ife))
      allocate(invA(maxRecCells,maxRecCells,ids:ide,jds:jde,ifs:ife))
      
      allocate(gstMatrix(maxGhostPointsOnCell,maxRecCells,ids:ide,jds:jde,ifs:ife))
      
      allocate(qC(nVar,                  ims:ime,jms:jme,ifs:ife))
      allocate(qL(nVar,nPointsOnEdge    ,ims:ime,jms:jme,ifs:ife))
      allocate(qR(nVar,nPointsOnEdge    ,ims:ime,jms:jme,ifs:ife))
      allocate(qB(nVar,nPointsOnEdge    ,ims:ime,jms:jme,ifs:ife))
      allocate(qT(nVar,nPointsOnEdge    ,ims:ime,jms:jme,ifs:ife))
      allocate(qQ(nVar,nQuadPointsOnCell,ims:ime,jms:jme,ifs:ife))
      
      allocate(Fe(nVar,ids:ide+1,jds:jde  ,ifs:ife))
      allocate(Ge(nVar,ids:ide  ,jds:jde+1,ifs:ife))
      
      allocate(FeP(nVar,nPointsOnEdge,ids:ide+1,jds:jde  ,ifs:ife))
      allocate(GeP(nVar,nPointsOnEdge,ids:ide  ,jds:jde+1,ifs:ife))
      
      allocate(src(nVar,ids:ide,jds:jde,ifs:ife))
      
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
      
      allocate(xq(nQuadPointsOncell))
      allocate(yq(nQuadPointsOncell))
      
      dV = dx * dy
      
      print*,''
      print*,'Set reconstruction cells on each stencil'
      print*,''
      
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            ! Reconstruction cells
            iCOS = 0
            do jRec = -recBdy,recBdy
              do iRec = -recBdy,recBdy
                if( .not.inCorner(i+iRec,j+jRec,iPatch) )then
                  iCOS = iCOS + 1
                  iRecCell(iCOS,i,j,iPatch) = i + iRec
                  jRecCell(iCOS,i,j,iPatch) = j + jRec
                endif
              enddo
            enddo
            nRecCells(i,j,iPatch) = iCOS
            
            ! Ghost interpolation cells
            iCOS = 0
            !do jRec = -recBdy,recBdy
            !  do iRec = -recBdy,recBdy
            !    if(inDomain(i+iRec,j+jRec,iPatch))then
            !      iCOS = iCOS + 1
            !      iGstCell(iCOS,i,j,iPatch) = i + iRec
            !      jGstCell(iCOS,i,j,iPatch) = j + jRec
            !    endif
            !  enddo
            !enddo
            
            irs = i - recBdy
            ire = i + recBdy
            jrs = j - recBdy
            jre = j + recBdy
            if(i-recBdy<ids)then
              irs = ids
              ire = ids + 2 * recBdy
            endif
            if(i+recBdy>ide)then
              irs = ide - 2 * recBdy
              ire = ide
            endif
            if(j-recBdy<jds)then
              jrs = jds
              jre = jds + 2 * recBdy
            endif
            if(j+recBdy>ide)then
              jrs = jde - 2 * recBdy
              jre = jde
            endif
            do jRec = jrs,jre
              do iRec = irs,ire
                iCOS = iCOS + 1
                iGstCell(iCOS,i,j,iPatch) = iRec
                jGstCell(iCOS,i,j,iPatch) = jRec
              enddo
            enddo
            
            nGstRecCells(i,j,iPatch) = iCOS
            
          enddo
        enddo
      enddo
      
      if(trim(reconstruct_scheme)=='WENO2D')then
        ! WENO 2D stencil
        do iPatch = ifs,ife
          do j = jds,jde
            do i = ids,ide
              iCOS = 2
              iWENOCell(1:nStencil1,iCOS,i,j,iPatch) = i
              jWENOCell(1:nStencil1,iCOS,i,j,iPatch) = j
              
              iStencil = 1 ! Left,Center,Top
              iCOS = 1
              iWENOCell(iStencil,iCOS,i,j,iPatch) = i - 1
              jWENOCell(iStencil,iCOS,i,j,iPatch) = j
              iCOS = 3
              iWENOCell(iStencil,iCOS,i,j,iPatch) = i
              jWENOCell(iStencil,iCOS,i,j,iPatch) = j + 1
              
              iStencil = 2 ! Bottom,Center,Left
              iCOS = 1
              iWENOCell(iStencil,iCOS,i,j,iPatch) = i
              jWENOCell(iStencil,iCOS,i,j,iPatch) = j - 1
              iCOS = 3
              iWENOCell(iStencil,iCOS,i,j,iPatch) = i - 1
              jWENOCell(iStencil,iCOS,i,j,iPatch) = j
              
              iStencil = 3 ! Right,Center,Bottom
              iCOS = 1
              iWENOCell(iStencil,iCOS,i,j,iPatch) = i + 1
              jWENOCell(iStencil,iCOS,i,j,iPatch) = j
              iCOS = 3
              iWENOCell(iStencil,iCOS,i,j,iPatch) = i
              jWENOCell(iStencil,iCOS,i,j,iPatch) = j - 1
              
              iStencil = 4 ! Top,Center,Right
              iCOS = 1
              iWENOCell(iStencil,iCOS,i,j,iPatch) = i
              jWENOCell(iStencil,iCOS,i,j,iPatch) = j + 1
              iCOS = 3
              iWENOCell(iStencil,iCOS,i,j,iPatch) = i + 1
              jWENOCell(iStencil,iCOS,i,j,iPatch) = j
              
              nWENOCells(1:nStencil1,i,j,iPatch) = 3
                
              do iStencil = nStencil1+1,nStencil_all
                iCOS = 0
                do jRec = -(iStencil-nStencil1)+1,(iStencil-nStencil1)-1
                  do iRec = -(iStencil-nStencil1)+1,(iStencil-nStencil1)-1
                    if( .not.inCorner(i+iRec,j+jRec,iPatch) )then
                      iCOS = iCOS + 1
                      iWENOCell(iStencil,iCOS,i,j,iPatch) = i + iRec
                      jWENOCell(iStencil,iCOS,i,j,iPatch) = j + jRec
                    endif
                  enddo
                enddo
                nWENOCells(iStencil,i,j,iPatch) = iCOS
              enddo
            enddo
          enddo
        enddo
      endif
      
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
      
      !$OMP PARALLEL DO PRIVATE(i,j,iCOS,iRec,jRec,nRC) COLLAPSE(3)
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            nRC = nRecCells(i,j,iPatch)
            do iCOS = 1,nRC
              iRec = iRecCell(iCOS,i,j,iPatch)
              jRec = jRecCell(iCOS,i,j,iPatch)
              
              xRel(:,iCOS,i,j,iPatch) = ( x(ccs:cce,iRec,jRec,iPatch) - x(cc,i,j,iPatch) ) * recdx
              yRel(:,iCOS,i,j,iPatch) = ( y(ccs:cce,iRec,jRec,iPatch) - y(cc,i,j,iPatch) ) * recdy
            enddo
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      if(trim(reconstruct_scheme)=='WENO2D')then
        !$OMP PARALLEL DO PRIVATE(i,j,iStencil,nRC,iCOS,iRec,jRec)
        do iPatch = ifs,ife
          do j = jds,jde
            do i = ids,ide
              do iStencil = 1,nStencil_all
                nRC = nWENOCells(iStencil,i,j,iPatch)
                do iCOS = 1,nRC
                  iRec = iWENOCell(iStencil,iCOS,i,j,iPatch)
                  jRec = jWENOCell(iStencil,iCOS,i,j,iPatch)
                  
                  xRelWENO(iStencil,:,iCOS,i,j,iPatch) = ( x(ccs:cce,iRec,jRec,iPatch) - x(cc,i,j,iPatch) ) * recdx
                  yRelWENO(iStencil,:,iCOS,i,j,iPatch) = ( y(ccs:cce,iRec,jRec,iPatch) - y(cc,i,j,iPatch) ) * recdy
                enddo
              enddo
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
        
        !$OMP PARALLEL DO PRIVATE(i,j,iStencil,nRC,iCOS,invstat,nxp,nyp,iidx,jidx,xdir,ydir,k,iR,jR,iRec,jRec,existWENOTerm,iPOC,xq,yq)
        do iPatch = ifs,ife
          do j = jds,jde
            do i = ids,ide
              ! Calculate WENO polynomial coefficient matrices
              do iStencil = 1,nStencil1
                iCenWENO(iStencil,i,j,iPatch) = 2
                
                nRC = nWENOCells(iStencil,i,j,iPatch)
                do iCOS = 1,nRC
                  call calc_polynomial_square_integration(1,xRelWENO(iStencil,1,iCOS,i,j,iPatch),xRelWENO(iStencil,2,iCOS,i,j,iPatch),&
                                                            yRelWENO(iStencil,1,iCOS,i,j,iPatch),yRelWENO(iStencil,4,iCOS,i,j,iPatch),&
                                                            WENOPoly(iStencil,iCOS,1:nRC,i,j,iPatch))
                enddo
                
                call BRINV(nRC,WENOPoly(iStencil,1:nRC,1:nRC,i,j,iPatch),invWENOPoly(iStencil,1:nRC,1:nRC,i,j,iPatch),invstat)
                if(invstat==0)then
                  print*,'Inverse Apoly dost not exist'
                  print*,'i,j,iPatch,nRC,nxp,nyp are'
                  print*,i,j,iPatch,nRC,nxp,nyp
                  stop 'Check BRINV for Special treamtment on boundary cells'
                endif
              enddo
              
              do iStencil = nStencil1+1,nStencil_all
                nRC = nWENOCells(iStencil,i,j,iPatch)
                nxp = maxval(iWENOCell(iStencil,1:nRC,i,j,iPatch)) - minval(iWENOCell(iStencil,1:nRC,i,j,iPatch)) + 1
                nyp = maxval(jWENOCell(iStencil,1:nRC,i,j,iPatch)) - minval(jWENOCell(iStencil,1:nRC,i,j,iPatch)) + 1
                
                ! Pick the cells that not in corner for reconstruction
                if(.not.noCorner(i,j,iPatch))then
                  iidx = minloc( abs( x(cc,iWENOCell(iStencil,1:nRC,i,j,iPatch),jWENOCell(iStencil,1:nRC,i,j,iPatch),iPatch) ) )
                  jidx = minloc( abs( y(cc,iWENOCell(iStencil,1:nRC,i,j,iPatch),jWENOCell(iStencil,1:nRC,i,j,iPatch),iPatch) ) )
                  
                  xdir = 1
                  ydir = 1
                  if(iWENOCell(iStencil,iidx(1),i,j,iPatch)>iWENOCell(iStencil,1,i,j,iPatch))xdir=-1
                  if(jWENOCell(iStencil,jidx(2),i,j,iPatch)>jWENOCell(iStencil,1,i,j,iPatch))ydir=-1
                  
                  k = 0
                  do jR = 0,nyp-1
                    do iR = 0,nxp-1
                      k = k + 1
                      iRec = iWENOCell(iStencil,iidx(1),i,j,iPatch) + xdir * iR
                      jRec = jWENOCell(iStencil,jidx(2),i,j,iPatch) + ydir * jR
                      if( .not.inCorner(iRec,jRec,iPatch) )then
                        existWENOTerm(iStencil,k) = 1
                      else
                        existWENOTerm(iStencil,k) = 0
                      endif
                    enddo
                  enddo
                else
                  existWENOTerm(iStencil,:) = 1
                endif
                
                k    = 0
                iCOS = 0
                do jR = 1,nyp
                  do iR = 1,nxp
                    k = k + 1
                    if(existWENOTerm(iStencil,k)>0)then
                      iCOS = iCOS + 1
                      
                      iRec = ( nxp + 1 ) / 2
                      jRec = ( nyp + 1 ) / 2
                      if(iRec==iR.and.jRec==jR)iCenWENO(iStencil,i,j,iPatch) = iCOS
                      
                      call calc_rectangle_poly_integration(nxp,nyp,&
                                                           xRelWENO(iStencil,1,iCOS,i,j,iPatch),xRelWENO(iStencil,2,iCOS,i,j,iPatch),&
                                                           yRelWENO(iStencil,1,iCOS,i,j,iPatch),yRelWENO(iStencil,4,iCOS,i,j,iPatch),&
                                                           WENOPoly(iStencil,iCOS,1:nRC,i,j,iPatch),existWENOTerm(iStencil,1:nxp*nyp))
                      
                      if(nxp==3)then
                        rematch_idx_3_to_3(iCOS,i,j,iPatch) = k
                      elseif(nxp==5)then
                        rematch_idx_5_to_5(iCOS,i,j,iPatch) = k
                      elseif(nxp==7)then
                        rematch_idx_7_to_7(iCOS,i,j,iPatch) = k
                      endif
                    endif
                  enddo
                enddo
                
                call BRINV(nRC,WENOPoly(iStencil,1:nRC,1:nRC,i,j,iPatch),invWENOPoly(iStencil,1:nRC,1:nRC,i,j,iPatch),invstat)
                if(invstat==0)then
                  print*,'Inverse Apoly dost not exist'
                  print*,'i,j,iPatch,nRC,nxp,nyp are'
                  print*,i,j,iPatch,nRC,nxp,nyp
                  stop 'Check BRINV for Special treamtment on boundary cells'
                endif
              enddo
              
              ! Calculate reconstruction matrix on edge
              nRC = stencil_width**2
              nxp = stencil_width
              nyp = stencil_width
              call calc_rectangle_poly_matrix(nxp,nyp,nPointsOnEdge,xL*recdx,yL*recdy,polyMatrixL(:,1:nRC,i,j,iPatch))
              call calc_rectangle_poly_matrix(nxp,nyp,nPointsOnEdge,xR*recdx,yR*recdy,polyMatrixR(:,1:nRC,i,j,iPatch))
              call calc_rectangle_poly_matrix(nxp,nyp,nPointsOnEdge,xB*recdx,yB*recdy,polyMatrixB(:,1:nRC,i,j,iPatch))
              call calc_rectangle_poly_matrix(nxp,nyp,nPointsOnEdge,xT*recdx,yT*recdy,polyMatrixT(:,1:nRC,i,j,iPatch))
              
              do iPOC = 1,nQuadPointsOncell
                xq(iPOC) = x(cqs+iPOC-1,i,j,iPatch) - x(cc,i,j,iPatch)
                yq(iPOC) = y(cqs+iPOC-1,i,j,iPatch) - y(cc,i,j,iPatch)
              enddo
              call calc_rectangle_poly_matrix(nxp,nyp,nQuadPointsOnCell,xq*recdx,yq*recdy,polyMatrixQ(:,1:nRC,i,j,iPatch))
              
              ! Calculate derivative matrix on quadrature points on cell
              !xq = 0
              !yq = 0
              call calc_rectangle_poly_deriv_matrix(nxp,nyp,nQuadPointsOnCell,xq,yq,&
                                                    recMatrixDx(:,1:nRC,i,j,iPatch),&
                                                    recMatrixDy(:,1:nRC,i,j,iPatch))
              recMatrixDx(:,1:nRC,i,j,iPatch) = recMatrixDx(:,1:nRC,i,j,iPatch) / dx
              recMatrixDy(:,1:nRC,i,j,iPatch) = recMatrixDy(:,1:nRC,i,j,iPatch) / dy
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
        
        k    = 0
        iCOS = 0
        do jRec = 1,3
          do iRec = 1,3
            k = k + 1
            if( iRec<=2 .and. jRec<=2 .and. (iRec+jRec)<=3 )then
              iCOS = iCOS + 1
              rematch_idx_2_to_3(iCOS) = k
            endif
          enddo
        enddo
        
        k    = 0
        iCOS = 0
        do jRec = 1,5
          do iRec = 1,5
            k = k + 1
            if(iRec<=3.and.jRec<=3)then
              iCOS = iCOS + 1
              rematch_idx_3_to_5(iCOS) = k
            endif
          enddo
        enddo
        
        k    = 0
        iCOS = 0
        do jRec = 1,7
          do iRec = 1,7
            k = k + 1
            if(iRec<=5.and.jRec<=5)then
              iCOS = iCOS + 1
              rematch_idx_5_to_7(iCOS) = k
            endif
          enddo
        enddo
        
        r = 0
        do jStencil = 1,nStencil
          do iStencil = 1,jStencil
            r(iStencil,jStencil) = 10**(2*(iStencil-1))
            !r(iStencil,jStencil) = 10**(1*iStencil-1)
          enddo
          r(:,jStencil) = r(:,jStencil) / sum(r(1:jStencil,jStencil))
        enddo
      endif
      
      if(trim(reconstruct_scheme)=='WLS-ENO')then
        !$OMP PARALLEL DO PRIVATE(i,j,iCOS,iRec,jRec,nRC,nRT,iPOC,xq,yq) COLLAPSE(3)
        do iPatch = ifs,ife
          do j = jds,jde
            do i = ids,ide
              nRC = nRecCells(i,j,iPatch)
              
              locPolyDegree(i,j,iPatch) = min( maxval(iRecCell(1:nRC,i,j,iPatch)) - minval(iRecCell(1:nRC,i,j,iPatch)), &
                                               maxval(jRecCell(1:nRC,i,j,iPatch)) - minval(jRecCell(1:nRC,i,j,iPatch)) )
              
              nRecTerms(i,j,iPatch) = ( locPolyDegree(i,j,iPatch) + 1 ) * ( locPolyDegree(i,j,iPatch) + 2 ) / 2
              
              nRT = nRecTerms(i,j,iPatch)
              
              do iCOS = 1,nRC
                iRec = iRecCell(iCOS,i,j,iPatch)
                jRec = jRecCell(iCOS,i,j,iPatch)
                
                if(iRec==i.and.jRec==j)iCenCell(i,j,iPatch) = iCOS
                
                call calc_polynomial_square_integration(locPolyDegree(i,j,iPatch),xRel(1,iCOS,i,j,iPatch),xRel(2,iCOS,i,j,iPatch),&
                                                                                  yRel(1,iCOS,i,j,iPatch),yRel(4,iCOS,i,j,iPatch),&
                                                                                  polyCoordCoef(iCOS,1:nRT,i,j,iPatch))
                ! Calculate distance between reconstruction cells and center cell
                dh(iCOS,i,j,iPatch) = sqrt( ( x(cc,iRec,jRec,iPatch) - x(cc,i,j,iPatch) )**2 + ( y(cc,iRec,jRec,iPatch) - y(cc,i,j,iPatch) )**2 )
                !dh(iCOS,i,j,iPatch) = spherical_distance(lat(cc,iRec,jRec,iPatch),lon(cc,iRec,jRec,iPatch),lat(cc,i,j,iPatch),lon(cc,i,j,iPatch),radius)
                !if(iRec==i.and.jRec==j)dh(iCOS,i,j,iPatch)=0
              enddo
              
              ! Calculate reconstruction matrix on edge
              call calc_polynomial_matrix(locPolyDegree(i,j,iPatch),nPointsOnEdge,nRT,xL*recdx,yL*recdy,recMatrixL(:,1:nRT,i,j,iPatch))
              call calc_polynomial_matrix(locPolyDegree(i,j,iPatch),nPointsOnEdge,nRT,xR*recdx,yR*recdy,recMatrixR(:,1:nRT,i,j,iPatch))
              call calc_polynomial_matrix(locPolyDegree(i,j,iPatch),nPointsOnEdge,nRT,xB*recdx,yB*recdy,recMatrixB(:,1:nRT,i,j,iPatch))
              call calc_polynomial_matrix(locPolyDegree(i,j,iPatch),nPointsOnEdge,nRT,xT*recdx,yT*recdy,recMatrixT(:,1:nRT,i,j,iPatch))
              
              do iPOC = 1,nQuadPointsOncell
                xq(iPOC) = x(cqs+iPOC-1,i,j,iPatch) - x(cc,i,j,iPatch)
                yq(iPOC) = y(cqs+iPOC-1,i,j,iPatch) - y(cc,i,j,iPatch)
              enddo
              
              call calc_polynomial_deriv_matrix(locPolyDegree(i,j,iPatch),nQuadPointsOnCell,nRT,&
                                                xq,yq,recMatrixDx(:,1:nRT,i,j,iPatch),recMatrixDy(:,1:nRT,i,j,iPatch))
              
              call calc_polynomial_matrix(locPolyDegree(i,j,iPatch),nQuadPointsOnCell,nRT,xq*recdx,yq*recdy,recMatrixQ(:,1:nRT,i,j,iPatch))
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
      endif
      
      if( trim(reconstruct_scheme)=='Polynomial' )then
        !$OMP PARALLEL DO PRIVATE(i,j,nRC,nxp,nyp,iidx,jidx,xdir,ydir,k,jR,iR,iRec,jRec,existPolyTerm,iCOS,invstat,iPOC,xq,yq) COLLAPSE(3)
        do iPatch = ifs,ife
          do j = jds,jde
            do i = ids,ide
              nRC = nRecCells(i,j,iPatch)
              nxp = maxval(iRecCell(1:nRC,i,j,iPatch)) - minval(iRecCell(1:nRC,i,j,iPatch)) + 1
              nyp = maxval(jRecCell(1:nRC,i,j,iPatch)) - minval(jRecCell(1:nRC,i,j,iPatch)) + 1
              
              ! Pick the cells that not in corner for reconstruction
              if(.not.noCorner(i,j,iPatch))then
                iidx = minloc( abs( x(cc,iRecCell(1:nRC,i,j,iPatch),jRecCell(1:nRC,i,j,iPatch),iPatch) ) )
                jidx = minloc( abs( y(cc,iRecCell(1:nRC,i,j,iPatch),jRecCell(1:nRC,i,j,iPatch),iPatch) ) )
                
                xdir = 1
                ydir = 1
                if(iRecCell(iidx(1),i,j,iPatch)>iRecCell(1,i,j,iPatch))xdir=-1
                if(jRecCell(jidx(2),i,j,iPatch)>jRecCell(1,i,j,iPatch))ydir=-1
                
                k = 0
                do jR = 0,nyp-1
                  do iR = 0,nxp-1
                    k = k + 1
                    iRec = iRecCell(iidx(1),i,j,iPatch) + xdir * iR
                    jRec = jRecCell(jidx(2),i,j,iPatch) + ydir * jR
                    if( .not.inCorner(iRec,jRec,iPatch) )then
                      existPolyTerm(k) = 1
                    else
                      existPolyTerm(k) = 0
                    endif
                  enddo
                enddo
              else
                existPolyTerm = 1
              endif
              
              k    = 0
              iCOS = 0
              do jR = 1,nyp
                do iR = 1,nxp
                  k = k + 1
                  if(existPolyTerm(k)>0)then
                    iCOS = iCOS + 1
                    call calc_rectangle_poly_integration(nxp,nyp,&
                                                         xRel(1,iCOS,i,j,iPatch),xRel(2,iCOS,i,j,iPatch),&
                                                         yRel(1,iCOS,i,j,iPatch),yRel(4,iCOS,i,j,iPatch),&
                                                         Apoly(iCOS,1:nRC,i,j,iPatch),existPolyTerm)
                  endif
                enddo
              enddo
              
              call BRINV(nRC,Apoly(1:nRC,1:nRC,i,j,iPatch),invApoly(1:nRC,1:nRC,i,j,iPatch),invstat)
              if(invstat==0)then
                print*,'Inverse Apoly dost not exist'
                print*,'i,j,iPatch,nRC,nxp,nyp are'
                print*,i,j,iPatch,nRC,nxp,nyp
                stop 'Check BRINV for Special treamtment on boundary cells'
              endif
              
              do iPOC = 1,nQuadPointsOncell
                xq(iPOC) = x(cqs+iPOC-1,i,j,iPatch) - x(cc,i,j,iPatch)
                yq(iPOC) = y(cqs+iPOC-1,i,j,iPatch) - y(cc,i,j,iPatch)
              enddo
            
              ! Calculate reconstruction matrix on edge
              call calc_rectangle_poly_matrix(nxp,nyp,nPointsOnEdge,xL*recdx,yL*recdy,polyMatrixL(:,1:nRC,i,j,iPatch),existPolyTerm)
              call calc_rectangle_poly_matrix(nxp,nyp,nPointsOnEdge,xR*recdx,yR*recdy,polyMatrixR(:,1:nRC,i,j,iPatch),existPolyTerm)
              call calc_rectangle_poly_matrix(nxp,nyp,nPointsOnEdge,xB*recdx,yB*recdy,polyMatrixB(:,1:nRC,i,j,iPatch),existPolyTerm)
              call calc_rectangle_poly_matrix(nxp,nyp,nPointsOnEdge,xT*recdx,yT*recdy,polyMatrixT(:,1:nRC,i,j,iPatch),existPolyTerm)
              
              call calc_rectangle_poly_matrix(nxp,nyp,nQuadPointsOnCell,xq*recdx,yq*recdy,polyMatrixQ(:,1:nRC,i,j,iPatch),existPolyTerm)
              
              polyMatrixL(:,1:nRC,i,j,iPatch) = matmul( polyMatrixL(:,1:nRC,i,j,iPatch), invApoly(1:nRC,1:nRC,i,j,iPatch) )
              polyMatrixR(:,1:nRC,i,j,iPatch) = matmul( polyMatrixR(:,1:nRC,i,j,iPatch), invApoly(1:nRC,1:nRC,i,j,iPatch) )
              polyMatrixB(:,1:nRC,i,j,iPatch) = matmul( polyMatrixB(:,1:nRC,i,j,iPatch), invApoly(1:nRC,1:nRC,i,j,iPatch) )
              polyMatrixT(:,1:nRC,i,j,iPatch) = matmul( polyMatrixT(:,1:nRC,i,j,iPatch), invApoly(1:nRC,1:nRC,i,j,iPatch) )
              
              polyMatrixQ(:,1:nRC,i,j,iPatch) = matmul( polyMatrixQ(:,1:nRC,i,j,iPatch), invApoly(1:nRC,1:nRC,i,j,iPatch) )
              
              ! Calculate derivative matrix
              !xq = 0
              !yq = 0
              call calc_rectangle_poly_deriv_matrix(nxp,nyp,nQuadPointsOnCell,xq,yq,recMatrixDx(:,1:nRC,i,j,iPatch),recMatrixDy(:,1:nRC,i,j,iPatch),existPolyTerm)
              
              recMatrixDx(:,1:nRC,i,j,iPatch) = matmul( recMatrixDx(:,1:nRC,i,j,iPatch), invApoly(1:nRC,1:nRC,i,j,iPatch) ) / dx
              recMatrixDy(:,1:nRC,i,j,iPatch) = matmul( recMatrixDy(:,1:nRC,i,j,iPatch), invApoly(1:nRC,1:nRC,i,j,iPatch) ) / dy
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
      endif
      
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
          
      do iPatch = ifs,ife
        do j = jms,jme
          do i = ims,ime
            ! Bottom Boundary Points
            ghsB(:,i,j,iPatch) = ghs(cbs+0*nPointsOnEdge:cbs+1*nPointsOnEdge-1,i,j,iPatch)
            ! Right Boundary Points
            ghsR(:,i,j,iPatch) = ghs(cbs+1*nPointsOnEdge:cbs+2*nPointsOnEdge-1,i,j,iPatch)
            ! Top Boundary Points
            ghsT(:,i,j,iPatch) = ghs(cbs+2*nPointsOnEdge:cbs+3*nPointsOnEdge-1,i,j,iPatch)
            ! Left Boundary Points
            ghsL(:,i,j,iPatch) = ghs(cbs+3*nPointsOnEdge:cbs+4*nPointsOnEdge-1,i,j,iPatch)
      
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
      
      integer(i_kind) :: iVar,i,j,iPatch,iPOE
      
      call fill_halo(stat%q,qQ(1,:,:,:,:))
      
      qC = stat%q
              
      if(case_num==5)then
        !$OMP PARALLEL DO PRIVATE(i,j) COLLAPSE(3)
        do iPatch = ifs,ife
          do j = jms,jme
            do i = ims,ime
              if(.not.inDomain(i,j,iPatch))then
                qQ(1,:,i,j,iPatch) = qQ(1,:,i,j,iPatch) + sqrtG(cqs:cqe,i,j,iPatch) * ghs(cqs:cqe,i,j,iPatch)
              endif
              qC(1,i,j,iPatch) = qC(1,i,j,iPatch) + sqrtGC(i,j,iPatch) * ghsC(i,j,iPatch)
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
      endif
      
      do iVar = 1,nVar
        call reconstruction(qC(iVar,:,:,:  ),&
                            qL(iVar,:,:,:,:),&
                            qR(iVar,:,:,:,:),&
                            qB(iVar,:,:,:,:),&
                            qT(iVar,:,:,:,:),&
                            qQ(iVar,:,:,:,:))
      enddo
      
      if(case_num==5)then
        !$OMP PARALLEL DO PRIVATE(i,j) COLLAPSE(3)
        do iPatch = ifs,ife
          do j = jdsm1,jdep1
            do i = idsm1,idep1
              qL(1,:,i,j,iPatch) = qL(1,:,i,j,iPatch) - sqrtGL(:,i,j,iPatch) * ghsL(:,i,j,iPatch)
              qR(1,:,i,j,iPatch) = qR(1,:,i,j,iPatch) - sqrtGR(:,i,j,iPatch) * ghsR(:,i,j,iPatch)
              qB(1,:,i,j,iPatch) = qB(1,:,i,j,iPatch) - sqrtGB(:,i,j,iPatch) * ghsB(:,i,j,iPatch)
              qT(1,:,i,j,iPatch) = qT(1,:,i,j,iPatch) - sqrtGT(:,i,j,iPatch) * ghsT(:,i,j,iPatch)
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
              
        qQ(1,:,:,:,:) = qQ(1,:,:,:,:) - sqrtG(cqs:cqe,:,:,:) * ghs(cqs:cqe,:,:,:)
        
        !$OMP PARALLEL DO PRIVATE(i,j) COLLAPSE(3)
        do iPatch = ifs,ife
          do j = jms,jme
            do i = ims,ime
              phit(:,i,j,iPatch) = qQ(1,:,i,j,iPatch) / sqrtG(cqs:cqe,i,j,iPatch) + ghs(cqs:cqe,i,j,iPatch)
              phitC(i,j,iPatch) = cell_quadrature(phit(:,i,j,iPatch))
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
        
        call reconstruction(phitC       ,&
                            dqdx=dphitdx,&
                            dqdy=dphitdy)
        
        !!$OMP PARALLEL DO PRIVATE(i,j) COLLAPSE(3)
        !do iPatch = ifs,ife
        !  do j = jds,jde
        !    do i = ids,ide
        !      !! 4th order
        !      !dphitdx(:,i,j,iPatch) = ( phitC(i-2,j,iPatch) - 8.*phitC(i-1,j,iPatch) + 8.*phitC(i+1,j,iPatch) - phitC(i+2,j,iPatch) )/(12.*dx)
        !      !dphitdy(:,i,j,iPatch) = ( phitC(i,j-2,iPatch) - 8.*phitC(i,j-1,iPatch) + 8.*phitC(i,j+1,iPatch) - phitC(i,j+2,iPatch) )/(12.*dy)
        !      ! 6th order
        !      dphitdx(:,i,j,iPatch) = (-phitC(i-3,j,iPatch) + 9.*phitC(i-2,j,iPatch) - 45.*phitC(i-1,j,iPatch) + 45.*phitC(i+1,j,iPatch) - 9.*phitC(i+2,j,iPatch) + phitC(i+3,j,iPatch) )/(60.*dx)
        !      dphitdy(:,i,j,iPatch) = (-phitC(i,j-3,iPatch) + 9.*phitC(i,j-2,iPatch) - 45.*phitC(i,j-1,iPatch) + 45.*phitC(i,j+1,iPatch) - 9.*phitC(i,j+2,iPatch) + phitC(i,j+3,iPatch) )/(60.*dy)
        !    enddo
        !  enddo
        !enddo
        !!$OMP END PARALLEL DO
      endif
      
      call fill_bdy_flux(qL,qR,qB,qT)
      
      !$OMP PARALLEL
      !$OMP DO PRIVATE(i,j,iPOE,iVar) COLLAPSE(3)
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,idep1
            do iPOE = 1,nPointsOnEdge
              FeP(:,iPOE,i,j,iPatch) = calc_F(sqrtGL(iPOE,i,j,iPatch),matrixIGL(:,:,iPOE,i,j,iPatch),qR(:,iPOE,i-1,j,iPatch),qL(:,iPOE,i,j,iPatch),ghsR(iPOE,i-1,j,iPatch),ghsL(iPOE,i,j,iPatch))
            enddo
            do iVar = 1,nVar
              Fe(iVar,i,j,iPatch) = Gaussian_quadrature_1d(FeP(iVar,:,i,j,iPatch))
            enddo
          enddo
        enddo
      enddo
      !$OMP END DO NOWAIT
      
      !$OMP DO PRIVATE(i,j,iPOE,iVar) COLLAPSE(3)
      do iPatch = ifs,ife
        do j = jds,jdep1
          do i = ids,ide
            do iPOE = 1,nPointsOnEdge
              GeP(:,iPOE,i,j,iPatch) = calc_G(sqrtGB(iPOE,i,j,iPatch),matrixIGB(:,:,iPOE,i,j,iPatch),qT(:,iPOE,i,j-1,iPatch),qB(:,iPOE,i,j,iPatch),ghsT(iPOE,i,j-1,iPatch),ghsB(iPOE,i,j,iPatch))
            enddo
            do iVar = 1,nVar
              Ge(iVar,i,j,iPatch) = Gaussian_quadrature_1d(GeP(iVar,:,i,j,iPatch))
            enddo
          enddo
        enddo
      enddo
      !$OMP END DO NOWAIT
      
      !$OMP DO PRIVATE(i,j) COLLAPSE(3)
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            src(:,i,j,iPatch) = calc_src(sqrtG(cqs:cqe,i,j,iPatch),matrixG(:,:,cqs:cqe,i,j,iPatch),matrixIG(:,:,cqs:cqe,i,j,iPatch),&
                                         qQ(:,:,i,j,iPatch),ghs(cqs:cqe,i,j,iPatch),dphitdx(:,i,j,iPatch),dphitdy(:,i,j,iPatch)    ,&
                                         tanx(cqs:cqe,i,j,iPatch),tany(cqs:cqe,i,j,iPatch)                                         ,&
                                         Coriolis(cqs:cqe,i,j,iPatch),delta(cqs:cqe,i,j,iPatch),iPatch)
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      
      !$OMP PARALLEL DO PRIVATE(i,j,iVar) COLLAPSE(4)
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            do iVar = 1,nVar
              tend%q(iVar,i,j,iPatch) = - ( Fe(iVar,i+1,j,iPatch) - Fe(iVar,i,j,iPatch) ) / dx &
                                        - ( Ge(iVar,i,j+1,iPatch) - Ge(iVar,i,j,iPatch) ) / dy &
                                        + src(iVar,i,j,iPatch)
            enddo
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      !print*,maxval(tend%q(:,ids:ide,jds:jde,ifs:ife)/stat%q(:,ids:ide,jds:jde,ifs:ife)),&
      !       minval(tend%q(:,ids:ide,jds:jde,ifs:ife)/stat%q(:,ids:ide,jds:jde,ifs:ife))
      !
      !call check_halo(stat%q)
      !call check_tend(tend%q)
      !
      !stop 'spatial_operator'
      
    end subroutine spatial_operator
    
    subroutine fill_halo(q,qQ)
      real(r_kind), dimension(             nVar,ims:ime,jms:jme,ifs:ife), intent(inout)          :: q
      real(r_kind), dimension(nQuadPointsOnCell,ims:ime,jms:jme,ifs:ife), intent(inout),optional :: qQ
      
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
              if(present(qQ)) qQ(:,i,j,iPatch) = tgq(1,:,i,j,iPatch)
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
      
      ! For WENO 2D
      real(r_kind), dimension(nStencil1) :: beta1     ! smooth indicator for 1st order stencil
      real(r_kind), dimension(nStencil1) :: sigma
      real(r_kind), dimension(nStencil ) :: beta      ! smooth indicator for high order stencil
      real(r_kind), dimension(nStencil ) :: alpha
      real(r_kind), dimension(nStencil ) :: w
      
      real(r_kind), dimension(1 ) :: a1         ! polynomial coefficients after rematch
      real(r_kind), dimension(9 ) :: a3         ! polynomial coefficients after rematch
      real(r_kind), dimension(25) :: a5         ! polynomial coefficients after rematch
      real(r_kind), dimension(49) :: a7         ! polynomial coefficients after rematch
      
      real(r_kind), dimension(            1 ) :: p1
      real(r_kind), dimension(0:nStencil1,3 ) :: p2
      real(r_kind), dimension(            9 ) :: p3
      real(r_kind), dimension(            25) :: p5
      real(r_kind), dimension(            49) :: p7
      
      real(r_kind), dimension(maxRecTerms) :: p
      
      real(r_kind), dimension(9 ) :: p1_on_3
      real(r_kind), dimension(9 ) :: p2_on_3
      real(r_kind), dimension(25) :: p1_on_5
      real(r_kind), dimension(25) :: p2_on_5
      real(r_kind), dimension(25) :: p3_on_5
      real(r_kind), dimension(49) :: p1_on_7
      real(r_kind), dimension(49) :: p2_on_7
      real(r_kind), dimension(49) :: p3_on_7
      real(r_kind), dimension(49) :: p5_on_7
      
      real(r_kind) :: tau
      real(r_kind) :: sigma_sum
      real(r_kind), parameter :: eps = 1.e-6
      
      
      integer(i_kind) :: iVar,i,j,iPatch,iCOS,iStencil
      integer(i_kind) :: iRec,jRec
      integer(i_kind) :: ic
      integer(i_kind) :: m,n
      
      if(trim(reconstruct_scheme)=='WLS-ENO')then
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
              
              polyCoef(1:n) = WLS_ENO(coordMtx(1:m,1:n),u(1:m),dh(1:m,i,j,iPatch),m,n,ic)
              
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
      elseif(trim(reconstruct_scheme)=='Polynomial')then
        !$OMP PARALLEL DO PRIVATE(j,i,m,iCOS,iRec,jRec,u) COLLAPSE(3)
        do iPatch = ifs,ife
          do j = jds,jde
            do i = ids,ide
              m = nRecCells(i,j,iPatch)
              do iCOS = 1,m
                iRec = iRecCell(iCOS,i,j,iPatch)
                jRec = jRecCell(iCOS,i,j,iPatch)
                
                u(iCOS) = q(iRec,jRec,iPatch)
              enddo
              
              if(present(qL  )) qL(:,i,j,iPatch) = matmul(polyMatrixL (:,1:m,i,j,iPatch),u(1:m))
              if(present(qR  )) qR(:,i,j,iPatch) = matmul(polyMatrixR (:,1:m,i,j,iPatch),u(1:m))
              if(present(qB  )) qB(:,i,j,iPatch) = matmul(polyMatrixB (:,1:m,i,j,iPatch),u(1:m))
              if(present(qT  )) qT(:,i,j,iPatch) = matmul(polyMatrixT (:,1:m,i,j,iPatch),u(1:m))
              
              if(present(qQ  )) qQ(:,i,j,iPatch) = matmul(polyMatrixQ (:,1:m,i,j,iPatch),u(1:m))
              
              if(present(dqdx)) dqdx(:,i,j,iPatch) = matmul(recMatrixDx(:,1:m,i,j,iPatch),u(1:m))
              if(present(dqdy)) dqdy(:,i,j,iPatch) = matmul(recMatrixDy(:,1:m,i,j,iPatch),u(1:m))
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
      elseif(trim(reconstruct_scheme)=='WENO')then
        !$OMP PARALLEL DO PRIVATE(j,i) COLLAPSE(3)
        do iPatch = ifs,ife
          do j = jds,jde
            do i = ids,ide
              if(present(qL  ))call WENO5(qL(1,i,j,iPatch),q(i-recBdy:i+recBdy,j,iPatch),-1)
              if(present(qR  ))call WENO5(qR(1,i,j,iPatch),q(i-recBdy:i+recBdy,j,iPatch), 1)
              if(present(qB  ))call WENO5(qB(1,i,j,iPatch),q(i,j-recBdy:j+recBdy,iPatch),-1)
              if(present(qT  ))call WENO5(qT(1,i,j,iPatch),q(i,j-recBdy:j+recBdy,iPatch), 1)
              
              if(present(qQ  )) qQ(:,i,j,iPatch) = q(i,j,iPatch)
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
      elseif(trim(reconstruct_scheme)=='WENO2D')then
        !$OMP PARALLEL DO PRIVATE(j,i,iStencil,m,iCOS,iRec,jRec,u,polyCoef,beta1,beta,sigma,sigma_sum,a1,a3,a5,a7,        &
        !$OMP                     p1,p2,p3,p5,p7,p1_on_3,p1_on_5,p1_on_7,p2_on_3,p2_on_5,p2_on_7,p3_on_5,p3_on_7,p5_on_7, &
        !$OMP                     tau,alpha,w,p) COLLAPSE(3) ! COLLAPSE must less than 3, iStencil cannot be parallelled
        do iPatch = ifs,ife
          do j = jds,jde
            do i = ids,ide
              do iStencil = 1,nStencil_all
                m = nWENOCells(iStencil,i,j,iPatch)
                do iCOS = 1,m
                  iRec = iWENOCell(iStencil,iCOS,i,j,iPatch)
                  jRec = jWENOCell(iStencil,iCOS,i,j,iPatch)
                  
                  u(iCOS) = q(iRec,jRec,iPatch)
                enddo
                !ic = iCenWENO(iStencil,i,j,iPatch)
                !u(1:m) = u(1:m) / u(ic)
                
                polyCoef(1:m) = matmul(invWENOPoly(iStencil,1:m,1:m,i,j,iPatch),u(1:m))
                
                !if( ( i==26.or.i==65 ) .and. j==45 .and. iPatch==3 .and. iStencil==nStencil1+2 )then
                !  print*,i,j,iPatch
                !  print*,u(1:m)
                !  print*,''
                !endif
                
                if(iStencil<=nStencil1)then
                  p2(iStencil,:) = polyCoef(1:m)
                  beta1(iStencil) = WENO_smooth_indicator_2(polyCoef(1:m))
                elseif(iStencil==nStencil1+1)then
                  do iCOS = 1,nStencil1
                    sigma(iCOS) = ( 1. + ( ( abs(beta1(1)-beta1(2)) + abs(beta1(1)-beta1(3)) &
                                           + abs(beta1(1)-beta1(4)) + abs(beta1(2)-beta1(3)) &
                                           + abs(beta1(2)-beta1(4)) + abs(beta1(3)-beta1(4)) ) / 6. )**2 / ( beta1(iCOS) + eps ) ) / nStencil1
                    !sigma(iCOS) = 0.25 / ( beta1(iStencil) + eps )**2
                  enddo
                  sigma_sum = sum(sigma)
                  sigma = sigma / sigma_sum
                  
                  a1 = polyCoef(1:m)
                  p1 = a1
                  
                  do iCOS = 1,3
                    p2(0,iCOS) = dot_product( sigma, p2(1:nStencil1,iCOS) )
                  enddo
                  beta(1) = WENO_smooth_indicator_2(p2(0,:))
                elseif(iStencil==nStencil1+2)then
                  ! Rematch array for calculating smooth indicator for 3rd order stencil
                  a3 = 0
                  do iCOS = 1,m
                    a3( rematch_idx_3_to_3(iCOS,i,j,iPatch) ) = polyCoef(iCOS)
                  enddo
                  
                  ! Rematch polynomial coefficients and calculate 3rd order polynomial
                  p1         = a1
                  p1_on_3    = 0
                  p1_on_3(1) = p1(1)
                  p2_on_3    = 0
                  do iCOS = 1,3
                    p2_on_3( rematch_idx_2_to_3(iCOS) ) = p2(0,iCOS)
                  enddo
                  p3 = ( a3 - p1_on_3 * r(1,2) ) / r(2,2)
                  !p3 = ( a3 - p2_on_3 * r(1,2) ) / r(2,2)
                  
                  beta(2) = WENO_smooth_indicator_3(p3)
                elseif(iStencil==nStencil1+3)then
                  ! Rematch array for calculating smooth indicator for 5th order stencil
                  a5 = 0
                  do iCOS = 1,m
                    a5( rematch_idx_5_to_5(iCOS,i,j,iPatch) ) =  polyCoef(iCOS)
                  enddo
                  
                  ! Rematch polynomial coefficients and calculate 5th order polynomial
                  p1_on_5    = 0
                  p1_on_5(1) = p1(1)
                  p2_on_5    = 0
                  do iCOS = 1,9
                    p2_on_5( rematch_idx_3_to_5(iCOS) ) = p2_on_3(iCOS)
                  enddo
                  p3_on_5    = 0
                  do iCOS = 1,9
                    p3_on_5( rematch_idx_3_to_5(iCOS) ) = p3(iCOS)
                  enddo
                  p5 = ( a5 - p1_on_5 * r(1,3) - p3_on_5 * r(2,3) ) / r(3,3)
                  !p5 = ( a5 - p2_on_5 * r(1,3) - p3_on_5 * r(2,3) ) / r(3,3)
                  
                  beta(3) = WENO_smooth_indicator_5(p5)
                elseif(iStencil==nStencil1+4)then
                  ! Rematch array for calculating smooth indicator for 7th order stencil
                  a7 = 0
                  do iCOS = 1,m
                    a7( rematch_idx_7_to_7(iCOS,i,j,iPatch) ) =  polyCoef(iCOS)
                  enddo
                  
                  ! Rematch polynomial coefficients and calculate 7th order polynomial
                  p1_on_7    = 0
                  p1_on_7(1) = p1(1)
                  p2_on_7    = 0
                  do iCOS = 1,25
                    p2_on_7( rematch_idx_5_to_7(iCOS) ) = p2_on_5(iCOS)
                  enddo
                  p3_on_7    = 0
                  do iCOS = 1,25
                    p3_on_7( rematch_idx_5_to_7(iCOS) ) = p3_on_5(iCOS)
                  enddo
                  p5_on_7 = 0
                  do iCOS = 1,25
                    p5_on_7( rematch_idx_5_to_7(iCOS) ) = p5(iCOS)
                  enddo
                  p7 = ( a7 - p1_on_7 * r(1,4) - p3_on_7 * r(2,4) - p5_on_7 * r(3,4) ) / r(4,4)
                  !p7 = ( a7 - p2_on_7 * r(1,4) - p3_on_7 * r(2,4) - p5_on_7 * r(3,4) ) / r(4,4)
                  
                  beta(4) = WENO_smooth_indicator_7(p7)
                endif
              enddo
              
              !tau = ( sum( abs( beta(nStencil) - beta(1:nStencil-1) ) ) / ( nStencil - 1. ) )**2
              !
              !do iStencil = 1,nStencil
              !  alpha(iStencil) = r(iStencil,nStencil) * ( 1. + tau / ( beta(iStencil) + eps ) )
              !enddo
              
              do iStencil = 1,nStencil
                alpha(iStencil) = r(iStencil,nStencil) / ( beta(iStencil) + eps )**2
              enddo
              
              w = alpha / sum(alpha)
              
              !if( ( i==26.or.i==65 ) .and. j==45 .and. iPatch==3 )then
              !  print*,i
              !  !print*,beta1
              !  print*,'beta   ',beta
              !  print*,'a3     ',a3
              !  !print*,'p1_on_3',p1_on_3
              !  !print*,'p3     ',p3
              !  !print*,alpha
              !  !print*,w
              !  !print*,r
              !  print*,''
              !endif
              
              if(nStencil==1)p = p1
              if(nStencil==2)p = w(1) * p2_on_3 + w(2) * p3
              if(nStencil==3)p = w(1) * p2_on_5 + w(2) * p3_on_5 + w(3) * p5
              if(nStencil==4)p = w(1) * p2_on_7 + w(2) * p3_on_7 + w(3) * p5_on_7 + w(4) * p7
              
              if(present(qL  )) qL(:,i,j,iPatch) = matmul(polyMatrixL (:,:,i,j,iPatch),p) !* u(ic)
              if(present(qR  )) qR(:,i,j,iPatch) = matmul(polyMatrixR (:,:,i,j,iPatch),p) !* u(ic)
              if(present(qB  )) qB(:,i,j,iPatch) = matmul(polyMatrixB (:,:,i,j,iPatch),p) !* u(ic)
              if(present(qT  )) qT(:,i,j,iPatch) = matmul(polyMatrixT (:,:,i,j,iPatch),p) !* u(ic)
              
              if(present(qQ  )) qQ(:,i,j,iPatch) = matmul(polyMatrixQ (:,:,i,j,iPatch),p) !* u(ic)
              
              if(present(dqdx)) dqdx(:,i,j,iPatch) = matmul(recMatrixDx(:,:,i,j,iPatch),p) !* u(ic)
              if(present(dqdy)) dqdy(:,i,j,iPatch) = matmul(recMatrixDy(:,:,i,j,iPatch),p) !* u(ic)
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
      endif
      
    end subroutine reconstruction
    
    function calc_F(sqrtG,matrixIG,qL,qR,ghsL,ghsR)
      real(r_kind), dimension(nVar) :: calc_F
      real(r_kind)                 , intent(in) :: sqrtG
      real(r_kind), dimension(2,2) , intent(in) :: matrixIG
      real(r_kind), dimension(nVar), intent(in) :: qL
      real(r_kind), dimension(nVar), intent(in) :: qR
      real(r_kind),                  intent(in) :: ghsL
      real(r_kind),                  intent(in) :: ghsR
      
      real(r_kind) :: m ! mach speed
      real(r_kind) :: p ! phi**2 / 2
      
      call AUSM_up(m,p,sqrtG,matrixIG,qL,qR,ghsL,ghsR,1)
      
      calc_F = 0.5 * m * ( qL + qR - sign(1._r_kind,m) * ( qR - qL ) )
      
      calc_F(2) = calc_F(2) + sqrtG * matrixIG(1,1) * p
      calc_F(3) = calc_F(3) + sqrtG * matrixIG(2,1) * p
    end function calc_F
    
    function calc_G(sqrtG,matrixIG,qL,qR,ghsL,ghsR)
      real(r_kind), dimension(nVar) :: calc_G
      real(r_kind)                 , intent(in) :: sqrtG
      real(r_kind), dimension(2,2) , intent(in) :: matrixIG
      real(r_kind), dimension(nVar), intent(in) :: qL
      real(r_kind), dimension(nVar), intent(in) :: qR
      real(r_kind),                  intent(in) :: ghsL
      real(r_kind),                  intent(in) :: ghsR
      
      real(r_kind) :: m ! mach speed
      real(r_kind) :: p ! phi**2 / 2
      
      call AUSM_up(m,p,sqrtG,matrixIG,qL,qR,ghsL,ghsR,2)
      
      calc_G = 0.5 * m * ( qL + qR - sign(1._r_kind,m) * ( qR - qL ) )
      
      calc_G(2) = calc_G(2) + sqrtG * matrixIG(1,2) * p
      calc_G(3) = calc_G(3) + sqrtG * matrixIG(2,2) * p
    end function calc_G
    
    function calc_src(sqrtG,matrixG,matrixIG,q,ghs,dphitdx,dphitdy,x,y,Coriolis,delta,iPatch)
      real(r_kind), dimension(nVar) :: calc_src
      real(r_kind), dimension(     nQuadPointsOnCell), intent(in) :: sqrtG
      real(r_kind), dimension(2, 2,nQuadPointsOnCell), intent(in) :: matrixG
      real(r_kind), dimension(2, 2,nQuadPointsOnCell), intent(in) :: matrixIG
      real(r_kind), dimension(nVar,nQuadPointsOnCell), intent(in) :: q
      real(r_kind), dimension(     nQuadPointsOnCell), intent(in) :: ghs
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
      psi_C(2,:) = Coriolis * (  G12 * phiu + G22 * phiv )
      psi_C(3,:) = Coriolis * ( -G11 * phiu - G12 * phiv )
      
      psi_B(1,:) = 0
      psi_B(2,:) = sqrtG * ghs * ( IG11 * dphitdx + IG12 * dphitdy )
      psi_B(3,:) = sqrtG * ghs * ( IG21 * dphitdx + IG22 * dphitdy )
      
      do iVar = 1,nVar
        calc_src(iVar) = cell_quadrature( psi_M(iVar,:) + psi_C(iVar,:) + psi_B(iVar,:) )
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
    
    subroutine AUSM_up(m,p,sqrtG,matrixIG,qL,qR,ghsL,ghsR,dir)
      real   (r_kind),                  intent(out) :: m
      real   (r_kind),                  intent(out) :: p
      real   (r_kind),                  intent(in ) :: sqrtG
      real   (r_kind), dimension(2,2) , intent(in ) :: matrixIG
      real   (r_kind), dimension(nVar), intent(in ) :: qL
      real   (r_kind), dimension(nVar), intent(in ) :: qR
      real   (r_kind),                  intent(in ) :: ghsL
      real   (r_kind),                  intent(in ) :: ghsR
      integer(i_kind)                 , intent(in ) :: dir ! 1 for x direction, 2 for y direction
      
      real(r_kind),parameter :: Ku    = 0.75
      real(r_kind),parameter :: Kp    = 0.25
      real(r_kind),parameter :: sigma = 1.
      real(r_kind),parameter :: sp    = 1.
      real(r_kind),parameter :: sn    = -1.
      
      real(r_kind) :: phi
      real(r_kind) :: u
      real(r_kind) :: v
      
      real(r_kind) :: phiL
      real(r_kind) :: uL
      real(r_kind) :: cL
      real(r_kind) :: pL
      real(r_kind) :: phiR
      real(r_kind) :: uR
      real(r_kind) :: cR
      real(r_kind) :: pR
      
      real(r_kind) :: a
      real(r_kind) :: ML
      real(r_kind) :: MR
      real(r_kind) :: Mbar2
      real(r_kind) :: Mh
      
      real(r_kind) :: P5MLsp
      real(r_kind) :: P5MRsn
      
      phiL = qL(1) / sqrtG + ghsL
      phiR = qR(1) / sqrtG + ghsR
      ! Convert wind to perpendicular to the edge on sphere
      if(dir==1)then
        uL = qL(2) / qL(1) / matrixIG(1,1) / radius
        uR = qR(2) / qR(1) / matrixIG(1,1) / radius
      elseif(dir==2)then
        uL = qL(3) / qL(1) / matrixIG(2,2) / radius
        uR = qR(3) / qR(1) / matrixIG(2,2) / radius
      endif
      cL   = sqrt( phiL )
      cR   = sqrt( phiR )
      pL   = 0.5 * phiL**2
      pR   = 0.5 * phiR**2
      
      phi = 0.5 * ( phiL + phiR )
      a   = 0.5 * ( cL + cR )
      
      ML = uL / a
      MR = uR / a
      
      Mbar2 = ( uL**2 + uR**2 ) / ( 2. * a**2 )
      
      Mh = M4( ML, sp ) + M4( MR, sn ) - Kp * max( 1. - sigma * Mbar2, 0. ) * ( phiR**2 - phiL**2 ) / ( 2. * phi * a**2 )
      m  = a * Mh
      
      ! Convert wind back to computational space
      if(dir==1)m = m * radius * matrixIG(1,1)
      if(dir==2)m = m * radius * matrixIG(2,2)
      
      P5MLsp = P5(ML,sp)
      P5MRsn = P5(MR,sn)
      
      p = P5MLsp * pL + P5MRsn * pR - Ku * P5MLsp * P5MRsn * ( phiL + phiR ) * a * ( uR - uL )
    end subroutine AUSM_up
    
    function M2(M,signal)
      real(r_kind) :: M2
      real(r_kind) :: M
      real(r_kind) :: signal ! must be 1 or -1
      
      M2 = signal * 0.25 * ( M + signal )**2
      
    end function M2
    
    function M4(M,signal)
      real(r_kind) :: M4
      real(r_kind) :: M
      real(r_kind) :: signal ! must be 1 or -1
      
      real(r_kind),parameter :: beta  = 0.125
      
      if(abs(M)>=1)then
        M4 = 0.5 * ( M + signal * abs(M) )
      else
        M4 = M2( M, signal ) * ( 1. - signal * 16. * beta * M2( M, -signal ) )
      endif
      
    end function M4
    
    function P5(M,signal)
      real(r_kind) :: P5
      real(r_kind) :: M
      real(r_kind) :: signal ! must be 1 or -1
      
      real(r_kind),parameter :: alpha = 0.1875
      
      if(abs(M)>=1)then
        P5 = 0.5 * ( 1. + signal * sign(1._r_kind,M) )
      else
        P5 = M2( M, signal ) * ( ( 2. * signal - M ) - signal * 16.*alpha * M * M2( M, -signal ) )
      endif
      
    end function P5
  
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
          write(123,'(1000e35.16)')(phi(i,j,iPatch),i=ims,ime)
        enddo
      enddo

      do iPatch = ifs,ife
        do j = jms,jme
          write(123,'(1000e35.16)')(us(i,j,iPatch),i=ims,ime)
        enddo
      enddo

      do iPatch = ifs,ife
        do j = jms,jme
          write(123,'(1000e35.16)')(vs(i,j,iPatch),i=ims,ime)
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
            write(123,'(1000e500.16)')(q(iVar,i,j,iPatch),i=ids,ide)
          enddo
        enddo
      enddo
      
      close(123)
      
    end subroutine check_tend
    
end module spatial_operators_mod

