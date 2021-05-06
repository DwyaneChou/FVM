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
  public trouble_cell_indicator,trouble_cell
      
  integer(i_kind),dimension(:,:,:  ), allocatable :: nRecCells ! number of cells for reconstruction
  integer(i_kind),dimension(:,:,:  ), allocatable :: nGstRecCells ! number of cells for ghost point reconstruction
  integer(i_kind),dimension(:,:,:  ), allocatable :: nRecTerms
  
  integer(i_kind),dimension(:,:,:), allocatable :: locPolyDegree ! degree of local reconstruction polynomial
  
  real   (r_kind), dimension(:,:,:,:,:), allocatable :: polyCoordCoef
  
  integer(i_kind), dimension(:,:,:,:), allocatable :: iRecCell ! x index of reconstruction cells
  integer(i_kind), dimension(:,:,:,:), allocatable :: jRecCell ! y index of reconstruction cells
  
  integer(i_kind), dimension(:,:,:,:), allocatable :: iGstCell ! x index of ghost reconstruction cells
  integer(i_kind), dimension(:,:,:,:), allocatable :: jGstCell ! y index of ghost reconstruction cells
  
  integer(i_kind), dimension(:,:,:,:), allocatable :: iWENOCell ! x index of WENO reconstruction cells
  integer(i_kind), dimension(:,:,:,:), allocatable :: jWENOCell ! y index of WENO reconstruction cells
  
  real   (r_kind), dimension(:,:,:,:,:,:), allocatable :: WENOPoly ! polynomial coefficients for WENO2D
  real   (r_kind), dimension(:,:,:,:,:,:), allocatable :: invWENOPoly ! polynomial coefficients for WENO2D
  
  logical, dimension(:,:,:,:), allocatable :: trouble_cell
  
  real   (r_kind), dimension(:,:), allocatable :: existWENOTerm
  
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
  
      real(r_kind), dimension(:,:,:,:,:), allocatable :: AGst
      real(r_kind), dimension(:,:,:,:,:), allocatable :: invAGst
  
      real(r_kind), dimension(:,:,:,:,:), allocatable :: Apoly
      real(r_kind), dimension(:,:,:,:,:), allocatable :: invApoly
  
      real(r_kind), dimension(:,:,:,:,:), allocatable :: xRel ! relative x coordinate of reconstruction cells
      real(r_kind), dimension(:,:,:,:,:), allocatable :: yRel ! relative y coordinate of reconstruction cells
      
      real(r_kind), dimension(:,:,:,:,:), allocatable :: xGst ! relative x coordinate of Ghost reconstruction cells
      real(r_kind), dimension(:,:,:,:,:), allocatable :: yGst ! relative y coordinate of Ghost reconstruction cells
      
      real(r_kind), dimension(:), allocatable :: existPolyTerm

      real(r_kind), dimension(:), allocatable :: xg
      real(r_kind), dimension(:), allocatable :: yg
      
      integer(i_kind), dimension(:), allocatable :: lack
      
      integer(i_kind) :: i,j,k,iPatch
      integer(i_kind) :: iCOS ! indices of Cells On Stencils
      integer(i_kind) :: iPOC ! indices of points on cell
      integer(i_kind) :: iType ! indices of WENO type
      integer(i_kind) :: ic
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
      
      allocate(iRecCell  (maxRecCells,ims:ime,jms:jme,ifs:ife))
      allocate(jRecCell  (maxRecCells,ims:ime,jms:jme,ifs:ife))
      
      allocate(iGstCell  (maxRecCells,ims:ime,jms:jme,ifs:ife))
      allocate(jGstCell  (maxRecCells,ims:ime,jms:jme,ifs:ife))
      
      allocate(xRel(4,maxRecCells,ids:ide,jds:jde,ifs:ife))
      allocate(yRel(4,maxRecCells,ids:ide,jds:jde,ifs:ife))
      
      allocate(xGst(4,maxRecCells,ids:ide,jds:jde,ifs:ife))
      allocate(yGst(4,maxRecCells,ids:ide,jds:jde,ifs:ife))

      allocate(xg(nPointsOnCell))
      allocate(yg(nPointsOnCell))
      
      if(trim(reconstruct_scheme)=='Polynomial')then
        allocate(polyMatrixL(nPointsOnEdge    ,maxRecTerms,ids:ide,jds:jde,ifs:ife))
        allocate(polyMatrixR(nPointsOnEdge    ,maxRecTerms,ids:ide,jds:jde,ifs:ife))
        allocate(polyMatrixB(nPointsOnEdge    ,maxRecTerms,ids:ide,jds:jde,ifs:ife))
        allocate(polyMatrixT(nPointsOnEdge    ,maxRecTerms,ids:ide,jds:jde,ifs:ife))
        allocate(polyMatrixQ(nQuadPointsOnCell,maxRecTerms,ids:ide,jds:jde,ifs:ife))
      
        allocate(Apoly   (maxRecCells,maxRecCells,ids:ide,jds:jde,ifs:ife))
        allocate(invApoly(maxRecCells,maxRecCells,ids:ide,jds:jde,ifs:ife))
      
        allocate(existPolyTerm(maxRecCells))
      elseif(trim(reconstruct_scheme)=='WENO')then
        allocate(wenoType(ids:ide,jds:jde,ifs:ife))
        
        allocate(iWENOCell(maxRecCells,ims:ime,jms:jme,ifs:ife))
        allocate(jWENOCell(maxRecCells,ims:ime,jms:jme,ifs:ife))
        
        allocate(lack(nWenoLack))
        
        allocate(xDirWENO(ids:ide,jds:jde))
        allocate(yDirWENO(ids:ide,jds:jde))
        
        if(use_trouble_cell_indicator)allocate(trouble_cell(nVar,ims:ime,jms:jme,ifs:ife))
      endif
      
      allocate(recMatrixDx(nQuadPointsOnCell,maxRecTerms,ids:ide,jds:jde,ifs:ife))
      allocate(recMatrixDy(nQuadPointsOnCell,maxRecTerms,ids:ide,jds:jde,ifs:ife))
      
      allocate(AGst   (maxRecCells,maxRecCells,ids:ide,jds:jde,ifs:ife))
      allocate(invAGst(maxRecCells,maxRecCells,ids:ide,jds:jde,ifs:ife))
      
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
      
      if(trim(reconstruct_scheme)=='WENO')then
        do iPatch = ifs,ife
          do j = jds,jde
            do i = ids,ide
              iCOS = 0
              iPOC = 0
              wenoType(i,j,iPatch) = 1
              
              if(iPatch==1)then
                xDirWENO(i,j) = 1
                yDirWENO(i,j) = 1
                if( i<ids+recBdy.and.j<jds+recBdy )then ! low left corner
                  xDirWENO(i,j) = -1
                  yDirWENO(i,j) = -1
                elseif(i>ide-recBdy.and.j<jds+recBdy)then ! low right corner
                  xDirWENO(i,j) = 1
                  yDirWENO(i,j) = -1
                !elseif(i>ide-recBdy.and.j>jde-recBdy)then ! up right corner
                !  xDirWENO(i,j) = 1
                !  yDirWENO(i,j) = 1
                elseif(i<ids+recBdy.and.j>jde-recBdy)then ! up left corner
                  xDirWENO(i,j) = -1
                  yDirWENO(i,j) = 1
                endif
              endif
              
              xdir = xDirWENO(i,j)
              ydir = yDirWENO(i,j)
              
              lack = 0
              do jRec = -recBdy,recBdy
                do iRec = -recBdy,recBdy
                  iCOS = iCOS + 1
                  iR = i + xdir * iRec
                  jR = j + ydir * jRec
                  
                  if( inCorner(iR,jR,iPatch) )then
                    iPOC = iPOC + 1
                    lack(iPOC) = iCOS
                  endif
                  
                  iWENOCell(iCOS,i,j,iPatch) = iR
                  jWENOCell(iCOS,i,j,iPatch) = jR
                enddo
              enddo
              
              do iType = 2,nWenoType
                if( sum( lack - wenoLack(iType,:) )==0 )then
                  wenoType(i,j,iPatch) = iType
                  exit
                endif
              enddo
              
            enddo
          enddo
        enddo
      endif
      
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
      
      if( trim(reconstruct_scheme)=='Polynomial' )then
        !$OMP PARALLEL DO PRIVATE(i,j,nRC,nxp,nyp,iidx,jidx,xdir,ydir,k,jR,iR,iRec,jRec,existPolyTerm,iCOS,invstat,iPOC) COLLAPSE(3)
        do iPatch = ifs,ife
          do j = jds,jde
            do i = ids,ide
              nRC = nRecCells(i,j,iPatch)
              nxp = maxval(iRecCell(1:nRC,i,j,iPatch)) - minval(iRecCell(1:nRC,i,j,iPatch)) + 1
              nyp = maxval(jRecCell(1:nRC,i,j,iPatch)) - minval(jRecCell(1:nRC,i,j,iPatch)) + 1
              
              ! Pick the cells that not in corner for reconstruction
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
            
              ! Calculate reconstruction matrix on edge
              call calc_rectangle_poly_matrix(nxp,nyp,nPointsOnEdge,xL,yL,polyMatrixL(:,1:nRC,i,j,iPatch),existPolyTerm)
              call calc_rectangle_poly_matrix(nxp,nyp,nPointsOnEdge,xR,yR,polyMatrixR(:,1:nRC,i,j,iPatch),existPolyTerm)
              call calc_rectangle_poly_matrix(nxp,nyp,nPointsOnEdge,xB,yB,polyMatrixB(:,1:nRC,i,j,iPatch),existPolyTerm)
              call calc_rectangle_poly_matrix(nxp,nyp,nPointsOnEdge,xT,yT,polyMatrixT(:,1:nRC,i,j,iPatch),existPolyTerm)
              
              call calc_rectangle_poly_matrix(nxp,nyp,nQuadPointsOnCell,xq,yq,polyMatrixQ(:,1:nRC,i,j,iPatch),existPolyTerm)
              
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
                                                               yGst(1,iCOS,i,j,iPatch),yGst(4,iCOS,i,j,iPatch),AGst(iCOS,1:nRC,i,j,iPatch))
                enddo
              enddo
              
              call BRINV(nRC,AGst(1:nRC,1:nRC,i,j,iPatch),invAGst(1:nRC,1:nRC,i,j,iPatch),invstat)
              if(invstat==0)then
                print*,'Inverse AGst dost not exist'
                print*,'i,j,iPatch,nRC,nxp,nyp are'
                print*,i,j,iPatch,nRC,nxp,nyp
                stop 'Check BRINV for Special treamtment on boundary cells'
              endif
              
              ! Calculate reconstruction matrix on edge
              xg(1:pg) = ( x(cgs:cgs+pg-1,i,j,iPatch) - x(cc,i,j,iPatch) ) * recdx
              yg(1:pg) = ( y(cgs:cgs+pg-1,i,j,iPatch) - y(cc,i,j,iPatch) ) * recdy
              call calc_rectangle_poly_matrix(nxp,nyp,pg,xg,yg,gstMatrix(1:pg,1:nRC,i,j,iPatch))
              
              gstMatrix(1:pg,1:nRC,i,j,iPatch) = matmul( gstMatrix(1:pg,1:nRC,i,j,iPatch), invAGst(1:nRC,1:nRC,i,j,iPatch) )
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
      
      if(reconstruct_scheme=='WENO'.and.use_trouble_cell_indicator)then
        do iVar = 1,nVar
          call reconstruction(qC(iVar,:,:,:  ),&
                              qL(iVar,:,:,:,:),&
                              qR(iVar,:,:,:,:),&
                              qB(iVar,:,:,:,:),&
                              qT(iVar,:,:,:,:),&
                              qQ(iVar,:,:,:,:),&
                              trouble_cell=trouble_cell(iVar,:,:,:))
        enddo
      else
        do iVar = 1,nVar
          call reconstruction(qC(iVar,:,:,:  ),&
                              qL(iVar,:,:,:,:),&
                              qR(iVar,:,:,:,:),&
                              qB(iVar,:,:,:,:),&
                              qT(iVar,:,:,:,:),&
                              qQ(iVar,:,:,:,:))
        enddo
      endif
      
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
        
        !call reconstruction(phitC       ,&
        !                    dqdx=dphitdx,&
        !                    dqdy=dphitdy)
        
        !$OMP PARALLEL DO PRIVATE(i,j) COLLAPSE(3)
        do iPatch = ifs,ife
          do j = jds,jde
            do i = ids,ide
              if(recBdy==1)then
                ! 2nd order
                dphitdx(:,i,j,iPatch) = ( phitC(i+1,j,iPatch) - phitC(i-1,j,iPatch) )/(2.*dx)
                dphitdy(:,i,j,iPatch) = ( phitC(i,j+1,iPatch) - phitC(i,j-1,iPatch) )/(2.*dy)
              elseif(recBdy==2)then
                ! 4th order
                dphitdx(:,i,j,iPatch) = ( phitC(i-2,j,iPatch) - 8.*phitC(i-1,j,iPatch) + 8.*phitC(i+1,j,iPatch) - phitC(i+2,j,iPatch) )/(12.*dx)
                dphitdy(:,i,j,iPatch) = ( phitC(i,j-2,iPatch) - 8.*phitC(i,j-1,iPatch) + 8.*phitC(i,j+1,iPatch) - phitC(i,j+2,iPatch) )/(12.*dy)
              elseif(recBdy>2)then
                ! 6th order
                dphitdx(:,i,j,iPatch) = (-phitC(i-3,j,iPatch) + 9.*phitC(i-2,j,iPatch) - 45.*phitC(i-1,j,iPatch) + 45.*phitC(i+1,j,iPatch) - 9.*phitC(i+2,j,iPatch) + phitC(i+3,j,iPatch) )/(60.*dx)
                dphitdy(:,i,j,iPatch) = (-phitC(i,j-3,iPatch) + 9.*phitC(i,j-2,iPatch) - 45.*phitC(i,j-1,iPatch) + 45.*phitC(i,j+1,iPatch) - 9.*phitC(i,j+2,iPatch) + phitC(i,j+3,iPatch) )/(60.*dy)
              endif
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
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
      
      !$OMP DO PRIVATE(i,j,iVar) COLLAPSE(4)
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
      !$OMP END DO
      !$OMP END PARALLEL
      
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
    
    subroutine reconstruction(q,qL,qR,qB,qT,qQ,dqdx,dqdy,trouble_cell)
      real(r_kind), dimension(                  ims:ime,jms:jme,ifs:ife), intent(in   )          :: q
      real(r_kind), dimension(nPointsOnEdge    ,ims:ime,jms:jme,ifs:ife), intent(  out),optional :: qL
      real(r_kind), dimension(nPointsOnEdge    ,ims:ime,jms:jme,ifs:ife), intent(  out),optional :: qR
      real(r_kind), dimension(nPointsOnEdge    ,ims:ime,jms:jme,ifs:ife), intent(  out),optional :: qB
      real(r_kind), dimension(nPointsOnEdge    ,ims:ime,jms:jme,ifs:ife), intent(  out),optional :: qT
      real(r_kind), dimension(nQuadPointsOnCell,ims:ime,jms:jme,ifs:ife), intent(inout),optional :: qQ
      real(r_kind), dimension(nQuadPointsOnCell,ims:ime,jms:jme,ifs:ife), intent(  out),optional :: dqdx ! x derivative on quadrature points
      real(r_kind), dimension(nQuadPointsOnCell,ims:ime,jms:jme,ifs:ife), intent(  out),optional :: dqdy ! y derivative on quadrature points
      logical     , dimension(                  ims:ime,jms:jme,ifs:ife), intent(in   ),optional :: trouble_cell
      
      real(r_kind), dimension(maxRecCells) :: u
      
      ! For WENO
      real   (r_kind), dimension(nWenoPoints) :: qrec
      real   (r_kind), dimension(maxRecCells) :: qH
      integer(i_kind) :: xdir,ydir
      logical         :: TC
      
      integer(i_kind) :: iVar,i,j,iPatch,iCOS,iCell,iStencil
      integer(i_kind) :: iRec,jRec
      integer(i_kind) :: ic
      integer(i_kind) :: m,n
      
      if(trim(reconstruct_scheme)=='WENO')then
        !$OMP PARALLEL DO PRIVATE(i,j,iCOS,iRec,jRec,u,qrec,TC,iCell,qH) COLLAPSE(3)
        do iPatch = ifs,ife
          do j = jds,jde
            do i = ids,ide
              do iCOS = 1,maxRecCells
                iRec = iWENOCell(iCOS,i,j,iPatch)
                jRec = jWENOCell(iCOS,i,j,iPatch)
                u(iCOS) = q(iRec,jRec,iPatch)
              enddo
              
              if(present(trouble_cell))then
                if( any( trouble_cell(iWENOCell(:,i,j,iPatch),jWENOCell(:,i,j,iPatch),iPatch) ) )then
                  TC = .true.
                else
                  TC = .false.
                endif
              else
                TC = .true.
              endif
              
              if(TC)then
                call WENO(qrec,u,wenoType(i,j,iPatch))
              else
               iCOS = 0
               do iCell = 1,maxRecCells
                 if(existPolyTerm(wenoType(i,j,iPatch),iCell)==1)then
                   iCOS = iCOS + 1
                   qH(iCOS) = u(iCell)
                 endif
               enddo
               
               qrec = matmul( iAPoly(wenoType(i,j,iPatch),:,1:iCOS), qH(1:iCOS) )
              endif
              
              if(present(qL)) qL(:,i,j,iPatch) = qrec(wenoLIdx(:,xDirWENO(i,j),yDirWENO(i,j)))
              if(present(qR)) qR(:,i,j,iPatch) = qrec(wenoRIdx(:,xDirWENO(i,j),yDirWENO(i,j)))
              if(present(qB)) qB(:,i,j,iPatch) = qrec(wenoBIdx(:,xDirWENO(i,j),yDirWENO(i,j)))
              if(present(qT)) qT(:,i,j,iPatch) = qrec(wenoTIdx(:,xDirWENO(i,j),yDirWENO(i,j)))
              if(present(qQ)) qQ(:,i,j,iPatch) = qrec(wenoQIdx(:,xDirWENO(i,j),yDirWENO(i,j)))
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
      elseif(trim(reconstruct_scheme)=='WENO3')then
        !$OMP PARALLEL DO PRIVATE(j,i) COLLAPSE(3)
        do iPatch = ifs,ife
          do j = jds,jde
            do i = ids,ide
              if(present(qL  ))call WENO3(qL(1,i,j,iPatch),q(i-recBdy:i+recBdy,j,iPatch),-1)
              if(present(qR  ))call WENO3(qR(1,i,j,iPatch),q(i-recBdy:i+recBdy,j,iPatch), 1)
              if(present(qB  ))call WENO3(qB(1,i,j,iPatch),q(i,j-recBdy:j+recBdy,iPatch),-1)
              if(present(qT  ))call WENO3(qT(1,i,j,iPatch),q(i,j-recBdy:j+recBdy,iPatch), 1)
              
              if(present(qQ  )) qQ(:,i,j,iPatch) = q(i,j,iPatch)
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
      elseif(trim(reconstruct_scheme)=='WENO5')then
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
    
    subroutine trouble_cell_indicator(trouble_cell,q)
      logical     , dimension(ims:ime,jms:jme,ifs:ife), intent(out) :: trouble_cell
      real(r_kind), dimension(ims:ime,jms:jme,ifs:ife), intent(in ) :: q
      
      integer :: i,j,iPatch
      
      integer :: imsp1,imem1,jmsp1,jmem1
      
      real(r_kind) :: c,c_up,c_down,c_left,c_right,ndir
      
      !$OMP PARALLEL DO PRIVATE(i,j,ndir,c_up,c_down,c_left,c_right,c) COLLAPSE(3)
      do iPatch = ifs,ife
        do j = jms,jme
          do i = ims,ime
            ndir = 0
            
            if(j<jme)then
              c_up = abs( q(i,j+1,iPatch) - q(i,j,iPatch) )
              ndir = ndir + 1
            else
              c_up = 0
            endif
            
            if(j>jms)then
              c_down = abs( q(i,j-1,iPatch) - q(i,j,iPatch) )
              ndir = ndir + 1
            else
              c_down = 0
            endif
            
            if(i<ime)then
              c_right = abs( q(i+1,j,iPatch) - q(i,j,iPatch) )
              ndir = ndir + 1
            else
              c_right = 0
            endif
            
            if(i>ims)then
              c_left = abs( q(i-1,j,iPatch) - q(i,j,iPatch) )
              ndir = ndir + 1
            else
              c_left = 0
            endif
            
            !c = ( c_up + c_down + c_right + c_left ) / ( ndir * ( sqrt(2.) * dx )**(dble(stencil_width)/2) * abs( q(i,j,iPatch) ) )
            c = ( c_up + c_down + c_right + c_left ) / ( ndir * ( sqrt(2.) * dx ) * abs( q(i,j,iPatch) ) )
            
            !print*,i,j,iPatch,c
            
            if(c>1)then
              trouble_cell(i,j,iPatch) = .true.
            else
              trouble_cell(i,j,iPatch) = .false.
            endif
            
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine trouble_cell_indicator
    
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

