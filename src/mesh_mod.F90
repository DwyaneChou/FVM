MODULE mesh_mod
  use constants_mod
  use parameters_mod
  use projection_mod
  use reconstruction_mod
  implicit none
  
  ! coordinate
  real(r_kind), dimension(:,:,:,:    ), allocatable :: x        ! central angle on x direction for cells on each patch, unit: radian
  real(r_kind), dimension(:,:,:,:    ), allocatable :: y        ! central angle on y direction for cells on each patch, unit: radian
  real(r_kind), dimension(:,:,:,:    ), allocatable :: lon      ! longitude on cells
  real(r_kind), dimension(:,:,:,:    ), allocatable :: lat      ! latitude on cells
  real(r_kind), dimension(:,:,:,:    ), allocatable :: sqrtG    ! jacobian of Transformation, sqrt(G)
  real(r_kind), dimension(:,:,:      ), allocatable :: sqrtGC   ! jacobian of Transformation, sqrt(G) on Cell
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixG  ! horizontal metric Tensor, which transform covariant vectors to contravariant vectors
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixIG ! horizontal metric Tensor, which transform contravariant vectors to covariant vectors
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixA  ! horizontal metric Tensor, which transform 
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixIA ! horizontal metric Tensor, which transform
  
  real(r_kind), dimension(:,:,:,:    ), allocatable :: sqrtGL   ! jacobian of Transformation, sqrt(G) on cell left edge
  real(r_kind), dimension(:,:,:,:    ), allocatable :: sqrtGR   ! jacobian of Transformation, sqrt(G) on cell right edge
  real(r_kind), dimension(:,:,:,:    ), allocatable :: sqrtGB   ! jacobian of Transformation, sqrt(G) on cell bottom edge
  real(r_kind), dimension(:,:,:,:    ), allocatable :: sqrtGT   ! jacobian of Transformation, sqrt(G) on cell top edge
  
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixIGL   ! matrixIG on cell left edge
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixIGR   ! matrixIG on cell right edge
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixIGB   ! matrixIG on cell bottom edge
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixIGT   ! matrixIG on cell top edge
  
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixAL   ! matrixA on cell left edge
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixAR   ! matrixA on cell right edge
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixAB   ! matrixA on cell bottom edge
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixAT   ! matrixA on cell top edge
  
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixIAL   ! matrixIA on cell left edge
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixIAR   ! matrixIA on cell right edge
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixIAB   ! matrixIA on cell bottom edge
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixIAT   ! matrixIA on cell top edge
  
  ! For cube boundary !
  real(r_kind), dimension(:,:,:,:    ), allocatable :: sqrtGL_adj   ! jacobian of Transformation, sqrt(G) on cell left edge, fix with adjacent panel
  real(r_kind), dimension(:,:,:,:    ), allocatable :: sqrtGR_adj   ! jacobian of Transformation, sqrt(G) on cell right edge, fix with adjacent panel
  real(r_kind), dimension(:,:,:,:    ), allocatable :: sqrtGB_adj   ! jacobian of Transformation, sqrt(G) on cell bottom edge, fix with adjacent panel
  real(r_kind), dimension(:,:,:,:    ), allocatable :: sqrtGT_adj   ! jacobian of Transformation, sqrt(G) on cell top edge, fix with adjacent panel
  
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixAL_adj   ! matrixA on cell left edge, fix with adjacent panel
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixAR_adj   ! matrixA on cell right edge, fix with adjacent panel
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixAB_adj   ! matrixA on cell bottom edge, fix with adjacent panel
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixAT_adj   ! matrixA on cell top edge, fix with adjacent panel
  
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixIAL_adj   ! matrixIA on cell left edge, fix with adjacent panel
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixIAR_adj   ! matrixIA on cell right edge, fix with adjacent panel
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixIAB_adj   ! matrixIA on cell bottom edge, fix with adjacent panel
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixIAT_adj   ! matrixIA on cell top edge, fix with adjacent panel
  
  ! For cube boundary !
  
  real(r_kind), dimension(:,:,:,:    ), allocatable :: Coriolis ! Coriolis parameter
  real(r_kind), dimension(:,:,:,:    ), allocatable :: delta    ! sqrt( 1 + tanx**2 + tany**2 )
  
  real(r_kind), dimension(:,:,:,:    ), allocatable :: sinlon  ! sin(longitude)
  real(r_kind), dimension(:,:,:,:    ), allocatable :: coslon  ! cos(longitude)
  real(r_kind), dimension(:,:,:,:    ), allocatable :: sinlat  ! sin(latitude)
  real(r_kind), dimension(:,:,:,:    ), allocatable :: coslat  ! cos(latitude)
  
  real(r_kind), dimension(:,:,:,:    ), allocatable :: sinx    ! trigonometric function
  real(r_kind), dimension(:,:,:,:    ), allocatable :: cosx    ! trigonometric function
  real(r_kind), dimension(:,:,:,:    ), allocatable :: tanx    ! trigonometric function
  real(r_kind), dimension(:,:,:,:    ), allocatable :: cotx    ! trigonometric function
  real(r_kind), dimension(:,:,:,:    ), allocatable :: secx    ! trigonometric function
  real(r_kind), dimension(:,:,:,:    ), allocatable :: cscx    ! trigonometric function
  
  real(r_kind), dimension(:,:,:,:    ), allocatable :: siny    ! trigonometric function
  real(r_kind), dimension(:,:,:,:    ), allocatable :: cosy    ! trigonometric function
  real(r_kind), dimension(:,:,:,:    ), allocatable :: tany    ! trigonometric function
  real(r_kind), dimension(:,:,:,:    ), allocatable :: coty    ! trigonometric function
  real(r_kind), dimension(:,:,:,:    ), allocatable :: secy    ! trigonometric function
  real(r_kind), dimension(:,:,:,:    ), allocatable :: cscy    ! trigonometric function
  
  ! Ghost points position
  real   (r_kind), dimension(:,:,:,:), allocatable :: ghost_x
  real   (r_kind), dimension(:,:,:,:), allocatable :: ghost_y
  integer(i_kind), dimension(:,:,:,:), allocatable :: ghost_i
  integer(i_kind), dimension(:,:,:,:), allocatable :: ghost_j
  integer(i_kind), dimension(:,:,:,:), allocatable :: ghost_p
  integer(i_kind), dimension(:,:,:,:), allocatable :: ghost_n
  
  real(r_kind), dimension(:,:,:,:), allocatable :: ghs    ! surface height
  real(r_kind), dimension(:,:,:,:), allocatable :: ghsL   ! surface height on left bondary on cell
  real(r_kind), dimension(:,:,:,:), allocatable :: ghsR   ! surface height on right bondary on cell
  real(r_kind), dimension(:,:,:,:), allocatable :: ghsB   ! surface height on bottom bondary on cell
  real(r_kind), dimension(:,:,:,:), allocatable :: ghsT   ! surface height on top bondary on cell
  real(r_kind), dimension(:,:,:  ), allocatable :: ghsC   ! surface height on cell
  
  real(r_kind), dimension(:,:,:    ), allocatable :: areaCell
  
  logical, dimension(:,:,:), allocatable :: inDomain
  logical, dimension(:,:,:), allocatable :: inCorner
  logical, dimension(:,:,:), allocatable :: noCorner
      
  contains
  
  subroutine init_mesh
    integer(i_kind) :: i, j, iPatch , iVar
    integer(i_kind) :: iQP,jQP,countQP
    integer(i_kind) :: iTOC,iVertex1,iVertex2
    integer(i_kind) :: iPOC,iPOE
    integer(i_kind) :: ig, jg, pg, ng
    
    real(r_kind) :: verticesCoord(3,2)
    
    real(r_kind), dimension(:,:), allocatable :: areaCell_temp
    
    ! Allocate arrays in structures
    allocate( x        (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( y        (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( lon      (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( lat      (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    
    allocate( sqrtG    (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( sqrtGC   (                     ims:ime, jms:jme, ifs:ife) )
    allocate( matrixG  (2, 2, nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( matrixIG (2, 2, nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( matrixA  (2, 2, nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( matrixIA (2, 2, nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    
    allocate( sqrtGL   (      nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( sqrtGR   (      nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( sqrtGB   (      nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( sqrtGT   (      nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    
    allocate( matrixIGL(2, 2, nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( matrixIGR(2, 2, nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( matrixIGB(2, 2, nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( matrixIGT(2, 2, nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    
    allocate( matrixAL (2, 2, nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( matrixAR (2, 2, nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( matrixAB (2, 2, nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( matrixAT (2, 2, nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    
    allocate( matrixIAL(2, 2, nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( matrixIAR(2, 2, nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( matrixIAB(2, 2, nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( matrixIAT(2, 2, nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    
    ! For cube boundary !
    allocate( sqrtGL_adj(nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( sqrtGR_adj(nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( sqrtGB_adj(nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( sqrtGT_adj(nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    
    allocate( matrixAL_adj(2, 2, nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( matrixAR_adj(2, 2, nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( matrixAB_adj(2, 2, nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( matrixAT_adj(2, 2, nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    
    allocate( matrixIAL_adj(2, 2, nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( matrixIAR_adj(2, 2, nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( matrixIAB_adj(2, 2, nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( matrixIAT_adj(2, 2, nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    ! For cube boundary !
  
    allocate( Coriolis (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( delta    (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    
    !allocate( sinlon   (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    !allocate( coslon   (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    !allocate( sinlat   (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    !allocate( coslat   (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    
    !allocate( sinx     (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    !allocate( cosx     (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( tanx     (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    !allocate( cotx     (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    !allocate( secx     (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    !allocate( cscx     (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    !allocate( siny     (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    !allocate( cosy     (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( tany     (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    !allocate( coty     (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    !allocate( secy     (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    !allocate( cscy     (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    
    allocate( ghost_x (nQuadPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( ghost_y (nQuadPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( ghost_i (nQuadPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( ghost_j (nQuadPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( ghost_p (nQuadPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( ghost_n (nQuadPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    
    allocate( ghs (nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( ghsL(nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( ghsR(nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( ghsB(nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( ghsT(nPointsOnEdge, ims:ime, jms:jme, ifs:ife) )
    allocate( ghsC(               ims:ime, jms:jme, ifs:ife) )
    
    allocate( areaCell (      ids:ide, jds:jde, ifs:ife) )
    
    allocate(inDomain  (ims:ime,jms:jme,ifs:ife))
    allocate(inCorner  (ims:ime,jms:jme,ifs:ife))
    allocate(noCorner  (ims:ime,jms:jme,ifs:ife))
      
    inDomain(ims:ime,jms:jme,ifs:ife) = .false. 
    inDomain(ids:ide,jds:jde,ifs:ife) = .true.
    
    inCorner(ims:ime,jms:jme,ifs:ife) = .false. 
    inCorner(ims:ids-1,jms:jds-1,ifs:ife) = .true. ! low left
    inCorner(ide+1:ime,jms:jds-1,ifs:ife) = .true. ! low right
    inCorner(ims:ids-1,jde+1:jme,ifs:ife) = .true. ! up left
    inCorner(ide+1:ime,jde+1:jme,ifs:ife) = .true. ! up right
    
    noCorner = .true.
    noCorner(ims:ids+recBdy-1,jms:jds+recBdy-1,ifs:ife) = .false.
    noCorner(ide-recBdy+1:ime,jms:jds+recBdy-1,ifs:ife) = .false.
    noCorner(ims:ids+recBdy-1,jde-recBdy+1:jme,ifs:ife) = .false.
    noCorner(ide-recBdy+1:ime,jde-recBdy+1:jme,ifs:ife) = .false.
    
    ! Init arrays
    ghost_i = 0
    ghost_j = 0
    ghost_p = 0
    ghost_n = 0
    
    nGhostPointsOnCell = 0
    
    !!$OMP PARALLEL DO PRIVATE(j,i,countQP,jQP,iQP,verticesCoord,iTOC,iVertex1,iVertex2,iPOC) COLLAPSE(3)
    do iPatch = ifs,ife
      do j = jms,jme
        do i = ims,ime
          ! Determine points on each cell
          ! cell center
          x(cc   ,i,j,iPatch) = (i - 0.5) * dx + x_min
          y(cc   ,i,j,iPatch) = (j - 0.5) * dy + y_min
          
          ! Low Left corner
          x(ccs+0,i,j,iPatch) = x(cc,i,j,iPatch) - 0.5 * dx
          y(ccs+0,i,j,iPatch) = y(cc,i,j,iPatch) - 0.5 * dy
          ! Low Right corner
          x(ccs+1,i,j,iPatch) = x(cc,i,j,iPatch) + 0.5 * dx
          y(ccs+1,i,j,iPatch) = y(cc,i,j,iPatch) - 0.5 * dy
          ! Up Right corner
          x(ccs+2,i,j,iPatch) = x(cc,i,j,iPatch) + 0.5 * dx
          y(ccs+2,i,j,iPatch) = y(cc,i,j,iPatch) + 0.5 * dy
          ! Up Left corner
          x(ccs+3,i,j,iPatch) = x(cc,i,j,iPatch) - 0.5 * dx
          y(ccs+3,i,j,iPatch) = y(cc,i,j,iPatch) + 0.5 * dy
          
          ! Bottom Boundary Points
          x(cbs+0*nPointsOnEdge:cbs+1*nPointsOnEdge-1,i,j,iPatch) = x(ccs+0,i,j,iPatch) + quad_pos_1d * dx
          y(cbs+0*nPointsOnEdge:cbs+1*nPointsOnEdge-1,i,j,iPatch) = y(ccs+0,i,j,iPatch)
          ! Right Boundary Points
          x(cbs+1*nPointsOnEdge:cbs+2*nPointsOnEdge-1,i,j,iPatch) = x(ccs+1,i,j,iPatch)
          y(cbs+1*nPointsOnEdge:cbs+2*nPointsOnEdge-1,i,j,iPatch) = y(ccs+1,i,j,iPatch) + quad_pos_1d * dy
          ! Top Boundary Points
          x(cbs+2*nPointsOnEdge:cbs+3*nPointsOnEdge-1,i,j,iPatch) = x(ccs+3,i,j,iPatch) + quad_pos_1d * dx
          y(cbs+2*nPointsOnEdge:cbs+3*nPointsOnEdge-1,i,j,iPatch) = y(ccs+3,i,j,iPatch)
          ! Left Boundary Points
          x(cbs+3*nPointsOnEdge:cbs+4*nPointsOnEdge-1,i,j,iPatch) = x(ccs+0,i,j,iPatch)
          y(cbs+3*nPointsOnEdge:cbs+4*nPointsOnEdge-1,i,j,iPatch) = y(ccs+0,i,j,iPatch) + quad_pos_1d * dy
          
          if(quad_opt==1)then
            ! Square Gaussian quadrature points
            countQP = 0
            do jQP = 1,nPointsOnEdge
              do iQP = 1,nPointsOnEdge
                x(cqs+countQP,i,j,iPatch) = x(ccs,i,j,iPatch) + dx * quad_pos_1d(iQP)
                y(cqs+countQP,i,j,iPatch) = y(ccs,i,j,iPatch) + dy * quad_pos_1d(jQP)
                !print*,i,j,iPatch,iQP,jQP,dx*R2D,quad_opt
                !print*,x(cqs+countQP,i,j,iPatch)*R2D,x(ccs,i,j,iPatch)*R2D
                !print*,y(cqs+countQP,i,j,iPatch)*R2D,y(ccs,i,j,iPatch)*R2D
                countQP = countQP + 1
              enddo
            enddo
          elseif(quad_opt==2)then
            ! Triangular quadrature points
            countQP = 0
            verticesCoord(1,:) = (/x(cc,i,j,iPatch), y(cc,i,j,iPatch)/)
            ! Loop 4 triangles on cell
            do iTOC = 1,nEdgesOnCell
              iVertex1 = ccs+iTOC-1
              iVertex2 = ccs+iTOC
              if(iTOC==nEdgesOnCell)iVertex2 = ccs
              verticesCoord(2,:) = (/x(iVertex1,i,j,iPatch), y(iVertex1,i,j,iPatch)/)
              verticesCoord(3,:) = (/x(iVertex2,i,j,iPatch), y(iVertex2,i,j,iPatch)/)
              do iQP = 1,nQuadOrder
                x(cqs+countQP,i,j,iPatch) = dot_product( verticesCoord(:,1), triQuad_pos(iQP,:) )
                y(cqs+countQP,i,j,iPatch) = dot_product( verticesCoord(:,2), triQuad_pos(iQP,:) )
                countQP = countQP + 1
              enddo
            enddo
          endif
          
          ! Calculate parameters on each points on this cell
          do iPOC = cc,cqe
            call pointProjPlane2Sphere(lon(iPOC,i,j,iPatch), lat(iPOC,i,j,iPatch), &
                                       x  (iPOC,i,j,iPatch), y  (iPOC,i,j,iPatch), iPatch)
            
            !sinlon(iPOC,i,j,iPatch) = sin(lon(iPOC,i,j,iPatch))
            !coslon(iPOC,i,j,iPatch) = cos(lon(iPOC,i,j,iPatch))
            
            !sinlat(iPOC,i,j,iPatch) = sin(lat(iPOC,i,j,iPatch))
            !coslat(iPOC,i,j,iPatch) = cos(lat(iPOC,i,j,iPatch))
            
            !sinx(iPOC,i,j,iPatch) = sin(x(iPOC,i,j,iPatch))
            !cosx(iPOC,i,j,iPatch) = cos(x(iPOC,i,j,iPatch))
            tanx(iPOC,i,j,iPatch) = tan(x(iPOC,i,j,iPatch))
            !cotx(iPOC,i,j,iPatch) = 1. / tanx(iPOC,i,j,iPatch)
            !secx(iPOC,i,j,iPatch) = 1. / cosx(iPOC,i,j,iPatch)
            !cscx(iPOC,i,j,iPatch) = 1. / sinx(iPOC,i,j,iPatch)
            !siny(iPOC,i,j,iPatch) = sin(y(iPOC,i,j,iPatch))
            !cosy(iPOC,i,j,iPatch) = cos(y(iPOC,i,j,iPatch))
            tany(iPOC,i,j,iPatch) = tan(y(iPOC,i,j,iPatch))
            !coty(iPOC,i,j,iPatch) = 1. / tany(iPOC,i,j,iPatch)
            !secy(iPOC,i,j,iPatch) = 1. / cosy(iPOC,i,j,iPatch)
            !cscy(iPOC,i,j,iPatch) = 1. / siny(iPOC,i,j,iPatch)
            
            call calc_matrixG (matrixG (:, :, iPOC,i,j,iPatch), x  (iPOC,i,j,iPatch), y  (iPOC,i,j,iPatch))
            call calc_matrixIG(matrixIG(:, :, iPOC,i,j,iPatch), x  (iPOC,i,j,iPatch), y  (iPOC,i,j,iPatch))
            call calc_matrixA (matrixA (:, :, iPOC,i,j,iPatch), lon(iPOC,i,j,iPatch), lat(iPOC,i,j,iPatch), iPatch)
            call calc_matrixIA(matrixIA(:, :, iPOC,i,j,iPatch), lon(iPOC,i,j,iPatch), lat(iPOC,i,j,iPatch), iPatch)
            call calc_Jacobian(sqrtG   (      iPOC,i,j,iPatch), x  (iPOC,i,j,iPatch), y  (iPOC,i,j,iPatch))
            
            Coriolis(iPOC,i,j,iPatch) = 2. * Omega * sin(lat(iPOC,i,j,iPatch))
          enddo
          sqrtGC(i,j,iPatch) = cell_quadrature(sqrtG(cqs:cqe,i,j,iPatch))
        enddo
      enddo
    enddo
    !!$OMP END PARALLEL DO
    
    ! Calculate ghost position
    !!$OMP PARALLEL DO PRIVATE(j,i,countQP,iQP,ig,jg,pg,ng)
    do iPatch = ifs,ife
      do j = jms,jme
        do i = ims,ime
          if( .not.inDomain(i,j,iPatch) )then
            countQP = 0
            do iQP = cqs,cqe
              countQP = countQP + 1
              call psp2ploc(ghost_x(countQP,i,j,iPatch),ghost_y(countQP,i,j,iPatch),ghost_p(countQP,i,j,iPatch),&
                            lon(iQP,i,j,iPatch),lat(iQP,i,j,iPatch))
              ghost_i(countQP,i,j,iPatch) = floor( ( ghost_x(countQP,i,j,iPatch) - x_min ) / dx ) + 1
              ghost_j(countQP,i,j,iPatch) = floor( ( ghost_y(countQP,i,j,iPatch) - y_min ) / dy ) + 1
              
              if( ghost_i(countQP,i,j,iPatch)==ids-1 ) ghost_i(countQP,i,j,iPatch) = ids
              if( ghost_j(countQP,i,j,iPatch)==jds-1 ) ghost_j(countQP,i,j,iPatch) = jds
              if( ghost_i(countQP,i,j,iPatch)==ide+1 ) ghost_i(countQP,i,j,iPatch) = ide
              if( ghost_j(countQP,i,j,iPatch)==jde+1 ) ghost_j(countQP,i,j,iPatch) = jde
              
              ! Set ghost points for interpolation
              ig = ghost_i(countQP,i,j,iPatch)
              jg = ghost_j(countQP,i,j,iPatch)
              pg = ghost_p(countQP,i,j,iPatch)
              nGhostPointsOnCell(ig,jg,pg) = nGhostPointsOnCell(ig,jg,pg) + 1
              ng = nGhostPointsOnCell(ig,jg,pg)
              x(cgs+ng-1,ig,jg,pg) = ghost_x(countQP,i,j,iPatch)
              y(cgs+ng-1,ig,jg,pg) = ghost_y(countQP,i,j,iPatch)
            
              ghost_n(countQP,i,j,iPatch) = ng
            enddo
          endif
        enddo
      enddo
    enddo
    !!$OMP END PARALLEL DO
    
    print*,''
    print*,'Actual max ghost points',maxval(nGhostPointsOnCell)
    
    ! Calculate mesh parameters on ghost points on in-domain cells
    do iPatch = ifs,ife
      do j = jds,jde
        do i = ids,ide
          ng = nGhostPointsOnCell(i,j,iPatch)
          if(ng>0)then
            ! Calculate parameters on each points on this cell
            do iPOC = cgs,cgs+ng-1
              call pointProjPlane2Sphere(lon(iPOC,i,j,iPatch), lat(iPOC,i,j,iPatch), &
                                         x  (iPOC,i,j,iPatch), y  (iPOC,i,j,iPatch), iPatch)
              
              !sinlon(iPOC,i,j,iPatch) = sin(lon(iPOC,i,j,iPatch))
              !coslon(iPOC,i,j,iPatch) = cos(lon(iPOC,i,j,iPatch))
              
              !sinlat(iPOC,i,j,iPatch) = sin(lat(iPOC,i,j,iPatch))
              !coslat(iPOC,i,j,iPatch) = cos(lat(iPOC,i,j,iPatch))
              
              !sinx(iPOC,i,j,iPatch) = sin(x(iPOC,i,j,iPatch))
              !cosx(iPOC,i,j,iPatch) = cos(x(iPOC,i,j,iPatch))
              tanx(iPOC,i,j,iPatch) = tan(x(iPOC,i,j,iPatch))
              !cotx(iPOC,i,j,iPatch) = 1. / tanx(iPOC,i,j,iPatch)
              !secx(iPOC,i,j,iPatch) = 1. / cosx(iPOC,i,j,iPatch)
              !cscx(iPOC,i,j,iPatch) = 1. / sinx(iPOC,i,j,iPatch)
              !siny(iPOC,i,j,iPatch) = sin(y(iPOC,i,j,iPatch))
              !cosy(iPOC,i,j,iPatch) = cos(y(iPOC,i,j,iPatch))
              tany(iPOC,i,j,iPatch) = tan(y(iPOC,i,j,iPatch))
              !coty(iPOC,i,j,iPatch) = 1. / tany(iPOC,i,j,iPatch)
              !secy(iPOC,i,j,iPatch) = 1. / cosy(iPOC,i,j,iPatch)
              !cscy(iPOC,i,j,iPatch) = 1. / siny(iPOC,i,j,iPatch)
              
              call calc_matrixG (matrixG (:, :, iPOC,i,j,iPatch), x  (iPOC,i,j,iPatch), y  (iPOC,i,j,iPatch))
              call calc_matrixIG(matrixIG(:, :, iPOC,i,j,iPatch), x  (iPOC,i,j,iPatch), y  (iPOC,i,j,iPatch))
              call calc_matrixA (matrixA (:, :, iPOC,i,j,iPatch), lon(iPOC,i,j,iPatch), lat(iPOC,i,j,iPatch), iPatch)
              call calc_matrixIA(matrixIA(:, :, iPOC,i,j,iPatch), lon(iPOC,i,j,iPatch), lat(iPOC,i,j,iPatch), iPatch)
              call calc_Jacobian(sqrtG   (      iPOC,i,j,iPatch), x  (iPOC,i,j,iPatch), y  (iPOC,i,j,iPatch))
              
              Coriolis(iPOC,i,j,iPatch) = 2. * Omega * sin(lat(iPOC,i,j,iPatch))
            enddo
          endif
        enddo
      enddo
    enddo
    
    do iPatch = ifs,ife
      do j = jms,jme
        do i = ims,ime
          do iPOE = 1,nPointsOnEdge
            sqrtGB(iPOE,i,j,iPatch) = sqrtG(cbs+0*nPointsOnEdge+iPOE-1,i,j,iPatch) 
            sqrtGR(iPOE,i,j,iPatch) = sqrtG(cbs+1*nPointsOnEdge+iPOE-1,i,j,iPatch) 
            sqrtGT(iPOE,i,j,iPatch) = sqrtG(cbs+2*nPointsOnEdge+iPOE-1,i,j,iPatch) 
            sqrtGL(iPOE,i,j,iPatch) = sqrtG(cbs+3*nPointsOnEdge+iPOE-1,i,j,iPatch) 
            
            matrixIGB(:,:,iPOE,i,j,iPatch) = matrixIG(:,:,cbs+0*nPointsOnEdge+iPOE-1,i,j,iPatch)
            matrixIGR(:,:,iPOE,i,j,iPatch) = matrixIG(:,:,cbs+1*nPointsOnEdge+iPOE-1,i,j,iPatch)
            matrixIGT(:,:,iPOE,i,j,iPatch) = matrixIG(:,:,cbs+2*nPointsOnEdge+iPOE-1,i,j,iPatch)
            matrixIGL(:,:,iPOE,i,j,iPatch) = matrixIG(:,:,cbs+3*nPointsOnEdge+iPOE-1,i,j,iPatch)
            
            matrixAB(:,:,iPOE,i,j,iPatch) = matrixA(:,:,cbs+0*nPointsOnEdge+iPOE-1,i,j,iPatch)
            matrixAR(:,:,iPOE,i,j,iPatch) = matrixA(:,:,cbs+1*nPointsOnEdge+iPOE-1,i,j,iPatch)
            matrixAT(:,:,iPOE,i,j,iPatch) = matrixA(:,:,cbs+2*nPointsOnEdge+iPOE-1,i,j,iPatch)
            matrixAL(:,:,iPOE,i,j,iPatch) = matrixA(:,:,cbs+3*nPointsOnEdge+iPOE-1,i,j,iPatch)
            
            matrixIAB(:,:,iPOE,i,j,iPatch) = matrixIA(:,:,cbs+0*nPointsOnEdge+iPOE-1,i,j,iPatch)
            matrixIAR(:,:,iPOE,i,j,iPatch) = matrixIA(:,:,cbs+1*nPointsOnEdge+iPOE-1,i,j,iPatch)
            matrixIAT(:,:,iPOE,i,j,iPatch) = matrixIA(:,:,cbs+2*nPointsOnEdge+iPOE-1,i,j,iPatch)
            matrixIAL(:,:,iPOE,i,j,iPatch) = matrixIA(:,:,cbs+3*nPointsOnEdge+iPOE-1,i,j,iPatch)
          enddo
        enddo
      enddo
    enddo
    
    do iPatch = ifs,ife
      do j = jms,jme
        do i = ims,ime
          do iPOC = 1,nPointsOnCell
            delta(iPOC,i,j,iPatch) = sqrt( 1. + tanx(iPOC,i,j,iPatch)**2 + tany(iPOC,i,j,iPatch)**2 )
          enddo
        enddo
      enddo
    enddo
    
    ! Calculate areaCell
    allocate( areaCell_temp (Nx,Ny) )
    areaCell = 0.
    call EquiangularAllAreas(Nx, areaCell_temp)
    areaCell_temp = areaCell_temp * radius**2
    
    do iPatch = ifs, ife
      areaCell(1:Nx,1:Ny,iPatch) = areaCell_temp
    enddo
    deallocate( areaCell_temp )
    
    print*,''
    print*,'max areaCell sqrtGC diff ratio',maxval(abs(sqrtGC(ids:ide,jds:jde,ifs:ife)*dx*dy-areaCell(ids:ide,jds:jde,ifs:ife))/areaCell(ids:ide,jds:jde,ifs:ife))
    
    deallocate( ghost_x )
    deallocate( ghost_y )
    !deallocate( sinx )
    !deallocate( cosx )
    !deallocate( tanx )
    !deallocate( cotx )
    !deallocate( secx )
    !deallocate( cscx )
    !deallocate( siny )
    !deallocate( cosy )
    !deallocate( tany )
    !deallocate( coty )
    !deallocate( secy )
    !deallocate( cscy )
    
  end subroutine init_mesh
  
  !------------------------------------------------------------------------------
  ! SUBROUTINE EquiangularAllAreas
  !
  ! Description:
  !   Compute the area of all cubed sphere grid cells, storing the results in
  !   a two dimensional array.
  !
  ! Parameters: 
  !   icube - Cell number (Nx or Ny) of the cubed sphere
  !   dA (OUT) - Output array containing the area of all cubed sphere grid cells
  !------------------------------------------------------------------------------
  SUBROUTINE EquiangularAllAreas(icube, dA)
    IMPLICIT NONE
    
    integer(i_kind),                         INTENT(IN)  :: icube
    real(r_kind)   , DIMENSION(icube,icube), INTENT(OUT) :: dA
    
    ! Local variables
    integer(i_kind)                           :: k, k1, k2
    real(r_kind)                              :: a1, a2, a3, a4
    real(r_kind) , DIMENSION(icube+1,icube+1) :: ang
    real(r_kind) , DIMENSION(icube+1)         :: gp
    
    !#ifdef DBG 
    real(r_kind)    :: dbg !DBG
    !#endif
    
    ! Recall that we are using equi-angular spherical gridding
    !   Compute the angle between equiangular cubed sphere projection grid lines.
    DO k = 1, icube+1
      gp(k) = -0.25 * pi + (pi/DBLE(2*(icube))) * DBLE(k-1)
    ENDDO
    
    DO k2=1,icube+1
      DO k1=1,icube+1
        ang(k1,k2) = ACOS(-SIN(gp(k1)) * SIN(gp(k2)))
      ENDDO
    ENDDO
    
    DO k2=1,icube
      DO k1=1,icube
        a1 =      ang(k1  , k2  )
        a2 = pi - ang(k1+1, k2  )
        a3 = pi - ang(k1  , k2+1)
        a4 =      ang(k1+1, k2+1)      
        ! area = r*r*(-2*pi+sum(interior angles))
        DA(k1,k2) = -2.0*pi+a1+a2+a3+a4
      ENDDO
    ENDDO
    
    ! Only for debugging - test consistency
    dbg = 0.0
    DO k2=1,icube
      DO k1=1,icube
        dbg = dbg + DA(k1,k2)
      ENDDO
    ENDDO
    
    !print*,''
    !print*,'total area error     : ', dbg - 4. * pi / 6. !DBG
  END SUBROUTINE EquiangularAllAreas
END MODULE mesh_mod

