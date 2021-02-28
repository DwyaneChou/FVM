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
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixG  ! horizontal metric Tensor, which transform covariant vectors to contravariant vectors
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixIG ! horizontal metric Tensor, which transform contravariant vectors to covariant vectors
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixA  ! horizontal metric Tensor, which transform 
  real(r_kind), dimension(:,:,:,:,:,:), allocatable :: matrixIA ! horizontal metric Tensor, which transform
  
  real(r_kind), dimension(:,:,:,:    ), allocatable :: Coriolis ! Coriolis parameter
  
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
  
  real(r_kind), dimension(:,:,:,:    ), allocatable :: zs    ! surface height
  
  real(r_kind), dimension(:,:,:    ), allocatable :: areaCell
  
  contains
  
  subroutine init_mesh
    integer(i_kind) :: i, j, iPatch , iVar
    integer(i_kind) :: iQP,jQP,countQP
    integer(i_kind) :: iTOC,iVertex1,iVertex2
    integer(i_kind) :: iPOC
    
    real(r_kind) :: verticesCoord(3,2)
    
    real(r_kind), dimension(:,:), allocatable :: areaCell_temp
    
    ! Allocate arrays in structures
    allocate( x        (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( y        (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( lon      (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( lat      (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    
    allocate( sqrtG    (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( matrixG  (2, 2, nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( matrixIG (2, 2, nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( matrixA  (2, 2, nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( matrixIA (2, 2, nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    
    allocate( Coriolis (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    
    allocate( sinlon   (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( coslon   (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( sinlat   (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( coslat   (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    
    allocate( sinx     (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( cosx     (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( tanx     (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( cotx     (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( secx     (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( cscx     (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( siny     (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( cosy     (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( tany     (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( coty     (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( secy     (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( cscy     (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    
    allocate( ghost_x (nTriQuadPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( ghost_y (nTriQuadPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( ghost_i (nTriQuadPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( ghost_j (nTriQuadPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    allocate( ghost_p (nTriQuadPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    
    allocate( zs       (      nPointsOnCell, ims:ime, jms:jme, ifs:ife) )
    
    allocate( areaCell (      ids:ide, jds:jde, ifs:ife) )
    
    !$OMP PARALLEL DO PRIVATE(j,i,countQP,jQP,iQP,verticesCoord,iTOC,iVertex1,iVertex2,iPOC) COLLAPSE(3)
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
          x(cbs+0*nPointsOnEdge:1*nPointsOnEdge,i,j,iPatch) = x(ccs+1,i,j,iPatch) + quad_pos_1d * dx
          y(cbs+0*nPointsOnEdge:1*nPointsOnEdge,i,j,iPatch) = y(ccs+1,i,j,iPatch)
          ! Right Boundary Points
          x(cbs+1*nPointsOnEdge:2*nPointsOnEdge,i,j,iPatch) = x(ccs+2,i,j,iPatch)
          y(cbs+1*nPointsOnEdge:2*nPointsOnEdge,i,j,iPatch) = y(ccs+2,i,j,iPatch) + quad_pos_1d * dy
          ! Top Boundary Points
          x(cbs+2*nPointsOnEdge:3*nPointsOnEdge,i,j,iPatch) = x(ccs+4,i,j,iPatch) + quad_pos_1d * dx
          y(cbs+2*nPointsOnEdge:3*nPointsOnEdge,i,j,iPatch) = y(ccs+4,i,j,iPatch)
          ! Left Boundary Points
          x(cbs+3*nPointsOnEdge:4*nPointsOnEdge,i,j,iPatch) = x(ccs+1,i,j,iPatch)
          y(cbs+3*nPointsOnEdge:4*nPointsOnEdge,i,j,iPatch) = y(ccs+1,i,j,iPatch) + quad_pos_1d * dy
          
          ! Gaussian quadrature points
          countQP = 0
          do jQP = 1,nQuadOrder
            do iQP = 1,nQuadOrder
              x(cqs+countQP,i,j,iPatch) = x(ccs,i,j,iPatch) + dx * quad_pos_1d(iQP)
              y(cqs+countQP,i,j,iPatch) = y(ccs,i,j,iPatch) + dy * quad_pos_1d(jQP)
              countQP = countQP + 1
            enddo
          enddo
          
          ! Triangle quadrature points
          countQP = 0
          verticesCoord(3,:) = (/x(cc,i,j,iPatch), y(cc,i,j,iPatch)/)
          ! Loop 4 triangles on cell
          do iTOC = 1,nEdgesOnCell
            iVertex1 = ccs+iTOC-1
            iVertex2 = ccs+iTOC
            if(iTOC==nEdgesOnCell)iVertex2 = ccs
            verticesCoord(1,:) = (/x(iVertex1,i,j,iPatch), y(iVertex1,i,j,iPatch)/)
            verticesCoord(2,:) = (/x(iVertex2,i,j,iPatch), y(iVertex2,i,j,iPatch)/)
            do iQP = 1,nTriQuadOrder
              x(cts+countQP,i,j,iPatch) = dot_product( verticesCoord(:,1), triQuad_pos(iQP,:) )
              y(cts+countQP,i,j,iPatch) = dot_product( verticesCoord(:,2), triQuad_pos(iQP,:) )
              countQP = countQP + 1
            enddo
          enddo
          
          ! Calculate parameters on each points on this cell
          do iPOC = 1,nPointsOnCell
            call pointProjPlane2Sphere(lon(iPOC,i,j,iPatch), lat(iPOC,i,j,iPatch), &
                                       x  (iPOC,i,j,iPatch), y  (iPOC,i,j,iPatch), iPatch)
            
            sinlon(iPOC,i,j,iPatch) = sin(lon(iPOC,i,j,iPatch))
            coslon(iPOC,i,j,iPatch) = cos(lon(iPOC,i,j,iPatch))
            
            sinlat(iPOC,i,j,iPatch) = sin(lat(iPOC,i,j,iPatch))
            coslat(iPOC,i,j,iPatch) = cos(lat(iPOC,i,j,iPatch))
            
            sinx(iPOC,i,j,iPatch) = sin(x(iPOC,i,j,iPatch))
            cosx(iPOC,i,j,iPatch) = cos(x(iPOC,i,j,iPatch))
            tanx(iPOC,i,j,iPatch) = tan(x(iPOC,i,j,iPatch))
            cotx(iPOC,i,j,iPatch) = 1. / tanx(iPOC,i,j,iPatch)
            secx(iPOC,i,j,iPatch) = 1. / cosx(iPOC,i,j,iPatch)
            cscx(iPOC,i,j,iPatch) = 1. / sinx(iPOC,i,j,iPatch)
            siny(iPOC,i,j,iPatch) = sin(y(iPOC,i,j,iPatch))
            cosy(iPOC,i,j,iPatch) = cos(y(iPOC,i,j,iPatch))
            tany(iPOC,i,j,iPatch) = tan(y(iPOC,i,j,iPatch))
            coty(iPOC,i,j,iPatch) = 1. / tany(iPOC,i,j,iPatch)
            secy(iPOC,i,j,iPatch) = 1. / cosy(iPOC,i,j,iPatch)
            cscy(iPOC,i,j,iPatch) = 1. / siny(iPOC,i,j,iPatch)
            
            call calc_matrixG (matrixG (:, :, iPOC,i,j,iPatch), x  (iPOC,i,j,iPatch), y  (iPOC,i,j,iPatch))
            call calc_matrixIG(matrixIG(:, :, iPOC,i,j,iPatch), x  (iPOC,i,j,iPatch), y  (iPOC,i,j,iPatch))
            call calc_matrixA (matrixA (:, :, iPOC,i,j,iPatch), lon(iPOC,i,j,iPatch), lat(iPOC,i,j,iPatch), iPatch)
            call calc_matrixIA(matrixIA(:, :, iPOC,i,j,iPatch), lon(iPOC,i,j,iPatch), lat(iPOC,i,j,iPatch), iPatch)
            call calc_Jacobian(sqrtG   (      iPOC,i,j,iPatch), x  (iPOC,i,j,iPatch), y  (iPOC,i,j,iPatch))
            
            Coriolis(iPOC,i,j,iPatch) = 2. * Omega * sinlat(iPOC,i,j,iPatch)
          enddo
          
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    
    ! Calculate ghost position
    !$OMP PARALLEL DO PRIVATE(j,i,countQP,iQP) COLLAPSE(3)
    do iPatch = ifs,ife
      ! Bottom
      do j = jms,jds-1
        do i = ims,ime
          countQP = 0
          do iQP = cts,cte
            countQP = countQP + 1
            call psp2ploc(ghost_x(countQP,i,j,iPatch),ghost_y(countQP,i,j,iPatch),ghost_p(countQP,i,j,iPatch),&
                          lon(iQP,i,j,iPatch),lat(iQP,i,j,iPatch))
            ghost_i(countQP,i,j,iPatch) = floor( ( ghost_x(countQP,i,j,iPatch) - x_min ) / dx ) + 1
            ghost_j(countQP,i,j,iPatch) = floor( ( ghost_y(countQP,i,j,iPatch) - y_min ) / dy ) + 1
          enddo
        enddo
      enddo
      
      ! Top
      do j = jde+1,jme
        do i = ims,ime
          countQP = 0
          do iQP = cts,cte
            countQP = countQP + 1
            call psp2ploc(ghost_x(countQP,i,j,iPatch),ghost_y(countQP,i,j,iPatch),ghost_p(countQP,i,j,iPatch),&
                          lon(iQP,i,j,iPatch),lat(iQP,i,j,iPatch))
            ghost_i(countQP,i,j,iPatch) = floor( ( ghost_x(countQP,i,j,iPatch) - x_min ) / dx ) + 1
            ghost_j(countQP,i,j,iPatch) = floor( ( ghost_y(countQP,i,j,iPatch) - y_min ) / dy ) + 1
          enddo
        enddo
      enddo
      
      ! Left
      do j = jds,jde
        do i = ims,ids-1
          countQP = 0
          do iQP = cts,cte
            countQP = countQP + 1
            call psp2ploc(ghost_x(countQP,i,j,iPatch),ghost_y(countQP,i,j,iPatch),ghost_p(countQP,i,j,iPatch),&
                          lon(iQP,i,j,iPatch),lat(iQP,i,j,iPatch))
            ghost_i(countQP,i,j,iPatch) = floor( ( ghost_x(countQP,i,j,iPatch) - x_min ) / dx ) + 1
            ghost_j(countQP,i,j,iPatch) = floor( ( ghost_y(countQP,i,j,iPatch) - y_min ) / dy ) + 1
          enddo
        enddo
      enddo
      
      ! Right
      do j = jds,jde
        do i = ide+1,ime
          countQP = 0
          do iQP = cts,cte
            countQP = countQP + 1
            call psp2ploc(ghost_x(countQP,i,j,iPatch),ghost_y(countQP,i,j,iPatch),ghost_p(countQP,i,j,iPatch),&
                          lon(iQP,i,j,iPatch),lat(iQP,i,j,iPatch))
            ghost_i(countQP,i,j,iPatch) = floor( ( ghost_x(countQP,i,j,iPatch) - x_min ) / dx ) + 1
            ghost_j(countQP,i,j,iPatch) = floor( ( ghost_y(countQP,i,j,iPatch) - y_min ) / dy ) + 1
          enddo
        enddo
      enddo
      
    enddo
    !$OMP END PARALLEL DO
    
    ! Calculate areaCell
    allocate( areaCell_temp (Nx,Ny) )
    areaCell = 0.
    call EquiangularAllAreas(Nx, areaCell_temp)
    areaCell_temp = areaCell_temp * radius**2
    
    do iPatch = ifs, ife
      areaCell(1:Nx,1:Ny,iPatch) = areaCell_temp
    enddo
    deallocate( areaCell_temp )
    
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
