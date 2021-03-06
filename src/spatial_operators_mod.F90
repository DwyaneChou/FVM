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
  
  real(r_kind), dimension(:,:,:,:), allocatable :: src   ! source term
  
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
      integer(i_kind) :: nRC
    
      real(r_kind) :: recCoef
      real(r_kind) :: recdx
      real(r_kind) :: recdy
      real(r_kind) :: recdV
      
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
      
      allocate(src(nVar,ids:ide,jds:jde,ifs:ife))
      
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
      
      !$OMP PARALLEL DO PRIVATE(i,j,iCOS,iRec,jRec,nRC) COLLAPSE(3)
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
            nRC = nRecTerms(i,j,iPatch)
            call calc_polynomial_matrix(locPolyDegree(i,j,iPatch),nPointsOnEdge,nRC,xL*recdx,yL*recdy,recMatrixL(:,1:nRC,i,j,iPatch))
            call calc_polynomial_matrix(locPolyDegree(i,j,iPatch),nPointsOnEdge,nRC,xR*recdx,yR*recdy,recMatrixR(:,1:nRC,i,j,iPatch))
            call calc_polynomial_matrix(locPolyDegree(i,j,iPatch),nPointsOnEdge,nRC,xB*recdx,yB*recdy,recMatrixB(:,1:nRC,i,j,iPatch))
            call calc_polynomial_matrix(locPolyDegree(i,j,iPatch),nPointsOnEdge,nRC,xT*recdx,yT*recdy,recMatrixT(:,1:nRC,i,j,iPatch))
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
      !$OMP PARALLEL DO PRIVATE(j,i,iPOC,xg,yg,nRC) COLLAPSE(4)
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            do iPOC = 1,nQuadPointsOncell
              !xg(iPOC) = ( x(cqs+iPOC-1,i,j,iPatch) - x(cc,i,j,iPatch) ) * recdx
              !yg(iPOC) = ( y(cqs+iPOC-1,i,j,iPatch) - y(cc,i,j,iPatch) ) * recdy
              xg(iPOC) = x(cqs+iPOC-1,i,j,iPatch) - x(cc,i,j,iPatch)
              yg(iPOC) = y(cqs+iPOC-1,i,j,iPatch) - y(cc,i,j,iPatch)
              nRC = nRecTerms(i,j,iPatch)
              call calc_polynomial_deriv_matrix(locPolyDegree(i,j,iPatch),nQuadPointsOnCell,nRC,&
                                                xg(1:nQuadPointsOnCell),yg(1:nQuadPointsOnCell),&
                                                recMatrixDx(:,1:nRC,i,j,iPatch)                ,&
                                                recMatrixDy(:,1:nRC,i,j,iPatch))
            enddo
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
            if(nRC/=0)then
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
              pg = nGhostPointsOnCell(i,j,iPatch)
              xg(1:pg) = ( x(cgs:cgs+pg-1,i,j,iPatch) - x(cc,i,j,iPatch) ) * recdx
              yg(1:pg) = ( y(cgs:cgs+pg-1,i,j,iPatch) - y(cc,i,j,iPatch) ) * recdy
              call calc_rectangle_poly_matrix(nxp,nyp,pg,xg,yg,gstMatrix(1:pg,1:nRC,i,j,iPatch))
              
              gstMatrix(1:pg,1:nRC,i,j,iPatch) = matmul( gstMatrix(1:pg,1:nRC,i,j,iPatch), invA(1:nRC,1:nRC,i,j,iPatch) )
            endif
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      
    end subroutine init_spatial_operator
    
    subroutine spatial_operator(stat,tend)
      type(stat_field), target, intent(inout) :: stat
      type(tend_field), target, intent(inout) :: tend
      
      integer(i_kind) :: iVar,i,j,iPatch
      
      call fill_halo(stat%q)
      
      
      
      !call check_halo(stat%q)
      
      !stop 'spatial_operator'
      
    end subroutine spatial_operator
    
    subroutine fill_halo(q)
      real(r_kind), dimension(nVar,ims:ime,jms:jme,ifs:ife),intent(inout) :: q
      
      real(r_kind), dimension(nVar,maxGhostPointsOnCell,ims:ime,jms:jme,ifs:ife) :: qg
      real(r_kind), dimension(nVar,nTriQuadPointsOnCell,ims:ime,jms:jme,ifs:ife) :: tgq ! value on triangle gaussian quadrature points
      
      integer(i_kind) :: iVar,i,j,iPatch,iPOC
      integer(i_kind) :: iRec,jRec
      integer(i_kind) :: ig,jg,pg,ng
      integer(i_kind) :: m,n
      integer(i_kind) :: igp,itp
      
      real(r_kind), dimension(maxRecCells) :: u
      
      real(r_kind) :: uc,vc,us,vs
      
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            m = nGstRecCells(i,j,iPatch)
            n = nGhostPointsOnCell(i,j,iPatch)
            if(n/=0)then
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
      
      do iPatch = ifs,ife
        do j = jms,jme
          do i = ims,ime
            !if( sum(ghost_n(:,i,j,iPatch)) /= 0 )then
            if( ghost_p(1,i,j,iPatch) /= 0 )then
              do iPOC = 1,nTriQuadPointsOnCell
                ig = ghost_i(iPOC,i,j,iPatch)
                jg = ghost_j(iPOC,i,j,iPatch)
                pg = ghost_p(iPOC,i,j,iPatch)
                ng = ghost_n(iPOC,i,j,iPatch)
                do iVar = 1,nVar
                  tgq(iVar,iPOC,i,j,iPatch) = qg(iVar,ng,ig,jg,pg)
                enddo
                
                igp = cgs + ng   - 1
                itp = cts + iPOC - 1
                
                uc = tgq(2,iPOC,i,j,iPatch) / tgq(1,iPOC,i,j,iPatch)
                vc = tgq(3,iPOC,i,j,iPatch) / tgq(1,iPOC,i,j,iPatch)
                call contravProjPlane2Sphere(us, vs, uc, vc, matrixA (:,:,igp,ig,jg,pg    ))
                call contravProjSphere2Plane(uc, vc, us, vs, matrixIA(:,:,itp,i ,j ,iPatch))
                
                tgq(1,iPOC,i,j,iPatch) = tgq(1,iPOC,i,j,iPatch) / sqrtG(igp,ig,jg,pg) * sqrtG(itp,i,j,iPatch)
                tgq(2,iPOC,i,j,iPatch) = tgq(1,iPOC,i,j,iPatch) * uc
                tgq(3,iPOC,i,j,iPatch) = tgq(1,iPOC,i,j,iPatch) * vc
                
                do iVar = 1,nVar
                  q(iVar,i,j,iPatch) = triangle_quadrature( tgq(iVar,:,i,j,iPatch) )
                enddo
              enddo
            endif
          enddo
        enddo
      enddo
      
    end subroutine fill_halo
    
    function calc_F(q,sqrtG,matrixG)
      real(r_kind), dimension(nVar) :: calc_F
      real(r_kind), dimension(nVar), intent(in) :: q
      real(r_kind)                 , intent(in) :: sqrtG
      real(r_kind), dimension(2,2) , intent(in) :: matrixG
      
      real(r_kind) :: phi,u,v
      real(r_kind) :: G11,G21
      
      G11 = matrixG(1,1)
      G21 = matrixG(2,1)
      
      phi = q(1) / sqrtG
      u   = q(2) / q(1)
      v   = q(3) / q(1)
      
      calc_F(1) = q(2)
      calc_F(2) = q(2) * u + 0.5 * G11 * sqrtG * phi**2
      calc_F(3) = q(2) * v + 0.5 * G21 * sqrtG * phi**2
    end function calc_F
    
    function calc_G(q,sqrtG,matrixG)
      real(r_kind), dimension(nVar) :: calc_G
      real(r_kind), dimension(nVar), intent(in) :: q
      real(r_kind)                 , intent(in) :: sqrtG
      real(r_kind), dimension(2,2) , intent(in) :: matrixG
      
      real(r_kind) :: phi,u,v
      real(r_kind) :: G12,G22
      
      G12 = matrixG(1,2)
      G22 = matrixG(2,2)
      
      phi = q(1) / sqrtG
      u   = q(2) / q(1)
      v   = q(3) / q(1)
      
      calc_G(1) = q(3)
      calc_G(2) = q(3) * u + 0.5 * G12 * sqrtG * phi**2
      calc_G(3) = q(3) * v + 0.5 * G22 * sqrtG * phi**2
    end function calc_G
    
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
          write(123,'(38f35.5)')(phi(i,j,iPatch),i=ims,ime)
        enddo
      enddo

      do iPatch = ifs,ife
        do j = jms,jme
          write(123,'(38f35.5)')(us(i,j,iPatch),i=ims,ime)
        enddo
      enddo

      do iPatch = ifs,ife
        do j = jms,jme
          write(123,'(38f35.5)')(vs(i,j,iPatch),i=ims,ime)
        enddo
      enddo
      
      close(123)
    
    end subroutine check_halo
    
end module spatial_operators_mod

