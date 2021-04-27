    module reconstruction_mod
      use constants_mod
      use parameters_mod
      use qr_solver_mod
      use math_mod
      implicit none
      
      real   (r_kind), dimension(:), allocatable :: quad_pos_1d
      real   (r_kind), dimension(:), allocatable :: quad_wts_1d
      real   (r_kind), dimension(:), allocatable :: quad_wts_2d
      
      real   (r_kind), dimension(:,:), allocatable :: triQuad_pos
      real   (r_kind), dimension(:  ), allocatable :: triQuad_wts
      
      integer(i_kind) :: maxRecCells
      integer(i_kind) :: maxRecTerms
  
      real   (r_kind) :: recCoef
      real   (r_kind) :: recdx
      real   (r_kind) :: recdy
      real   (r_kind) :: recdV

      real(r_kind), dimension(:), allocatable :: xL
      real(r_kind), dimension(:), allocatable :: xR
      real(r_kind), dimension(:), allocatable :: xB
      real(r_kind), dimension(:), allocatable :: xT
      
      real(r_kind), dimension(:), allocatable :: yL
      real(r_kind), dimension(:), allocatable :: yR
      real(r_kind), dimension(:), allocatable :: yB
      real(r_kind), dimension(:), allocatable :: yT
      
      real(r_kind), dimension(:), allocatable :: xq
      real(r_kind), dimension(:), allocatable :: yq
      
      ! For WENO 2D
      integer(i_kind), dimension(:,:,:,:), allocatable :: iCenWENO           ! center cell index on reconstruction stencil for WENO2D
      real   (r_kind), dimension(:,:    ), allocatable :: r                  ! optimal coefficients for WENO 2D
      integer(i_kind), dimension(:,:,:,:,:), allocatable :: rematch_idx_3_to_3_SI ! Just for calculating smooth indicator(SI) for 1st order stencil
      integer(i_kind), dimension(:,:,:,:), allocatable :: rematch_idx_3_to_3
      integer(i_kind), dimension(:,:,:,:), allocatable :: rematch_idx_5_to_5
      integer(i_kind), dimension(:,:,:,:), allocatable :: rematch_idx_7_to_7
      integer(i_kind), dimension(:      ), allocatable :: rematch_idx_2_to_3
      integer(i_kind), dimension(:      ), allocatable :: rematch_idx_3_to_5
      integer(i_kind), dimension(:      ), allocatable :: rematch_idx_5_to_7
      ! For WENO 2D
      
      ! For WENO
      ! WENO stencil type, record the indices of cells in corner
      ! Type 1: []
      ! Type 2: [19,20,24,25]
      ! Type 3: [20,25]
      ! Type 4: [24,25]
      ! Type 5: [25]
      integer(i_kind), parameter :: nWenoStencil= 9
      integer(i_kind), parameter :: nWenoCells  = 9           ! nCells in each WENO reconstruction stencil
      integer(i_kind), parameter :: nWenoTerms  = nWenoCells  ! nTerms in each WENO reconstruction stencil
      integer(i_kind), parameter :: nWenoType   = 5
      integer(i_kind), parameter :: nWenoPoints = 12 ! nPoints on cell those are reconstructed by WENO
      
      integer(i_kind), dimension(:,:,:), allocatable :: wenoType
        
      logical,dimension(nWenoStencil,nWenoType) :: availableStencil = .true.
      
      real   (r_kind), dimension(nWenoPoints,nWenoStencil,nWenoType ) :: wenoCoef
      real   (r_kind), dimension(nWenoStencil,nWenoTerms ,nWenoTerms) :: AWENO      ! (nStencil,nWenoTerms,nWenoTerms)
      real   (r_kind), dimension(nWenoStencil,nWenoTerms ,nWenoTerms) :: invAWENO   ! (nStencil,nWenoTerms,nWenoTerms)
      integer(i_kind), dimension(nWenoStencil,nWenoCells            ) :: wenoidx    ! (nStencil,nWenoCells)
      
      integer(i_kind), dimension(2,-1:1,-1:1) :: wenoLIdx
      integer(i_kind), dimension(2,-1:1,-1:1) :: wenoRIdx
      integer(i_kind), dimension(2,-1:1,-1:1) :: wenoBIdx
      integer(i_kind), dimension(2,-1:1,-1:1) :: wenoTIdx
      integer(i_kind), dimension(4,-1:1,-1:1) :: wenoQIdx
      
      real(r_kind), dimension(nWenoPoints) :: xw
      real(r_kind), dimension(nWenoPoints) :: yw
      
      real(r_kind), dimension(nWenoPoints,nWenoTerms) :: ps ! position polynomial of reconstructed points
      ! For WENO
      
    contains
      subroutine init_reconstruction
        integer(i_kind) :: i,j,k,iStencil,iCOS
        integer(i_kind) :: iQP,jQP,iPOC
        
        integer(i_kind) :: n
        
        real(r_kind) :: ic(nWenoCells),jc(nWenoCells)
        real(r_kind) :: iP,jP
        real(r_kind) :: x_min,x_max,y_min,y_max
        
        real(r_kind), dimension(:,:), allocatable :: quad_wts_tmp_1d
        real(r_kind), dimension(:,:), allocatable :: quad_wts_tmp_2d
        
        real(r_kind), dimension(:,:), allocatable :: quad_wts_tri_tmp_1d
        real(r_kind), dimension(:,:), allocatable :: quad_wts_tri_tmp_2d
        
        real   (r_kind), dimension(:  ), allocatable :: quad_pos_tri_1d
        real   (r_kind), dimension(:  ), allocatable :: quad_wts_tri_1d
        real   (r_kind), dimension(:,:), allocatable :: quad_pos_tri_2d
        real   (r_kind), dimension(:  ), allocatable :: quad_wts_tri_2d
      
        maxRecCells = stencil_width**2
        maxRecTerms = maxRecCells
        
        n = nint( sqrt(real(nQuadOrder,r_kind)) )
        
        allocate(quad_pos_1d(nPointsOnEdge   ))
        allocate(quad_wts_1d(nPointsOnEdge   ))
        allocate(quad_wts_2d(nPointsOnEdge**2))
        
        allocate(triQuad_pos(nQuadOrder,3))
        allocate(triQuad_wts(nQuadOrder  ))
        
        allocate(quad_pos_tri_1d(n   ))
        allocate(quad_wts_tri_1d(n   ))
        allocate(quad_pos_tri_2d(2,nQuadOrder))
        allocate(quad_wts_tri_2d(  nQuadOrder))
        
        allocate(quad_wts_tmp_1d(nPointsOnEdge,1            ))
        allocate(quad_wts_tmp_2d(nPointsOnEdge,nPointsOnEdge))
        
        allocate(quad_wts_tri_tmp_1d(n,1))
        allocate(quad_wts_tri_tmp_2d(n,n))
    
        allocate(xL(nPointsOnEdge))
        allocate(xR(nPointsOnEdge))
        allocate(xB(nPointsOnEdge))
        allocate(xT(nPointsOnEdge))
        
        allocate(yL(nPointsOnEdge))
        allocate(yR(nPointsOnEdge))
        allocate(yB(nPointsOnEdge))
        allocate(yT(nPointsOnEdge))
        
        allocate(xq(nQuadPointsOncell))
        allocate(yq(nQuadPointsOncell))
        
        if(nPointsOnEdge>1)then
          call Gaussian_Legendre(nPointsOnEdge, quad_pos_1d, quad_wts_1d)
        elseif(nPointsOnEdge==1)then
          quad_pos_1d = 0
          quad_wts_1d = 2
        endif
        
        ! Square quadrature
        quad_pos_1d = ( quad_pos_1d + 1. ) / 2.
        quad_wts_1d = quad_wts_1d / 2.
        
        quad_wts_tmp_1d(:,1) = quad_wts_1d
        
        quad_wts_tmp_2d = matmul(quad_wts_tmp_1d,transpose(quad_wts_tmp_1d))
        
        k = 0
        do j = 1,nPointsOnEdge
          do i = 1,nPointsOnEdge
            k = k + 1
            quad_wts_2d(k) = quad_wts_tmp_2d(i,j)
          enddo
        enddo
        
        ! Triangle quadrature from square quadratrue
        if(nQuadOrder==1)then
          ! Order 1
          triQuad_pos(1,:) = (/ 1./3., 1./3., 1./3. /)
          
          triQuad_wts      = (/ 1 /)
        else
          call Gaussian_Legendre(n, quad_pos_tri_1d, quad_wts_tri_1d)
          
          quad_pos_tri_1d = ( quad_pos_tri_1d + 1. ) / 2.
          quad_wts_tri_1d = quad_wts_tri_1d / 2.
          
          quad_wts_tri_tmp_1d(:,1) = quad_wts_tri_1d
          quad_wts_tri_tmp_2d = matmul(quad_wts_tri_tmp_1d,transpose(quad_wts_tri_tmp_1d))
          
          k = 0
          do j = 1,n
            do i = 1,n
              k = k + 1
              quad_pos_tri_2d(1,k) = quad_pos_tri_1d(i)
              quad_pos_tri_2d(2,k) = quad_pos_tri_1d(j)
              
              quad_wts_tri_2d(  k) = quad_wts_tri_tmp_2d(i,j) * ( 1. - quad_pos_tri_2d(1,k) )
              
              triQuad_pos(k,2) = quad_pos_tri_2d(1,k)
              triQuad_pos(k,3) = ( 1. - quad_pos_tri_2d(1,k) ) * quad_pos_tri_2d(2,k)
              triQuad_pos(k,1) = 1. - quad_pos_tri_2d(1,k) - ( 1. - quad_pos_tri_2d(1,k) ) * quad_pos_tri_2d(2,k)
              
              triQuad_wts(k) = quad_wts_tri_2d(  k) * 2.
            enddo
          enddo
        endif
    
        recCoef = 1.
        recdx   = 1. / ( dx * recCoef )
        recdy   = 1. / ( dy * recCoef )
        recdV   = 1. / ( recCoef**2 )
        
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
        
        xL = xL * recdx
        xR = xR * recdx
        xB = xB * recdx
        xT = xT * recdx
        yL = yL * recdy
        yR = yR * recdy
        yB = yB * recdy
        yT = yT * recdy
        
        ! Square Gaussian quadrature points
        iPOC = 0
        do jQP = 1,nPointsOnEdge
          do iQP = 1,nPointsOnEdge
            iPOC = iPOC + 1
            xq(iPOC) = quad_pos_1d(iQP) - 0.5
            yq(iPOC) = quad_pos_1d(jQP) - 0.5
            !print*,iQP,jQP
            !print*,xq(iPOC)
            !print*,yq(iPOC)
            !print*,''
          enddo
        enddo
        
        if(trim(reconstruct_scheme)=='WENO')then
          open(10,file='WENO_COEF.txt',status='old')
          do k = 1,nWenoType
            do j = 1,nWenoPoints
              read(10,*)( wenoCoef( j, i, k ), i = 1,nStencil )
            enddo
            
            !do j = 1,nPointsOnEdge
            !  write(*,'(9e)')wenoCoefL( j, :, k )
            !  write(*,'(9e)')wenoCoefR( j, :, k )
            !  write(*,'(9e)')wenoCoefB( j, :, k )
            !  write(*,'(9e)')wenoCoefT( j, :, k )
            !enddo
            !do j = 1,nQuadPointsOnCell
            !  write(*,'(9e)')wenoCoefQ( j, :, k )
            !enddo
            if( abs(wenoCoef(1,1,k)) < 1.e-14 ) availableStencil(:,k) = .false.
          enddo
          close(10)
          
          ! Calculate invA
          ! middle middle
          ic(1) = 0
          jc(1) = 0
          ! middle left
          ic(2) = -1
          jc(2) = 0
          ! low left
          ic(3) = -1
          jc(3) = -1
          ! low middle
          ic(4) = 0
          jc(4) = -1
          ! low right
          ic(5) = 1
          jc(5) = -1
          ! middle right
          ic(6) = 1
          jc(6) = 0
          ! up right
          ic(7) = 1
          jc(7) = 1
          ! up middle
          ic(8) = 0
          jc(8) = 1
          ! up left
          ic(9) = -1
          jc(9) = 1
          do iStencil = 1,nWenoStencil
            iCOS = 0
            do j = -1,1
              do i = -1,1
                iCOS = iCOS + 1
                iP = i+ic(iStencil);
                jP = j+jc(iStencil);
                x_min = iP-0.5;
                x_max = iP+0.5;
                y_min = jP-0.5;
                y_max = jP+0.5;
                call calc_rectangle_poly_integration(3,3,x_min,x_max,y_min,y_max,AWENO(iStencil,iCOS,:))
              enddo
            enddo
            call BRINV(9,AWENO(iStencil,:,:),invAWENO(iStencil,:,:))
            where(abs(invAWENO)<1.e-15)invAWENO=0
          enddo
          
          wenoidx(1,1:3) = (/ 7, 8, 9/); wenoidx(1,4:6) = (/12,13,14/); wenoidx(1,7:9) = (/17,18,19/);
          wenoidx(2,1:3) = (/ 6, 7, 8/); wenoidx(2,4:6) = (/11,12,13/); wenoidx(2,7:9) = (/16,17,18/);
          wenoidx(3,1:3) = (/ 1, 2, 3/); wenoidx(3,4:6) = (/ 6, 7, 8/); wenoidx(3,7:9) = (/11,12,13/);
          wenoidx(4,1:3) = (/ 2, 3, 4/); wenoidx(4,4:6) = (/ 7, 8, 9/); wenoidx(4,7:9) = (/12,13,14/);
          wenoidx(5,1:3) = (/ 3, 4, 5/); wenoidx(5,4:6) = (/ 8, 9,10/); wenoidx(5,7:9) = (/13,14,15/);
          wenoidx(6,1:3) = (/ 8, 9,10/); wenoidx(6,4:6) = (/13,14,15/); wenoidx(6,7:9) = (/18,19,20/);
          wenoidx(7,1:3) = (/13,14,15/); wenoidx(7,4:6) = (/18,19,20/); wenoidx(7,7:9) = (/23,24,25/);
          wenoidx(8,1:3) = (/12,13,14/); wenoidx(8,4:6) = (/17,18,19/); wenoidx(8,7:9) = (/22,23,24/);
          wenoidx(9,1:3) = (/11,12,13/); wenoidx(9,4:6) = (/16,17,18/); wenoidx(9,7:9) = (/21,22,23/);
          
          xw(1:2 ) = xL
          yw(1:2 ) = yL
          xw(3:4 ) = xR
          yw(3:4 ) = yR
          xw(5:6 ) = xB
          yw(5:6 ) = yB
          xw(7:8 ) = xT
          yw(7:8 ) = yT
          xw(9:12) = xq
          yw(9:12) = yq
          
          call calc_rectangle_poly_matrix(3,3,nWENOPoints,xw,yw,ps)
          
          wenoLIdx(1,1,1) =  1; wenoLIdx(2,1,1) =  2;
          wenoRIdx(1,1,1) =  3; wenoRIdx(2,1,1) =  4;
          wenoBIdx(1,1,1) =  5; wenoBIdx(2,1,1) =  6;
          wenoTIdx(1,1,1) =  7; wenoTIdx(2,1,1) =  8;
          wenoQIdx(1,1,1) =  9; wenoQIdx(2,1,1) = 10;
          wenoQIdx(3,1,1) = 11; wenoQIdx(4,1,1) = 12;
          
          wenoLIdx(1,-1,1) =  3; wenoLIdx(2,-1,1) =  4;
          wenoRIdx(1,-1,1) =  1; wenoRIdx(2,-1,1) =  2;
          wenoBIdx(1,-1,1) =  6; wenoBIdx(2,-1,1) =  5;
          wenoTIdx(1,-1,1) =  8; wenoTIdx(2,-1,1) =  7;
          wenoQIdx(1,-1,1) = 10; wenoQIdx(2,-1,1) =  9;
          wenoQIdx(3,-1,1) = 12; wenoQIdx(4,-1,1) = 11;
          
          wenoLIdx(1,-1,-1) =  4; wenoLIdx(2,-1,-1) =  3;
          wenoRIdx(1,-1,-1) =  2; wenoRIdx(2,-1,-1) =  1;
          wenoBIdx(1,-1,-1) =  8; wenoBIdx(2,-1,-1) =  7;
          wenoTIdx(1,-1,-1) =  6; wenoTIdx(2,-1,-1) =  5;
          wenoQIdx(1,-1,-1) = 12; wenoQIdx(2,-1,-1) = 11;
          wenoQIdx(3,-1,-1) = 10; wenoQIdx(4,-1,-1) =  9;
          
          wenoLIdx(1,1,-1) =  2; wenoLIdx(2,1,-1) =  1;
          wenoRIdx(1,1,-1) =  4; wenoRIdx(2,1,-1) =  3;
          wenoBIdx(1,1,-1) =  7; wenoBIdx(2,1,-1) =  8;
          wenoTIdx(1,1,-1) =  5; wenoTIdx(2,1,-1) =  6;
          wenoQIdx(1,1,-1) = 11; wenoQIdx(2,1,-1) = 12;
          wenoQIdx(3,1,-1) =  9; wenoQIdx(4,1,-1) = 10;
        endif
        
      end subroutine init_reconstruction
      
      function Gaussian_quadrature_1d(q)
        real(r_kind) :: Gaussian_quadrature_1d
        real(r_kind) :: q(nPointsOnEdge)
        
        Gaussian_quadrature_1d = dot_product(quad_wts_1d,q)
        
      end function Gaussian_quadrature_1d
      
      function Gaussian_quadrature_2d(q)
        real(r_kind) :: Gaussian_quadrature_2d
        real(r_kind) :: q(nPointsOnEdge**2)
        
        Gaussian_quadrature_2d = dot_product(quad_wts_2d,q)
        
      end function Gaussian_quadrature_2d
      
      ! Integration field on cell by 4 triangle quadrature
      function cell_quadrature(q)
        real(r_kind) :: cell_quadrature
        real(r_kind) :: q(nQuadPointsOnCell)
        
        integer(i_kind) :: iT
        integer(i_kind) :: is,ie
        
        if(quad_opt==1)then
          cell_quadrature = Gaussian_quadrature_2d(q)
        elseif(quad_opt==2)then
          cell_quadrature = 0
          do iT = 1,nEdgesOnCell
            is = 1 + ( it - 1 ) * nQuadOrder
            ie = it * nQuadOrder
            cell_quadrature = cell_quadrature + dot_product( triQuad_wts,q(is:ie) )
          enddo
          cell_quadrature = cell_quadrature / real(nEdgesOnCell,r_kind)
        endif
    
      end function cell_quadrature
    
      subroutine WENO(qrec,q,wenoType)
        real   (r_kind), intent(out) :: qrec(nWenoPoints)
        real   (r_kind), intent(in ) :: q   (25)
        integer(i_kind), intent(in ) :: wenoType
        
        real   (r_kind), parameter :: eps   = 1.e-15
        real   (r_kind), parameter :: theta = 3
        
        integer(i_kind) :: iStencil,jStencil,iCell,iPoint,iTerm,iCount
        
        real(r_kind), dimension(nWenoCells ) :: qC
        
        real(r_kind), dimension(nWenoStencil,nWenoTerms ) :: a  ! coef of reconstruction polynomial
        real(r_kind), dimension(nWenoStencil            ) :: SI
        real(r_kind), dimension(nWenoStencil,nWenoPoints) :: w
        real(r_kind), dimension(nWenoPoints ,nWenoTerms ) :: wa ! weighted polynomial coef
        
        real(r_kind), dimension(nWenoStencil) :: rp
        real(r_kind), dimension(nWenoStencil) :: rn
        
        real(r_kind), dimension(nWenoStencil,nWenoPoints) :: wp
        real(r_kind), dimension(nWenoStencil,nWenoPoints) :: wn
        
        real(r_kind), dimension(nWenoPoints ,nWenoTerms ) :: wap ! weighted polynomial coef
        real(r_kind), dimension(nWenoPoints ,nWenoTerms ) :: wan ! weighted polynomial coef
        
        real(r_kind), dimension(nWenoPoints) :: qrecp
        real(r_kind), dimension(nWenoPoints) :: qrecn
        
        real(r_kind) :: sigmap
        real(r_kind) :: sigman
        
        real(r_kind) :: tau
        
        ! Rematch cells on each stencil
        do iStencil = 1,nWenoStencil
          do iCell = 1,nWenoCells
            qC(iCell) = q(wenoidx(iStencil,iCell))
          enddo
          
          a(iStencil,:) = matmul( invAWENO(iStencil,:,:),qC(:) )
          
          if( availableStencil(iStencil,wenoType) )then
            SI(iStencil) = WENO_smooth_indicator_3(a(iStencil,:))
          else
            SI(iStencil) = Inf
          endif
        enddo
        
        ! For WENO-Z
        tau = 0
        iCount = 0
        do jStencil = 1,nWENOStencil-1
          if( availableStencil(jStencil,wenoType) )then
            do iStencil = jStencil+1,nWENOStencil
              if( availableStencil(iStencil,wenoType) )then
                iCount = iCount + 1
                tau = tau + abs( SI(iStencil) - SI(jStencil) )
              endif
            enddo
          endif
        enddo
        tau = tau / iCount
        ! For WENO-Z
        
        do iPoint = 1,nWenoPoints
          if( .not.any(wenoCoef(iPoint,:,wenoType)<0) )then
            !w(:,iPoint) = wenoCoef(iPoint,:,wenoType) / ( SI + eps )**2 ! Origin
            w(:,iPoint) = wenoCoef(iPoint,:,wenoType) * ( 1. + tau / ( SI + eps ) ) ! WENO-Z
            w(:,iPoint) = w(:,iPoint) / sum(w(:,iPoint))
            
            do iTerm = 1,nWenoTerms
              wa(iPoint,iTerm) = dot_product( w(:,iPoint), a(:,iTerm) )
            enddo
            
            qrec(iPoint) = dot_product( wa(iPoint,:), ps(iPoint,:) )
          else
            ! Deal with negative linear weights, accroding to Shi J., Hu C. and Shu C.,2002
            rp = 0.5 * ( wenoCoef(iPoint,:,wenoType) + theta * abs( wenoCoef(iPoint,:,wenoType) ) )
            rn = rp - wenoCoef(iPoint,:,wenoType)
            
            sigmap = sum(rp)
            sigman = sum(rn)
            
            rp = rp / sigmap
            rn = rn / sigman
            
            !wp(:,iPoint) = rp / ( SI + eps )**2 ! Origin
            wp(:,iPoint) = rp * ( 1. + tau / ( SI + eps ) ) ! WENO-Z
            wp(:,iPoint) = wp(:,iPoint) / sum(wp(:,iPoint))
            
            !wn(:,iPoint) = rn / ( SI + eps )**2 ! Origin
            wn(:,iPoint) = rn * ( 1. + tau / ( SI + eps ) ) ! WENO-Z
            wn(:,iPoint) = wn(:,iPoint) / sum(wn(:,iPoint))
            
            do iTerm = 1,nWenoTerms
              wap(iPoint,iTerm) = dot_product( wp(:,iPoint), a(:,iTerm) )
              wan(iPoint,iTerm) = dot_product( wn(:,iPoint), a(:,iTerm) )
            enddo
            
            qrecp(iPoint) = dot_product( wap(iPoint,:), ps(iPoint,:) )
            qrecn(iPoint) = dot_product( wan(iPoint,:), ps(iPoint,:) )
            
            qrec(iPoint) = sigmap * qrecp(iPoint) - sigman * qrecn(iPoint)
          endif
        enddo
        
      end subroutine WENO
      
      function WLS_ENO(A,u,h,m,n,ic,x0)
        ! WLS_ENO
        ! Ax = b -> WAx=Wb -> min|| WAx - Wb || -> x
        ! where A->A, x->WLS_ENO, b->u
        ! WLS_ENO : unknown vector x (the vector of reconstruction polynomial coefficients)
        ! A       : matrix of coordinate coefficient
        ! u       : vector of known values
        ! h       : average edge length of cell
        ! m       : number of known values (volumn integration value on cell, or in other words, number of cells in a stencil)
        ! n       : number of coefficients of reconstruction polynomial
        ! ic      : index of center cell on stencil
        ! x0      : initial value for Conjugate Gradient Method
        real   (r_kind), dimension(n  )             :: WLS_ENO
        integer(i_kind)                , intent(in) :: m
        integer(i_kind)                , intent(in) :: n
        real   (r_kind), dimension(m,n), intent(in) :: A
        real   (r_kind), dimension(m  ), intent(in) :: u
        real   (r_kind), dimension(m  ), intent(in) :: h
        integer(i_kind)                , intent(in) :: ic
        
        real   (r_kind), dimension(n  ), intent(in),optional :: x0
        
        real(r_kind),parameter :: alpha   = 1.5
        !real(r_kind),parameter :: epsilon = 1.e-2 ! 1.e-2 for dx>0.25 degree, 1.e+2 for dx<=0.25 degree
        
        real(r_kind), dimension(m,n) :: WA
        real(r_kind), dimension(m  ) :: Wu
        real(r_kind), dimension(m  ) :: W ! weights on each cells
        real(r_kind), dimension(m  ) :: beta

        real(r_kind), dimension(m  ) :: u_bar
        real(r_kind)                 :: u_avg

        !  For LAPACK only
        real   (r_kind), dimension(m+n) :: work
        integer(i_kind) :: INFO
        
        integer(i_kind) :: i,j,k
        
        u_bar = u
        u_avg = sum( abs(u) ) / m
        if(u_avg/=0) u_bar = ( u - u_avg ) / u_avg
        
        !u_bar = std(u,m)

        do j = 1,m
          beta(j) = ( u_bar(j) - u_bar(ic) )**2 + epsilon * h(j)**2
        enddo
        beta(ic) = minval(beta,beta>0)
        
        W = 1./beta
        W(ic) = alpha * W(ic)
        W = W / sum(W)
        
        !print*,'u_bar'
        !write(*,'(5e)')u_bar
        !print*,'W'
        !write(*,'(5e)')W
        !print*,'h'
        !write(*,'(5e)')h
        !print*,'diff u'
        !write(*,'(5e)')( u_bar - u_bar(ic) )**2
        !print*,'eps*h**2'
        !write(*,'(5e)')epsilon * h**2
        !print*,''
        
        do j = 1,m
          WA(j,:) = W(j) * A(j,:)
          Wu(j  ) = W(j) * u(j  )
        enddo
        
        ! Solver by Tsinghua
        call qr_solver(M,N,WA,Wu,WLS_ENO)
        
        !if(present(x0))then
        !  call qr_solver(M,N,WA,Wu,WLS_ENO,x0)
        !else
        !  call qr_solver(M,N,WA,Wu,WLS_ENO)
        !endif
        
        !! Solver by LAPACK DGELS
        !call DGELS( 'N', M, N, 1, WA, M, Wu, M, WORK, M+N, INFO )
        !WLS_ENO = Wu
        
      end function WLS_ENO
      
      function std(q,m)
        integer(i_kind) :: m
        real(r_kind), dimension(m) :: std
        real(r_kind), dimension(m) :: q
        
        real   (r_kind) :: avg
        real   (r_kind) :: n
        
        n = real(m,r_kind)
        
        avg = sum(q) / n
        
        std = sqrt( sum( ( q - avg )**2 ) / n )
        
        std = ( q - avg ) / std
        
      end function std
      
      subroutine WENO_limiter(Qrec,Q,dir)
        real   (r_kind)              , intent(out) :: Qrec
        real   (r_kind), dimension(5), intent(in ) :: Q
        integer(i_kind)              , intent(in ) :: dir
        
        integer(i_kind), parameter :: nStencil = 3
        real   (r_kind), parameter :: weno_coef(3)  = [0.1, 0.6, 0.3]
        real   (r_kind), parameter :: eps           = 1.E-16
        
        real(r_kind) Qim(nStencil-1)
        real(r_kind) Qip(nStencil-1)
        real(r_kind) Qi
        
        integer(i_kind) iStencil
        
        Qim(2) = Q(1)
        Qim(1) = Q(2)
        Qi     = Q(3)
        Qip(1) = Q(4)
        Qip(2) = Q(5)
        
        if( .not. any(Q==FillValue) )then
          call WENO5(Qrec,Q,dir)
        elseif( any(Q==FillValue) .and. Qi/=FillValue )then
          call WENO3(Qrec,Q(2:4),dir)
        else
          print*,'Center point is FillValue, WENO cannot be implemented.'
        endif
        
      end subroutine WENO_limiter
      
      ! 1D WENO 5th order slope limiter, according to Sun,2015
      ! "A Slope Constrained 4th OrderMulti-Moment Finite Volume Method with WENO Limiter"
      ! and Jiang and Shu, 1996
      subroutine WENO5(Qrec,Q,dir)
        real   (r_kind)              , intent(out) :: Qrec
        real   (r_kind), dimension(5), intent(in ) :: Q
        integer(i_kind)              , intent(in ) :: dir
        
        integer(i_kind), parameter :: nStencil = 3
        real   (r_kind), parameter :: weno_coef(3)  = [0.1, 0.6, 0.3]
        real   (r_kind), parameter :: eps           = 1.E-16
        
        real(r_kind) Qim(nStencil-1)
        real(r_kind) Qip(nStencil-1)
        real(r_kind) Qi
        
        real(r_kind), dimension(nStencil) :: stencil
        real(r_kind), dimension(nStencil) :: coefA
        real(r_kind), dimension(nStencil) :: coefB
        real(r_kind), dimension(nStencil) :: alpha
        real(r_kind), dimension(nStencil) :: beta
        real(r_kind), dimension(nStencil) :: omega
        
        real(r_kind) tau40
        real(r_kind) tau41
        real(r_kind) tau5
        
        integer(i_kind) iStencil
        
        Qim(2) = Q(1)
        Qim(1) = Q(2)
        Qi     = Q(3)
        Qip(1) = Q(4)
        Qip(2) = Q(5)
        
        if(dir>0)then
          stencil (1) =  Qim(2)/3. - 7./6. * Qim(1) + 11./6. * Qi     
          stencil (2) = -Qim(1)/6. + 5./6. * Qi     +  1./3. * Qip(1) 
          stencil (3) =  Qi    /3. + 5./6. * Qip(1) -  1./6. * Qip(2)
          
          coefA(1) = Qim(2) - 2. * Qim(1) + Qi
          coefA(2) = Qim(1) - 2. * Qi     + Qip(1)
          coefA(3) = Qi     - 2. * Qip(1) + Qip(2)
          
          coefB(1) =      Qim(2) - 4. * Qim(1) + 3. * Qi
          coefB(2) =      Qim(1) -      Qip(1)
          coefB(3) = 3. * Qi     - 4. * Qip(1) +      Qip(2)
        elseif(dir<0)then
          stencil (1) =  Qip(2)/3. - 7./6. * Qip(1) + 11./6. * Qi     
          stencil (2) = -Qip(1)/6. + 5./6. * Qi     +  1./3. * Qim(1) 
          stencil (3) =  Qi    /3. + 5./6. * Qim(1) -  1./6. * Qim(2)
          
          coefA(1) = Qip(2) - 2. * Qip(1) + Qi
          coefA(2) = Qip(1) - 2. * Qi     + Qim(1)
          coefA(3) = Qi     - 2. * Qim(1) + Qim(2)
          
          coefB(1) =      Qip(2) - 4. * Qip(1) + 3. * Qi
          coefB(2) =      Qip(1) -      Qim(1)
          coefB(3) = 3. * Qi     - 4. * Qim(1) +      Qim(2)
        endif
        
        beta = coefA**2 * 13. / 12. + coefB**2 * 0.25
        
        !! WENO-Z
        !tau40 = abs( beta(1) - beta(2) )
        !tau41 = abs( beta(2) - beta(3) )
        !tau5  = abs( beta(3) - beta(1) )
        !
        !if( tau40<=minval(beta) .and. tau41>minval(beta) )then
        !  omega = [1./3.,2./3.,0.]
        !elseif( tau40>minval(beta) .and. tau41<minval(beta) )then
        !  omega = [0.,2./3.,1./3.]
        !else
        !  alpha = weno_coef * ( 1. + tau5 / ( eps + beta ) )
        !  omega = alpha / sum(alpha)
        !endif
        !! WENO-Z
        !
        !!! origin WENO
        !!alpha = weno_coef / ( eps + beta )**2
        !!omega = alpha / sum(alpha)
        !!! origin WENO
        
        if(.not.any(Q==FillValue))then
          ! WENO-Z
          tau40 = abs( beta(1) - beta(2) )
          tau41 = abs( beta(2) - beta(3) )
          tau5  = abs( beta(3) - beta(1) )
          
          if( tau40<=minval(beta) .and. tau41>minval(beta) )then
            omega = [1./3.,2./3.,0.]
          elseif( tau40>minval(beta) .and. tau41<minval(beta) )then
            omega = [0.,2./3.,1./3.]
          else
            alpha = weno_coef * ( 1. + tau5 / ( eps + beta ) )
            omega = alpha / sum(alpha)
          endif
          ! WENO-Z
        else
          ! origin WENO
          alpha = weno_coef / ( eps + beta )**2
          omega = alpha / sum(alpha)
          ! origin WENO
        endif
        
        Qrec = dot_product( stencil, omega )
        
      end subroutine WENO5
      
      ! 1D 3th order WENO limiter, accroding to 
      ! T. Lunnet et al., 2017, Monthly Weather Review
      subroutine WENO3(Qrec,Q,dir)
        real   (r_kind)              , intent(out) :: Qrec
        real   (r_kind), dimension(3), intent(in ) :: Q
        integer(i_kind)              , intent(in ) :: dir
        
        integer(i_kind), parameter :: nStencil = 2
        real   (r_kind), parameter :: weno_coef(2)  = [1./3., 2./3.]
        real   (r_kind), parameter :: eps           = 1.E-16
        
        real(r_kind) Qim(nStencil-1)
        real(r_kind) Qip(nStencil-1)
        real(r_kind) Qi
        
        real(r_kind), dimension(nStencil) :: stencil
        real(r_kind), dimension(nStencil) :: coefA
        real(r_kind), dimension(nStencil) :: coefB
        real(r_kind), dimension(nStencil) :: alpha
        real(r_kind), dimension(nStencil) :: beta
        real(r_kind), dimension(nStencil) :: omega
        
        real(r_kind) :: tau
        
        integer(i_kind) iStencil
        
        Qim(1) = Q(1)
        Qi     = Q(2)
        Qip(1) = Q(3)
        
        if(dir>0)then
          stencil (1) = -Qim(1)/2. + 3./2. * Qi
          stencil (2) =  Qi / 2. + Qip(1) / 2.
          
          coefA(1) =-Qim(1) + Qi
          coefA(2) = Qip(1) - Qi
        elseif(dir<0)then
          stencil (1) = -Qip(1)/2. + 3./2. * Qi
          stencil (2) = Qi / 2. + Qim(1) / 2.
          
          coefA(1) =-Qip(1) + Qi
          coefA(2) = Qim(1) - Qi
        endif
        
        beta = coefA**2
        
        !! origin WENO
        !alpha = weno_coef / ( eps + beta )**2
        !! origin WENO
        
        ! WENO-Z
        tau = abs( beta(1) -beta(2) )
        alpha = weno_coef * ( 1. + tau / ( eps + beta ) )
        ! WENO-Z
        
        omega = alpha / sum(alpha)
        
        Qrec = dot_product( stencil, omega )
        
      end subroutine WENO3
      
      subroutine WENO2D(polyCoef,nWenoCells,rematch_idx_3_to_3_SI,rematch_idx_3_to_3,rematch_idx_5_to_5,rematch_idx_7_to_7,p,i,j,iPatch)
        real   (r_kind),dimension(:,:),intent(in ) :: polyCoef
        integer(i_kind),dimension(:  ),intent(in ) :: nWenoCells    ! nWenoCells
        integer(i_kind),dimension(:,:),intent(in ) :: rematch_idx_3_to_3_SI
        integer(i_kind),dimension(:  ),intent(in ) :: rematch_idx_3_to_3
        integer(i_kind),dimension(:  ),intent(in ) :: rematch_idx_5_to_5
        integer(i_kind),dimension(:  ),intent(in ) :: rematch_idx_7_to_7
        real   (r_kind),dimension(:  ),intent(out) :: p
        integer(i_kind), intent(in) :: i,j,iPatch
        
        real(r_kind), dimension(nStencil1) :: beta1     ! smooth indicator for 1st order stencil
        real(r_kind), dimension(nStencil1) :: sigma
        real(r_kind), dimension(nStencil ) :: beta      ! smooth indicator for high order stencil
        real(r_kind), dimension(nStencil ) :: alpha
        real(r_kind), dimension(nStencil ) :: w
        
        real(r_kind), dimension(1 ) :: a1         ! polynomial coefficients after rematch
        real(r_kind), dimension(9 ) :: a3         ! polynomial coefficients after rematch
        real(r_kind), dimension(25) :: a5         ! polynomial coefficients after rematch
        real(r_kind), dimension(49) :: a7         ! polynomial coefficients after rematch
        
        real(r_kind), dimension(          1 ) :: p1
        real(r_kind), dimension(nStencil1,9 ) :: p3_SI
        real(r_kind), dimension(          9 ) :: p3
        real(r_kind), dimension(          25) :: p5
        real(r_kind), dimension(          49) :: p7
        
        real(r_kind), dimension(9 ) :: p1_on_3
        real(r_kind), dimension(25) :: p1_on_5
        real(r_kind), dimension(25) :: p3_on_5
        real(r_kind), dimension(49) :: p1_on_7
        real(r_kind), dimension(49) :: p3_on_7
        real(r_kind), dimension(49) :: p5_on_7
      
        real(r_kind) :: tau
        real(r_kind) :: sigma_sum
        real(r_kind), parameter :: eps = 1.e-15
        
        integer(i_kind) :: iCOS,iStencil
        integer(i_kind) :: m
        
        do iStencil = 1,nStencil_all
          m = nWenoCells(iStencil)
          if(iStencil<=nStencil1)then
            if(m>0)then
              do iCOS = 1,m
                p3_SI( iStencil,rematch_idx_3_to_3_SI(iStencil,iCOS) ) = polyCoef(iStencil,iCOS)
              enddo
              beta1(iStencil) = WENO_smooth_indicator_3(p3_SI(iStencil,:))
            else
              beta1(iStencil) = abs(Inf)
            endif
          elseif(iStencil==nStencil1+1)then
            a1 = polyCoef(iStencil,1:m)
            p1 = a1
            beta(1) = minval(beta1) / r(2,2)
          elseif(iStencil==nStencil1+2)then
            ! Rematch array for calculating smooth indicator for 3rd order stencil
            a3 = 0
            do iCOS = 1,m
              a3( rematch_idx_3_to_3(iCOS) ) = polyCoef(iStencil,iCOS)
            enddo
            
            ! Rematch polynomial coefficients and calculate 3rd order polynomial
            p1_on_3    = 0
            p1_on_3(1) = p1(1)
            p3 = ( a3 - p1_on_3 * r(1,2) ) / r(2,2)
            
            beta(2) = WENO_smooth_indicator_3(p3)
            beta(1) = min( beta(1), beta(2) )
          elseif(iStencil==nStencil1+3)then
            ! Rematch array for calculating smooth indicator for 5th order stencil
            a5 = 0
            do iCOS = 1,m
              a5( rematch_idx_5_to_5(iCOS) ) = polyCoef(iStencil,iCOS)
            enddo
            
            ! Rematch polynomial coefficients and calculate 5th order polynomial
            p1_on_5    = 0
            p1_on_5(1) = p1(1)
            p3_on_5    = 0
            do iCOS = 1,9
              p3_on_5( rematch_idx_3_to_5(iCOS) ) = p3(iCOS)
            enddo
            p5 = ( a5 - p1_on_5 * r(1,3) - p3_on_5 * r(2,3) ) / r(3,3)
            
            beta(3) = WENO_smooth_indicator_5(p5)
          elseif(iStencil==nStencil1+4)then
            ! Rematch array for calculating smooth indicator for 7th order stencil
            a7 = 0
            do iCOS = 1,m
              a7( rematch_idx_7_to_7(iCOS) ) = polyCoef(iStencil,iCOS)
            enddo
            
            ! Rematch polynomial coefficients and calculate 7th order polynomial
            p1_on_7    = 0
            p1_on_7(1) = p1(1)
            p3_on_7    = 0
            do iCOS = 1,25
              p3_on_7( rematch_idx_5_to_7(iCOS) ) = p3_on_5(iCOS)
            enddo
            p5_on_7 = 0
            do iCOS = 1,25
              p5_on_7( rematch_idx_5_to_7(iCOS) ) = p5(iCOS)
            enddo
            p7 = ( a7 - p1_on_7 * r(1,4) - p3_on_7 * r(2,4) - p5_on_7 * r(3,4) ) / r(4,4)
            
            beta(4) = WENO_smooth_indicator_7(p7)
          endif
        enddo
        
        !! Engineer solution of low order smooth indicator
        !if(nStencil>=3)then
        !  if( abs( beta(3)-beta(2) ) / ( beta(2) + eps ) < 0.1 ) then
        !    do iStencil = 2,nStencil
        !      beta(iStencil) = beta(1)
        !    enddo
        !  endif
        !endif
        
        !! In smooth region, the Smooth Indicator(SI) of high order stencil is smaller than SI for low order stencil
        !! once the highest order stencil SI is small, we directly use high order stencil
        !if(beta(nStencil)<=beta(nStencil-1))then
        !  do iStencil = 1,nStencil-1
        !    beta(iStencil) = beta(nStencil)
        !  enddo
        !endif
        
        !if(i==23.and.j==23.and.iPatch==1)then
        !  print*,'beta1'
        !  print*,beta1
        !  print*,'beta'
        !  print*,beta
        !  print*,''
        !  do iStencil = 1,nStencil1
        !   print*,iStencil
        !   print*,p3_SI(iStencil,:)
        !  enddo
        !  print*,9
        !  do iCOS = 1,3
        !    write(*,'(3e)'),a3(3*(iCOS-1)+1:3*iCOS)
        !  enddo
        !  print*,10
        !  do iCOS = 1,5
        !    write(*,'(5e)'),a5(5*(iCOS-1)+1:5*iCOS)
        !  enddo
        !  print*,11
        !  do iCOS = 1,7
        !    write(*,'(7e)'),a7(7*(iCOS-1)+1:7*iCOS)
        !  enddo
        !endif
        
        if(maxval(beta)/=0)beta = beta / maxval(beta)
        
        tau = ( sum( abs( beta(nStencil) - beta(1:nStencil-1) ) ) / ( nStencil - 1. ) )**2
        do iStencil = 1,nStencil
          alpha(iStencil) = r(iStencil,nStencil) * ( 1. + tau / ( beta(iStencil) + eps ) )
        enddo
        
        !do iStencil = 1,nStencil
        !  alpha(iStencil) = r(iStencil,nStencil) / ( beta(iStencil) + eps )**2
        !enddo
        
        w = alpha / sum(alpha)
        
        !if(i==23.and.j==23.and.iPatch==1)then
        !  print*,beta1
        !  print*,beta
        !  print*,alpha
        !  print*,w
        !  print*,r(:,nStencil)
        !  print*,''
        !endif
        
        if(nStencil==1)p = p1
        if(nStencil==2)p = w(1) * p1_on_3 + w(2) * p3
        if(nStencil==3)p = w(1) * p1_on_5 + w(2) * p3_on_5 + w(3) * p5
        if(nStencil==4)p = w(1) * p1_on_7 + w(2) * p3_on_7 + w(3) * p5_on_7 + w(4) * p7
        
      end subroutine WENO2D
      
      ! Reconstruct on the left most boundary
      function left_side_recon3(q)
        real(r_kind) :: left_side_recon3
        real(r_kind),dimension(3),intent(in) :: q
        
        real(r_kind) Qip(2)
        real(r_kind) Qi
        
        Qi     = Q(1)
        Qip(1) = Q(2)
        Qip(2) = Q(3)
        
        left_side_recon3 = Qip(2)/3. - 7./6. * Qip(1) + 11./6. * Qi
        !left_side_recon3 = 1.5 * Qi - 0.5 * Qip(1)
        !left_side_recon3 = 3./8. * Qip(2) - 5./4. * Qip(1) + 15./8. * Qi
      
      end function left_side_recon3
      
      ! Reconstruct on the right most boundary
      function right_side_recon3(q)
        real(r_kind) :: right_side_recon3
        real(r_kind),dimension(3),intent(in) :: q
        
        real(r_kind) Qim(2)
        real(r_kind) Qi
        
        Qim(2) = Q(1)
        Qim(1) = Q(2)
        Qi     = Q(3)
        
        right_side_recon3 = Qim(2)/3. - 7./6. * Qim(1) + 11./6. * Qi
        !right_side_recon3 = 1.5 * Qi - 0.5 * Qim(1)
        !right_side_recon3 = 3./8. * Qim(2) - 5./4. * Qim(1) + 15./8. * Qi
      
      end function right_side_recon3
      
      function WENO_smooth_indicator_2(a_in)
        real(r_kind) :: WENO_smooth_indicator_2
        real(r_kind) :: a_in(:)
        real(r16) :: a(lbound(a_in,1):ubound(a_in,1))
        
        a = a_in
        
        WENO_smooth_indicator_2 = a(2)**2 + a(3)**2
        
      end function WENO_smooth_indicator_2
      
      function WENO_smooth_indicator_3(a) !(a_in)
        real(r_kind) :: WENO_smooth_indicator_3
        real(r_kind) :: a(:)
        !real(r_kind) :: a_in(:)
        !real(r16) :: a(lbound(a_in,1):ubound(a_in,1))
        
        !a = a_in
        
        WENO_smooth_indicator_3 = ( 720  * a(2) * a(2) + 3120 * a(3) * a(3) + 720  * a(4) * a(4) + 840   * a(5) * a(5)  &
                                  + 120  * a(4) * a(6) + 3389 * a(6) * a(6) + 3120 * a(7) * a(7) + 120   * a(2) * a(8)  &
                                  + 3389 * a(8) * a(8) + 520  * a(3) * a(9) + 520  * a(7) * a(9) + 13598 * a(9) * a(9) ) / 720
      end function WENO_smooth_indicator_3
      
      function WENO_smooth_indicator_5(a_in)
        real(r_kind) :: WENO_smooth_indicator_5
        real(r_kind) :: a_in(:)
        real(r16) :: a(lbound(a_in,1):ubound(a_in,1))
        
        a = a_in
        
        WENO_smooth_indicator_5 = a(2)**2. + (13.*a(3)**2)/3. + (3129.*a(4)**2)/80. + (87617.*a(5)**2)/140. + a(6)**2. + (7.*a(7)**2)/6. + (1./6)*a(6)*a(8) + (3389.*a(8)**2)/720. + (17./30)*a(7)*a(9) + (47459.*a(9)**2)/1120. + (1./40)*a(6)*a(10) +    &
                                  (5101.*a(8)*a(10))/1120. + (54673043.*a(10)**2)/80640. + (13.*a(11)**2)/3. + (1./24)*a(4)*a(12) + (3389.*a(12)**2)/720. + (7./20)*a(5)*a(13) + (13./18)*a(11)*a(13) + (6799.*a(13)**2)/360. + (1043./160)*a(4)*a(14) +   &
                                  (73./32)*a(12)*a(14) + (22846129.*a(14)**2)/134400. + (87617./840)*a(5)*a(15) + (13./120)*a(11)*a(15) + (306967.*a(13)*a(15))/16800. + (469977913.*a(15)**2)/172800. + (1./2)*a(6)*a(16) + (1./24)*a(8)*a(16) +          &
                                  (1./160)*a(10)*a(16) + (3129.*a(16)**2)/80. + (17./30)*a(7)*a(17) + (11./80)*a(9)*a(17) + (47459.*a(17)**2)/1120. + (1./24)*a(6)*a(18) + (73./32)*a(8)*a(18) + (24721.*a(10)*a(18))/22400. +                             &
                                  (1043./160)*a(16)*a(18) + (22846129.*a(18)**2)/134400. + (11./80)*a(7)*a(19) + (114997.*a(9)*a(19))/5600. + (114997.*a(17)*a(19))/5600. + (19583517.*a(19)**2)/12800. + (1./160)*a(6)*a(20) +                            &
                                  (24721.*a(8)*a(20))/22400. + (37850569.*a(10)*a(20))/115200. + (3129.*a(16)*a(20))/3200. + (2105043.*a(18)*a(20))/12800. + (368483712607.*a(20)**2)/15052800. + (21./5)*a(11)*a(21) + (7./20)*a(13)*a(21) +              &
                                  (21./400)*a(15)*a(21) + (87617.*a(21)**2)/140. + (1./160)*a(4)*a(22) + (5101.*a(12)*a(22))/1120. + (24721.*a(14)*a(22))/22400. + (54673043.*a(22)**2)/80640. + (21./400)*a(5)*a(23) + (7./20)*a(11)*a(23) +              &
                                  (306967.*a(13)*a(23))/16800. + (7071./800)*a(15)*a(23) + (87617./840)*a(21)*a(23) + (469977913.*a(23)**2)/172800. + (3129.*a(4)*a(24))/3200. + (24721.*a(12)*a(24))/22400. + (2105043.*a(14)*a(24))/12800. +             &
                                  (37850569.*a(22)*a(24))/115200. + (368483712607.*a(24)**2)/15052800. + (1./480)*a(2)*(240.*a(4) + 80.*a(12) + 20.*a(14) + 12.*a(22) + 3.*a(24)) + (87617.*a(5)*a(25))/5600. + (21./400)*a(11)*a(25) +                    &
                                  (7071./800)*a(13)*a(25) + (2475532433.*a(15)*a(25))/940800. + (87617.*a(21)*a(25))/5600. + (2475532433.*a(23)*a(25))/940800. + (2210903809027.*a(25)**2)/5644800. +                                                      &
                                  (a(3)*(15120.*a(5) + 2600.*a(13) + 1260.*a(15) + 390.*a(23) + 189.*a(25)))/3600
                              
      end function WENO_smooth_indicator_5
      
      function WENO_smooth_indicator_7(a_in)
        real(r_kind) :: WENO_smooth_indicator_7
        real(r_kind) :: a_in(:)
        real(r16) :: a(lbound(a_in,1):ubound(a_in,1))
        
        a = a_in
        
        WENO_smooth_indicator_7 = a(2)**2. + (13.*a(3)**2)/3. + (3129.*a(4)**2)/80. + (87617.*a(5)**2)/140. + (14127./224)*a(4)*a(6) + (252337135.*a(6)**2)/16128. + (508579./336)*a(5)*a(7) + (11102834003.*a(7)**2)/19712. +  &
                                  a(8)**2. + (7.*a(9)**2)/6. + (1./6)*a(8)*a(10) + (3389.*a(10)**2)/720. + (17./30)*a(9)*a(11) + (47459.*a(11)**2)/1120. + (1./40)*a(8)*a(12) + (5101.*a(10)*a(12))/1120. +                     &
                                  (54673043.*a(12)**2)/80640. + (47./336)*a(9)*a(13) + (34435./504)*a(11)*a(13) + (18042105247.*a(13)**2)/1064448. + (1./224)*a(8)*a(14) + (13579.*a(10)*a(14))/8064. +                         &
                                  (581814439.*a(12)*a(14))/354816. + (7505515786259.*a(14)**2)/12300288. + (13.*a(15)**2)/3. + (1./24)*a(4)*a(16) + (1./96)*a(6)*a(16) + (3389.*a(16)**2)/720. + (7./20)*a(5)*a(17) +           &
                                  (29./224)*a(7)*a(17) + (13./18)*a(15)*a(17) + (6799.*a(17)**2)/360. + (1043./160)*a(4)*a(18) + (4709./896)*a(6)*a(18) + (73./32)*a(16)*a(18) + (22846129.*a(18)**2)/134400. +                 &
                                  (87617./840)*a(5)*a(19) + (508579.*a(7)*a(19))/4032. + (13./120)*a(15)*a(19) + (306967.*a(17)*a(19))/16800. + (469977913.*a(19)**2)/172800. + (4709./896)*a(4)*a(20) +                        &
                                  (252337135.*a(6)*a(20))/96768. + (7561.*a(16)*a(20))/13440. + (18944567.*a(18)*a(20))/69120. + (82716113321.*a(20)**2)/1216512. + (508579.*a(5)*a(21))/4032. +                                &
                                  (11102834003.*a(7)*a(21))/118272. + (13./672)*a(15)*a(21) + (29183.*a(17)*a(21))/4320. + (1667122157.*a(19)*a(21))/253440. + (827161133251.*a(21)**2)/337920. +                               &
                                  (1./2)*a(8)*a(22) + (1./24)*a(10)*a(22) + (1./160)*a(12)*a(22) + (1./896)*a(14)*a(22) + (3129.*a(22)**2)/80. + (17./30)*a(9)*a(23) + (11./80)*a(11)*a(23) +                                   &
                                  (19./560)*a(13)*a(23) + (47459.*a(23)**2)/1120. + (1./24)*a(8)*a(24) + (73./32)*a(10)*a(24) + (24721.*a(12)*a(24))/22400. + (9401.*a(14)*a(24))/23040. + (1043./160)*a(22)*a(24) +            &
                                  (22846129.*a(24)**2)/134400. + (11./80)*a(9)*a(25) + (114997.*a(11)*a(25))/5600. + (190717.*a(13)*a(25))/11520. + (114997.*a(23)*a(25))/5600. + (19583517.*a(25)**2)/12800. +                 &
                                  (1./160)*a(8)*a(26) + (24721.*a(10)*a(26))/22400. + (37850569.*a(12)*a(26))/115200. + (44754957.*a(14)*a(26))/112640. + (3129.*a(22)*a(26))/3200. +                                           &
                                  (2105043.*a(24)*a(26))/12800. + (368483712607.*a(26)**2)/15052800. + (19./560)*a(9)*a(27) + (190717.*a(11)*a(27))/11520. + (138785425.*a(13)*a(27))/16896. +                                  &
                                  (45373.*a(23)*a(27))/8960. + (232084199.*a(25)*a(27))/94080. + (121599625855067.*a(27)**2)/198696960. + (1./896)*a(8)*a(28) + (9401.*a(10)*a(28))/23040. +                                    &
                                  (44754957.*a(12)*a(28))/112640. + (1732042104523.*a(14)*a(28))/5857280. + (447.*a(22)*a(28))/2560. + (18304651.*a(24)*a(28))/301056. + (3921295089347.*a(26)*a(28))/66232320. +               &
                                  (50585444357410783.*a(28)**2)/2296053760. + (21./5)*a(15)*a(29) + (7./20)*a(17)*a(29) + (21./400)*a(19)*a(29) + (3./320)*a(21)*a(29) + (87617.*a(29)**2)/140. +                               &
                                  (1./160)*a(4)*a(30) + (1./640)*a(6)*a(30) + (5101.*a(16)*a(30))/1120. + (24721.*a(18)*a(30))/22400. + (4877.*a(20)*a(30))/17920. + (54673043.*a(30)**2)/80640. +                              &
                                  (21./400)*a(5)*a(31) + (87.*a(7)*a(31))/4480. + (7./20)*a(15)*a(31) + (306967.*a(17)*a(31))/16800. + (7071./800)*a(19)*a(31) + (245947.*a(21)*a(31))/75264. +                                 &
                                  (87617./840)*a(29)*a(31) + (469977913.*a(31)**2)/172800. + (3129.*a(4)*a(32))/3200. + (14127.*a(6)*a(32))/17920. + (24721.*a(16)*a(32))/22400. + (2105043.*a(18)*a(32))/12800. +              &
                                  (199574873.*a(20)*a(32))/1505280. + (37850569.*a(30)*a(32))/115200. + (368483712607.*a(32)**2)/15052800. + (87617.*a(5)*a(33))/5600. + (508579.*a(7)*a(33))/26880. +                          &
                                  (21./400)*a(15)*a(33) + (7071./800)*a(17)*a(33) + (2475532433.*a(19)*a(33))/940800. + (6585971677.*a(21)*a(33))/2069760. + (87617.*a(29)*a(33))/5600. +                                       &
                                  (2475532433.*a(31)*a(33))/940800. + (2210903809027.*a(33)**2)/5644800. + (14127.*a(4)*a(34))/17920. + (50467427.*a(6)*a(34))/129024. + (4877.*a(16)*a(34))/17920. +                           &
                                  (199574873.*a(18)*a(34))/1505280. + (13070811329953.*a(20)*a(34))/198696960. + (365888837.*a(30)*a(34))/4515840. + (712963154333.*a(32)*a(34))/18063360. +                                    &
                                  (140082866134879519.*a(34)**2)/14306181120. + (508579.*a(5)*a(35))/26880. + (11102834003.*a(7)*a(35))/788480. + (3./320)*a(15)*a(35) + (245947.*a(17)*a(35))/75264. +                         &
                                  (6585971677.*a(19)*a(35))/2069760. + (679682189184289.*a(21)*a(35))/287006720. + (87617.*a(29)*a(35))/31360. + (878623885.*a(31)*a(35))/903168. +                                             &
                                  (282333442237789.*a(33)*a(35))/298045440. + (7284309039256637977.*a(35)**2)/20664483840. + (1./8)*a(8)*a(36) + (1./96)*a(10)*a(36) + (1./640)*a(12)*a(36) +                                   &
                                  (a(14)*a(36))/3584. + (14127./224)*a(22)*a(36) + (4709./896)*a(24)*a(36) + (14127.*a(26)*a(36))/17920. + (14127.*a(28)*a(36))/100352. + (252337135.*a(36)**2)/16128. +                        &
                                  (47./336)*a(9)*a(37) + (19./560)*a(11)*a(37) + (15.*a(13)*a(37))/1792. + (34435./504)*a(23)*a(37) + (190717.*a(25)*a(37))/11520. + (921799.*a(27)*a(37))/225792. +                            &
                                  (18042105247.*a(37)**2)/1064448. + (1./96)*a(8)*a(38) + (7561.*a(10)*a(38))/13440. + (4877.*a(12)*a(38))/17920. + (90877.*a(14)*a(38))/903168. + (4709./896)*a(22)*a(38) +                    &
                                  (18944567.*a(24)*a(38))/69120. + (199574873.*a(26)*a(38))/1505280. + (59028127.*a(28)*a(38))/1204224. + (252337135.*a(36)*a(38))/96768. + (82716113321.*a(38)**2)/1216512. +                  &
                                  (19./560)*a(9)*a(39) + (45373.*a(11)*a(39))/8960. + (921799.*a(13)*a(39))/225792. + (190717.*a(23)*a(39))/11520. + (232084199.*a(25)*a(39))/94080. +                                          &
                                  (1197465737.*a(27)*a(39))/602112. + (138785425.*a(37)*a(39))/16896. + (121599625855067.*a(39)**2)/198696960. + (1./640)*a(8)*a(40) + (4877.*a(10)*a(40))/17920. +                             &
                                  (365888837.*a(12)*a(40))/4515840. + (1297893755.*a(14)*a(40))/13246464. + (14127.*a(22)*a(40))/17920. + (199574873.*a(24)*a(40))/1505280. +                                                   &
                                  (712963154333.*a(26)*a(40))/18063360. + (22761431403211.*a(28)*a(40))/476872704. + (50467427.*a(36)*a(40))/129024. + (13070811329953.*a(38)*a(40))/198696960. +                               &
                                  (140082866134879519.*a(40)**2)/14306181120. + (15.*a(9)*a(41))/1792. + (921799.*a(11)*a(41))/225792. + (40247773253.*a(13)*a(41))/19869696. + (921799.*a(23)*a(41))/225792. +                 &
                                  (1197465737.*a(25)*a(41))/602112. + (176458381700353.*a(27)*a(41))/178827264. + (40247773253.*a(37)*a(41))/19869696. + (176458381700353.*a(39)*a(41))/178827264. +                            &
                                  (1400828669297420455.*a(41)**2)/5722472448. + (a(8)*a(42))/3584. + (90877.*a(10)*a(42))/903168. + (1297893755.*a(12)*a(42))/13246464. +                                                       &
                                  (16743073677063.*a(14)*a(42))/229605376. + (14127.*a(22)*a(42))/100352. + (59028127.*a(24)*a(42))/1204224. + (22761431403211.*a(26)*a(42))/476872704. +                                       &
                                  (293626747159200335.*a(28)*a(42))/8265793536. + (252337135.*a(36)*a(42))/3612672. + (34793506056791.*a(38)*a(42))/1430618112. + (45173351014373731.*a(40)*a(42))/1907490816. +                &
                                  (6410191990918725761813.*a(42)**2)/727389831168. + (87./56)*a(15)*a(43) + (29./224)*a(17)*a(43) + (87.*a(19)*a(43))/4480. + (87.*a(21)*a(43))/25088. +                                        &
                                  (508579./336)*a(29)*a(43) + (508579.*a(31)*a(43))/4032. + (508579.*a(33)*a(43))/26880. + (508579.*a(35)*a(43))/150528. + (11102834003.*a(43)**2)/19712. + (1./896)*a(4)*a(44) +               &
                                  (a(6)*a(44))/3584. + (13579.*a(16)*a(44))/8064. + (9401.*a(18)*a(44))/23040. + (90877.*a(20)*a(44))/903168. + (581814439.*a(30)*a(44))/354816. + (44754957.*a(32)*a(44))/112640. +            &
                                  (1297893755.*a(34)*a(44))/13246464. + (7505515786259.*a(44)**2)/12300288. + (3./320)*a(5)*a(45) + (87.*a(7)*a(45))/25088. + (29./224)*a(15)*a(45) + (29183.*a(17)*a(45))/4320. +              &
                                  (245947.*a(19)*a(45))/75264. + (181859.*a(21)*a(45))/150528. + (508579.*a(29)*a(45))/4032. + (1667122157.*a(31)*a(45))/253440. + (6585971677.*a(33)*a(45))/2069760. +                         &
                                  (140250847273.*a(35)*a(45))/119218176. + (11102834003.*a(43)*a(45))/118272. + (827161133251.*a(45)**2)/337920. + (447.*a(4)*a(46))/2560. + (14127.*a(6)*a(46))/100352. +                      &
                                  (9401.*a(16)*a(46))/23040. + (18304651.*a(18)*a(46))/301056. + (59028127.*a(20)*a(46))/1204224. + (44754957.*a(30)*a(46))/112640. + (3921295089347.*a(32)*a(46))/66232320. +                  &
                                  (22761431403211.*a(34)*a(46))/476872704. + (1732042104523.*a(44)*a(46))/5857280. + (50585444357410783.*a(46)**2)/2296053760. + (87617.*a(5)*a(47))/31360. +                                   &
                                  (508579.*a(7)*a(47))/150528. + (87.*a(15)*a(47))/4480. + (245947.*a(17)*a(47))/75264. + (878623885.*a(19)*a(47))/903168. + (140250847273.*a(21)*a(47))/119218176. +                           &
                                  (508579.*a(29)*a(47))/26880. + (6585971677.*a(31)*a(47))/2069760. + (282333442237789.*a(33)*a(47))/298045440. + (45522868146575.*a(35)*a(47))/39739392. +                                     &
                                  (11102834003.*a(43)*a(47))/788480. + (679682189184289.*a(45)*a(47))/287006720. + (7284309039256637977.*a(47)**2)/20664483840. + (14127.*a(4)*a(48))/100352. +                                 &
                                  (252337135.*a(6)*a(48))/3612672. + (90877.*a(16)*a(48))/903168. + (59028127.*a(18)*a(48))/1204224. + (34793506056791.*a(20)*a(48))/1430618112. +                                              &
                                  (1297893755.*a(30)*a(48))/13246464. + (22761431403211.*a(32)*a(48))/476872704. + (45173351014373731.*a(34)*a(48))/1907490816. + (16743073677063.*a(44)*a(48))/229605376. +                    &
                                  (293626747159200335.*a(46)*a(48))/8265793536. + (6410191990918725761813.*a(48)**2)/727389831168. +                                                                                            &
                                  (a(2)*(26880.*a(4) + 6720.*a(6) + 8960.*a(16) + 2240.*a(18) + 560.*a(20) + 1344.*a(30) + 336.*a(32) + 84.*a(34) + 240.*a(44) + 60.*a(46) + 15.*a(48)))/53760. +                               &
                                  a(3)*((21.*a(5))/5. + (87.*a(7))/56. + (13.*a(17))/18. + (7.*a(19))/20. + (29.*a(21))/224. + (13.*a(31))/120. + (21.*a(33))/400. + (87.*a(35))/4480. + (13.*a(45))/672. + (3.*a(47))/320. +   &
                                    (87.*a(49))/25088) + (508579.*a(5)*a(49))/150528. + (11102834003.*a(7)*a(49))/4415488. + (87.*a(15)*a(49))/25088. + (181859.*a(17)*a(49))/150528. +                                         &
                                  (140250847273.*a(19)*a(49))/119218176. + (452315578754789.*a(21)*a(49))/516612096. + (508579.*a(29)*a(49))/150528. + (140250847273.*a(31)*a(49))/119218176. +                                 &
                                  (45522868146575.*a(33)*a(49))/39739392. + (25839156781083324157.*a(35)*a(49))/30307909632. + (11102834003.*a(43)*a(49))/4415488. + (452315578754789.*a(45)*a(49))/516612096. +                &
                                  (25839156781083324157.*a(47)*a(49))/30307909632. + (12820383982264910635167.*a(49)**2)/40410546176
      end function WENO_smooth_indicator_7
      
    end module reconstruction_mod
    