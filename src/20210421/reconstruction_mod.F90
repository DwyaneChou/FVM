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
      
      ! For WENO 2D
      integer(i_kind), dimension(:,:,:,:), allocatable :: iCenWENO           ! center cell index on reconstruction stencil for WENO2D
      real   (r_kind), dimension(:,:    ), allocatable :: r                  ! optimal coefficients for WENO 2D
      integer(i_kind), dimension(:,:,:,:), allocatable :: rematch_idx_3_to_3
      integer(i_kind), dimension(:,:,:,:), allocatable :: rematch_idx_5_to_5
      integer(i_kind), dimension(:,:,:,:), allocatable :: rematch_idx_7_to_7
      integer(i_kind), dimension(:      ), allocatable :: rematch_idx_2_to_3
      integer(i_kind), dimension(:      ), allocatable :: rematch_idx_3_to_5
      integer(i_kind), dimension(:      ), allocatable :: rematch_idx_5_to_7
      
    contains
      subroutine init_reconstruction
        integer(i_kind) :: i,j,k
        
        integer(i_kind) :: n
        
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
        
        ! origin WENO
        alpha = weno_coef / ( eps + beta )**2
        omega = alpha / sum(alpha)
        ! origin WENO
        
        Qrec = dot_product( stencil, omega )
        
      end subroutine WENO3
      
      subroutine WENO2D(polyCoef,nCellsOnStencil,rematch_idx_3_to_3,rematch_idx_5_to_5,rematch_idx_7_to_7,p,i,j,iPatch)
        real   (r_kind),dimension(:,:),intent(in ) :: polyCoef
        integer(i_kind),dimension(:  ),intent(in ) :: nCellsOnStencil    ! nWENOCells
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
        
        real(r_kind), dimension(            1 ) :: p1
        real(r_kind), dimension(0:nStencil1,3 ) :: p2
        real(r_kind), dimension(            9 ) :: p3
        real(r_kind), dimension(            25) :: p5
        real(r_kind), dimension(            49) :: p7
        
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
        real(r_kind), parameter :: eps = 1.e-15
        
        integer(i_kind) :: iCOS,iStencil
        integer(i_kind) :: m
        
        do iStencil = 1,nStencil_all
          m = nCellsOnStencil(iStencil)
          if(iStencil<=nStencil1)then
            p2(iStencil,:) = polyCoef(iStencil,1:m)
            beta1(iStencil) = WENO_smooth_indicator_2(polyCoef(iStencil,1:m))
          elseif(iStencil==nStencil1+1)then
            tau = ( ( abs(beta1(1)-beta1(2)) + abs(beta1(1)-beta1(3)) &
                    + abs(beta1(1)-beta1(4)) + abs(beta1(2)-beta1(3)) &
                    + abs(beta1(2)-beta1(4)) + abs(beta1(3)-beta1(4)) ) / 6. )**2
            do iCOS = 1,nStencil1
              sigma(iCOS) = ( 1. + tau / ( beta1(iCOS) + eps ) ) / nStencil1
            enddo
            sigma_sum = sum(sigma)
            sigma = sigma / sigma_sum
            
            a1 = polyCoef(iStencil,1:m)
            p1 = a1
            
            do iCOS = 1,3
              p2(0,iCOS) = dot_product( sigma, p2(1:nStencil1,iCOS) )
            enddo
            beta(1) = WENO_smooth_indicator_2(p2(0,:))
          elseif(iStencil==nStencil1+2)then
            ! Rematch array for calculating smooth indicator for 3rd order stencil
            a3 = 0
            do iCOS = 1,m
              a3( rematch_idx_3_to_3(iCOS) ) = polyCoef(iStencil,iCOS)
            enddo
            
            ! Rematch polynomial coefficients and calculate 3rd order polynomial
            p1_on_3    = 0
            p1_on_3(1) = p1(1)
            p2_on_3    = 0
            do iCOS = 1,3
              p2_on_3( rematch_idx_2_to_3(iCOS) ) = p2(0,iCOS)
            enddo
            p3 = ( a3 - p1_on_3 * r(1,2) ) / r(2,2)
            
            beta(2) = WENO_smooth_indicator_3(p3)
          elseif(iStencil==nStencil1+3)then
            ! Rematch array for calculating smooth indicator for 5th order stencil
            a5 = 0
            do iCOS = 1,m
              a5( rematch_idx_5_to_5(iCOS) ) =  polyCoef(iStencil,iCOS)
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
            
            beta(3) = WENO_smooth_indicator_5(p5)
          elseif(iStencil==nStencil1+4)then
            ! Rematch array for calculating smooth indicator for 7th order stencil
            a7 = 0
            do iCOS = 1,m
              a7( rematch_idx_7_to_7(iCOS) ) =  polyCoef(iStencil,iCOS)
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
            
            beta(4) = WENO_smooth_indicator_7(p7)
          endif
        enddo
        
        ! Engineer solution of low order smooth indicator
        if(nStencil>=3)then
          if( abs( beta(3)-beta(2) ) / ( beta(2) + eps ) < 0.1 ) then
            do iStencil = 2,nStencil
              beta(iStencil) = beta(1)
            enddo
          endif
        endif
        
        !if(i==66.and.j==45.and.iPatch==3)then
        !  print*,a1
        !  print*,''
        !  print*,a3
        !  print*,''
        !  print*,a5
        !  print*,''
        !  print*,a7
        !  print*,''
        !  print*,beta
        !endif
        
        !tau = ( sum( abs( beta(nStencil) - beta(1:nStencil-1) ) ) / ( nStencil - 1. ) )**2
        !
        !do iStencil = 1,nStencil
        !  alpha(iStencil) = r(iStencil,nStencil) * ( 1. + tau / ( beta(iStencil) + eps ) )
        !enddo
        
        do iStencil = 1,nStencil
          alpha(iStencil) = r(iStencil,nStencil) / ( beta(iStencil) + eps )**2
        enddo
        
        w = alpha / sum(alpha)
        
        !if(i==23.and.j==23.and.iPatch==1)then
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
      
      function WENO_smooth_indicator_3(a_in)
        real(r_kind) :: WENO_smooth_indicator_3
        real(r_kind) :: a_in(:)
        real(r16) :: a(lbound(a_in,1):ubound(a_in,1))
        
        a = a_in
        
        WENO_smooth_indicator_3 = (1./720)*(720*a(2)**2 + 3120*a(3)**2 + 720*a(4)**2 + 840*a(5)**2 + &
                                  120*a(4)*a(6) + 3389*a(6)**2 + 3120*a(7)**2 + 120*a(2)*a(8) +      &
                                  3389*a(8)**2 + 520*a(3)*a(9) + 520*a(7)*a(9) + 13598*a(9)**2)
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
        
        WENO_smooth_indicator_7 = a(1)**2. + (13.*a(2)**2)/12. + (1043.*a(3)**2)/240. + (87617.*a(4)**2)/2240. + (4709.*a(3)*a(5))/1120. + (50467427.*a(5)**2)/80640. + (508579.*a(4)*a(6))/8064. + (11102834003.*a(6)**2)/709632. +               &
                                  (12535.*a(3)*a(7))/8064. + (179019829.*a(5)*a(7))/118272. + (2309389472701.*a(7)**2)/4100096. + (13.*a(8)**2)/12. + (169.*a(9)**2)/144. + (13./72)*a(8)*a(10) + (2039.*a(10)**2)/2880. +                         &
                                  (91./160)*a(9)*a(11) + (90701.*a(11)**2)/26880. + (13./480)*a(8)*a(12) + (7457.*a(10)*a(12))/13440. + (50502407.*a(12)**2)/967680. + (377.*a(9)*a(13))/2688. + (515143.*a(11)*a(13))/96768. +                    &
                                  (11102927159.*a(13)**2)/8515584. + (13.*a(8)*a(14))/2688. + (17803.*a(10)*a(14))/96768. + (179049649.*a(12)*a(14))/1419264. + (177645356221.*a(14)**2)/3784704. + (1./72)*a(3)*a(15) +                           &
                                  (1./480)*a(5)*a(15) + (a(7)*a(15))/2688. + (83.*a(15)**2)/240. + (7./160)*a(4)*a(16) + (29.*a(6)*a(16))/2688. + (1079.*a(16)**2)/2880. + (1043.*a(3)*a(17))/1440. + (4709.*a(5)*a(17))/13440. +                  &
                                  (12535.*a(7)*a(17))/96768. + (83.*a(15)*a(17))/1440. + (9769.*a(17)**2)/57600. + (87617.*a(4)*a(18))/13440. + (508579.*a(6)*a(18))/96768. + (581.*a(16)*a(18))/3200. + (283411.*a(18)**2)/537600. +              &
                                  (4709.*a(3)*a(19))/13440. + (50467427.*a(5)*a(19))/483840. + (179019829.*a(7)*a(19))/1419264. + (83.*a(15)*a(19))/9600. + (32447.*a(17)*a(19))/268800. + (151635481.*a(19)**2)/19353600. +                       &
                                  (508579.*a(4)*a(20))/96768. + (11102834003.*a(6)*a(20))/4257792. + (2407.*a(16)*a(20))/53760. + (1569497.*a(18)*a(20))/1935360. + (33309123049.*a(20)**2)/170311680. + (12535.*a(3)*a(21))/96768. +              &
                                  (179019829.*a(5)*a(21))/1419264. + (2309389472701.*a(7)*a(21))/24600576. + (83.*a(15)*a(21))/53760. + (14545.*a(17)*a(21))/387072. + (537258287.*a(19)*a(21))/28385280. +                                        &
                                  (6928169472583.*a(21)**2)/984023040. + (21./40)*a(8)*a(22) + (7./160)*a(10)*a(22) + (21.*a(12)*a(22))/3200. + (3.*a(14)*a(22))/2560. + (257.*a(22)**2)/2240. + (91./160)*a(9)*a(23) +                            &
                                  (441.*a(11)*a(23))/3200. + (87.*a(13)*a(23))/2560. + (3341.*a(23)**2)/26880. + (7./160)*a(8)*a(24) + (901.*a(10)*a(24))/3200. + (1327.*a(12)*a(24))/12800. + (203.*a(14)*a(24))/6144. +                          &
                                  (257.*a(22)*a(24))/13440. + (3733.*a(24)**2)/76800. + (441.*a(9)*a(25))/3200. + (13251.*a(11)*a(25))/12800. + (24739.*a(13)*a(25))/30720. + (771.*a(23)*a(25))/12800. + (502849.*a(25)**2)/5017600. +            &
                                  (21.*a(8)*a(26))/3200. + (1327.*a(10)*a(26))/12800. + (2405987.*a(12)*a(26))/153600. + (59689843.*a(14)*a(26))/3153920. + (257.*a(22)*a(26))/89600. + (81253.*a(24)*a(26))/2508800. +                            &
                                  (50614343.*a(26)**2)/36126720. + (87.*a(9)*a(27))/2560. + (24739.*a(11)*a(27))/30720. + (3700996421.*a(13)*a(27))/9461760. + (7453.*a(23)*a(27))/501760. + (2680739.*a(25)*a(27))/18063360. +                    &
                                  (55516126291.*a(27)**2)/1589575680. + (3.*a(8)*a(28))/2560. + (203.*a(10)*a(28))/6144. + (59689843.*a(12)*a(28))/3153920. + (2309389736321.*a(14)*a(28))/164003840. + (257.*a(22)*a(28))/501760. +               &
                                  (173303.*a(24)*a(28))/18063360. + (179145073.*a(26)*a(28))/52985856. + (11546950685117.*a(28)**2)/9184215040. + (1./480)*a(3)*a(29) + (a(5)*a(29))/3200. + (a(7)*a(29))/17920. +                                 &
                                  (229.*a(15)*a(29))/1120. + (229.*a(17)*a(29))/13440. + (229.*a(19)*a(29))/89600. + (229.*a(21)*a(29))/501760. + (583.*a(29)**2)/16128. + (21.*a(4)*a(30))/3200. + (29.*a(6)*a(30))/17920. +                      &
                                  (2977.*a(16)*a(30))/13440. + (687.*a(18)*a(30))/12800. + (6641.*a(20)*a(30))/501760. + (7579.*a(30)**2)/193536. + (1043.*a(3)*a(31))/9600. + (4709.*a(5)*a(31))/89600. + (2507.*a(7)*a(31))/129024. +            &
                                  (229.*a(15)*a(31))/13440. + (3401.*a(17)*a(31))/38400. + (74841.*a(19)*a(31))/2508800. + (161011.*a(21)*a(31))/18063360. + (583.*a(29)*a(31))/96768. + (55109.*a(31)**2)/3870720. +                              &
                                  (87617.*a(4)*a(32))/89600. + (508579.*a(6)*a(32))/645120. + (687.*a(16)*a(32))/12800. + (495653.*a(18)*a(32))/2508800. + (2665423.*a(20)*a(32))/18063360. + (583.*a(30)*a(32))/30720. +                          &
                                  (761351.*a(32)**2)/36126720. + (4709.*a(3)*a(33))/89600. + (50467427.*a(5)*a(33))/3225600. + (179019829.*a(7)*a(33))/9461760. + (229.*a(15)*a(33))/89600. + (74841.*a(17)*a(33))/2508800. +                      &
                                  (50598019.*a(19)*a(33))/18063360. + (179131157.*a(21)*a(33))/52985856. + (583.*a(29)*a(33))/645120. + (164867.*a(31)*a(33))/18063360. + (354951029.*a(33)**2)/1300561920. +                                      &
                                  (508579.*a(4)*a(34))/645120. + (11102834003.*a(6)*a(34))/28385280. + (6641.*a(16)*a(34))/501760. + (2665423.*a(18)*a(34))/18063360. + (55515908927.*a(20)*a(34))/794787840. +                                    &
                                  (16907.*a(30)*a(34))/3612672. + (3875125.*a(32)*a(34))/130056192. + (1009406617.*a(34)**2)/148635648. + (2507.*a(3)*a(35))/129024. + (179019829.*a(5)*a(35))/9461760. +                                          &
                                  (2309389472701.*a(7)*a(35))/164003840. + (229.*a(15)*a(35))/501760. + (161011.*a(17)*a(35))/18063360. + (179131157.*a(19)*a(35))/52985856. + (11546950316049.*a(21)*a(35))/4592107520. +                         &
                                  (583.*a(29)*a(35))/3612672. + (340609.*a(31)*a(35))/130056192. + (16293119.*a(33)*a(35))/24772608. + (209944596119.*a(35)**2)/858783744. + (29./224)*a(8)*a(36) + (29.*a(10)*a(36))/2688. +                      &
                                  (29.*a(12)*a(36))/17920. + (29.*a(14)*a(36))/100352. + (547.*a(22)*a(36))/8064. + (547.*a(24)*a(36))/96768. + (547.*a(26)*a(36))/645120. + (547.*a(28)*a(36))/3612672. + (1109.*a(36)**2)/101376. +              &
                                  (377.*a(9)*a(37))/2688. + (87.*a(11)*a(37))/2560. + (841.*a(13)*a(37))/100352. + (7111.*a(23)*a(37))/96768. + (547.*a(25)*a(37))/30720. + (15863.*a(27)*a(37))/3612672. + (14417.*a(37)**2)/1216512. +           &
                                  (29.*a(8)*a(38))/2688. + (481.*a(10)*a(38))/7680. + (11121.*a(12)*a(38))/501760. + (24827.*a(14)*a(38))/3612672. + (547.*a(22)*a(38))/96768. + (52121.*a(24)*a(38))/1935360. +                                   &
                                  (156623.*a(26)*a(38))/18063360. + (324805.*a(28)*a(38))/130056192. + (1109.*a(36)*a(38))/608256. + (100687.*a(38)**2)/24330240. + (87.*a(9)*a(39))/2560. + (94813.*a(11)*a(39))/501760. +                        &
                                  (523895.*a(13)*a(39))/3612672. + (547.*a(23)*a(39))/30720. + (752099.*a(25)*a(39))/18063360. + (3855433.*a(27)*a(39))/130056192. + (7763.*a(37)*a(39))/1351680. + (1071253.*a(39)**2)/227082240. +               &
                                  (29.*a(8)*a(40))/17920. + (11121.*a(10)*a(40))/501760. + (50549047.*a(12)*a(40))/18063360. + (179089409.*a(14)*a(40))/52985856. + (547.*a(22)*a(40))/645120. + (156623.*a(24)*a(40))/18063360. +                 &
                                  (354846089.*a(26)*a(40))/650280960. + (25601647.*a(28)*a(40))/38928384. + (1109.*a(36)*a(40))/4055040. + (294281.*a(38)*a(40))/113541120. + (457413343.*a(40)**2)/8174960640. +                                  &
                                  (841.*a(9)*a(41))/100352. + (523895.*a(11)*a(41))/3612672. + (11103051367.*a(13)*a(41))/158957568. + (15863.*a(23)*a(41))/3612672. + (3855433.*a(25)*a(41))/130056192. +                                         &
                                  (11103432863.*a(27)*a(41))/817496064. + (32161.*a(37)*a(41))/22708224. + (5178911.*a(39)*a(41))/817496064. + (99934045327.*a(41)**2)/71939653632. + (29.*a(8)*a(42))/100352. +                                   &
                                  (24827.*a(10)*a(42))/3612672. + (179089409.*a(12)*a(42))/52985856. + (2309389841769.*a(14)*a(42))/918421504. + (547.*a(22)*a(42))/3612672. + (324805.*a(24)*a(42))/130056192. +                                  &
                                  (25601647.*a(26)*a(42))/38928384. + (2309390489521.*a(28)*a(42))/4723310592. + (1109.*a(36)*a(42))/22708224. + (595715.*a(38)*a(42))/817496064. + (1613911961.*a(40)*a(42))/11989942272. +                       &
                                  (20784519753409.*a(42)**2)/415651332096. + (a(3)*a(43))/2688. + (a(5)*a(43))/17920. + (a(7)*a(43))/100352. + (439.*a(15)*a(43))/8064. + (439.*a(17)*a(43))/96768. + (439.*a(19)*a(43))/645120. +                 &
                                  (439.*a(21)*a(43))/3612672. + (355.*a(29)*a(43))/16896. + (355.*a(31)*a(43))/202752. + (71.*a(33)*a(43))/270336. + (355.*a(35)*a(43))/7569408. + (1883.*a(43)**2)/585728. + (3.*a(4)*a(44))/2560. +              &
                                  (29.*a(6)*a(44))/100352. + (5707.*a(16)*a(44))/96768. + (439.*a(18)*a(44))/30720. + (12731.*a(20)*a(44))/3612672. + (4615.*a(30)*a(44))/202752. + (497.*a(32)*a(44))/90112. +                                    &
                                  (10295.*a(34)*a(44))/7569408. + (1883.*a(44)**2)/540672. + (149.*a(3)*a(45))/7680. + (4709.*a(5)*a(45))/501760. + (12535.*a(7)*a(45))/3612672. + (439.*a(15)*a(45))/96768. +                                     &
                                  (43157.*a(17)*a(45))/1935360. + (131891.*a(19)*a(45))/18063360. + (277393.*a(21)*a(45))/130056192. + (355.*a(29)*a(45))/202752. + (6469.*a(31)*a(45))/811008. + (18947.*a(33)*a(45))/7569408. +                  &
                                  (192133.*a(35)*a(45))/272498688. + (1883.*a(43)*a(45))/3514368. + (166849.*a(45)**2)/140574720. + (87617.*a(4)*a(46))/501760. + (508579.*a(6)*a(46))/3612672. + (439.*a(16)*a(46))/30720. +                      &
                                  (724343.*a(18)*a(46))/18063360. + (3796357.*a(20)*a(46))/130056192. + (497.*a(30)*a(46))/90112. + (70663.*a(32)*a(46))/7569408. + (1718281.*a(34)*a(46))/272498688. + (39543.*a(44)*a(46))/23429120. +           &
                                  (206413.*a(46)**2)/187432960. + (4709.*a(3)*a(47))/501760. + (50467427.*a(5)*a(47))/18063360. + (179019829.*a(7)*a(47))/52985856. + (439.*a(15)*a(47))/645120. + (131891.*a(17)*a(47))/18063360. +               &
                                  (354531269.*a(19)*a(47))/650280960. + (179173189.*a(21)*a(47))/272498688. + (71.*a(29)*a(47))/270336. + (18947.*a(31)*a(47))/7569408. + (152428361.*a(33)*a(47))/1362493440. +                                   &
                                  (537934207.*a(35)*a(47))/3996647424. + (1883.*a(43)*a(47))/23429120. + (68641.*a(45)*a(47))/93716480. + (7280501.*a(47)**2)/613416960. +                                                                         &
                                  a(2)*((21.*a(4))/40. + (29.*a(6))/224. + (13.*a(16))/72. + (7.*a(18))/160. + (29.*a(20))/2688. + (13.*a(30))/480. + (21.*a(32))/3200. + (29.*a(34))/17920. + (13.*a(44))/2688. + (3.*a(46))/2560. +              &
                                    (29.*a(48))/100352) + (508579.*a(4)*a(48))/3612672. + (11102834003.*a(6)*a(48))/158957568. + (12731.*a(16)*a(48))/3612672. + (3796357.*a(18)*a(48))/130056192. +                                               &
                                  (11103313091.*a(20)*a(48))/817496064. + (10295.*a(30)*a(48))/7569408. + (1718281.*a(32)*a(48))/272498688. + (33311234585.*a(34)*a(48))/11989942272. + (7801.*a(44)*a(48))/18743296. +                            &
                                  (945479.*a(46)*a(48))/674758656. + (122145706369.*a(48)**2)/415651332096. + a(1)*(a(3)/6. + a(5)/40. + a(7)/224. + a(15)/6. + a(17)/72. + a(19)/480. + a(21)/2688. + a(29)/40. + a(31)/480. + a(33)/3200. +      &
                                    a(35)/17920. + a(43)/224. + a(45)/2688. + a(47)/17920. + a(49)/100352) + (12535.*a(3)*a(49))/3612672. + (179019829.*a(5)*a(49))/52985856. + (2309389472701.*a(7)*a(49))/918421504. +                           &
                                  (439.*a(15)*a(49))/3612672. + (277393.*a(17)*a(49))/130056192. + (179173189.*a(19)*a(49))/272498688. + (2309390286157.*a(21)*a(49))/4723310592. + (355.*a(29)*a(49))/7569408. +                                  &
                                  (192133.*a(31)*a(49))/272498688. + (537934207.*a(33)*a(49))/3996647424. + (6928173057815.*a(35)*a(49))/69275222016. + (269.*a(43)*a(49))/18743296. + (137099.*a(45)*a(49))/674758656. +                          &
                                  (1973870039.*a(47)*a(49))/69275222016. + (25403308874543.*a(49)**2)/2401541029888
      end function WENO_smooth_indicator_7
      
    end module reconstruction_mod
    