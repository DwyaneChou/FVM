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
      
      integer(i_kind),dimension(:,:,:), allocatable :: nRecCells ! number of cells for reconstruction
      integer(i_kind),dimension(:,:,:), allocatable :: nGstRecCells ! number of cells for ghost point reconstruction
      integer(i_kind),dimension(:,:,:), allocatable :: nRecTerms
      
      integer(i_kind),dimension(:,:,:), allocatable :: locPolyDegree ! degree of local reconstruction polynomial
      
      real   (r_kind), dimension(:,:,:,:,:), allocatable :: polyCoordCoef
      
      real   (r_kind) :: recCoef
      real   (r_kind) :: recdx
      real   (r_kind) :: recdy
      real   (r_kind) :: recdV
      
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
              
              !print*,k
              !print*,quad_wts_tri_2d(  k)
              !print*,triQuad_pos(k,:)
              !print*,''
            enddo
          enddo
          !stop
        endif
        
        allocate(nRecCells   (ids:ide,jds:jde,ifs:ife))
        allocate(nGstRecCells(ids:ide,jds:jde,ifs:ife))
        allocate(nRecTerms   (ids:ide,jds:jde,ifs:ife))
        
        allocate(locPolyDegree(ids:ide,jds:jde,ifs:ife))
        
        allocate(polyCoordCoef(maxRecCells,maxRecTerms,ids:ide,jds:jde,ifs:ife))
        
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
    end module reconstruction_mod
    