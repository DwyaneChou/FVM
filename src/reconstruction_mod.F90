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
        
        WENO_smooth_indicator_3 = a(2)**2._r16 + a(2)*(2._r16*a(3) + a(5) + a(6) + (2._r16*a(8))/3._r16 + (2._r16*a(9))/3) +                                        &
                                  (1._r16/90)*(480._r16*a(3)**2._r16 + 90._r16*a(4)**2._r16 + 150._r16*a(5)**2._r16 + 285._r16*a(5)*a(6) + 658._r16*a(6)**2._r16 +  &
                                    90._r16*a(5)*a(7) + 60._r16*a(6)*a(7) + 480._r16*a(7)**2._r16 + 285._r16*a(5)*a(8) + 270._r16*a(6)*a(8) +                       &
                                    480._r16*a(7)*a(8) + 658._r16*a(8)**2._r16 + 270._r16*a(5)*a(9) + 1236._r16*a(6)*a(9) + 320._r16*a(7)*a(9) +                    &
                                    1236._r16*a(8)*a(9) + 2752._r16*a(9)**2._r16 + 30._r16*a(4)*(3._r16*a(5) + 2._r16*a(6) + 6._r16*a(7) + 3._r16*a(8) +            &
                                      2._r16*a(9)) + 10._r16*a(3)*(9._r16*a(5) + 48._r16*a(6) + 6._r16*a(8) + 32._r16*a(9)))
      
      end function WENO_smooth_indicator_3
      
      function WENO_smooth_indicator_5(a_in)
        real(r_kind) :: WENO_smooth_indicator_5
        real(r_kind) :: a_in(:)
        real(r16) :: a(lbound(a_in,1):ubound(a_in,1))
        
        a = a_in
        
        WENO_smooth_indicator_5 = a(2)**2._r16 + (16._r16*a(3)**2)/3._r16 + (249._r16*a(4)**2)/5._r16 + 184._r16*a(4)*a(5) + (27968._r16*a(5)**2)/35._r16 + a(6)**2._r16 + a(4)*a(7) + a(5)*a(7) + a(6)*a(7) + (5._r16*a(7)**2)/3._r16 + (15._r16/2)*a(4)*a(8) + (48._r16/5)*a(5)*a(8) + (2._r16/3)*a(6)*a(8) + (19._r16/6)*a(7)*a(8) +                                  &
                                  (329._r16*a(8)**2)/45._r16 + (249._r16/5)*a(4)*a(9) + 92._r16*a(5)*a(9) + (1._r16/2)*a(6)*a(9) + (46._r16/15)*a(7)*a(9) + (61._r16/3)*a(8)*a(9) + (2329._r16*a(9)**2)/35._r16 + 92._r16*a(4)*a(10) + (27968._r16/35)*a(5)*a(10) + (2._r16/5)*a(6)*a(10) + 3._r16*a(7)*a(10) +                                                          &
                                  (906._r16/35)*a(8)*a(10) + (2947._r16/12)*a(9)*a(10) + (335651._r16*a(10)**2)/315._r16 + 2._r16*a(6)*a(11) + a(7)*a(11) + (2._r16/3)*a(8)*a(11) + (1._r16/2)*a(9)*a(11) + (2._r16/5)*a(10)*a(11) + (16._r16*a(11)**2)/3._r16 + (2._r16/3)*a(4)*a(12) + (2._r16/3)*a(5)*a(12) +                                                         &
                                  a(6)*a(12) + (19._r16/6)*a(7)*a(12) + 3._r16*a(8)*a(12) + (29._r16/10)*a(9)*a(12) + (17._r16/6)*a(10)*a(12) + (16._r16/3)*a(11)*a(12) + (329._r16*a(12)**2)/45._r16 + 5._r16*a(4)*a(13) + (32._r16/5)*a(5)*a(13) + (2._r16/3)*a(6)*a(13) + 3._r16*a(7)*a(13) +                                                                         &
                                  (206._r16/15)*a(8)*a(13) + (229._r16/12)*a(9)*a(13) + (170._r16/7)*a(10)*a(13) + (32._r16/9)*a(11)*a(13) + (206._r16/15)*a(12)*a(13) + (1376._r16*a(13)**2)/45._r16 + (166._r16/5)*a(4)*a(14) + (184._r16/3)*a(5)*a(14) + (1._r16/2)*a(6)*a(14) + (29._r16/10)*a(7)*a(14) +                                                            &
                                  (229._r16/12)*a(8)*a(14) + (1747._r16/14)*a(9)*a(14) + (921._r16/4)*a(10)*a(14) + (8._r16/3)*a(11)*a(14) + (66._r16/5)*a(12)*a(14) + (763._r16/9)*a(13)*a(14) + (145069._r16*a(14)**2)/525._r16 + (184._r16/3)*a(4)*a(15) + (55936._r16/105)*a(5)*a(15) +                                                                              &
                                  (2._r16/5)*a(6)*a(15) + (17._r16/6)*a(7)*a(15) + (170._r16/7)*a(8)*a(15) + (921._r16/4)*a(9)*a(15) + (125870._r16/63)*a(10)*a(15) + (32._r16/15)*a(11)*a(15) + (578._r16/45)*a(12)*a(15) + (56576._r16/525)*a(13)*a(15) + (15292._r16/15)*a(14)*a(15) +                                                                                &
                                  (20894896._r16*a(15)**2)/4725._r16 + 2._r16*a(6)*a(16) + a(7)*a(16) + (2._r16/3)*a(8)*a(16) + (1._r16/2)*a(9)*a(16) + (2._r16/5)*a(10)*a(16) + 15._r16*a(11)*a(16) + (15._r16/2)*a(12)*a(16) + 5._r16*a(13)*a(16) + (15._r16/4)*a(14)*a(16) + 3._r16*a(15)*a(16) +                                                                     &
                                  (249._r16*a(16)**2)/5._r16 + (1._r16/2)*a(4)*a(17) + (1._r16/2)*a(5)*a(17) + a(6)*a(17) + (46._r16/15)*a(7)*a(17) + (29._r16/10)*a(8)*a(17) + (14._r16/5)*a(9)*a(17) + (41._r16/15)*a(10)*a(17) + (15._r16/2)*a(11)*a(17) + (61._r16/3)*a(12)*a(17) +                                                                                  &
                                  (229._r16/12)*a(13)*a(17) + (55._r16/3)*a(14)*a(17) + (107._r16/6)*a(15)*a(17) + (249._r16/5)*a(16)*a(17) + (2329._r16*a(17)**2)/35._r16 + (15._r16/4)*a(4)*a(18) + (24._r16/5)*a(5)*a(18) + (2._r16/3)*a(6)*a(18) + (29._r16/10)*a(7)*a(18) + (66._r16/5)*a(8)*a(18) +                                                                &
                                  (55._r16/3)*a(9)*a(18) + (4082._r16/175)*a(10)*a(18) + 5._r16*a(11)*a(18) + (229._r16/12)*a(12)*a(18) + (763._r16/9)*a(13)*a(18) + (235._r16/2)*a(14)*a(18) + (5227._r16/35)*a(15)*a(18) + (166._r16/5)*a(16)*a(18) + (1747._r16/14)*a(17)*a(18) +                                                                                     &
                                  (145069._r16*a(18)**2)/525._r16 + (249._r16/10)*a(4)*a(19) + 46._r16*a(5)*a(19) + (1._r16/2)*a(6)*a(19) + (14._r16/5)*a(7)*a(19) + (55._r16/3)*a(8)*a(19) + (20966._r16/175)*a(9)*a(19) + (4421._r16/20)*a(10)*a(19) + (15._r16/4)*a(11)*a(19) + (55._r16/3)*a(12)*a(19) +                                                             &
                                  (235._r16/2)*a(13)*a(19) + (26801._r16/35)*a(14)*a(19) + (33901._r16/24)*a(15)*a(19) + (249._r16/10)*a(16)*a(19) + (20966._r16/175)*a(17)*a(19) + (26801._r16/35)*a(18)*a(19) + (436497._r16*a(19)**2)/175._r16 + 46._r16*a(4)*a(20) + (13984._r16/35)*a(5)*a(20) +                                                                    &
                                  (2._r16/5)*a(6)*a(20) + (41._r16/15)*a(7)*a(20) + (4082._r16/175)*a(8)*a(20) + (4421._r16/20)*a(9)*a(20) + (3020894._r16*a(10)*a(20))/1575._r16 + 3._r16*a(11)*a(20) + (107._r16/6)*a(12)*a(20) + (5227._r16/35)*a(13)*a(20) + (33901._r16/24)*a(14)*a(20) +                                                                           &
                                  (428901._r16/35)*a(15)*a(20) + (498._r16/25)*a(16)*a(20) + (4077._r16/35)*a(17)*a(20) + (170298._r16/175)*a(18)*a(20) + (1288271._r16/140)*a(19)*a(20) + (146684527._r16*a(20)**2)/3675._r16 + 2._r16*a(6)*a(21) + a(7)*a(21) + (2._r16/3)*a(8)*a(21) +                                                                                &
                                  (1._r16/2)*a(9)*a(21) + (2._r16/5)*a(10)*a(21) + (96._r16/5)*a(11)*a(21) + (48._r16/5)*a(12)*a(21) + (32._r16/5)*a(13)*a(21) + (24._r16/5)*a(14)*a(21) + (96._r16/25)*a(15)*a(21) + 184._r16*a(16)*a(21) + 92._r16*a(17)*a(21) + (184._r16/3)*a(18)*a(21) +                                                                            &
                                  46._r16*a(19)*a(21) + (184._r16/5)*a(20)*a(21) + (27968._r16*a(21)**2)/35._r16 + (2._r16/5)*a(4)*a(22) + (2._r16/5)*a(5)*a(22) + a(6)*a(22) + 3._r16*a(7)*a(22) + (17._r16/6)*a(8)*a(22) + (41._r16/15)*a(9)*a(22) + (8._r16/3)*a(10)*a(22) + (48._r16/5)*a(11)*a(22) +                                                                &
                                  (906._r16/35)*a(12)*a(22) + (170._r16/7)*a(13)*a(22) + (4082._r16/175)*a(14)*a(22) + (794._r16/35)*a(15)*a(22) + 92._r16*a(16)*a(22) + (2947._r16/12)*a(17)*a(22) + (921._r16/4)*a(18)*a(22) + (4421._r16/20)*a(19)*a(22) + (2579._r16/12)*a(20)*a(22) +                                                                               &
                                  (27968._r16/35)*a(21)*a(22) + (335651._r16*a(22)**2)/315._r16 + 3._r16*a(4)*a(23) + (96._r16/25)*a(5)*a(23) + (2._r16/3)*a(6)*a(23) + (17._r16/6)*a(7)*a(23) + (578._r16/45)*a(8)*a(23) + (107._r16/6)*a(9)*a(23) + (794._r16/35)*a(10)*a(23) + (32._r16/5)*a(11)*a(23) +                                                              &
                                  (170._r16/7)*a(12)*a(23) + (56576._r16/525)*a(13)*a(23) + (5227._r16/35)*a(14)*a(23) + (33216._r16/175)*a(15)*a(23) + (184._r16/3)*a(16)*a(23) + (921._r16/4)*a(17)*a(23) + (15292._r16/15)*a(18)*a(23) + (33901._r16/24)*a(19)*a(23) +                                                                                                &
                                  (62828._r16/35)*a(20)*a(23) + (55936._r16/105)*a(21)*a(23) + (125870._r16/63)*a(22)*a(23) + (20894896._r16*a(23)**2)/4725._r16 + (498._r16/25)*a(4)*a(24) + (184._r16/5)*a(5)*a(24) + (1._r16/2)*a(6)*a(24) + (41._r16/15)*a(7)*a(24) + (107._r16/6)*a(8)*a(24) +                                                                      &
                                  (4077._r16/35)*a(9)*a(24) + (2579._r16/12)*a(10)*a(24) + (24._r16/5)*a(11)*a(24) + (4082._r16/175)*a(12)*a(24) + (5227._r16/35)*a(13)*a(24) + (170298._r16/175)*a(14)*a(24) + (62828._r16/35)*a(15)*a(24) + 46._r16*a(16)*a(24) + (4421._r16/20)*a(17)*a(24) +                                                                         &
                                  (33901._r16/24)*a(18)*a(24) + (1288271._r16/140)*a(19)*a(24) + 16974._r16*a(20)*a(24) + (13984._r16/35)*a(21)*a(24) + (3020894._r16*a(22)*a(24))/1575._r16 + (428901._r16/35)*a(23)*a(24) + (146684527._r16*a(24)**2)/3675._r16 +                                                                                                      &
                                  a(2)*(2._r16*a(3) + 2._r16*a(4) + 2._r16*a(5) + a(7) + a(8) + a(9) + a(10) + (2._r16*a(12))/3._r16 + (2._r16*a(13))/3._r16 + (2._r16*a(14))/3._r16 + (2._r16*a(15))/3._r16 + a(17)/2._r16 + a(18)/2._r16 + a(19)/2._r16 + a(20)/2._r16 + (2._r16*a(22))/5._r16 + (2._r16*a(23))/5._r16 + (2._r16*a(24))/5._r16 + (2._r16*a(25))/5) +   &
                                  (184._r16/5)*a(4)*a(25) + (55936._r16/175)*a(5)*a(25) + (2._r16/5)*a(6)*a(25) + (8._r16/3)*a(7)*a(25) + (794._r16/35)*a(8)*a(25) + (2579._r16/12)*a(9)*a(25) + (83914._r16/45)*a(10)*a(25) + (96._r16/25)*a(11)*a(25) + (794._r16/35)*a(12)*a(25) +                                                                                    &
                                  (33216._r16/175)*a(13)*a(25) + (62828._r16/35)*a(14)*a(25) + (57230368._r16*a(15)*a(25))/3675._r16 + (184._r16/5)*a(16)*a(25) + (2579._r16/12)*a(17)*a(25) + (62828._r16/35)*a(18)*a(25) + 16974._r16*a(19)*a(25) + (46384376._r16/315)*a(20)*a(25) +                                                                                  &
                                  (55936._r16/175)*a(21)*a(25) + (83914._r16/45)*a(22)*a(25) + (57230368._r16*a(23)*a(25))/3675._r16 + (46384376._r16/315)*a(24)*a(25) + (7041838976._r16*a(25)**2)/11025._r16 +                                                                                                                                                         &
                                  a(3)*(15._r16*a(4) + (96._r16*a(5))/5._r16 + a(7) + (16._r16*a(8))/3._r16 + (15._r16*a(9))/2._r16 + (48._r16*a(10))/5._r16 + (2._r16*a(12))/3._r16 + (32._r16*a(13))/9._r16 + 5._r16*a(14) + (32._r16*a(15))/5._r16 + a(17)/2._r16 + (8._r16*a(18))/3._r16 + (15._r16*a(19))/4._r16 + (24._r16*a(20))/5._r16 + (2._r16*a(22))/5._r16 + &
                                    (32._r16*a(23))/15._r16 + 3._r16*a(24) + (96._r16*a(25))/25)
                              
      end function WENO_smooth_indicator_5
      
      function WENO_smooth_indicator_7(a_in)
        real(r_kind) :: WENO_smooth_indicator_7
        real(r_kind) :: a_in(:)
        real(r16) :: a(lbound(a_in,1):ubound(a_in,1))
        
        a = a_in
        
        WENO_smooth_indicator_7 = a(2)**2._r16 + (16._r16*a(3)**2)/3._r16 + (249._r16*a(4)**2)/5._r16 + 184._r16*a(4)*a(5) + (27968._r16*a(5)**2)/35._r16 + (2046._r16/7)*a(4)*a(6) + 3685._r16*a(5)*a(6) + (1258735._r16*a(6)**2)/63._r16 + (849._r16/2)*a(4)*a(7) + (147424._r16/21)*a(5)*a(7) + 110556._r16*a(6)*a(7) + (55384592._r16*a(7)**2)/77._r16 + a(8)**2._r16 +                                            &
                                  a(4)*a(9) + a(5)*a(9) + a(6)*a(9) + a(7)*a(9) + a(8)*a(9) + (5._r16*a(9)**2)/3._r16 + (15._r16/2)*a(4)*a(10) + (48._r16/5)*a(5)*a(10) + (35._r16/3)*a(6)*a(10) + (96._r16/7)*a(7)*a(10) + (2._r16/3)*a(8)*a(10) + (19._r16/6)*a(9)*a(10) + (329._r16*a(10)**2)/45._r16 +                                                                                                             &
                                  (249._r16/5)*a(4)*a(11) + 92._r16*a(5)*a(11) + (1023._r16/7)*a(6)*a(11) + (849._r16/4)*a(7)*a(11) + (1._r16/2)*a(8)*a(11) + (46._r16/15)*a(9)*a(11) + (61._r16/3)*a(10)*a(11) + (2329._r16*a(11)**2)/35._r16 + 92._r16*a(4)*a(12) + (27968._r16/35)*a(5)*a(12) + (3685._r16/2)*a(6)*a(12) +                                                                                          &
                                  (73712._r16/21)*a(7)*a(12) + (2._r16/5)*a(8)*a(12) + 3._r16*a(9)*a(12) + (906._r16/35)*a(10)*a(12) + (2947._r16/12)*a(11)*a(12) + (335651._r16*a(12)**2)/315._r16 + (1023._r16/7)*a(4)*a(13) + (3685._r16/2)*a(5)*a(13) + (1258735._r16/63)*a(6)*a(13) + 55278._r16*a(7)*a(13) +                                                                                                     &
                                  (1._r16/3)*a(8)*a(13) + (62._r16/21)*a(9)*a(13) + (1129._r16/36)*a(10)*a(13) + (24566._r16/63)*a(11)*a(13) + (73703._r16/15)*a(12)*a(13) + (55384529._r16*a(13)**2)/2079._r16 + (849._r16/4)*a(4)*a(14) + (73712._r16/21)*a(5)*a(14) + 55278._r16*a(6)*a(14) + (55384592._r16/77)*a(7)*a(14) +                                                                                       &
                                  (2._r16/7)*a(8)*a(14) + (35._r16/12)*a(9)*a(14) + (2318._r16/63)*a(10)*a(14) + (2831._r16/5)*a(11)*a(14) + (6486782._r16/693)*a(12)*a(14) + (884449._r16/6)*a(13)*a(14) + (2879999015._r16*a(14)**2)/3003._r16 + 2._r16*a(8)*a(15) + a(9)*a(15) + (2._r16/3)*a(10)*a(15) +                                                                                                           &
                                  (1._r16/2)*a(11)*a(15) + (2._r16/5)*a(12)*a(15) + (1._r16/3)*a(13)*a(15) + (2._r16/7)*a(14)*a(15) + (16._r16*a(15)**2)/3._r16 + (2._r16/3)*a(4)*a(16) + (2._r16/3)*a(5)*a(16) + (2._r16/3)*a(6)*a(16) + (2._r16/3)*a(7)*a(16) + a(8)*a(16) + (19._r16/6)*a(9)*a(16) + 3._r16*a(10)*a(16) +                                                                                           &
                                  (29._r16/10)*a(11)*a(16) + (17._r16/6)*a(12)*a(16) + (39._r16/14)*a(13)*a(16) + (11._r16/4)*a(14)*a(16) + (16._r16/3)*a(15)*a(16) + (329._r16*a(16)**2)/45._r16 + 5._r16*a(4)*a(17) + (32._r16/5)*a(5)*a(17) + (70._r16/9)*a(6)*a(17) + (64._r16/7)*a(7)*a(17) + (2._r16/3)*a(8)*a(17) +                                                                                             &
                                  3._r16*a(9)*a(17) + (206._r16/15)*a(10)*a(17) + (229._r16/12)*a(11)*a(17) + (170._r16/7)*a(12)*a(17) + (353._r16/12)*a(13)*a(17) + (2174._r16/63)*a(14)*a(17) + (32._r16/9)*a(15)*a(17) + (206._r16/15)*a(16)*a(17) + (1376._r16*a(17)**2)/45._r16 + (166._r16/5)*a(4)*a(18) +                                                                                                       &
                                  (184._r16/3)*a(5)*a(18) + (682._r16/7)*a(6)*a(18) + (283._r16/2)*a(7)*a(18) + (1._r16/2)*a(8)*a(18) + (29._r16/10)*a(9)*a(18) + (229._r16/12)*a(10)*a(18) + (1747._r16/14)*a(11)*a(18) + (921._r16/4)*a(12)*a(18) + (46063._r16/126)*a(13)*a(18) + (21233._r16/40)*a(14)*a(18) +                                                                                                     &
                                  (8._r16/3)*a(15)*a(18) + (66._r16/5)*a(16)*a(18) + (763._r16/9)*a(17)*a(18) + (145069._r16*a(18)**2)/525._r16 + (184._r16/3)*a(4)*a(19) + (55936._r16/105)*a(5)*a(19) + (3685._r16/3)*a(6)*a(19) + (147424._r16/63)*a(7)*a(19) + (2._r16/5)*a(8)*a(19) + (17._r16/6)*a(9)*a(19) +                                                                                                    &
                                  (170._r16/7)*a(10)*a(19) + (921._r16/4)*a(11)*a(19) + (125870._r16/63)*a(12)*a(19) + (92129._r16/20)*a(13)*a(19) + (2027122._r16/231)*a(14)*a(19) + (32._r16/15)*a(15)*a(19) + (578._r16/45)*a(16)*a(19) + (56576._r16/525)*a(17)*a(19) + (15292._r16/15)*a(18)*a(19) +                                                                                                              &
                                  (20894896._r16*a(19)**2)/4725._r16 + (682._r16/7)*a(4)*a(20) + (3685._r16/3)*a(5)*a(20) + (2517470._r16/189)*a(6)*a(20) + 36852._r16*a(7)*a(20) + (1._r16/3)*a(8)*a(20) + (39._r16/14)*a(9)*a(20) + (353._r16/12)*a(10)*a(20) + (46063._r16/126)*a(11)*a(20) + (92129._r16/20)*a(12)*a(20) +                                                                                         &
                                  (69230677._r16*a(13)*a(20))/1386._r16 + (829171._r16/6)*a(14)*a(20) + (16._r16/9)*a(15)*a(20) + (1322._r16/105)*a(16)*a(20) + (1174._r16/9)*a(17)*a(20) + (1529482._r16/945)*a(18)*a(20) + (101957._r16/5)*a(19)*a(20) + (229846019._r16*a(20)**2)/2079._r16 + (283._r16/2)*a(4)*a(21) +                                                                                             &
                                  (147424._r16/63)*a(5)*a(21) + 36852._r16*a(6)*a(21) + (110769184._r16/231)*a(7)*a(21) + (2._r16/7)*a(8)*a(21) + (11._r16/4)*a(9)*a(21) + (2174._r16/63)*a(10)*a(21) + (21233._r16/40)*a(11)*a(21) + (2027122._r16/231)*a(12)*a(21) + (829171._r16/6)*a(13)*a(21) +                                                                                                                   &
                                  (1799999394._r16*a(14)*a(21))/1001._r16 + (32._r16/21)*a(15)*a(21) + (62._r16/5)*a(16)*a(21) + (144544._r16/945)*a(17)*a(21) + (70499._r16/30)*a(18)*a(21) + (134601472._r16*a(19)*a(21))/3465._r16 + (27528484._r16/45)*a(20)*a(21) + (59759980928._r16*a(21)**2)/15015._r16 + 2._r16*a(8)*a(22) +                                                                                  &
                                  a(9)*a(22) + (2._r16/3)*a(10)*a(22) + (1._r16/2)*a(11)*a(22) + (2._r16/5)*a(12)*a(22) + (1._r16/3)*a(13)*a(22) + (2._r16/7)*a(14)*a(22) + 15._r16*a(15)*a(22) + (15._r16/2)*a(16)*a(22) + 5._r16*a(17)*a(22) + (15._r16/4)*a(18)*a(22) + 3._r16*a(19)*a(22) + (5._r16/2)*a(20)*a(22) +                                                                                               &
                                  (15._r16/7)*a(21)*a(22) + (249._r16*a(22)**2)/5._r16 + (1._r16/2)*a(4)*a(23) + (1._r16/2)*a(5)*a(23) + (1._r16/2)*a(6)*a(23) + (1._r16/2)*a(7)*a(23) + a(8)*a(23) + (46._r16/15)*a(9)*a(23) + (29._r16/10)*a(10)*a(23) + (14._r16/5)*a(11)*a(23) + (41._r16/15)*a(12)*a(23) +                                                                                                        &
                                  (94._r16/35)*a(13)*a(23) + (53._r16/20)*a(14)*a(23) + (15._r16/2)*a(15)*a(23) + (61._r16/3)*a(16)*a(23) + (229._r16/12)*a(17)*a(23) + (55._r16/3)*a(18)*a(23) + (107._r16/6)*a(19)*a(23) + (367._r16/21)*a(20)*a(23) + (413._r16/24)*a(21)*a(23) + (249._r16/5)*a(22)*a(23) +                                                                                                        &
                                  (2329._r16*a(23)**2)/35._r16 + (15._r16/4)*a(4)*a(24) + (24._r16/5)*a(5)*a(24) + (35._r16/6)*a(6)*a(24) + (48._r16/7)*a(7)*a(24) + (2._r16/3)*a(8)*a(24) + (29._r16/10)*a(9)*a(24) + (66._r16/5)*a(10)*a(24) + (55._r16/3)*a(11)*a(24) + (4082._r16/175)*a(12)*a(24) + (113._r16/4)*a(13)*a(24) +                                                                                    &
                                  (10438._r16/315)*a(14)*a(24) + 5._r16*a(15)*a(24) + (229._r16/12)*a(16)*a(24) + (763._r16/9)*a(17)*a(24) + (235._r16/2)*a(18)*a(24) + (5227._r16/35)*a(19)*a(24) + (13015._r16/72)*a(20)*a(24) + (4451._r16/21)*a(21)*a(24) + (166._r16/5)*a(22)*a(24) + (1747._r16/14)*a(23)*a(24) +                                                                                                &
                                  (145069._r16*a(24)**2)/525._r16 + (249._r16/10)*a(4)*a(25) + 46._r16*a(5)*a(25) + (1023._r16/14)*a(6)*a(25) + (849._r16/8)*a(7)*a(25) + (1._r16/2)*a(8)*a(25) + (14._r16/5)*a(9)*a(25) + (55._r16/3)*a(10)*a(25) + (20966._r16/175)*a(11)*a(25) + (4421._r16/20)*a(12)*a(25) +                                                                                                       &
                                  (110554._r16/315)*a(13)*a(25) + (2548._r16/5)*a(14)*a(25) + (15._r16/4)*a(15)*a(25) + (55._r16/3)*a(16)*a(25) + (235._r16/2)*a(17)*a(25) + (26801._r16/35)*a(18)*a(25) + (33901._r16/24)*a(19)*a(25) + (47093._r16/21)*a(20)*a(25) + 3256._r16*a(21)*a(25) + (249._r16/10)*a(22)*a(25) +                                                                                             &
                                  (20966._r16/175)*a(23)*a(25) + (26801._r16/35)*a(24)*a(25) + (436497._r16*a(25)**2)/175._r16 + 46._r16*a(4)*a(26) + (13984._r16/35)*a(5)*a(26) + (3685._r16/4)*a(6)*a(26) + (36856._r16/21)*a(7)*a(26) + (2._r16/5)*a(8)*a(26) + (41._r16/15)*a(9)*a(26) + (4082._r16/175)*a(10)*a(26) +                                                                                             &
                                  (4421._r16/20)*a(11)*a(26) + (3020894._r16*a(12)*a(26))/1575._r16 + (22111._r16/5)*a(13)*a(26) + (3243398._r16/385)*a(14)*a(26) + 3._r16*a(15)*a(26) + (107._r16/6)*a(16)*a(26) + (5227._r16/35)*a(17)*a(26) + (33901._r16/24)*a(18)*a(26) + (428901._r16/35)*a(19)*a(26) +                                                                                                          &
                                  (169519._r16/6)*a(20)*a(26) + (37299217._r16/693)*a(21)*a(26) + (498._r16/25)*a(22)*a(26) + (4077._r16/35)*a(23)*a(26) + (170298._r16/175)*a(24)*a(26) + (1288271._r16/140)*a(25)*a(26) + (146684527._r16*a(26)**2)/3675._r16 + (1023._r16/14)*a(4)*a(27) + (3685._r16/4)*a(5)*a(27) +                                                                                               &
                                  (1258735._r16/126)*a(6)*a(27) + 27639._r16*a(7)*a(27) + (1._r16/3)*a(8)*a(27) + (94._r16/35)*a(9)*a(27) + (113._r16/4)*a(10)*a(27) + (110554._r16/315)*a(11)*a(27) + (22111._r16/5)*a(12)*a(27) + (11076910._r16/231)*a(13)*a(27) + (3980021._r16/30)*a(14)*a(27) +                                                                                                                  &
                                  (5._r16/2)*a(15)*a(27) + (367._r16/21)*a(16)*a(27) + (13015._r16/72)*a(17)*a(27) + (47093._r16/21)*a(18)*a(27) + (169519._r16/6)*a(19)*a(27) + (636922745._r16*a(20)*a(27))/2079._r16 + (3390389._r16/4)*a(21)*a(27) + (83._r16/5)*a(22)*a(27) + (3994._r16/35)*a(23)*a(27) +                                                                                                        &
                                  (70667._r16/60)*a(24)*a(27) + (10737358._r16/735)*a(25)*a(27) + (32208643._r16/175)*a(26)*a(27) + (24203066389._r16*a(27)**2)/24255._r16 + (849._r16/8)*a(4)*a(28) + (36856._r16/21)*a(5)*a(28) + 27639._r16*a(6)*a(28) + (27692296._r16/77)*a(7)*a(28) + (2._r16/7)*a(8)*a(28) +                                                                                                    &
                                  (53._r16/20)*a(9)*a(28) + (10438._r16/315)*a(10)*a(28) + (2548._r16/5)*a(11)*a(28) + (3243398._r16/385)*a(12)*a(28) + (3980021._r16/30)*a(13)*a(28) + (8639997122._r16*a(14)*a(28))/5005._r16 + (15._r16/7)*a(15)*a(28) + (413._r16/24)*a(16)*a(28) + (4451._r16/21)*a(17)*a(28) +                                                                                                   &
                                  3256._r16*a(18)*a(28) + (37299217._r16/693)*a(19)*a(28) + (3390389._r16/4)*a(20)*a(28) + (33119989481._r16*a(21)*a(28))/3003._r16 + (498._r16/35)*a(22)*a(28) + (15727._r16/140)*a(23)*a(28) + (1014982._r16/735)*a(24)*a(28) + (3711873._r16/175)*a(25)*a(28) +                                                                                                                     &
                                  (2834741878._r16*a(26)*a(28))/8085._r16 + (386504357._r16/70)*a(27)*a(28) + (1258559602819._r16*a(28)**2)/35035._r16 + 2._r16*a(8)*a(29) + a(9)*a(29) + (2._r16/3)*a(10)*a(29) + (1._r16/2)*a(11)*a(29) + (2._r16/5)*a(12)*a(29) + (1._r16/3)*a(13)*a(29) + (2._r16/7)*a(14)*a(29) +                                                                                                 &
                                  (96._r16/5)*a(15)*a(29) + (48._r16/5)*a(16)*a(29) + (32._r16/5)*a(17)*a(29) + (24._r16/5)*a(18)*a(29) + (96._r16/25)*a(19)*a(29) + (16._r16/5)*a(20)*a(29) + (96._r16/35)*a(21)*a(29) + 184._r16*a(22)*a(29) + 92._r16*a(23)*a(29) + (184._r16/3)*a(24)*a(29) + 46._r16*a(25)*a(29) +                                                                                                &
                                  (184._r16/5)*a(26)*a(29) + (92._r16/3)*a(27)*a(29) + (184._r16/7)*a(28)*a(29) + (27968._r16*a(29)**2)/35._r16 + (2._r16/5)*a(4)*a(30) + (2._r16/5)*a(5)*a(30) + (2._r16/5)*a(6)*a(30) + (2._r16/5)*a(7)*a(30) + a(8)*a(30) + 3._r16*a(9)*a(30) + (17._r16/6)*a(10)*a(30) +                                                                                                           &
                                  (41._r16/15)*a(11)*a(30) + (8._r16/3)*a(12)*a(30) + (55._r16/21)*a(13)*a(30) + (31._r16/12)*a(14)*a(30) + (48._r16/5)*a(15)*a(30) + (906._r16/35)*a(16)*a(30) + (170._r16/7)*a(17)*a(30) + (4082._r16/175)*a(18)*a(30) + (794._r16/35)*a(19)*a(30) + (778._r16/35)*a(20)*a(30) +                                                                                                     &
                                  (766._r16/35)*a(21)*a(30) + 92._r16*a(22)*a(30) + (2947._r16/12)*a(23)*a(30) + (921._r16/4)*a(24)*a(30) + (4421._r16/20)*a(25)*a(30) + (2579._r16/12)*a(26)*a(30) + (5895._r16/28)*a(27)*a(30) + (829._r16/4)*a(28)*a(30) + (27968._r16/35)*a(29)*a(30) + (335651._r16*a(30)**2)/315._r16 +                                                                                          &
                                  3._r16*a(4)*a(31) + (96._r16/25)*a(5)*a(31) + (14._r16/3)*a(6)*a(31) + (192._r16/35)*a(7)*a(31) + (2._r16/3)*a(8)*a(31) + (17._r16/6)*a(9)*a(31) + (578._r16/45)*a(10)*a(31) + (107._r16/6)*a(11)*a(31) + (794._r16/35)*a(12)*a(31) + (989._r16/36)*a(13)*a(31) + (290._r16/9)*a(14)*a(31) +                                                                                         &
                                  (32._r16/5)*a(15)*a(31) + (170._r16/7)*a(16)*a(31) + (56576._r16/525)*a(17)*a(31) + (5227._r16/35)*a(18)*a(31) + (33216._r16/175)*a(19)*a(31) + (3446._r16/15)*a(20)*a(31) + (197984._r16/735)*a(21)*a(31) + (184._r16/3)*a(22)*a(31) + (921._r16/4)*a(23)*a(31) +                                                                                                                   &
                                  (15292._r16/15)*a(24)*a(31) + (33901._r16/24)*a(25)*a(31) + (62828._r16/35)*a(26)*a(31) + (26071._r16/12)*a(27)*a(31) + (160480._r16/63)*a(28)*a(31) + (55936._r16/105)*a(29)*a(31) + (125870._r16/63)*a(30)*a(31) + (20894896._r16*a(31)**2)/4725._r16 + (498._r16/25)*a(4)*a(32) +                                                                                                 &
                                  (184._r16/5)*a(5)*a(32) + (2046._r16/35)*a(6)*a(32) + (849._r16/10)*a(7)*a(32) + (1._r16/2)*a(8)*a(32) + (41._r16/15)*a(9)*a(32) + (107._r16/6)*a(10)*a(32) + (4077._r16/35)*a(11)*a(32) + (2579._r16/12)*a(12)*a(32) + (3071._r16/9)*a(13)*a(32) + (9909._r16/20)*a(14)*a(32) +                                                                                                     &
                                  (24._r16/5)*a(15)*a(32) + (4082._r16/175)*a(16)*a(32) + (5227._r16/35)*a(17)*a(32) + (170298._r16/175)*a(18)*a(32) + (62828._r16/35)*a(19)*a(32) + (2094626._r16/735)*a(20)*a(32) + (1448217._r16/350)*a(21)*a(32) + 46._r16*a(22)*a(32) + (4421._r16/20)*a(23)*a(32) +                                                                                                              &
                                  (33901._r16/24)*a(24)*a(32) + (1288271._r16/140)*a(25)*a(32) + 16974._r16*a(26)*a(32) + (6790711._r16/252)*a(27)*a(32) + (3130037._r16/80)*a(28)*a(32) + (13984._r16/35)*a(29)*a(32) + (3020894._r16*a(30)*a(32))/1575._r16 + (428901._r16/35)*a(31)*a(32) +                                                                                                                         &
                                  (146684527._r16*a(32)**2)/3675._r16 + (184._r16/5)*a(4)*a(33) + (55936._r16/175)*a(5)*a(33) + 737._r16*a(6)*a(33) + (147424._r16/105)*a(7)*a(33) + (2._r16/5)*a(8)*a(33) + (8._r16/3)*a(9)*a(33) + (794._r16/35)*a(10)*a(33) + (2579._r16/12)*a(11)*a(33) + (83914._r16/45)*a(12)*a(33) +                                                                                            &
                                  (128981._r16/30)*a(13)*a(33) + (810850._r16/99)*a(14)*a(33) + (96._r16/25)*a(15)*a(33) + (794._r16/35)*a(16)*a(33) + (33216._r16/175)*a(17)*a(33) + (62828._r16/35)*a(18)*a(33) + (57230368._r16*a(19)*a(33))/3675._r16 + (6283261._r16/175)*a(20)*a(33) +                                                                                                                           &
                                  (553001536._r16*a(21)*a(33))/8085._r16 + (184._r16/5)*a(22)*a(33) + (2579._r16/12)*a(23)*a(33) + (62828._r16/35)*a(24)*a(33) + 16974._r16*a(25)*a(33) + (46384376._r16/315)*a(26)*a(33) + (13579961._r16/40)*a(27)*a(33) + (149399660._r16/231)*a(28)*a(33) +                                                                                                                        &
                                  (55936._r16/175)*a(29)*a(33) + (83914._r16/45)*a(30)*a(33) + (57230368._r16*a(31)*a(33))/3675._r16 + (46384376._r16/315)*a(32)*a(33) + (7041838976._r16*a(33)**2)/11025._r16 + (2046._r16/35)*a(4)*a(34) + 737._r16*a(5)*a(34) + (503494._r16/63)*a(6)*a(34) + (110556._r16/5)*a(7)*a(34) +                                                                                          &
                                  (1._r16/3)*a(8)*a(34) + (55._r16/21)*a(9)*a(34) + (989._r16/36)*a(10)*a(34) + (3071._r16/9)*a(11)*a(34) + (128981._r16/30)*a(12)*a(34) + (13846139._r16/297)*a(13)*a(34) + (773893._r16/6)*a(14)*a(34) + (16._r16/5)*a(15)*a(34) + (778._r16/35)*a(16)*a(34) +                                                                                                                       &
                                  (3446._r16/15)*a(17)*a(34) + (2094626._r16/735)*a(18)*a(34) + (6283261._r16/175)*a(19)*a(34) + (9443072306._r16*a(20)*a(34))/24255._r16 + (37699652._r16/35)*a(21)*a(34) + (92._r16/3)*a(22)*a(34) + (5895._r16/28)*a(23)*a(34) + (26071._r16/12)*a(24)*a(34) +                                                                                                                      &
                                  (6790711._r16/252)*a(25)*a(34) + (13579961._r16/40)*a(26)*a(34) + (10204611013._r16*a(27)*a(34))/2772._r16 + (61109921._r16/6)*a(28)*a(34) + (27968._r16/105)*a(29)*a(34) + (4027882._r16*a(30)*a(34))/2205._r16 + (17811074._r16/945)*a(31)*a(34) + (515465914._r16*a(32)*a(34))/2205._r16 +                                                                                        &
                                  (4638690187._r16*a(33)*a(34))/1575._r16 + (3485726218351._r16*a(34)**2)/218295._r16 + (849._r16/10)*a(4)*a(35) + (147424._r16/105)*a(5)*a(35) + (110556._r16/5)*a(6)*a(35) + (110769184._r16/385)*a(7)*a(35) + (2._r16/7)*a(8)*a(35) + (31._r16/12)*a(9)*a(35) + (290._r16/9)*a(10)*a(35) +                                                                                          &
                                  (9909._r16/20)*a(11)*a(35) + (810850._r16/99)*a(12)*a(35) + (773893._r16/6)*a(13)*a(35) + (719999762._r16/429)*a(14)*a(35) + (96._r16/35)*a(15)*a(35) + (766._r16/35)*a(16)*a(35) + (197984._r16/735)*a(17)*a(35) + (1448217._r16/350)*a(18)*a(35) +                                                                                                                                 &
                                  (553001536._r16*a(19)*a(35))/8085._r16 + (37699652._r16/35)*a(20)*a(35) + (44639985856._r16*a(21)*a(35))/3185._r16 + (184._r16/7)*a(22)*a(35) + (829._r16/4)*a(23)*a(35) + (160480._r16/63)*a(24)*a(35) + (3130037._r16/80)*a(25)*a(35) + (149399660._r16/231)*a(26)*a(35) +                                                                                                         &
                                  (61109921._r16/6)*a(27)*a(35) + (12059996196._r16/91)*a(28)*a(35) + (55936._r16/245)*a(29)*a(35) + (566422._r16/315)*a(30)*a(35) + (48726976._r16*a(31)*a(35))/2205._r16 + (356389813._r16*a(32)*a(35))/1050._r16 + (408259618912._r16*a(33)*a(35))/72765._r16 +                                                                                                                     &
                                  (9277394428._r16/105)*a(34)*a(35) + (181257782850736._r16*a(35)**2)/315315._r16 + 2._r16*a(8)*a(36) + a(9)*a(36) + (2._r16/3)*a(10)*a(36) + (1._r16/2)*a(11)*a(36) + (2._r16/5)*a(12)*a(36) + (1._r16/3)*a(13)*a(36) + (2._r16/7)*a(14)*a(36) + (70._r16/3)*a(15)*a(36) +                                                                                                            &
                                  (35._r16/3)*a(16)*a(36) + (70._r16/9)*a(17)*a(36) + (35._r16/6)*a(18)*a(36) + (14._r16/3)*a(19)*a(36) + (35._r16/9)*a(20)*a(36) + (10._r16/3)*a(21)*a(36) + (2046._r16/7)*a(22)*a(36) + (1023._r16/7)*a(23)*a(36) + (682._r16/7)*a(24)*a(36) + (1023._r16/14)*a(25)*a(36) +                                                                                                          &
                                  (2046._r16/35)*a(26)*a(36) + (341._r16/7)*a(27)*a(36) + (2046._r16/49)*a(28)*a(36) + 3685._r16*a(29)*a(36) + (3685._r16/2)*a(30)*a(36) + (3685._r16/3)*a(31)*a(36) + (3685._r16/4)*a(32)*a(36) + 737._r16*a(33)*a(36) + (3685._r16/6)*a(34)*a(36) + (3685._r16/7)*a(35)*a(36) +                                                                                                      &
                                  (1258735._r16*a(36)**2)/63._r16 + (1._r16/3)*a(4)*a(37) + (1._r16/3)*a(5)*a(37) + (1._r16/3)*a(6)*a(37) + (1._r16/3)*a(7)*a(37) + a(8)*a(37) + (62._r16/21)*a(9)*a(37) + (39._r16/14)*a(10)*a(37) + (94._r16/35)*a(11)*a(37) + (55._r16/21)*a(12)*a(37) + (18._r16/7)*a(13)*a(37) +                                                                                                  &
                                  (71._r16/28)*a(14)*a(37) + (35._r16/3)*a(15)*a(37) + (1129._r16/36)*a(16)*a(37) + (353._r16/12)*a(17)*a(37) + (113._r16/4)*a(18)*a(37) + (989._r16/36)*a(19)*a(37) + (323._r16/12)*a(20)*a(37) + (53._r16/2)*a(21)*a(37) + (1023._r16/7)*a(22)*a(37) + (24566._r16/63)*a(23)*a(37) +                                                                                                 &
                                  (46063._r16/126)*a(24)*a(37) + (110554._r16/315)*a(25)*a(37) + (3071._r16/9)*a(26)*a(37) + (147410._r16/441)*a(27)*a(37) + (82919._r16/252)*a(28)*a(37) + (3685._r16/2)*a(29)*a(37) + (73703._r16/15)*a(30)*a(37) + (92129._r16/20)*a(31)*a(37) + (22111._r16/5)*a(32)*a(37) +                                                                                                       &
                                  (128981._r16/30)*a(33)*a(37) + (147407._r16/35)*a(34)*a(37) + (165833._r16/40)*a(35)*a(37) + (1258735._r16/63)*a(36)*a(37) + (55384529._r16*a(37)**2)/2079._r16 + (5._r16/2)*a(4)*a(38) + (16._r16/5)*a(5)*a(38) + (35._r16/9)*a(6)*a(38) + (32._r16/7)*a(7)*a(38) + (2._r16/3)*a(8)*a(38) +                                                                                         &
                                  (39._r16/14)*a(9)*a(38) + (1322._r16/105)*a(10)*a(38) + (367._r16/21)*a(11)*a(38) + (778._r16/35)*a(12)*a(38) + (323._r16/12)*a(13)*a(38) + (13922._r16/441)*a(14)*a(38) + (70._r16/9)*a(15)*a(38) + (353._r16/12)*a(16)*a(38) + (1174._r16/9)*a(17)*a(38) +                                                                                                                         &
                                  (13015._r16/72)*a(18)*a(38) + (3446._r16/15)*a(19)*a(38) + (5005._r16/18)*a(20)*a(38) + (61618._r16/189)*a(21)*a(38) + (682._r16/7)*a(22)*a(38) + (46063._r16/126)*a(23)*a(38) + (1529482._r16/945)*a(24)*a(38) + (47093._r16/21)*a(25)*a(38) + (2094626._r16/735)*a(26)*a(38) +                                                                                                     &
                                  (2607541._r16/756)*a(27)*a(38) + (594470._r16/147)*a(28)*a(38) + (3685._r16/3)*a(29)*a(38) + (92129._r16/20)*a(30)*a(38) + (101957._r16/5)*a(31)*a(38) + (169519._r16/6)*a(32)*a(38) + (6283261._r16/175)*a(33)*a(38) + (347637._r16/8)*a(34)*a(38) +                                                                                                                                &
                                  (16049039._r16/315)*a(35)*a(38) + (2517470._r16/189)*a(36)*a(38) + (69230677._r16*a(37)*a(38))/1386._r16 + (229846019._r16*a(38)**2)/2079._r16 + (83._r16/5)*a(4)*a(39) + (92._r16/3)*a(5)*a(39) + (341._r16/7)*a(6)*a(39) + (283._r16/4)*a(7)*a(39) + (1._r16/2)*a(8)*a(39) +                                                                                                       &
                                  (94._r16/35)*a(9)*a(39) + (367._r16/21)*a(10)*a(39) + (3994._r16/35)*a(11)*a(39) + (5895._r16/28)*a(12)*a(39) + (147410._r16/441)*a(13)*a(39) + (16987._r16/35)*a(14)*a(39) + (35._r16/6)*a(15)*a(39) + (113._r16/4)*a(16)*a(39) + (13015._r16/72)*a(17)*a(39) +                                                                                                                     &
                                  (70667._r16/60)*a(18)*a(39) + (26071._r16/12)*a(19)*a(39) + (2607541._r16/756)*a(20)*a(39) + (240379._r16/48)*a(21)*a(39) + (1023._r16/14)*a(22)*a(39) + (110554._r16/315)*a(23)*a(39) + (47093._r16/21)*a(24)*a(39) + (10737358._r16/735)*a(25)*a(39) +                                                                                                                             &
                                  (6790711._r16/252)*a(26)*a(39) + (6288722._r16/147)*a(27)*a(39) + (6521974._r16/105)*a(28)*a(39) + (3685._r16/4)*a(29)*a(39) + (22111._r16/5)*a(30)*a(39) + (169519._r16/6)*a(31)*a(39) + (32208643._r16/175)*a(32)*a(39) + (13579961._r16/40)*a(33)*a(39) +                                                                                                                         &
                                  (169777157._r16/315)*a(34)*a(39) + (3912761._r16/5)*a(35)*a(39) + (1258735._r16/126)*a(36)*a(39) + (11076910._r16/231)*a(37)*a(39) + (636922745._r16*a(38)*a(39))/2079._r16 + (24203066389._r16*a(39)**2)/24255._r16 + (92._r16/3)*a(4)*a(40) + (27968._r16/105)*a(5)*a(40) +                                                                                                        &
                                  (3685._r16/6)*a(6)*a(40) + (73712._r16/63)*a(7)*a(40) + (2._r16/5)*a(8)*a(40) + (55._r16/21)*a(9)*a(40) + (778._r16/35)*a(10)*a(40) + (5895._r16/28)*a(11)*a(40) + (4027882._r16*a(12)*a(40))/2205._r16 + (147407._r16/35)*a(13)*a(40) + (12973606._r16*a(14)*a(40))/1617._r16 +                                                                                                     &
                                  (14._r16/3)*a(15)*a(40) + (989._r16/36)*a(16)*a(40) + (3446._r16/15)*a(17)*a(40) + (26071._r16/12)*a(18)*a(40) + (17811074._r16/945)*a(19)*a(40) + (347637._r16/8)*a(20)*a(40) + (57367834._r16/693)*a(21)*a(40) + (2046._r16/35)*a(22)*a(40) + (3071._r16/9)*a(23)*a(40) +                                                                                                          &
                                  (2094626._r16/735)*a(24)*a(40) + (6790711._r16/252)*a(25)*a(40) + (515465914._r16*a(26)*a(40))/2205._r16 + (169777157._r16/315)*a(27)*a(40) + (1358399890._r16*a(28)*a(40))/1323._r16 + 737._r16*a(29)*a(40) + (128981._r16/30)*a(30)*a(40) + (6283261._r16/175)*a(31)*a(40) +                                                                                                       &
                                  (13579961._r16/40)*a(32)*a(40) + (4638690187._r16*a(33)*a(40))/1575._r16 + (13580699._r16/2)*a(34)*a(40) + (452750829._r16/35)*a(35)*a(40) + (503494._r16/63)*a(36)*a(40) + (13846139._r16/297)*a(37)*a(40) + (9443072306._r16*a(38)*a(40))/24255._r16 +                                                                                                                             &
                                  (10204611013._r16*a(39)*a(40))/2772._r16 + (3485726218351._r16*a(40)**2)/218295._r16 + (341._r16/7)*a(4)*a(41) + (3685._r16/6)*a(5)*a(41) + (1258735._r16/189)*a(6)*a(41) + 18426._r16*a(7)*a(41) + (1._r16/3)*a(8)*a(41) + (18._r16/7)*a(9)*a(41) + (323._r16/12)*a(10)*a(41) +                                                                                                     &
                                  (147410._r16/441)*a(11)*a(41) + (147407._r16/35)*a(12)*a(41) + (221538242._r16*a(13)*a(41))/4851._r16 + (5306695._r16/42)*a(14)*a(41) + (35._r16/9)*a(15)*a(41) + (323._r16/12)*a(16)*a(41) + (5005._r16/18)*a(17)*a(41) + (2607541._r16/756)*a(18)*a(41) +                                                                                                                          &
                                  (347637._r16/8)*a(19)*a(41) + (3918459695._r16*a(20)*a(41))/8316._r16 + (11732773._r16/9)*a(21)*a(41) + (341._r16/7)*a(22)*a(41) + (147410._r16/441)*a(23)*a(41) + (2607541._r16/756)*a(24)*a(41) + (6288722._r16/147)*a(25)*a(41) + (169777157._r16/315)*a(26)*a(41) +                                                                                                              &
                                  (23196074042._r16*a(27)*a(41))/3969._r16 + (226369893._r16/14)*a(28)*a(41) + (3685._r16/6)*a(29)*a(41) + (147407._r16/35)*a(30)*a(41) + (347637._r16/8)*a(31)*a(41) + (169777157._r16/315)*a(32)*a(41) + (13580699._r16/2)*a(33)*a(41) + (1546237109._r16/21)*a(34)*a(41) +                                                                                                          &
                                  (12222647561._r16/60)*a(35)*a(41) + (1258735._r16/189)*a(36)*a(41) + (221538242._r16*a(37)*a(41))/4851._r16 + (3918459695._r16*a(38)*a(41))/8316._r16 + (23196074042._r16*a(39)*a(41))/3969._r16 + (1546237109._r16/21)*a(40)*a(41) + (17428710403085._r16*a(41)**2)/43659._r16 +                                                                                                    &
                                  (283._r16/4)*a(4)*a(42) + (73712._r16/63)*a(5)*a(42) + 18426._r16*a(6)*a(42) + (55384592._r16/231)*a(7)*a(42) + (2._r16/7)*a(8)*a(42) + (71._r16/28)*a(9)*a(42) + (13922._r16/441)*a(10)*a(42) + (16987._r16/35)*a(11)*a(42) + (12973606._r16*a(12)*a(42))/1617._r16 +                                                                                                               &
                                  (5306695._r16/42)*a(13)*a(42) + (11519996214._r16*a(14)*a(42))/7007._r16 + (10._r16/3)*a(15)*a(42) + (53._r16/2)*a(16)*a(42) + (61618._r16/189)*a(17)*a(42) + (240379._r16/48)*a(18)*a(42) + (57367834._r16/693)*a(19)*a(42) + (11732773._r16/9)*a(20)*a(42) +                                                                                                                       &
                                  (50939983882._r16*a(21)*a(42))/3003._r16 + (2046._r16/49)*a(22)*a(42) + (82919._r16/252)*a(23)*a(42) + (594470._r16/147)*a(24)*a(42) + (6521974._r16/105)*a(25)*a(42) + (1358399890._r16*a(26)*a(42))/1323._r16 + (226369893._r16/14)*a(27)*a(42) +                                                                                                                                  &
                                  (13268155815766._r16*a(28)*a(42))/63063._r16 + (3685._r16/7)*a(29)*a(42) + (165833._r16/40)*a(30)*a(42) + (16049039._r16/315)*a(31)*a(42) + (3912761._r16/5)*a(32)*a(42) + (452750829._r16/35)*a(33)*a(42) + (12222647561._r16/60)*a(34)*a(42) +                                                                                                                                     &
                                  (13266715817221._r16*a(35)*a(42))/5005._r16 + (2517470._r16/441)*a(36)*a(42) + (13846141._r16/308)*a(37)*a(42) + (24119988934._r16*a(38)*a(42))/43659._r16 + (5880459043._r16/693)*a(39)*a(42) + (2041307389622._r16*a(40)*a(42))/14553._r16 + (9184662275813._r16*a(41)*a(42))/4158._r16 +                                                                                          &
                                  (9969223416919553._r16*a(42)**2)/693693._r16 + 2._r16*a(8)*a(43) + a(9)*a(43) + (2._r16/3)*a(10)*a(43) + (1._r16/2)*a(11)*a(43) + (2._r16/5)*a(12)*a(43) + (1._r16/3)*a(13)*a(43) + (2._r16/7)*a(14)*a(43) + (192._r16/7)*a(15)*a(43) + (96._r16/7)*a(16)*a(43) + (64._r16/7)*a(17)*a(43) +                                                                                          &
                                  (48._r16/7)*a(18)*a(43) + (192._r16/35)*a(19)*a(43) + (32._r16/7)*a(20)*a(43) + (192._r16/49)*a(21)*a(43) + (849._r16/2)*a(22)*a(43) + (849._r16/4)*a(23)*a(43) + (283._r16/2)*a(24)*a(43) + (849._r16/8)*a(25)*a(43) + (849._r16/10)*a(26)*a(43) + (283._r16/4)*a(27)*a(43) +                                                                                                       &
                                  (849._r16/14)*a(28)*a(43) + (147424._r16/21)*a(29)*a(43) + (73712._r16/21)*a(30)*a(43) + (147424._r16/63)*a(31)*a(43) + (36856._r16/21)*a(32)*a(43) + (147424._r16/105)*a(33)*a(43) + (73712._r16/63)*a(34)*a(43) + (147424._r16/147)*a(35)*a(43) + 110556._r16*a(36)*a(43) +                                                                                                        &
                                  55278._r16*a(37)*a(43) + 36852._r16*a(38)*a(43) + 27639._r16*a(39)*a(43) + (110556._r16/5)*a(40)*a(43) + 18426._r16*a(41)*a(43) + (110556._r16/7)*a(42)*a(43) + (55384592._r16*a(43)**2)/77._r16 + (2._r16/7)*a(4)*a(44) + (2._r16/7)*a(5)*a(44) + (2._r16/7)*a(6)*a(44) + (2._r16/7)*a(7)*a(44) +                                                                                   &
                                  a(8)*a(44) + (35._r16/12)*a(9)*a(44) + (11._r16/4)*a(10)*a(44) + (53._r16/20)*a(11)*a(44) + (31._r16/12)*a(12)*a(44) + (71._r16/28)*a(13)*a(44) + (5._r16/2)*a(14)*a(44) + (96._r16/7)*a(15)*a(44) + (2318._r16/63)*a(16)*a(44) + (2174._r16/63)*a(17)*a(44) +                                                                                                                       &
                                  (10438._r16/315)*a(18)*a(44) + (290._r16/9)*a(19)*a(44) + (13922._r16/441)*a(20)*a(44) + (1958._r16/63)*a(21)*a(44) + (849._r16/4)*a(22)*a(44) + (2831._r16/5)*a(23)*a(44) + (21233._r16/40)*a(24)*a(44) + (2548._r16/5)*a(25)*a(44) + (9909._r16/20)*a(26)*a(44) +                                                                                                                  &
                                  (16987._r16/35)*a(27)*a(44) + (38221._r16/80)*a(28)*a(44) + (73712._r16/21)*a(29)*a(44) + (6486782._r16/693)*a(30)*a(44) + (2027122._r16/231)*a(31)*a(44) + (3243398._r16/385)*a(32)*a(44) + (810850._r16/99)*a(33)*a(44) + (12973606._r16*a(34)*a(44))/1617._r16 +                                                                                                                  &
                                  (608138._r16/77)*a(35)*a(44) + 55278._r16*a(36)*a(44) + (884449._r16/6)*a(37)*a(44) + (829171._r16/6)*a(38)*a(44) + (3980021._r16/30)*a(39)*a(44) + (773893._r16/6)*a(40)*a(44) + (5306695._r16/42)*a(41)*a(44) + (373127._r16/3)*a(42)*a(44) + (55384592._r16/77)*a(43)*a(44) +                                                                                                     &
                                  (2879999015._r16*a(44)**2)/3003._r16 + (15._r16/7)*a(4)*a(45) + (96._r16/35)*a(5)*a(45) + (10._r16/3)*a(6)*a(45) + (192._r16/49)*a(7)*a(45) + (2._r16/3)*a(8)*a(45) + (11._r16/4)*a(9)*a(45) + (62._r16/5)*a(10)*a(45) + (413._r16/24)*a(11)*a(45) + (766._r16/35)*a(12)*a(45) +                                                                                                     &
                                  (53._r16/2)*a(13)*a(45) + (1958._r16/63)*a(14)*a(45) + (64._r16/7)*a(15)*a(45) + (2174._r16/63)*a(16)*a(45) + (144544._r16/945)*a(17)*a(45) + (4451._r16/21)*a(18)*a(45) + (197984._r16/735)*a(19)*a(45) + (61618._r16/189)*a(20)*a(45) + (56192._r16/147)*a(21)*a(45) +                                                                                                             &
                                  (283._r16/2)*a(22)*a(45) + (21233._r16/40)*a(23)*a(45) + (70499._r16/30)*a(24)*a(45) + 3256._r16*a(25)*a(45) + (1448217._r16/350)*a(26)*a(45) + (240379._r16/48)*a(27)*a(45) + (1233041._r16/210)*a(28)*a(45) + (147424._r16/63)*a(29)*a(45) + (2027122._r16/231)*a(30)*a(45) +                                                                                                      &
                                  (134601472._r16*a(31)*a(45))/3465._r16 + (37299217._r16/693)*a(32)*a(45) + (553001536._r16*a(33)*a(45))/8085._r16 + (57367834._r16/693)*a(34)*a(45) + (1412505632._r16*a(35)*a(45))/14553._r16 + 36852._r16*a(36)*a(45) + (829171._r16/6)*a(37)*a(45) + (27528484._r16/45)*a(38)*a(45) +                                                                                             &
                                  (3390389._r16/4)*a(39)*a(45) + (37699652._r16/35)*a(40)*a(45) + (11732773._r16/9)*a(41)*a(45) + (10699380._r16/7)*a(42)*a(45) + (110769184._r16/231)*a(43)*a(45) + (1799999394._r16*a(44)*a(45))/1001._r16 + (59759980928._r16*a(45)**2)/15015._r16 + (498._r16/35)*a(4)*a(46) +                                                                                                     &
                                  (184._r16/7)*a(5)*a(46) + (2046._r16/49)*a(6)*a(46) + (849._r16/14)*a(7)*a(46) + (1._r16/2)*a(8)*a(46) + (53._r16/20)*a(9)*a(46) + (413._r16/24)*a(10)*a(46) + (15727._r16/140)*a(11)*a(46) + (829._r16/4)*a(12)*a(46) + (82919._r16/252)*a(13)*a(46) + (38221._r16/80)*a(14)*a(46) +                                                                                                &
                                  (48._r16/7)*a(15)*a(46) + (10438._r16/315)*a(16)*a(46) + (4451._r16/21)*a(17)*a(46) + (1014982._r16/735)*a(18)*a(46) + (160480._r16/63)*a(19)*a(46) + (594470._r16/147)*a(20)*a(46) + (1233041._r16/210)*a(21)*a(46) + (849._r16/8)*a(22)*a(46) + (2548._r16/5)*a(23)*a(46) +                                                                                                        &
                                  3256._r16*a(24)*a(46) + (3711873._r16/175)*a(25)*a(46) + (3130037._r16/80)*a(26)*a(46) + (6521974._r16/105)*a(27)*a(46) + (3607401._r16/40)*a(28)*a(46) + (36856._r16/21)*a(29)*a(46) + (3243398._r16/385)*a(30)*a(46) + (37299217._r16/693)*a(31)*a(46) +                                                                                                                           &
                                  (2834741878._r16*a(32)*a(46))/8085._r16 + (149399660._r16/231)*a(33)*a(46) + (1358399890._r16*a(34)*a(46))/1323._r16 + (3443692649._r16*a(35)*a(46))/2310._r16 + 27639._r16*a(36)*a(46) + (3980021._r16/30)*a(37)*a(46) + (3390389._r16/4)*a(38)*a(46) + (386504357._r16/70)*a(39)*a(46) +                                                                                           &
                                  (61109921._r16/6)*a(40)*a(46) + (226369893._r16/14)*a(41)*a(46) + (939064079._r16/40)*a(42)*a(46) + (27692296._r16/77)*a(43)*a(46) + (8639997122._r16*a(44)*a(46))/5005._r16 + (33119989481._r16*a(45)*a(46))/3003._r16 + (1258559602819._r16*a(46)**2)/35035._r16 + (184._r16/7)*a(4)*a(47) +                                                                                       &
                                  (55936._r16/245)*a(5)*a(47) + (3685._r16/7)*a(6)*a(47) + (147424._r16/147)*a(7)*a(47) + (2._r16/5)*a(8)*a(47) + (31._r16/12)*a(9)*a(47) + (766._r16/35)*a(10)*a(47) + (829._r16/4)*a(11)*a(47) + (566422._r16/315)*a(12)*a(47) + (165833._r16/40)*a(13)*a(47) +                                                                                                                      &
                                  (608138._r16/77)*a(14)*a(47) + (192._r16/35)*a(15)*a(47) + (290._r16/9)*a(16)*a(47) + (197984._r16/735)*a(17)*a(47) + (160480._r16/63)*a(18)*a(47) + (48726976._r16*a(19)*a(47))/2205._r16 + (16049039._r16/315)*a(20)*a(47) + (1412505632._r16*a(21)*a(47))/14553._r16 +                                                                                                            &
                                  (849._r16/10)*a(22)*a(47) + (9909._r16/20)*a(23)*a(47) + (1448217._r16/350)*a(24)*a(47) + (3130037._r16/80)*a(25)*a(47) + (356389813._r16*a(26)*a(47))/1050._r16 + (3912761._r16/5)*a(27)*a(47) + (3443692649._r16*a(28)*a(47))/2310._r16 + (147424._r16/105)*a(29)*a(47) +                                                                                                          &
                                  (810850._r16/99)*a(30)*a(47) + (553001536._r16*a(31)*a(47))/8085._r16 + (149399660._r16/231)*a(32)*a(47) + (408259618912._r16*a(33)*a(47))/72765._r16 + (452750829._r16/35)*a(34)*a(47) + (119542288576._r16*a(35)*a(47))/4851._r16 + (110556._r16/5)*a(36)*a(47) +                                                                                                                  &
                                  (773893._r16/6)*a(37)*a(47) + (37699652._r16/35)*a(38)*a(47) + (61109921._r16/6)*a(39)*a(47) + (9277394428._r16/105)*a(40)*a(47) + (12222647561._r16/60)*a(41)*a(47) + (268934398220._r16/693)*a(42)*a(47) + (110769184._r16/385)*a(43)*a(47) +                                                                                                                                      &
                                  (719999762._r16/429)*a(44)*a(47) + (44639985856._r16*a(45)*a(47))/3185._r16 + (12059996196._r16/91)*a(46)*a(47) + (181257782850736._r16*a(47)**2)/315315._r16 + (2046._r16/49)*a(4)*a(48) + (3685._r16/7)*a(5)*a(48) + (2517470._r16/441)*a(6)*a(48) + (110556._r16/7)*a(7)*a(48) +                                                                                                  &
                                  (1._r16/3)*a(8)*a(48) + (71._r16/28)*a(9)*a(48) + (53._r16/2)*a(10)*a(48) + (82919._r16/252)*a(11)*a(48) + (165833._r16/40)*a(12)*a(48) + (13846141._r16/308)*a(13)*a(48) + (373127._r16/3)*a(14)*a(48) + (32._r16/7)*a(15)*a(48) + (13922._r16/441)*a(16)*a(48) +                                                                                                                   &
                                  (61618._r16/189)*a(17)*a(48) + (594470._r16/147)*a(18)*a(48) + (16049039._r16/315)*a(19)*a(48) + (24119988934._r16*a(20)*a(48))/43659._r16 + (10699380._r16/7)*a(21)*a(48) + (283._r16/4)*a(22)*a(48) + (16987._r16/35)*a(23)*a(48) + (240379._r16/48)*a(24)*a(48) +                                                                                                                 &
                                  (6521974._r16/105)*a(25)*a(48) + (3912761._r16/5)*a(26)*a(48) + (5880459043._r16/693)*a(27)*a(48) + (939064079._r16/40)*a(28)*a(48) + (73712._r16/63)*a(29)*a(48) + (12973606._r16*a(30)*a(48))/1617._r16 + (57367834._r16/693)*a(31)*a(48) + (1358399890._r16*a(32)*a(48))/1323._r16 +                                                                                              &
                                  (452750829._r16/35)*a(33)*a(48) + (2041307389622._r16*a(34)*a(48))/14553._r16 + (268934398220._r16/693)*a(35)*a(48) + 18426._r16*a(36)*a(48) + (5306695._r16/42)*a(37)*a(48) + (11732773._r16/9)*a(38)*a(48) + (226369893._r16/14)*a(39)*a(48) +                                                                                                                                     &
                                  (12222647561._r16/60)*a(40)*a(48) + (9184662275813._r16*a(41)*a(48))/4158._r16 + 6111332994._r16*a(42)*a(48) + (55384592._r16/231)*a(43)*a(48) + (11519996214._r16*a(44)*a(48))/7007._r16 + (50939983882._r16*a(45)*a(48))/3003._r16 + (13268155815766._r16*a(46)*a(48))/63063._r16 +                                                                                                &
                                  (13266715817221._r16*a(47)*a(48))/5005._r16 + (9969223416919553._r16*a(48)**2)/693693._r16 + a(2)*(2._r16*a(3) + 2._r16*a(4) + 2._r16*a(5) + 2._r16*a(6) + 2._r16*a(7) + a(9) + a(10) + a(11) + a(12) + a(13) + a(14) + (2._r16*a(16))/3._r16 + (2._r16*a(17))/3._r16 + (2._r16*a(18))/3._r16 + (2._r16*a(19))/3._r16 +                                                              &
                                    (2._r16*a(20))/3._r16 + (2._r16*a(21))/3._r16 + a(23)/2._r16 + a(24)/2._r16 + a(25)/2._r16 + a(26)/2._r16 + a(27)/2._r16 + a(28)/2._r16 + (2._r16*a(30))/5._r16 + (2._r16*a(31))/5._r16 + (2._r16*a(32))/5._r16 + (2._r16*a(33))/5._r16 + (2._r16*a(34))/5._r16 + (2._r16*a(35))/5._r16 + a(37)/3._r16 + a(38)/3._r16 + a(39)/3._r16 + a(40)/3._r16 + a(41)/3._r16 +               &
                                    a(42)/3._r16 + (2._r16*a(44))/7._r16 + (2._r16*a(45))/7._r16 + (2._r16*a(46))/7._r16 + (2._r16*a(47))/7._r16 + (2._r16*a(48))/7._r16 + (2._r16*a(49))/7) + (849._r16/14)*a(4)*a(49) + (147424._r16/147)*a(5)*a(49) + (110556._r16/7)*a(6)*a(49) + (110769184._r16/539)*a(7)*a(49) + (2._r16/7)*a(8)*a(49) +                                                                        &
                                  (5._r16/2)*a(9)*a(49) + (1958._r16/63)*a(10)*a(49) + (38221._r16/80)*a(11)*a(49) + (608138._r16/77)*a(12)*a(49) + (373127._r16/3)*a(13)*a(49) + (1619999470._r16*a(14)*a(49))/1001._r16 + (192._r16/49)*a(15)*a(49) + (1958._r16/63)*a(16)*a(49) + (56192._r16/147)*a(17)*a(49) +                                                                                                    &
                                  (1233041._r16/210)*a(18)*a(49) + (1412505632._r16*a(19)*a(49))/14553._r16 + (10699380._r16/7)*a(20)*a(49) + (1254239603488._r16*a(21)*a(49))/63063._r16 + (849._r16/14)*a(22)*a(49) + (38221._r16/80)*a(23)*a(49) + (1233041._r16/210)*a(24)*a(49) + (3607401._r16/40)*a(25)*a(49) +                                                                                                 &
                                  (3443692649._r16*a(26)*a(49))/2310._r16 + (939064079._r16/40)*a(27)*a(49) + (3057839035777._r16*a(28)*a(49))/10010._r16 + (147424._r16/147)*a(29)*a(49) + (608138._r16/77)*a(30)*a(49) + (1412505632._r16*a(31)*a(49))/14553._r16 + (3443692649._r16*a(32)*a(49))/2310._r16 +                                                                                                        &
                                  (119542288576._r16*a(33)*a(49))/4851._r16 + (268934398220._r16/693)*a(34)*a(49) + (1167627951869504._r16*a(35)*a(49))/231231._r16 + (110556._r16/7)*a(36)*a(49) + (373127._r16/3)*a(37)*a(49) + (10699380._r16/7)*a(38)*a(49) + (939064079._r16/40)*a(39)*a(49) +                                                                                                                    &
                                  (268934398220._r16/693)*a(40)*a(49) + 6111332994._r16*a(41)*a(49) + (238801244711212._r16*a(42)*a(49))/3003._r16 + (110769184._r16/539)*a(43)*a(49) + (1619999470._r16*a(44)*a(49))/1001._r16 + (1254239603488._r16*a(45)*a(49))/63063._r16 + (3057839035777._r16*a(46)*a(49))/10010._r16 +                                                                                          &
                                  (1167627951869504._r16*a(47)*a(49))/231231._r16 + (238801244711212._r16*a(48)*a(49))/3003._r16 + (39876897932311200._r16*a(49)**2)/77077._r16 +                                                                                                                                                                                                                                      &
                                  a(3)*(15._r16*a(4) + (96._r16*a(5))/5._r16 + (70._r16*a(6))/3._r16 + (192._r16*a(7))/7._r16 + a(9) + (16._r16*a(10))/3._r16 + (15._r16*a(11))/2._r16 + (48._r16*a(12))/5._r16 + (35._r16*a(13))/3._r16 + (96._r16*a(14))/7._r16 + (2._r16*a(16))/3._r16 + (32._r16*a(17))/9._r16 + 5._r16*a(18) + (32._r16*a(19))/5._r16 + (70._r16*a(20))/9._r16 + (64._r16*a(21))/7._r16 +         &
                                    a(23)/2._r16 + (8._r16*a(24))/3._r16 + (15._r16*a(25))/4._r16 + (24._r16*a(26))/5._r16 + (35._r16*a(27))/6._r16 + (48._r16*a(28))/7._r16 + (2._r16*a(30))/5._r16 + (32._r16*a(31))/15._r16 + 3._r16*a(32) + (96._r16*a(33))/25._r16 + (14._r16*a(34))/3._r16 + (192._r16*a(35))/35._r16 + a(37)/3._r16 + (16._r16*a(38))/9._r16 + (5._r16*a(39))/2._r16 + (16._r16*a(40))/5._r16 + &
                                    (35._r16*a(41))/9._r16 + (32._r16*a(42))/7._r16 + (2._r16*a(44))/7._r16 + (32._r16*a(45))/21._r16 + (15._r16*a(46))/7._r16 + (96._r16*a(47))/35._r16 + (10._r16*a(48))/3._r16 + (192._r16*a(49))/49)
      end function WENO_smooth_indicator_7
      
    end module reconstruction_mod
    