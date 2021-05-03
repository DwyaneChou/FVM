    module weno_coef_mod
      use constants_mod
      implicit none
      private
      public calc_weno_coef
      
    contains
      subroutine calc_weno_coef(coef,availableStencil,iAWENO,wenoIdx,xw,yw,stencil_width_sub,lack_pos,stencil_width,nPointsOnEdge)
        real   (r_kind), dimension(:,:,:), allocatable, intent(out) :: coef
        logical        , dimension(:,:  ), allocatable, intent(out) :: availableStencil
        real   (r_kind), dimension(:,:,:), allocatable, intent(out) :: iAWENO
        integer(i_kind), dimension(:,:  ), allocatable, intent(out) :: wenoIdx
        real   (r_kind), dimension(:    ), allocatable, intent(out) :: xw
        real   (r_kind), dimension(:    ), allocatable, intent(out) :: yw
        integer(i_kind)                               , intent(out) :: stencil_width_sub
        integer(i_kind), dimension(:,:  ), allocatable, intent(out) :: lack_pos
        integer(i_kind)                               , intent(in ) :: stencil_width
        integer(i_kind)                               , intent(in ) :: nPointsOnEdge
        
        real   (r16), dimension(:,:), allocatable :: c
        
        real   (r16), dimension(:), allocatable :: quad_pos_1d
        real   (r16), dimension(:), allocatable :: quad_wts_1d
        
        real   (r16), dimension(:), allocatable :: xC
        real   (r16), dimension(:), allocatable :: yC
        
        real   (r16), dimension(:), allocatable :: xP
        real   (r16), dimension(:), allocatable :: yP
        
        integer(i_kind), dimension(:  ), allocatable :: iSC  ! center cell index for each stencil
        integer(i_kind), dimension(:,:), allocatable :: POS  ! index of cells on full stencil
        integer(i_kind), dimension(:,:), allocatable :: COSL ! cell indices in substencil
        integer(i_kind), dimension(:  ), allocatable :: COSH ! cell indices in full stencil
        real   (r16   ), dimension(:,:), allocatable :: existTermL
        real   (r16   ), dimension(:  ), allocatable :: existTermH
        integer(i_kind), dimension(:  ), allocatable :: nExistTermL
        integer(i_kind)                              :: nExistTermH
        
        integer(i_kind), dimension(:  ), allocatable :: nAvailableStencil
        
        real   (r16), dimension(:,:,:), allocatable :: AL
        real   (r16), dimension(:,:,:), allocatable :: iAL
        real   (r16), dimension(:,:  ), allocatable :: AH
        real   (r16), dimension(:,:  ), allocatable :: iAH
        
        real   (r16), dimension(:,:,:), allocatable :: PL
        real   (r16), dimension(  :,:), allocatable :: PH
        
        real   (r16), dimension(:,:,:), allocatable :: RL
        real   (r16), dimension(  :,:), allocatable :: RH
        real   (r16), dimension(:,:,:), allocatable :: RL_on_RH
        real   (r16), dimension(  :,:), allocatable :: RL_LS
        real   (r16), dimension(  :,:), allocatable :: iRL_LS
        
        integer(i_kind), dimension(:  ), allocatable :: lack
        
        real   (r16), dimension(:,:), allocatable :: post
        real   (r16), dimension(  :), allocatable :: post_check
        
        integer(i_kind) :: nStencil, nPoints, nWenoType, recBdy, nWenoLack
        
        integer(i_kind) :: recBdy_sub
        
        integer(i_kind) :: nCOSL ! Cell number On Stencil for Low order stencil(substencil)
        integer(i_kind) :: nCOSH ! Cell number On Stencil for High order stencil(full stencil)

        integer(i_kind) :: ix_min,ix_max
        integer(i_kind) :: iy_min,iy_max
        
        integer(i_kind) :: nx,ny
        real   (r16) :: x_min, x_max
        real   (r16) :: y_min, y_max
        
        integer(i_kind) :: i,j,k,iCOS,iRCOS,iType,iPOS,iStencil,jStencil,iPoint
        integer(i_kind) :: iC,jC
        integer(i_kind) :: nRCL,nRCH
        integer(i_kind) :: nAS
        integer(i_kind) :: invStat
        
        recBdy = ( stencil_width - 1 ) / 2
        
        if(stencil_width==3)then
          nStencil          = 4
          stencil_width_sub = 2
          nCOSL             = 6
        else
          nStencil          = ( stencil_width - 2 )**2
          stencil_width_sub = 3
          nCOSL             = stencil_width_sub**2
        endif
        recBdy_sub = 1
        nCOSH      = stencil_width**2
        
        nWenoType = recBdy**2 + 1
        nWenoLack = recBdy**2
        nPoints   = 4 * nPointsOnEdge + nPointsOnEdge**2
        
        ix_min = -recBdy
        ix_max = recBdy
        iy_min = -recBdy
        iy_max = recBdy
        
        allocate( coef            (nWenoType,nPoints,nStencil) )
        allocate( availableStencil(nWenoType,nStencil) )
        allocate( iAWENO          (nStencil,nCOSL,nCOSL) )
        allocate( wenoIdx         (nStencil,nCOSL) )
        allocate( xw              (nPoints) )
        allocate( yw              (nPoints) )
        allocate( lack_pos        (nWenoType,nWenoLack) )
        
        allocate( quad_pos_1d(nPointsOnEdge)            )
        allocate( quad_wts_1d(nPointsOnEdge)            )
        allocate( xC      (nCOSH)                       )
        allocate( yC      (nCOSH)                       )
        allocate( xP      (nPoints)                     )
        allocate( yP      (nPoints)                     )
        allocate( lack    (nWenoLack)                   )
        allocate( c       (nPoints,nStencil)            )
        allocate( iSC     (nStencil)                    )
        allocate( POS     (ix_min:ix_max,iy_min:iy_max) )
        allocate( COSL    (nStencil,nCOSL)              )
        allocate( COSH    (         nCOSH)              )
        allocate( AL      (nStencil,nCOSL,nCOSL)        )
        allocate( iAL     (nStencil,nCOSL,nCOSL)        )
        allocate( AH      (nCOSH,nCOSH)                 )
        allocate( iAH     (nCOSH,nCOSH)                 )
        allocate( PL      (nStencil,nPoints,nCOSL)      )
        allocate( PH      (         nPoints,nCOSH)      )
        allocate( RL      (nPoints,nStencil,nCOSL)      )
        allocate( RH      (nPoints,         nCOSH)      )
        allocate( RL_on_RH(nPoints,nStencil,nCOSH)      )
        allocate( RL_LS   (nStencil,nStencil)           )
        allocate( iRL_LS  (nStencil,nStencil)           )
        
        allocate( existTermL(nStencil,nCOSL) )
        allocate( existTermH(         nCOSH) )
        
        allocate( nAvailableStencil(nWenoType) )
        
        allocate( nExistTermL(nStencil) )
        
        allocate( post      (nStencil, nCOSH) )
        allocate( post_check(          nCOSH) )
        
        call Gaussian_Legendre(nPointsOnEdge, quad_pos_1d, quad_wts_1d)
        
        quad_pos_1d = quad_pos_1d / 2._r16
        
        ! Left
        xP(nPointsOnEdge*0+1:nPointsOnEdge*1) = -0.5_r16
        yP(nPointsOnEdge*0+1:nPointsOnEdge*1) = quad_pos_1d
        ! Right
        xP(nPointsOnEdge*1+1:nPointsOnEdge*2) = 0.5_r16
        yP(nPointsOnEdge*1+1:nPointsOnEdge*2) = quad_pos_1d
        ! Bottom
        xP(nPointsOnEdge*2+1:nPointsOnEdge*3) = quad_pos_1d
        yP(nPointsOnEdge*2+1:nPointsOnEdge*3) = -0.5_r16
        ! Top
        xP(nPointsOnEdge*3+1:nPointsOnEdge*4) = quad_pos_1d
        yP(nPointsOnEdge*3+1:nPointsOnEdge*4) = 0.5_r16
        ! Quadrature Points
        do j = 1,nPointsOnEdge
          xP(nPointsOnEdge*(4+j-1)+1:nPointsOnEdge*(4+j)) = quad_pos_1d
          yP(nPointsOnEdge*(4+j-1)+1:nPointsOnEdge*(4+j)) = quad_pos_1d(j)
        enddo
        
        xw = xP
        yw = yP
        
        ! lack_pos number must be ascending
        lack_pos = 0
        if(stencil_width==3)then
          lack_pos(2,1) = 9
        elseif(stencil_width==5)then
          lack_pos(2,:) = (/19,20,24,25/)
          lack_pos(3,:) = (/20,25, 0, 0/)
          lack_pos(4,:) = (/24,25, 0, 0/)
          lack_pos(5,:) = (/25, 0, 0, 0/)
        elseif(stencil_width==7)then
          lack_pos( 2,:) = (/33,34,35,40,41,42,47,48,49/)
          lack_pos( 3,:) = (/34,35,41,42,48,49, 0, 0, 0/)
          lack_pos( 4,:) = (/35,42,49, 0, 0, 0, 0, 0, 0/)
          lack_pos( 5,:) = (/40,41,42,47,48,49, 0, 0, 0/)
          lack_pos( 6,:) = (/41,42,48,49, 0, 0, 0, 0, 0/)
          lack_pos( 7,:) = (/42,49, 0, 0, 0, 0, 0, 0, 0/)
          lack_pos( 8,:) = (/47,48,49, 0, 0, 0, 0, 0, 0/)
          lack_pos( 9,:) = (/48,49, 0, 0, 0, 0, 0, 0, 0/)
          lack_pos(10,:) = (/49, 0, 0, 0, 0, 0, 0, 0, 0/)
        endif
        
        iCOS = 0
        do j = -recBdy,recBdy
          do i = -recBdy,recBdy
            iCOS = iCOS + 1
            POS(i,j) = iCOS
            xC(iCOS) = i
            yC(iCOS) = j
          enddo
        enddo
        
        if(stencil_width==3)then
          iSC(1) = 1
          iSC(2) = 3
          iSC(3) = 7
          iSC(4) = 9
        elseif(stencil_width>3)then
          iStencil = 0
          do j = -recBdy+recBdy_sub,recBdy-recBdy_sub
            do i = -recBdy+recBdy_sub,recBdy-recBdy_sub
              iStencil = iStencil + 1
              iPOS = POS(i,j)
              iSC(iStencil) = iPOS
            enddo
          enddo
        endif
        
        availableStencil  = .false.
        nAvailableStencil = 0
        
        do iType = 1,nWenoType
          existTermL = 0
          existTermH = 0
          lack       = lack_pos( iType, : )
          COSL       = 0
          COSH       = 0
          AL         = 0
          iAL        = 0
          AH         = 0
          iAH        = 0
          PL         = 0
          PH         = 0
          RL         = 0
          RH         = 0
          RL_on_RH   = 0
          
          if(stencil_width==3)then
            if(iType==1)then
              availableStencil(iType,:) = .true.
              nAvailableStencil(iType) = 4
            elseif(iType==2)then
              availableStencil(iType,1) = .true.
              nAvailableStencil(iType) = 1
            endif
          elseif(stencil_width>3)then
            do jStencil = 1,nStencil
              if( .not.any(lack==iSC(jStencil)) )then
                availableStencil(iType,jStencil) = .true.
                nAvailableStencil(iType) = nAvailableStencil(iType) + 1
              endif
            enddo
          endif
          nAS = nAvailableStencil(iType)
          
          ! Calculate inverse matrices on substencil
          if(stencil_width==3)then
            existTermL = 1
            COSL(1,:) = (/1,2,3,4,5,7/)
            COSL(2,:) = (/1,2,3,5,6,9/)
            COSL(3,:) = (/1,4,5,7,8,9/)
            COSL(4,:) = (/3,5,6,7,8,9/)
            
            iStencil = 0
            do jStencil = 1,nStencil
              if(availableStencil(iType,jStencil))then
                iStencil = iStencil + 1
                
                nRCL = nCOSL
                do iCOS = 1,nRCL
                  x_min = xC(COSL(iStencil,iCOS)) - 0.5_r16
                  x_max = xC(COSL(iStencil,iCOS)) + 0.5_r16
                  y_min = yC(COSL(iStencil,iCOS)) - 0.5_r16
                  y_max = yC(COSL(iStencil,iCOS)) + 0.5_r16
                  
                  call calc_polynomial_square_integration(2,x_min,x_max,y_min,y_max,AL(iStencil,iCOS,1:nRCL))
                enddo
                
                call BRINV(nRCL,AL(iStencil,1:nRCL,1:nRCL),iAL(iStencil,1:nRCL,1:nRCL),invStat)
                if(invStat==0)then
                  print*,'WENO substencil coef calculation failed, iType, iStencil', iType, iStencil
                  stop
                endif
                
                call calc_polynomial_matrix(2,nPoints,nRCL,xP,yP,PL(iStencil,1:nPoints,1:nRCL))
                
                RL(1:nPoints,iStencil,1:nRCL) = matmul( PL(iStencil,1:nPoints,1:nRCL), iAL(iStencil,1:nRCL,1:nRCL) )
              endif
            enddo
          elseif(stencil_width>3)then
            iStencil = 0
            do jStencil = 1,nStencil
              if( availableStencil(iType,jStencil) )then
                iStencil = iStencil + 1
                iCOS = 0
                iRCOS = 0
                do j = -recBdy_sub,recBdy_sub
                  do i = -recBdy_sub,recBdy_sub
                    iC = xC(iSC(jStencil)) + i
                    jC = yC(iSC(jStencil)) + j
                    if( iC>=ix_min .and. iC<=ix_max .and. jC>=iy_min .and. jC<=iy_max )then
                      iCOS = iCOS + 1
                      if( .not.any(lack==POS(iC,jC)) )then
                        iRCOS = iRCOS + 1
                        existTermL(iStencil,iCOS) = 1
                        COSL(iStencil,iRCOS) = POS(iC,jC)
                      endif
                    endif
                  enddo ! i
                enddo ! j
                nExistTermL(iStencil) = sum(existTermL(iStencil,:))
                
                nRCL = nExistTermL(iStencil)
                nx   = stencil_width_sub
                ny   = stencil_width_sub
                do iCOS = 1,nRCL
                  x_min = xC(COSL(iStencil,iCOS)) - 0.5_r16
                  x_max = xC(COSL(iStencil,iCOS)) + 0.5_r16
                  y_min = yC(COSL(iStencil,iCOS)) - 0.5_r16
                  y_max = yC(COSL(iStencil,iCOS)) + 0.5_r16
                  
                  call calc_rectangle_poly_integration( nx,ny,x_min,x_max,y_min,y_max,AL(iStencil,iCOS,1:nRCL),existTermL(iStencil,:))
                enddo
                
                call BRINV(nRCL,AL(iStencil,1:nRCL,1:nRCL),iAL(iStencil,1:nRCL,1:nRCL),invStat)
                if(invStat==0)then
                  print*,'WENO substencil coef calculation failed, iType, iStencil', iType, iStencil
                  stop
                endif
                
                call calc_rectangle_poly_matrix(nx,ny,nPoints,xP,yP,PL(iStencil,1:nPoints,1:nRCL),existTermL(iStencil,:))
                
                RL(1:nPoints,iStencil,1:nRCL) = matmul( PL(iStencil,1:nPoints,1:nRCL), iAL(iStencil,1:nRCL,1:nRCL) )
              endif ! availableStencil
            enddo ! jStencil
          endif
          
          if(iType==1)then
            iAWENO  = iAL ! iAWENO in Full Type
            wenoIdx = COSL
          endif
          
          ! Calculate inverse matrices on full stencil
          iCOS = 0
          iRCOS = 0
          do j = -recBdy,recBdy
            do i = -recBdy,recBdy
              iCOS = iCOS + 1
              if( .not.any(lack==POS(i,j)) )then
                iRCOS = iRCOS + 1
                COSH(iCOS) = iRCOS
                existTermH(iCOS) = 1
              endif
            enddo
          enddo
          nExistTermH = sum(existTermH)
          
          nRCH = nExistTermH
          nx   = stencil_width
          ny   = stencil_width
          iRCOS = 0
          do j = -recBdy,recBdy
            do i = -recBdy,recBdy
              if( .not.any(lack==POS(i,j)) )then
                iRCOS = iRCOS + 1
                x_min = real(i,r16) - 0.5_r16
                x_max = real(i,r16) + 0.5_r16
                y_min = real(j,r16) - 0.5_r16
                y_max = real(j,r16) + 0.5_r16
                call calc_rectangle_poly_integration(nx,ny,x_min,x_max,y_min,y_max,AH(iRCOS,1:nRCH),existTermH)
              endif
            enddo
          enddo
          
          call BRINV(nRCH,AH(1:nRCH,1:nRCH),iAH(1:nRCH,1:nRCH),invStat)
          if(invStat==0)then
            print*,'WENO full stencil coef calculation failed'
            stop
          endif
          
          call calc_rectangle_poly_matrix(nx,ny,nPoints,xP,yP,PH(1:nPoints,1:nRCH),existTermH)
          
          RH(1:nPoints,1:nRCH) = matmul( PH(1:nPoints,1:nRCH), iAH(1:nRCH,1:nRCH) )
          
          if(stencil_width==3)then
            iStencil = 0
            do jStencil = 1,nStencil
              if( availableStencil(iType,jStencil) )then
                iStencil = iStencil + 1
                do iRCOS = 1,nCOSL
                  RL_on_RH(1:nPoints,iStencil,COSH(COSL(iStencil,iRCOS))) = RL(1:nPoints,iStencil,iRCOS)
                enddo
              endif
            enddo
          elseif(stencil_width>3)then
            iStencil = 0
            do jStencil = 1,nStencil
              if( availableStencil(iType,jStencil) )then
                iStencil = iStencil + 1
                iCOS     = 0
                iRCOS    = 0
                do j = -recBdy_sub,recBdy_sub
                  do i = -recBdy_sub,recBdy_sub
                    iC = xC(iSC(jStencil)) + i
                    jC = yC(iSC(jStencil)) + j
                    if( iC>=ix_min .and. iC<=ix_max .and. jC>=iy_min .and. jC<=iy_max )then
                      iCOS = iCOS + 1
                      if( .not.any(lack==POS(iC,jC)) )then
                        iRCOS = iRCOS + 1
                        RL_on_RH(1:nPoints,iStencil,COSH(COSL(iStencil,iRCOS))) = RL(1:nPoints,iStencil,iRCOS)
                      endif
                    endif
                  enddo
                enddo
              endif
            enddo
          endif
          
          do iPoint = 1,nPoints
            RL_LS(1:nAS,1:nAS) = matmul( RL_on_RH(iPoint,1:nAS,1:nRCH), transpose(RL_on_RH(iPoint,1:nAS,1:nRCH)) )
            call BRINV(nAS, RL_LS(1:nAS,1:nAS), iRL_LS(1:nAS,1:nAS), invStat)
            if(invStat==0)then
              print*,'WENO optimal coef calculatio fail'
              stop
            endif
            
            c(iPoint,1:nAS) = matmul( matmul( iRL_LS(1:nAS,1:nAS), RL_on_RH(iPoint,1:nAS,1:nRCH) ), RH(iPoint,1:nRCH) )
            
            coef(iType,iPoint,:) = 0
            iStencil = 0
            do jStencil = 1,nStencil
              if( availableStencil(iType,jStencil) )then
                iStencil = iStencil + 1
                coef(iType,iPoint,jStencil) = c(iPoint,iStencil)
              endif
            enddo
            
            if(stencil_width==3.and.iType==2)then
              coef(iType,iPoint,:) = 0
              coef(iType,iPoint,1) = 1
            endif
            
            !print*,coef(iType,iPoint,:)
            !print*,''
          enddo ! iPoint
          
          ! Post Check
          do iPoint = 1,nPoints
            do iStencil = 1,nAS
              post(iStencil,1:nRCH) = 0
              post(iStencil,1:nRCH) = c(iPoint,iStencil) * RL_on_RH(iPoint,iStencil,1:nRCH)
            enddo
            post_check = sum(post(1:nAS,:),1)
            !print*,'iType ',iType,' iPoint',iPoint,' diff=',sum(post_check(1:nRCH)-RH(iPoint,1:nRCH))
          enddo
        enddo ! iType
        
        !stop 'calc_weno_coef'
      end subroutine calc_weno_coef
    
      subroutine calc_rectangle_poly_integration(nx,ny,x_min,x_max,y_min,y_max,c,existPolyTerm)
        integer(i_kind), intent(in   ) :: nx  ! number of points on x direction
        integer(i_kind), intent(in   ) :: ny  ! number of points on y direction
        real   (r16), intent(in   ) :: x_min
        real   (r16), intent(in   ) :: x_max
        real   (r16), intent(in   ) :: y_min
        real   (r16), intent(in   ) :: y_max
        real   (r16), intent(inout) :: c(:)
        real   (r16), intent(in   ),optional :: existPolyTerm(nx*ny)
        
        real   (r16) :: ext(nx*ny)
        integer(i_kind) :: i,j,k,iCOS
        
        ext = 1
        if(present(existPolyTerm))ext = existPolyTerm
                
        k    = 0
        c    = 0
        iCOS = 0
        do j = 0,ny-1
          do i = 0,nx-1
            k = k + 1
            if(ext(k)>0)then
              iCOS = iCOS + 1
              c(iCOS) = ( x_max**(i+1) - x_min**(i+1) ) * ( y_max**(j+1) - y_min**(j+1) ) / real( ( i + 1 ) * ( j + 1 ), r16 )
            endif
          enddo
        enddo
        
      end subroutine  calc_rectangle_poly_integration
      
      subroutine calc_rectangle_poly_matrix(nx,ny,m,xi,eta,A,existPolyTerm)
        integer(i_kind), intent(in   ) :: nx ! number of points on x direction for reconstruction
        integer(i_kind), intent(in   ) :: ny ! number of points on y direction for reconstruction
        integer(i_kind), intent(in   ) :: m  ! number of unkonwn point values
        real   (r16), intent(in   ) :: xi (m)
        real   (r16), intent(in   ) :: eta(m)
        real   (r16), intent(inout) :: A  (:,:)
        real   (r16), intent(in   ),optional :: existPolyTerm(nx*ny)
        
        real   (r16) :: ext(nx*ny)
        
        real   (r16) :: x
        real   (r16) :: y
        integer(i_kind) :: iPOC
        integer(i_kind)  :: i,j,k,iCOS
        
        ext = 1
        if(present(existPolyTerm))ext = existPolyTerm
        
        do iPOC = 1,m
          x = xi (iPOC)
          y = eta(iPOC)
          
          k = 0
          iCOS = 0
          do j = 0,ny-1
            do i = 0,nx-1
              k = k + 1
              if(ext(k)>0)then
                iCOS = iCOS + 1
                A(iPOC,iCOS) = x**real(i,r16) * y**real(j,r16)
              endif
            enddo
          enddo
        enddo
      
      end subroutine calc_rectangle_poly_matrix
      
      subroutine  calc_polynomial_square_integration(d,x_min,x_max,y_min,y_max,c)
        integer(i_kind), intent(in ) :: d ! degree of polynomial
        real   (r16), intent(in ) :: x_min
        real   (r16), intent(in ) :: x_max
        real   (r16), intent(in ) :: y_min
        real   (r16), intent(in ) :: y_max
        real   (r16), intent(out) :: c(:)
        
        integer(i_kind) :: i,j,k
        
        k = 0
        c = 0
        do j = 0,d
          do i = 0,j
            k = k + 1
            c(k) = ( x_max**(j-i+1) - x_min**(j-i+1) ) * ( y_max**(i+1) - y_min**(i+1) ) / real( ( i + 1 ) * ( j - i + 1 ), r16 )
          enddo
        enddo
        
      end subroutine  calc_polynomial_square_integration
      
      subroutine calc_polynomial_matrix(d,m,n,xi,eta,A)
        integer(i_kind), intent(in ) :: d ! polynomial degree
        integer(i_kind), intent(in ) :: m ! number of points
        integer(i_kind), intent(in ) :: n ! number of terms of polynomial, n = (d+1)*(d+2)/2
        real   (r16), intent(in ) :: xi (m)
        real   (r16), intent(in ) :: eta(m)
        real   (r16), intent(out) :: A  (m,n)
        
        real   (r16) :: x
        real   (r16) :: y
        integer(i_kind) :: iPOC
        integer(i_kind)  :: i,j,k
        
        do iPOC = 1,m
          x = xi (iPOC)
          y = eta(iPOC)
          
          k = 0
          do j = 0,d
            do i = 0,j
              k = k + 1
              A(iPOC,k) = x**real(j-i,r16)*y**real(i,r16)
            enddo
          enddo
        enddo
      
      end subroutine calc_polynomial_matrix
      
      ! calculate inverse matrix of A_input
      ! N is the order of matrix A_input and A
      ! A is inverse A_input
      ! L is status info
      SUBROUTINE BRINV(N,A_input,A,L)
      implicit none
      integer(i_kind),intent(in )           :: N
      real   (r16),intent(in )           :: A_input(N,N)
      real   (r16),intent(out)           :: A      (N,N)
      integer(i_kind),intent(out), optional :: L
      
      real(r16):: T,D
      integer :: IS(N),JS(N)
      integer :: i,j,k
      
      A = A_input
      
      if(present(L))L=1
      do K=1,N
        D=0._r16
        do I=K,N
          do J=K,N
            IF (ABS(A(I,J)).GT.D) THEN
              D=ABS(A(I,J))
              IS(K)=I
              JS(K)=J
            END IF
          enddo
        enddo
      
        IF (D+1._r16.EQ.1._r16) THEN
          if(present(L))L=0
          WRITE(*,*)'ERR**NOT INV'
          RETURN
        END IF
        
        do J=1,N
          T=A(K,J)
          A(K,J)=A(IS(K),J)
          A(IS(K),J)=T
        enddo
        
        do I=1,N
          T=A(I,K)
          A(I,K)=A(I,JS(K))
          A(I,JS(K))=T
        enddo
        
        A(K,K)=1._r16/A(K,K)
        do J=1,N
          IF (J.NE.K) THEN
            A(K,J)=A(K,J)*A(K,K)
          END IF
        enddo
        
        do I=1,N
          IF (I.NE.K) THEN
            do J=1,N
              IF (J.NE.K) THEN
                A(I,J)=A(I,J)-A(I,K)*A(K,J)
              END IF
            enddo
          END IF
        enddo
        
        do I=1,N
          IF (I.NE.K) THEN
            A(I,K)=-A(I,K)*A(K,K)
          END IF
        enddo
      enddo
      
      do K=N,1,-1
        do J=1,N
          T=A(K,J)
          A(K,J)=A(JS(K),J)
          A(JS(K),J)=T
        enddo
        do I=1,N
          T=A(I,K)
          A(I,K)=A(I,IS(K))
          A(I,IS(K))=T
        enddo
      enddo
      RETURN
      END SUBROUTINE BRINV
      
      ! Gaussian quadrature from http://bbs.fcode.cn/thread-219-1-1.html, http://fcode.cn/algorithm-73-1.html
      Subroutine Gaussian_Legendre(n, p, w)
        Integer(i_kind),intent(in ) :: n
        Real   (r16),intent(out) :: p(n), w(n)
        
        Real   (r16   ) :: fn(n), ak(n)
        Real   (r16   ) :: m
        !定义数组,大小n由module开始声明。    
        Integer(i_kind) :: i, j
        j = 0 !赋值控制循环变量的初值           
        m = -1.000001 !设置计算域[-1，1] 的下限，即代替-1 
        Do i = 1, 200000 !这个循环次数应该是由步长0.00001决 定,计算方法：200000=2/0.000001     
          If (legendreP(m,n)*legendreP(m+0.00001,n)<0) Then !从下限处开始往上逐步累加，
            !由步长0.00001说明最多求解10^5个解
            j = j + 1 !记录这是第几个解
            fn(j) = bis(m, m+0.00001, n)
            !调用二分法求解程序在分好的一小段上求解，
            !将解存储在fn（j）
            ak(j) = 2.0_r16/(n*dLegendreP1(fn(j),n)*dLegendrePn(fn(j),n)) !高斯点的权重
            !Write (*, *) '高斯点序号', j
            !Write (*, *) '高斯点', fn(j)
            !Write (*, *) '高斯点权重', ak(j)
            p(j) = fn(j)
            w(j) = ak(j)
          End If
          m = m + 0.00001 !执行完一次判断m向前推进一步
        End Do
      End Subroutine Gaussian_Legendre
      
      Real(r16) Function legendreP(x,n) !定义Legendre函数
        Real   (r16   ),intent(in ) :: x
        Integer(i_kind),intent(in ) :: n
        
        Real   (r16   ) :: a(n) !a(n)代表n阶勒让德多项式
        Integer(i_kind) :: i
        a(1) = x !1阶勒让德多项式
        a(2) = 1.5_r16*(x**2) - 0.5_r16 !2阶勒让德多项式
        Do i = 3, n
          a(i) = (2*i-1)*x*a(i-1)/i - (i-1)*a(i-2)/i
          !利用递推关系产生n阶勒让德多项式
        End Do
        legendreP = a(n) !生成的n阶勒让德多项式   
      End Function legendreP
      
      Real(r16) Function dLegendreP1(x,n)!生成的（n-1）阶勒让德多项式  
        Real   (r16   ),intent(in ) :: x
        Integer(i_kind),intent(in ) :: n   
        
        Real   (r16   ) :: a(n) !a(n-1)代表（n-1）阶勒让德多项式
        Integer :: i
        a(1) = x
        a(2) = 1.5_r16*x**2 - 0.5_r16
        Do i = 3, n - 1
          a(i) = (2*i-1)*x*a(i-1)/i - (i-1)*a(i-2)/i
        End Do
        dLegendreP1 = a(n-1) 
      End Function dLegendreP1
      
      Real (r16) Function dLegendrePn(x,n)!生成n阶勒让德多项式的导数表达式
        Real   (r16   ),intent(in ) :: x
        Integer(i_kind),intent(in ) :: n
        
        Real   (r16   ) :: a(n)
        Integer :: i
        a(1) = x
        a(2) = 1.5_r16*x**2 - 0.5_r16
        Do i = 3, n
          a(i) = (2*i-1)*x*a(i-1)/i - (i-1)*a(i-2)/i
        End Do
        dLegendrePn = n*a(n-1)/(1-x**2) - n*x*a(n)/(1-x**2)
        
      End Function dLegendrePn
      
      Real (r16) Function bis(a_in, b_in,n) !二分法求解函数的解
        Real   (r16   ), intent(in ) :: a_in, b_in
        integer(i_kind), intent(in ) :: n
        
        Real   (r16   )  :: a, b
        Real   (r16   ) :: c
        !a,b是传递进来的划分好的有一个解存在的区间
        a = a_in
        b = b_in
        
        Do
          c = (a+b)/2.0_r16
          If (legendreP(c,n)*legendreP(a,n)<0) Then
            b = c
          Else
            a = c
          End If
          If ((b-a)<1e-34) exit 
        End Do
        bis = c!bis即是利用二分法求得的解
      End Function bis
    end module weno_coef_mod