    module weno_coef_mod
      use constants_mod
      implicit none
      private
      public calc_weno_coef
      
    contains
      subroutine calc_weno_coef(coef,stencil_width,nPointsOnEdge)
        real   (r_kind), dimension(:,:), allocatable, intent(out) :: coef
        integer(i_kind)                             , intent(in ) :: stencil_width
        integer(i_kind)                             , intent(in ) :: nPointsOnEdge
        
        real   (r16), dimension(nPointsOnEdge) :: quad_pos_1d
        real   (r16), dimension(nPointsOnEdge) :: quad_wts_1d
        
        real   (r16), dimension(stencil_width**2) :: xC
        real   (r16), dimension(stencil_width**2) :: yC
        
        real   (r16), dimension(:), allocatable :: x
        real   (r16), dimension(:), allocatable :: y
        
        integer(i_kind), dimension(:), allocatable :: lack_pos
        
        integer(i_kind) :: nStencil, nPoints, nWenoType, recBdy
        
        integer(i_kind) :: nCOSL ! Cell number On Stencil for Low order stencil
        integer(i_kind) :: nCOSH ! Cell number On Stencil for High order stencil
        
        integer(i_kind) :: i,j,k,iCOS,iRCOS
        
        recBdy = ( stencil_width - 1 ) / 2
        
        if(stencil_width==3)then
          nStencil  = 4
          nCOSL     = 4
        else
          nStencil  = ( stencil_width - 2 )**2
          nCOSL     = 9
        endif
        nCOSH = stencil_width**2
        
        nWenoType = recBdy**2 + 1
        
        nPoints = 4*nPointsOnEdge+nPointsOnEdge**2
        
        allocate( x(nPoints) )
        allocate( y(nPoints) )
        allocate( lack_pos(recBdy**2) )
        allocate( coef(nPoints*nWenoType,nStencil) )
        
        call Gaussian_Legendre(nPointsOnEdge, quad_pos_1d, quad_wts_1d)
        
        quad_pos_1d = quad_pos_1d / 2._r16
        
        ! Left
        x(nPointsOnEdge*0+1:nPointsOnEdge*1) = -0.5
        y(nPointsOnEdge*0+1:nPointsOnEdge*1) = quad_pos_1d
        ! Right
        x(nPointsOnEdge*1+1:nPointsOnEdge*2) = 0.5
        y(nPointsOnEdge*1+1:nPointsOnEdge*2) = quad_pos_1d
        ! Bottom
        x(nPointsOnEdge*2+1:nPointsOnEdge*3) = quad_pos_1d
        y(nPointsOnEdge*2+1:nPointsOnEdge*3) = -0.5
        ! Top
        x(nPointsOnEdge*3+1:nPointsOnEdge*4) = quad_pos_1d
        y(nPointsOnEdge*3+1:nPointsOnEdge*4) = 0.5
        ! Quadrature Points
        do j = 1,nPointsOnEdge
          x(nPointsOnEdge*(4+j-1)+1:nPointsOnEdge*(4+j)) = quad_pos_1d
          y(nPointsOnEdge*(4+j-1)+1:nPointsOnEdge*(4+j)) = quad_pos_1d(j)
        enddo
        
        coef = 0
      end subroutine calc_weno_coef
    
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
      
      real    :: T,D
      integer :: IS(N),JS(N)
      integer :: i,j,k
      
      A = A_input
      
      if(present(L))L=1
      do K=1,N
        D=0.
        do I=K,N
          do J=K,N
            IF (ABS(A(I,J)).GT.D) THEN
              D=ABS(A(I,J))
              IS(K)=I
              JS(K)=J
            END IF
          enddo
        enddo
      
        IF (D+1.0.EQ.1.0) THEN
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
        
        A(K,K)=1/A(K,K)
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
        !��������,��Сn��module��ʼ������    
        Integer(i_kind) :: i, j
        j = 0 !��ֵ����ѭ�������ĳ�ֵ           
        m = -1.000001 !���ü�����[-1��1] �����ޣ�������-1 
        Do i = 1, 200000 !���ѭ������Ӧ�����ɲ���0.00001�� ��,���㷽����200000=2/0.000001     
          If (legendreP(m,n)*legendreP(m+0.00001,n)<0) Then !�����޴���ʼ�������ۼӣ�
            !�ɲ���0.00001˵��������10^5����
            j = j + 1 !��¼���ǵڼ�����
            fn(j) = bis(m, m+0.00001, n)
            !���ö��ַ��������ڷֺõ�һС������⣬
            !����洢��fn��j��
            ak(j) = 2.0_r16/(n*dLegendreP1(fn(j),n)*dLegendrePn(fn(j),n)) !��˹���Ȩ��
            !Write (*, *) '��˹�����', j
            !Write (*, *) '��˹��', fn(j)
            !Write (*, *) '��˹��Ȩ��', ak(j)
            p(j) = fn(j)
            w(j) = ak(j)
          End If
          m = m + 0.00001 !ִ����һ���ж�m��ǰ�ƽ�һ��
        End Do
      End Subroutine Gaussian_Legendre
      
      Real(r16) Function legendreP(x,n) !����Legendre����
        Real   (r16   ),intent(in ) :: x
        Integer(i_kind),intent(in ) :: n
        
        Real   (r16   ) :: a(n) !a(n)����n�����õ¶���ʽ
        Integer(i_kind) :: i
        a(1) = x !1�����õ¶���ʽ
        a(2) = 1.5_r16*(x**2) - 0.5_r16 !2�����õ¶���ʽ
        Do i = 3, n
          a(i) = (2*i-1)*x*a(i-1)/i - (i-1)*a(i-2)/i
          !���õ��ƹ�ϵ����n�����õ¶���ʽ
        End Do
        legendreP = a(n) !���ɵ�n�����õ¶���ʽ   
      End Function legendreP
      
      Real(r16) Function dLegendreP1(x,n)!���ɵģ�n-1�������õ¶���ʽ  
        Real   (r16   ),intent(in ) :: x
        Integer(i_kind),intent(in ) :: n   
        
        Real   (r16   ) :: a(n) !a(n-1)������n-1�������õ¶���ʽ
        Integer :: i
        a(1) = x
        a(2) = 1.5_r16*x**2 - 0.5_r16
        Do i = 3, n - 1
          a(i) = (2*i-1)*x*a(i-1)/i - (i-1)*a(i-2)/i
        End Do
        dLegendreP1 = a(n-1) 
      End Function dLegendreP1
      
      Real (r16) Function dLegendrePn(x,n)!����n�����õ¶���ʽ�ĵ�������ʽ
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
      
      Real (r16) Function bis(a_in, b_in,n) !���ַ���⺯���Ľ�
        Real   (r16   ), intent(in ) :: a_in, b_in
        integer(i_kind), intent(in ) :: n
        
        Real   (r16   )  :: a, b
        Real   (r16   ) :: c
        !a,b�Ǵ��ݽ����Ļ��ֺõ���һ������ڵ�����
        a = a_in
        b = b_in
        
        Do
          c = (a+b)/2.0_r16
          If (legendreP(c,n)*legendreP(a,n)<0) Then
            b = c
          Else
            a = c
          End If
          If ((b-a)<1.e-16) exit 
        End Do
        bis = c!bis�������ö��ַ���õĽ�
      End Function bis
    end module weno_coef_mod