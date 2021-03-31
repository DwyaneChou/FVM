    module qr_solver_mod
    use constants_mod
    implicit none
    !----------------------------------------module coment
    !  Version     :  V1.0    
    !  Coded by    :  syz 
    !  Date        :  
    !-----------------------------------------------------
    !  Description :   ������Gram-Schimdt����������С��
    !                    ������ģ��
    !    
    !-----------------------------------------------------
    !  Parameters  :
    !      1.    solve     �ⳬ������ ��������
    !      2.    gram_dec  G-S ,QR�ֽ�
    !      3.    uptri     �����Ƿ��̻ش�����
    !      4.
    !-----------------------------------------------------
    !  Post Script :
    !      1.      �����Ե������� QR�ֽ⺯��
    !      2.      Ҳ���Ե��ýⷽ�̺���
    !-----------------------------------------------------
    
    contains
    
      subroutine qr_solver(M,N,A,b,x,x0)
        !---------------------------------subroutine  comment
        !  Version   :  V1.0    
        !  Coded by  :  syz 
        !  Date      :  
        !-----------------------------------------------------
        !  Purpose   :  ͨ��������Gram-Schmidt����������С������
        !               ��������
        !-----------------------------------------------------
        !  Method    :
        !               �Գ�������  A����QR�ֽ��   ���̱�Ϊ
        !                   QR x=b
        !                => Rx=Q'b   RΪ��������
        !                => �ش������������С���������µĽ�
        !-----------------------------------------------------
        !  Post Script :
        !       1.       ����� ����������  Ax=b    ���� A(M,N)  M>N    
        !       2.
        !----------------------------------------------------
        integer(i_kind)           :: M,N
        real   (r_kind)           :: A (M,N)
        real   (r_kind)           :: b (M)
        real   (r_kind)           :: x (N)
        real   (r_kind) ,optional :: x0(N)
        
        real   (r_kind) ::Q(M,N)
        real   (r_kind) ::R(N,N)
        real   (r_kind) ::QT(N,M)  !Q��ת�þ���
        real   (r_kind) ::QTb(N)   !Q'b
        real   (r_kind) ::AT (N,M)
        real   (r_kind) ::ATA(N,N)
        real   (r_kind) ::ATb(N)
        
        real   (r_kind) ::x0CG(N)
        
        !  For LAPACK only
        real   (r_kind), dimension(m+n) :: work
        real   (r_kind), dimension(n  ) :: tau
        integer(i_kind) :: INFO
        
        integer i,j
        
        ! Modified Gram-Schmidt
        call gram_dec(A,Q,R,M,N)
        
        QT  = transpose(Q)
        QTb = matmul(QT,b)  !  Rx=Q'b
        
        call uptri(R,QTb,x,N) !�ش�
        
        !! Cholesky QR(factorization)
        !AT  = transpose(A)
        !ATA = matmul( AT, A )
        !ATb = matmul( AT, b )
        !call chol_eq(ATA,ATb,x,N)
        
        !! Conjugate Gradient Method
        !AT  = transpose(A)
        !ATA = matmul( AT, A )
        !ATb = matmul( AT, b )
        !if(present(x0))then
        !  x0CG = x0
        !else
        !  x0CG = 1.
        !endif
        !call CG(ATA,ATb,x,x0CG,N)
        
      end subroutine qr_solver
      
      subroutine gram_dec(A,Q,R,M,N)
        !---------------------------------subroutine  comment
        !  Version   :  V1.0    
        !  Coded by  :  syz 
        !  Date      :  
        !-----------------------------------------------------
        !  Purpose   :   ���������� Gram-Schmidt�ֽ�������QR�ֽ�
        !    
        !-----------------------------------------------------
        !  Input  parameters  :
        !       1.    Aԭʼ����
        !       2.    A(M,N)
        !  Output parameters  :
        !       1.    �ֽ���Ϊ   Q(M,N):ע�� Q���Ƿ���Q������Ϊ��׼������
        !       2.                 R(N,N)��R�Ƿ���
        !       3.   
        !----------------------------------------------------
        !  Post Script :
        !       1.  ע������ά�����ֽ��Q��������������
        !       2.  ���ڱ�̷������Բο������������Ӧ�á����ʹ����
        !       3.  ��ϸ����ѧ���ͣ����Բο� ��ʡ����ѧԺ��
        !           ���Դ����̲ġ�Linear Algebra with Application��
        !----------------------------------------------------
        
        integer(i_kind)                , intent(in ) :: M,N
        real   (r_kind), dimension(M,N), intent(in ) :: A
        real   (r_kind), dimension(M,N), intent(out) :: Q
        real   (r_kind), dimension(N,N), intent(out) :: R
        
        real   (r_kind), dimension(M  ) :: T(M)
        
        integer(i_kind) :: i,j,k
        
        R(1,1) = sqrt( dot_product( A(:,1), A(:,1) ) )
        Q(:,1) = A(:,1) / R(1,1)
        
        do k = 2,N
          do j = 1,k-1
            R(j,k) = dot_product( Q(:,j), A(:,k) )
          end do
          
          T = A(:,k)
          
          do j = 1,k-1
            T = T - Q(:,j) * R(j,k)
          end do
        
          R(k,k) = sqrt( dot_product( T, T ) )
          Q(:,k) = T / R(k,k)
        end do
        
      end subroutine gram_dec
      
      subroutine chol_eq(A,b,x,N)
        !---------------------------------subroutine  comment
        !  Version   :  V1.0    
        !  Coded by  :  syz 
        !  Date      :  
        !-----------------------------------------------------
        !  Purpose   :  ��cholesky�ֽⷽ���ⷽ��
        !    
        !-----------------------------------------------------
        integer(i_kind)N
        real   (r_kind)::A(N,N),b(N),x(N)
        real   (r_kind)::L(N,N),y(N),LT(N,N)!LT ΪL��ת�þ���
        
        integer(i_kind)i,j
        
        call chol(A,L,N)
        
        call downtri(L,b,y,N)
        
        LT = transpose(L)
        
        call uptri(LT,y,x,N)  !��һ���Ѿ������x
      
      end subroutine chol_eq
      
      subroutine  chol(A,L,N)
        !---------------------------------subroutine  comment
        !  Version   :  V1.0    
        !  Coded by  :  syz 
        !  Date      :  
        !-----------------------------------------------------
        !  Purpose   :  Cholesky�ֽ��ӳ���
        !-----------------------------------------------------  
        integer(i_kind):: N
        real   (r_kind):: A(N,N),L(N,N)
        real   (r_kind):: s
        integer(i_kind):: i,j,k
        
        L = 0
        
        L(1,1) = sqrt(A(1,1))
        
        L(2:,1) = A(2:,1) / L(1,1)
        
        do j=2,N
          s = dot_product(L(j,1:j-1),L(j,1:j-1))
          
          L(j,j) = sqrt( A(j,j) - s )
          
          !ע��i��Χ
          do i=j+1,N
            s = dot_product(L(i,1:j-1),L(j,1:j-1))
            
            L(i,j) = ( A(i,j) - s ) / L(j,j)
          end do
        end do 
      
      end subroutine chol
      
      subroutine uptri(A,b,x,N)
        !---------------------------------subroutine  comment
        !  Version   :  V1.0    
        !  Coded by  :  syz 
        !  Date      :  2010-4-8
        !-----------------------------------------------------
        !  Purpose   :  �����Ƿ�����Ļش�����
        !                 Ax=b
        !-----------------------------------------------------
        !  Input  parameters  :
        !       1.   A(N,N)ϵ������
        !       2.   b(N)������
        !       3.   N����ά��
        !  Output parameters  :
        !       1.  x  ���̵ĸ�
        !       2.
        !  Common parameters  :
        !
        !----------------------------------------------------
        integer(i_kind):: N
        real   (r_kind):: A(N,N),b(N),x(N)
        integer(i_kind):: i,j,k
        
        x(N) = b(N) / A(N,N)
        
        !�ش�����
        do i = n-1,1,-1
          x(i) = ( b(i) - dot_product( A(i,i+1:N), x(i+1:N)) ) / A(i,i)
        end do
      end subroutine uptri
      
      subroutine downtri(A,b,x,N)
        !---------------------------------subroutine  comment
        !  Version   :  V1.0    
        !  Coded by  :  syz 
        !  Date      :  2010-4-9
        !-----------------------------------------------------
        !  Purpose   :  �����Ƿ�����Ļش�����
        !                 Ax=b
        !-----------------------------------------------------
        integer(i_kind):: N
        real   (r_kind):: A(N,N),b(N),x(N)
        integer(i_kind):: i,j,k
        
        x(1) = b(1) / a(1,1)
        
        do k = 2,N
          x(k) = ( b(k) - dot_product( A(k,1:k-1), x(1:k-1)) ) / A(k,k)
        end do
      
      end subroutine downtri
      
      subroutine CG(A,b,x,x0,N)
        !---------------------------------subroutine  comment
        !  Version   :  V1.0    
        !  Coded by  :  syz 
        !  Date      :  
        !-----------------------------------------------------
        !  Purpose   :  �����ݶȷ�(Conjugate Gradient Method)
        !               ���ڼ��㷽�� AX=b
        !-----------------------------------------------------
        !  Input  parameters  :
        !       1.  A,b ���弴  AX=b
        !       2.  x0������ֵ
        !       3.  N ���̵�ά��
        !  Output parameters  :
        !       1. x ���̵Ľ�
        !       2.
        !  Common parameters  :
        !
        !----------------------------------------------------
        
        integer(i_kind)::N
        real   (r_kind)::A(N,N),b(N),x(N),x0(N)
        
        real   (r_kind)::r0(N),r1(N),p0(N),p1(N),x1(N),x2(N)
        
        integer::i,j,k
        
        integer(i_kind), parameter :: IMAX = 200    ! max iter number
        real   (r_kind), parameter :: tol  = 1.e-14 ! error tolerence
        
        real(r_kind) alpha, beta
        real(r_kind) r0r0
        real(r_kind) Ap0 (N)
        
        r0 = b - matmul(A,x0)
        
        p0 = r0
        
        x1 = x0
        
        do k = 1,IMAX
          r0r0  = dot_product( r0, r0 )
          Ap0   = matmul(A,p0)
          alpha = r0r0 / dot_product( p0, Ap0 )
          x2    = x1 + alpha * p0
          
          if( r0r0 < tol )exit
          
          r1    = r0 - alpha * Ap0
          beta  = dot_product(r1,r1) / r0r0
          p1    = r1 + beta * p0
          
          r0    = r1
          p0    = p1
          x1    = x2
        enddo
        
        x = x1
      end subroutine CG
      
    end module qr_solver_mod
    