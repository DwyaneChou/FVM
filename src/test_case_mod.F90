module test_case_mod
  use parameters_mod
  use mesh_mod
  use stat_mod
  use projection_mod
  implicit none
  
  contains
    
  subroutine init_testcase
    integer i,j,iPatch,iPOC

    real(r_kind), dimension(:,:,:), allocatable :: phi
    real(r_kind), dimension(:,:,:), allocatable :: u
    real(r_kind), dimension(:,:,:), allocatable :: v
    real(r_kind), dimension(:,:,:), allocatable :: uc
    real(r_kind), dimension(:,:,:), allocatable :: vc
    
    allocate(phi(ims:ime,jms:jme,ifs:ife))
    allocate(u  (ims:ime,jms:jme,ifs:ife))
    allocate(v  (ims:ime,jms:jme,ifs:ife))
    allocate(uc (ims:ime,jms:jme,ifs:ife))
    allocate(vc (ims:ime,jms:jme,ifs:ife))
    
    if(case_num == 2)then
      print*,''
      print*,'test case 2 is selected'
      call gaussian_hill(stat(0))
    endif
    
    ! Check fields and calculate zs on cell
    do iPatch = ifs,ife
      do j = jms,jme
        do i = ims,ime
          phi(i,j,iPatch) = stat(0)%q(1,i,j,iPatch) / sqrtGC(i,j,iPatch)
          uc (i,j,iPatch) = stat(0)%q(2,i,j,iPatch) / stat(0)%q(1,i,j,iPatch)
          vc (i,j,iPatch) = stat(0)%q(3,i,j,iPatch) / stat(0)%q(1,i,j,iPatch)
          call contravProjPlane2Sphere(u(i,j,iPatch), v(i,j,iPatch), uc(i,j,iPatch), vc(i,j,iPatch), matrixA(:,:,cc,i,j,iPatch))
          zsc(i,j,iPatch) = cell_quadrature(zs(cqs:cqe,i,j,iPatch))
        enddo
      enddo
    enddo
    
    print*,''
    print*,'max/min value of phi: ',maxval(phi(ids:ide,jds:jde,ifs:ife)),minval(phi(ids:ide,jds:jde,ifs:ife))
    print*,'max/min value of u  : ',maxval(u  (ids:ide,jds:jde,ifs:ife)),minval(u  (ids:ide,jds:jde,ifs:ife))
    print*,'max/min value of v  : ',maxval(v  (ids:ide,jds:jde,ifs:ife)),minval(v  (ids:ide,jds:jde,ifs:ife))
    print*,'max/min value of zsc: ',maxval(zsc(ids:ide,jds:jde,ifs:ife)),minval(zsc(ids:ide,jds:jde,ifs:ife))
    
    deallocate(phi)
    deallocate(u  )
    deallocate(v  )
    deallocate(uc )
    deallocate(vc )
  end subroutine init_testcase
  
  ! Global steady flow
  subroutine gaussian_hill(stat)
    type(stat_field), intent(inout) :: stat
    
    real(r_kind),dimension(:,:,:,:), allocatable :: phi
    real(r_kind),dimension(:,:,:,:), allocatable :: u
    real(r_kind),dimension(:,:,:,:), allocatable :: v
    real(r_kind),dimension(:,:,:,:), allocatable :: uc
    real(r_kind),dimension(:,:,:,:), allocatable :: vc
    
    real    :: u0
    real    :: alpha0
    real    :: b0
    real    :: lambda_c
    real    :: theta_c
    real    :: gh0
    
    real    :: X,Y,Z
    real    :: Xc,Yc,Zc
    
    integer :: i,j,iPatch,iPOC
    
    allocate(phi (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(u   (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(v   (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(uc  (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(vc  (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    
    gh0      = 100. * gravity
    alpha0   = 0.
    b0       = 0.5 * 10.e-14
    lambda_c = 0. !3. * pi / 2.
    theta_c  = 0.
    u0       = 2. * pi * radius / (12. * 86400.)
    Xc       = radius * cos(lambda_c) * cos(theta_c)
    Yc       = radius * sin(lambda_c) * cos(theta_c)
    Zc       = radius * sin(theta_c)
    
    do iPatch = ifs, ife
      do j = jms, jme
        do i = ims, ime
          do iPOC = 1,nPointsOnCell
            X = radius * cos(lon(iPOC,i,j,iPatch)) * cos(lat(iPOC,i,j,iPatch))
            Y = radius * sin(lon(iPOC,i,j,iPatch)) * cos(lat(iPOC,i,j,iPatch))
            Z = radius * sin(lat(iPOC,i,j,iPatch))
            
            phi(iPOC,i,j,iPatch) = gh0 * exp(-b0 * ((X-Xc)**2 + (Y-Yc)**2 + (Z-Zc)**2))
            u  (iPOC,i,j,iPatch) = u0 * (cos(lat(iPOC,i,j,iPatch)) * cos(alpha0) + sin(alpha0) * cos(lon(iPOC,i,j,iPatch)) * sin(lat(iPOC,i,j,iPatch)))
            v  (iPOC,i,j,iPatch) = -u0 * sin(alpha0) * sin(lon(iPOC,i,j,iPatch))
          enddo
        enddo
      enddo
    enddo
    
    do iPatch = ifs, ife
      do j = jms, jme
        do i = ims, ime
          do iPOC = 1,nPointsOnCell
            call contravProjSphere2Plane(uc(iPOC,i,j,iPatch), vc(iPOC,i,j,iPatch), u(iPOC,i,j,iPatch), v(iPOC,i,j,iPatch), matrixIA(:,:,iPOC,i,j,iPatch))
          enddo
        enddo
      enddo
    enddo
    
    phi = sqrtG * phi
    uc  = phi * uc
    vc  = phi * vc
    
    do iPatch = ifs, ife
      do j = jms, jme
        do i = ims, ime
          stat%q(1,i,j,iPatch) = cell_quadrature(phi(cqs:cqe,i,j,iPatch))
          stat%q(2,i,j,iPatch) = cell_quadrature(uc (cqs:cqe,i,j,iPatch))
          stat%q(3,i,j,iPatch) = cell_quadrature(vc (cqs:cqe,i,j,iPatch))
        enddo
      enddo
    enddo
    
    deallocate(phi )
    deallocate(u   )
    deallocate(v   )
    deallocate(uc  )
    deallocate(vc  )
    
    zs = 0
  end subroutine gaussian_hill
  
end module test_case_mod
    