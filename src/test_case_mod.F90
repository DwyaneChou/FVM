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
    
    allocate(phi(ims:ime,jms:jme,ifs:ife))
    allocate(u  (ims:ime,jms:jme,ifs:ife))
    allocate(v  (ims:ime,jms:jme,ifs:ife))
    
    if(case_num == 2)then
      print*,''
      print*,'test case 2 is selected'
      call case2(stat(0))
    elseif(case_num==5)then
      print*,''
      print*,'test case 5 is selected'
      call case5(stat(0))
    elseif(case_num == 6)then
      print*,''
      print*,'test case 6 is selected'
      call case6(stat(0))
    endif
    
    ! Check fields and calculate zs on cell
    do iPatch = ifs,ife
      do j = jms,jme
        do i = ims,ime
          phi(i,j,iPatch) = stat(0)%q(1,i,j,iPatch) / sqrtGC(i,j,iPatch)
          u  (i,j,iPatch) = stat(0)%q(2,i,j,iPatch) / stat(0)%q(1,i,j,iPatch)
          v  (i,j,iPatch) = stat(0)%q(3,i,j,iPatch) / stat(0)%q(1,i,j,iPatch)
          call contravProjPlane2Sphere(u(i,j,iPatch), v(i,j,iPatch), u(i,j,iPatch), v(i,j,iPatch), matrixA(:,:,cc,i,j,iPatch))
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
  end subroutine init_testcase

  ! Global steady flow
  subroutine case2(stat)
    type(stat_field), intent(inout) :: stat
    
    real(r_kind),dimension(:,:,:,:), allocatable :: phi
    real(r_kind),dimension(:,:,:,:), allocatable :: u
    real(r_kind),dimension(:,:,:,:), allocatable :: v
    real(r_kind),dimension(:,:,:,:), allocatable :: uc
    real(r_kind),dimension(:,:,:,:), allocatable :: vc
    
    real(r_kind)    :: u0
    real(r_kind)    :: gh0 = 29400.
    real(r_kind)    :: gh
    
    integer :: i,j,iPatch,iPOC
    
    allocate(phi (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(u   (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(v   (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(uc  (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(vc  (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    
    u0 = 2. * pi * radius / (12. * 86400.)
    
    do iPatch = ifs, ife
      do j = jms, jme
        do i = ims, ime
          do iPOC = 1,nPointsOnCell
            phi(iPOC,i,j,iPatch) = gh0 - (radius * Omega * u0 + u0**2 / 2.) * sinlat(iPOC,i,j,iPatch)**2
            u  (iPOC,i,j,iPatch) = u0 * coslat(iPOC,i,j,iPatch)
          enddo
        enddo
      enddo
    enddo
    
    v = 0.
    
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
    
  end subroutine case2
  
  ! Isolated mountain
  subroutine case5(stat)
    type(stat_field), intent(inout) :: stat
    
    real(r_kind),dimension(:,:,:,:), allocatable :: phi
    real(r_kind),dimension(:,:,:,:), allocatable :: u
    real(r_kind),dimension(:,:,:,:), allocatable :: v
    real(r_kind),dimension(:,:,:,:), allocatable :: uc
    real(r_kind),dimension(:,:,:,:), allocatable :: vc
    
    real(r_kind),dimension(nPointsOnCell,ims:ime,jms:jme,ifs:ife) :: ghs
    real(r_kind),dimension(nPointsOnCell,ims:ime,jms:jme,ifs:ife) :: longitude
    real(r_kind),dimension(nPointsOnCell,ims:ime,jms:jme,ifs:ife) :: latitude
    real(r_kind),dimension(nPointsOnCell,ims:ime,jms:jme,ifs:ife) :: r
    
    real(r_kind),parameter :: hs0      = 2000.
    real(r_kind),parameter :: u0       = 20.
    real(r_kind),parameter :: alpha    = 0.
    real(r_kind),parameter :: gh0      = 5960.*gravity
                                    
    real(r_kind) :: rr
    real(r_kind) :: labmda_c
    real(r_kind) :: theta_c
    
    integer :: i,j,iPatch,iPOC
    
    allocate(phi (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(u   (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(v   (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(uc  (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(vc  (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    
    rr       = pi/9.
    labmda_c = 1.5*pi
    theta_c  = pi/6.
    
    longitude = lon
    latitude  = lat
    
    where(longitude<0) longitude = 2. * pi + longitude
    
    ghs = 0.
    
    do iPatch = ifs, ife
      do j = jms, jme
        do i = ims, ime
          do iPOC = 1,nPointsOnCell
            r  (iPOC,i,j,iPatch) = sqrt(min(rr**2,(longitude(iPOC,i,j,iPatch)-labmda_c)**2+(latitude(iPOC,i,j,iPatch)-theta_c)**2))
            ghs(iPOC,i,j,iPatch) = gravity*hs0*(1.-r(iPOC,i,j,iPatch)/rr)
            u  (iPOC,i,j,iPatch) = u0*(cos(latitude(iPOC,i,j,iPatch))*cos(alpha)+cos(longitude(iPOC,i,j,iPatch))*sin(latitude(iPOC,i,j,iPatch))*sin(alpha))
            v  (iPOC,i,j,iPatch) = -u0*sin(longitude(iPOC,i,j,iPatch))*sin(alpha)
            phi(iPOC,i,j,iPatch) = gh0 - (radius*Omega*u0 + u0**2/2.d0)*(-cos(longitude(iPOC,i,j,iPatch))*cos(latitude(iPOC,i,j,iPatch))*sin(alpha) + sin(latitude(iPOC,i,j,iPatch))*cos(alpha))**2 - ghs(iPOC,i,j,iPatch)
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
    
    zs = ghs / gravity
    
  end subroutine case5
  
  ! Rossby-Haurwitz wave with wavenumber 4
  subroutine case6(stat)
    type(stat_field), intent(inout) :: stat
    
    real(r_kind),dimension(:,:,:,:), allocatable :: phi
    real(r_kind),dimension(:,:,:,:), allocatable :: u
    real(r_kind),dimension(:,:,:,:), allocatable :: v
    real(r_kind),dimension(:,:,:,:), allocatable :: uc
    real(r_kind),dimension(:,:,:,:), allocatable :: vc
    
    real(r_kind),dimension(nPointsOnCell,ims:ime,jms:jme,ifs:ife) :: u1,u2,u3                ! working array
    real(r_kind),dimension(nPointsOnCell,ims:ime,jms:jme,ifs:ife) :: AA1,Ac,A21,A22,A23,Ah   ! working array
    real(r_kind),dimension(nPointsOnCell,ims:ime,jms:jme,ifs:ife) :: Bc,BB1,BB2,Bh           ! working array
    real(r_kind),dimension(nPointsOnCell,ims:ime,jms:jme,ifs:ife) :: CC,CC1,CC2,Ch           ! working array
    !real(r_kind),dimension(nPointsOnCell,ims:ime,jms:jme,ifs:ife) :: sinlat                  ! working array
    !real(r_kind),dimension(nPointsOnCell,ims:ime,jms:jme,ifs:ife) :: coslat                  ! working array
    real(r_kind),dimension(nPointsOnCell,ims:ime,jms:jme,ifs:ife) :: longitude               ! working array
    
    real(r_kind),parameter :: omg  = 7.848d-6       ! angular velocity of RH wave
    real(r_kind),parameter :: R    = 4              ! wave number of RH wave
    real(r_kind),parameter :: h0   = 8000.          ! wave number of RH wave
    
    integer :: i,j,iPatch,iPOC
    
    allocate(phi (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(u   (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(v   (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(uc  (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(vc  (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    
    !sinlat    = sin(lat)
    !coslat    = cos(lat)
    longitude = lon
    
    do iPatch = ifs, ife
      do j = jms, jme
        do i = ims, ime
          do iPOC = 1,nPointsOnCell
            u1(iPOC,i,j,iPatch) = coslat(iPOC,i,j,iPatch)
            u2(iPOC,i,j,iPatch) = R*coslat(iPOC,i,j,iPatch)**(R-1)*sinlat(iPOC,i,j,iPatch)**2*cos(R*longitude(iPOC,i,j,iPatch))
            u3(iPOC,i,j,iPatch) = coslat(iPOC,i,j,iPatch)**(R+1)*cos(R*longitude(iPOC,i,j,iPatch))
            u (iPOC,i,j,iPatch) = radius*omg*(u1(iPOC,i,j,iPatch)+u2(iPOC,i,j,iPatch)-u3(iPOC,i,j,iPatch))
            
            v (iPOC,i,j,iPatch) = -radius*omg*R*coslat(iPOC,i,j,iPatch)**(R-1)*sinlat(iPOC,i,j,iPatch)*sin(R*longitude(iPOC,i,j,iPatch))
            
            AA1 (iPOC,i,j,iPatch) = omg*0.5*(2.*Omega+omg)*coslat(iPOC,i,j,iPatch)**2
            Ac  (iPOC,i,j,iPatch) = 0.25*omg**2
            A21 (iPOC,i,j,iPatch) = (R+1.)*coslat(iPOC,i,j,iPatch)**(2.*R+2.)
            A22 (iPOC,i,j,iPatch) = (2.*R**2-R-2.)*coslat(iPOC,i,j,iPatch)**(2.*R)
            A23 (iPOC,i,j,iPatch) = 2.*R**2*coslat(iPOC,i,j,iPatch)**(2.*R-2)
            Ah  (iPOC,i,j,iPatch) = AA1(iPOC,i,j,iPatch)+Ac(iPOC,i,j,iPatch)*(A21(iPOC,i,j,iPatch)+A22(iPOC,i,j,iPatch)-A23(iPOC,i,j,iPatch))
            
            Bc  (iPOC,i,j,iPatch) = 2.*(Omega+omg)*omg/((R+1)*(R+2))*coslat(iPOC,i,j,iPatch)**R
            BB1 (iPOC,i,j,iPatch) = R**2+2.*R+2.
            BB2 (iPOC,i,j,iPatch) = (R+1.)**2.*coslat(iPOC,i,j,iPatch)**2.
            Bh  (iPOC,i,j,iPatch) = Bc(iPOC,i,j,iPatch)*(BB1(iPOC,i,j,iPatch)-BB2(iPOC,i,j,iPatch))
            
            CC  (iPOC,i,j,iPatch) = 0.25*omg**2*coslat(iPOC,i,j,iPatch)**(2.*R)
            CC1 (iPOC,i,j,iPatch) = (R+1.)*coslat(iPOC,i,j,iPatch)**2;
            CC2 (iPOC,i,j,iPatch) = R+2.
            Ch  (iPOC,i,j,iPatch) = CC(iPOC,i,j,iPatch)*(CC1(iPOC,i,j,iPatch)-CC2(iPOC,i,j,iPatch))
            
            phi (iPOC,i,j,iPatch) = gravity*h0+radius**2*(Ah(iPOC,i,j,iPatch) + Bh(iPOC,i,j,iPatch)*cos(R*longitude(iPOC,i,j,iPatch)) + Ch(iPOC,i,j,iPatch)*cos(2.0*R*longitude(iPOC,i,j,iPatch)))
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
  end subroutine case6
  
end module test_case_mod
    