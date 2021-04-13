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
      print*,'test case 2 Steady-state geostrophically balanced flow is selected'
      call case2(stat(0))
    elseif(case_num==5)then
      print*,''
      print*,'test case 5 Zonal flow over an isolated mountain is selected'
      call case5(stat(0))
    elseif(case_num == 6)then
      print*,''
      print*,'test case 6 Rossby¨CHaurwitz wave is selected'
      call case6(stat(0))
    elseif(case_num == 8)then
      print*,''
      print*,'test case 8 is Barotropic instability selected'
      call case8(stat(0))
    elseif(case_num == 9)then
      print*,''
      print*,'test case 9 is Shock wave selected'
      call case9(stat(0))
    endif
    
    ! Check fields and calculate ghs on cell
    do iPatch = ifs,ife
      do j = jms,jme
        do i = ims,ime
          phi(i,j,iPatch) = stat(0)%q(1,i,j,iPatch) / sqrtGC(i,j,iPatch)
          uc (i,j,iPatch) = stat(0)%q(2,i,j,iPatch) / stat(0)%q(1,i,j,iPatch)
          vc (i,j,iPatch) = stat(0)%q(3,i,j,iPatch) / stat(0)%q(1,i,j,iPatch)
          call contravProjPlane2Sphere(u(i,j,iPatch), v(i,j,iPatch), uc(i,j,iPatch), vc(i,j,iPatch), matrixA(:,:,cc,i,j,iPatch))
        enddo
      enddo
    enddo
    
    print*,''
    print*,'max/min value of phi : ',maxval(phi (ids:ide,jds:jde,ifs:ife)),minval(phi (ids:ide,jds:jde,ifs:ife))
    print*,'max/min value of u   : ',maxval(u   (ids:ide,jds:jde,ifs:ife)),minval(u   (ids:ide,jds:jde,ifs:ife))
    print*,'max/min value of v   : ',maxval(v   (ids:ide,jds:jde,ifs:ife)),minval(v   (ids:ide,jds:jde,ifs:ife))
    print*,'max/min value of ghsC: ',maxval(ghsC(ids:ide,jds:jde,ifs:ife)),minval(ghsC(ids:ide,jds:jde,ifs:ife))
    
    deallocate(phi)
    deallocate(u  )
    deallocate(v  )
    deallocate(uc )
    deallocate(vc )
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
    real(r_kind)    :: alpha = 0! pi/4.
    
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
            phi(iPOC,i,j,iPatch) = gh0 - (radius * Omega * u0 + u0**2 / 2.) * ( -cos(lon(iPOC,i,j,iPatch))*cos(lat(iPOC,i,j,iPatch))*sin(alpha) + sin(lat(iPOC,i,j,iPatch))*cos(alpha) )**2
            u  (iPOC,i,j,iPatch) = u0 * ( cos(lat(iPOC,i,j,iPatch))*cos(alpha) + cos(lon(iPOC,i,j,iPatch))*sin(lat(iPOC,i,j,iPatch))*sin(alpha) )
            v  (iPOC,i,j,iPatch) = -u0 * sin(lon(iPOC,i,j,iPatch)) * sin(alpha)
            
            Coriolis(iPOC,i,j,iPatch) = 2. * Omega * ( -cos(lon(iPOC,i,j,iPatch))*cos(lat(iPOC,i,j,iPatch))*sin(alpha) + sin(lat(iPOC,i,j,iPatch))*cos(alpha) )
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
    
    ghs  = 0
    ghsC = 0
  end subroutine case2
  
  ! Isolated mountain
  subroutine case5(stat)
    type(stat_field), intent(inout) :: stat
    
    real(r_kind),dimension(:,:,:,:), allocatable :: phi
    real(r_kind),dimension(:,:,:,:), allocatable :: u
    real(r_kind),dimension(:,:,:,:), allocatable :: v
    real(r_kind),dimension(:,:,:,:), allocatable :: uc
    real(r_kind),dimension(:,:,:,:), allocatable :: vc
    
    real(r_kind),dimension(nPointsOnCell,ims:ime,jms:jme,ifs:ife) :: longitude
    real(r_kind),dimension(nPointsOnCell,ims:ime,jms:jme,ifs:ife) :: latitude
    real(r_kind),dimension(nPointsOnCell,ims:ime,jms:jme,ifs:ife) :: r
    
    real(r_kind) :: dzsdlon
    real(r_kind) :: dzsdlat
    
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
            
            !u  (iPOC,i,j,iPatch) = 0
            !v  (iPOC,i,j,iPatch) = 0
            !phi(iPOC,i,j,iPatch) = gh0 - ghs(iPOC,i,j,iPatch)
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
          
          ghsC(i,j,iPatch) = cell_quadrature(ghs(cqs:cqe,i,j,iPatch))
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
            u1(iPOC,i,j,iPatch) = cos(lat(iPOC,i,j,iPatch))
            u2(iPOC,i,j,iPatch) = R*cos(lat(iPOC,i,j,iPatch))**(R-1)*sin(lat(iPOC,i,j,iPatch))**2*cos(R*longitude(iPOC,i,j,iPatch))
            u3(iPOC,i,j,iPatch) = cos(lat(iPOC,i,j,iPatch))**(R+1)*cos(R*longitude(iPOC,i,j,iPatch))
            u (iPOC,i,j,iPatch) = radius*omg*(u1(iPOC,i,j,iPatch)+u2(iPOC,i,j,iPatch)-u3(iPOC,i,j,iPatch))
            
            v (iPOC,i,j,iPatch) = -radius*omg*R*cos(lat(iPOC,i,j,iPatch))**(R-1)*sin(lat(iPOC,i,j,iPatch))*sin(R*longitude(iPOC,i,j,iPatch))
            
            AA1 (iPOC,i,j,iPatch) = omg*0.5*(2.*Omega+omg)*cos(lat(iPOC,i,j,iPatch))**2
            Ac  (iPOC,i,j,iPatch) = 0.25*omg**2
            A21 (iPOC,i,j,iPatch) = (R+1.)*cos(lat(iPOC,i,j,iPatch))**(2.*R+2.)
            A22 (iPOC,i,j,iPatch) = (2.*R**2-R-2.)*cos(lat(iPOC,i,j,iPatch))**(2.*R)
            A23 (iPOC,i,j,iPatch) = 2.*R**2*cos(lat(iPOC,i,j,iPatch))**(2.*R-2)
            Ah  (iPOC,i,j,iPatch) = AA1(iPOC,i,j,iPatch)+Ac(iPOC,i,j,iPatch)*(A21(iPOC,i,j,iPatch)+A22(iPOC,i,j,iPatch)-A23(iPOC,i,j,iPatch))
            
            Bc  (iPOC,i,j,iPatch) = 2.*(Omega+omg)*omg/((R+1)*(R+2))*cos(lat(iPOC,i,j,iPatch))**R
            BB1 (iPOC,i,j,iPatch) = R**2+2.*R+2.
            BB2 (iPOC,i,j,iPatch) = (R+1.)**2.*cos(lat(iPOC,i,j,iPatch))**2.
            Bh  (iPOC,i,j,iPatch) = Bc(iPOC,i,j,iPatch)*(BB1(iPOC,i,j,iPatch)-BB2(iPOC,i,j,iPatch))
            
            CC  (iPOC,i,j,iPatch) = 0.25*omg**2*cos(lat(iPOC,i,j,iPatch))**(2.*R)
            CC1 (iPOC,i,j,iPatch) = (R+1.)*cos(lat(iPOC,i,j,iPatch))**2;
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
    
    ghs  = 0
    ghsC = 0
  end subroutine case6
  
  ! Barotropic instability
  subroutine case8(stat)
    type(stat_field), intent(inout) :: stat
    
    real(r_kind),dimension(nPointsOnCell,ims:ime,jms:jme,ifs:ife) :: longitude
    
    real(r_kind),dimension(:,:,:,:), allocatable :: phi
    real(r_kind),dimension(:,:,:,:), allocatable :: u
    real(r_kind),dimension(:,:,:,:), allocatable :: v
    real(r_kind),dimension(:,:,:,:), allocatable :: uc
    real(r_kind),dimension(:,:,:,:), allocatable :: vc
    
    integer(i_kind) :: neval, ierr
    real   (r_kind) :: abserr
    
    real(r_kind), parameter :: gh0 = gravity * 1.0e4
    real(r_kind), parameter :: ghd = gravity * 120
    real(r_kind), parameter :: lat2 = pi / 4.0
    real(r_kind), parameter :: alpha = 1.0 / 3.0
    real(r_kind), parameter :: beta = 1.0 / 15.0
    
    integer :: i,j,iPatch,iPOC
    
    allocate(phi (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(u   (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(v   (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(uc  (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(vc  (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    
    longitude = lon
    where(longitude> pi) longitude = longitude - 2. * pi
    where(longitude<-pi) longitude = 2. * pi + longitude
    
    do iPatch = ifs,ife
      do j = jms,jme
        do i = ims,ime
          do iPOC = 1,nPointsOnCell
            u(iPOC,i,j,iPatch) = u_function(lat(iPOC,i,j,iPatch))
            v(iPOC,i,j,iPatch) = 0.
            call qags(gh_integrand, -0.5*pi, lat(iPOC,i,j,iPatch), 1.0e-15, 1.0e-10, phi(iPOC,i,j,iPatch), abserr, neval, ierr)
            phi(iPOC,i,j,iPatch) = gh0 - phi(iPOC,i,j,iPatch)
            phi(iPOC,i,j,iPatch) = phi(iPOC,i,j,iPatch) + ghd * cos(lat(iPOC,i,j,iPatch)) * exp(-(longitude(iPOC,i,j,iPatch) / alpha)**2) * exp(-((lat2 - lat(iPOC,i,j,iPatch)) / beta)**2)
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
    
    ghs  = 0
    ghsC = 0
  end subroutine case8
  
  ! Shock wave
  subroutine case9(stat)
    type(stat_field), intent(inout) :: stat
    
    real(r_kind),dimension(:,:,:,:), allocatable :: phi
    real(r_kind),dimension(:,:,:,:), allocatable :: u
    real(r_kind),dimension(:,:,:,:), allocatable :: v
    real(r_kind),dimension(:,:,:,:), allocatable :: uc
    real(r_kind),dimension(:,:,:,:), allocatable :: vc

    real(r_kind),dimension(:,:,:,:), allocatable :: longitude
    real(r_kind),dimension(:,:,:,:), allocatable :: latitude
    real(r_kind),dimension(:,:,:,:), allocatable :: r
         
    real(r_kind) :: rr
    real(r_kind) :: labmda_c
    real(r_kind) :: theta_c
    real(r_kind) :: gh0 = 30000.
    
    integer :: i,j,iPatch,iPOC
    
    allocate(phi (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(u   (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(v   (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(uc  (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(vc  (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    
    allocate(longitude (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(latitude  (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    allocate(r         (nPointsOnCell,ims:ime,jms:jme,ifs:ife))
    
    rr       = pi/9.
    labmda_c = 180 * D2R
    theta_c  = 0
    
    longitude = lon
    latitude  = lat
    
    where(longitude<0) longitude = 2. * pi + longitude
    
    do iPatch = ifs, ife
      do j = jms, jme
        do i = ims, ime
          do iPOC = 1,nPointsOnCell
            r(iPOC,i,j,iPatch) = sqrt( (longitude(iPOC,i,j,iPatch)-labmda_c)**2 + (latitude(iPOC,i,j,iPatch)-theta_c)**2 )
            if(r(iPOC,i,j,iPatch)<=rr)phi(iPOC,i,j,iPatch) = 2. * gh0
            if(r(iPOC,i,j,iPatch)> rr)phi(iPOC,i,j,iPatch) = gh0
            !if(     longitude(iPOC,i,j,iPatch) > (180-30) * D2R &
            !  .and. longitude(iPOC,i,j,iPatch) < (180+30) * D2R &
            !  .and. lat      (iPOC,i,j,iPatch) > (0  -30) * D2R & 
            !  .and. lat      (iPOC,i,j,iPatch) < (0  +30) * D2R )then
            !  phi(iPOC,i,j,iPatch) = phi(iPOC,i,j,iPatch) + 30000
            !endif
            
            u  (iPOC,i,j,iPatch) = 0
            v  (iPOC,i,j,iPatch) = 0
            
            Coriolis(iPOC,i,j,iPatch) = 0
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
    
    ghs  = 0
    ghsC = 0
  end subroutine case9
  
  real(r_kind) function gh_integrand(lat) result(res)

    real(r_kind), intent(in) :: lat

    real(r_kind) u, f

    u = u_function(lat)
    f = 2. * omega * sin(lat)
    res = radius * u * (f + tan(lat) / radius * u)

  end function gh_integrand

  real(r_kind) function u_function(lat) result(res)
    real(r_kind), intent(in) :: lat

    real(r_kind), parameter :: lat0  = pi / 7.0
    real(r_kind), parameter :: lat1  = pi / 2.0 - lat0
    real(r_kind), parameter :: u_max = 80.0
    real(r_kind), parameter :: en    = exp(-4.0 / (lat1 - lat0)**2)

    if (lat <= lat0 .or. lat >= lat1) then
      res = 0.
    else
      res = u_max / en * exp(1.0 / (lat - lat0) / (lat - lat1))
    end if

  end function u_function
  
end module test_case_mod
    