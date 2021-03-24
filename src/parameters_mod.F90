module parameters_mod
  use constants_mod
  implicit none
  
  integer(i_kind) :: nVar = 3 ! Just for Intel fortran 2019, this compiler cannot get correct nVar with parameter property
  
  ! Namelist parameters
  ! time_settings
  integer(i_kind) :: run_days
  integer(i_kind) :: run_hours
  integer(i_kind) :: run_minutes
  integer(i_kind) :: run_seconds
  real(r_kind)    :: dt               ! time step
  real(i_kind)    :: history_interval ! output interval in seconds
  
  character*200 :: integral_scheme
  
  ! Case select
  integer(i_kind) :: case_num
  
  ! Domain
  real(r_kind)    :: dx     !  grid-spacing in the x-direction
  real(r_kind)    :: dy     !  grid-spacing in the y-direction
  
  integer(i_kind) :: xhalo  !  halo number of x-diretion
  integer(i_kind) :: yhalo  !  halo number of y-diretion
  
  integer(i_kind) :: nEdgesOnCell
  integer(i_kind) :: nPointsOnEdge
  integer(i_kind) :: nQuadPointsOnCell
  integer(i_kind) :: nPointsOnCell
  integer(i_kind) :: maxGhostPointsOnCell
  
  integer(i_kind),dimension(:,:,:),allocatable :: nGhostPointsOnCell
  
  ! dynamic options
  integer(i_kind) :: nQuadOrder
  integer(i_kind) :: quad_opt
  character*200   :: reconstruct_scheme = 'WENO'
  integer(i_kind) :: stencil_width = 5
  integer(i_kind) :: recBdy        = 2
  real   (r_kind) :: epsilon       = 1.e-2 ! Coefficient for avioding divide 0 in WLS-ENO, 1.e-2 for dx>0.25 degree, 1.e+2 for dx<=0.25 degree
                                           ! or in other words, increasing by finer grid
  ! Index parameter
  integer(i_kind) :: ids      ! The starting index in the x-direction (Physical domain)
  integer(i_kind) :: ide      ! The ending index in the x-direction  (Physical domain)
  integer(i_kind) :: jds      ! The starting index in the y-direction  (Physical domain)
  integer(i_kind) :: jde      ! The ending index in the y-direction  (Physical domain)
  
  integer(i_kind) :: idsm1
  integer(i_kind) :: idep1
  integer(i_kind) :: jdsm1
  integer(i_kind) :: jdep1
  
  integer(i_kind) :: ims      ! The starting index in the x-direction (Physical cell/element domain)
  integer(i_kind) :: ime      ! The ending index in the x-direction  (Physical cell/element domain)
  integer(i_kind) :: jms      ! The starting index in the y-direction  (Physical cell/element domain)
  integer(i_kind) :: jme      ! The ending index in the y-direction  (Physical cell/element domain)
  
  integer(i_kind) :: cc       ! center point index on cell
  integer(i_kind) :: ccs      ! corner points indices start
  integer(i_kind) :: cce      ! corner points indices end
  integer(i_kind) :: cbs      ! boundary points indices start
  integer(i_kind) :: cbe      ! boundary points indices end
  integer(i_kind) :: cqs      ! Gaussian quadrature points indices start
  integer(i_kind) :: cqe      ! Gaussian quadrature points indices end
  integer(i_kind) :: cgs      ! ghost points indices start
  integer(i_kind) :: cge      ! ghost points indices end
  
  integer(i_kind) :: ifs      ! The starting index of patch(face)
  integer(i_kind) :: ife      ! The ending index of patch(face)
                      
  integer(i_kind) :: Nx       ! Element numbers in the x-direction
  integer(i_kind) :: Ny       ! Element numbers in the y-direction
  
  integer(i_kind) :: Nx_halo  ! Element numbers in the x-direction with halo
  integer(i_kind) :: Ny_halo  ! Element numbers in the y-directionwith halo
  
  integer(i_kind) :: nIntegralSubSteps ! number of integral substeps in temporal integration scheme
  integer(i_kind) :: nsteps            ! total integral steps
  
  integer(i_kind), parameter :: Nf = 6           ! Number of cube faces
  
  real(r_kind)    :: x_min = -45.   !  start location of x-direction
  real(r_kind)    :: x_max =  45.   !  end location of x-direction
  real(r_kind)    :: y_min = -45.   !  start location of y-direction
  real(r_kind)    :: y_max =  45.   !  end location of y-direction
  
  ! Model run time control variables
  integer(i8    ) :: total_run_time   ! total run time for this model in seconds, this variable is determined by run_days, run_hours ...
  integer(i_kind) :: total_run_steps  ! total run steps for this model in seconds, this variable is determined by total_run_time and dt
  
  namelist /time_settings/ dt               ,&
                           run_days         ,&
                           run_hours        ,&
                           run_minutes      ,&
                           run_seconds      ,&
                           history_interval ,&
                           integral_scheme
  
  namelist /case_select/   case_num
  
  namelist /domain/        dx                  ,&
                           dy                  ,&
                           nPointsOnEdge
  namelist /dynamic_opt/   reconstruct_scheme,&
                           quad_opt          ,&
                           nQuadOrder        ,&
                           stencil_width     ,&
                           epsilon
  
  contains
  
  subroutine readNamelist
    
    open(1, file = 'namelist.input',status='old')
    read(1, nml  = time_settings)
    read(1, nml  = case_select  )
    read(1, nml  = domain       )
    read(1, nml  = dynamic_opt  )
    close(1)
    
  end subroutine readNamelist
  
  subroutine init_Parameters
    
    ! Setting default values
    dt            = 300.
    dx            = 2.
    dy            = 2.
    nPointsOnEdge = 3
    
    run_days         = 1
    run_hours        = 0
    run_minutes      = 0
    run_seconds      = 0
    history_interval = 360
    integral_scheme  = 'RK4'
    
    ! Read namelist
    call readNamelist
    
    if(trim(reconstruct_scheme)=='WENO')then
      if(stencil_width/=5)then
        print*,'stencil_width is not 5, during using WENO, stencil_width has been reset to 5'
        stencil_width = 5
      endif
      if(nPointsOnEdge/=1)then
        print*,'nPointsOnEdge is not 1, during using WENO, nPointsOnEdge has been reset to 1'
        nPointsOnEdge = 1
      endif
    endif
    
    ! Calculate total run time in seconds
    total_run_time = run_days * 86400 + run_hours * 3600 + run_minutes * 60 + run_seconds
    
    ! Calculate total run steps
    total_run_steps = ceiling(total_run_time/dt)
    
    ! Check if dx and dy are avaliable to achieve the integral element number
    if( (x_max - x_min)/dx - int((x_max - x_min)/dx) /= 0. )then
      stop '90 divide dx must be integer(i_kind), choose another dx'
    end if
    
    if( (y_max - y_min)/dy - int((y_max - y_min)/dy) /= 0. )then
        stop '90 divide dy must be integer(i_kind), choose another dy'
    end if
    
    ! Calculate element numbers on x/y direction
    Nx = nint((x_max - x_min)/dx)
    Ny = nint((y_max - y_min)/dy)
    
    recBdy = ( stencil_width - 1 ) / 2
    xhalo  = recBdy
    yhalo  = recBdy
    
    ! Calculate starting and ending index for physical domain
    ids  = 1
    ide  = nx
    jds  = 1
    jde  = ny
    
    idsm1 = ids - 1
    idep1 = ide + 1
    jdsm1 = jds - 1
    jdep1 = jde + 1
    
    ims  = 1  - xhalo
    ime  = Nx + xhalo
    jms  = 1  - yhalo
    jme  = Ny + yhalo
    
    ! Setting the starting patch index and ending patch index
    ifs = 1
    ife = Nf
    
    ! Convert Degree to Radian
    dx    = dx    * D2R
    dy    = dy    * D2R
    x_min = x_min * D2R
    x_max = x_max * D2R
    y_min = y_min * D2R
    y_max = y_max * D2R
    
    ! Setting the number of substeps in temporal integration scheme
    if(trim(adjustl(integral_scheme)) == 'RK3_TVD')then
      nIntegralSubSteps = 3
    elseif(trim(adjustl(integral_scheme)) == 'RK4')then
      nIntegralSubSteps = 4
    else
      stop 'Unknown integral scheme, please select from RK3_TVD or RK4 ...'
    endif
    
    nsteps = total_run_time / dt
    
    nEdgesOnCell         = 4
    
    if(quad_opt==1)then
      nQuadPointsOnCell = nPointsOnEdge**2
    elseif(quad_opt==2)then
      nQuadPointsOnCell = nQuadOrder * nEdgesOnCell
    endif
    
    maxGhostPointsOnCell = 4 * nQuadPointsOnCell
    ! nPointsOnCell = center point + corner points + boundary points + triangle quadrature points + max ghost points
    nPointsOnCell        = 1 + nEdgesOnCell + nEdgesOnCell * nPointsOnEdge + nQuadPointsOnCell + maxGhostPointsOnCell
    
    cc  = 1
    ccs = cc + 1
    cce = ccs + nEdgesOnCell - 1
    cbs = cce + 1
    cbe = cbs + nEdgesOnCell * nPointsOnEdge - 1
    cqs = cbe + 1
    cqe = cqs + nQuadPointsOnCell- 1
    cgs = cqe + 1
    cge = cgs + maxGhostPointsOnCell - 1
    
    allocate( nGhostPointsOnCell(ids:ide,jds:jde,ifs:ife) )
    nGhostPointsOnCell = 0
    
  end subroutine init_Parameters
  
end module parameters_mod
    