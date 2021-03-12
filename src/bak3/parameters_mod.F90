module parameters_mod
  use constants_mod
  implicit none
  
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
  
  integer(i_kind) :: nPointsOnEdge
  integer(i_kind) :: nQuadOrder
  integer(i_kind) :: nQuadPointsOnCell
  integer(i_kind) :: nEdgesOnCell
  integer(i_kind) :: nPointsOnCell
  integer(i_kind) :: maxGhostPointsOnCell
  
  ! dynamic options
  character*200   :: reconstruct_scheme = 'WENO'
  integer(i_kind) :: stencil_width = 5
  integer(i_kind) :: recBdy        = 2
  
  ! Index parameter
  integer(i_kind) :: ids      ! The starting index in the x-direction (Physical domain)
  integer(i_kind) :: ide      ! The ending index in the x-direction  (Physical domain)
  integer(i_kind) :: jds      ! The starting index in the y-direction  (Physical domain)
  integer(i_kind) :: jde      ! The ending index in the y-direction  (Physical domain)
  
  integer(i_kind) :: ims      ! The starting index in the x-direction (Physical cell/element domain)
  integer(i_kind) :: ime      ! The ending index in the x-direction  (Physical cell/element domain)
  integer(i_kind) :: jms      ! The starting index in the y-direction  (Physical cell/element domain)
  integer(i_kind) :: jme      ! The ending index in the y-direction  (Physical cell/element domain)
  
  integer(i_kind) :: ifs      ! The starting index of patch(face)
  integer(i_kind) :: ife      ! The ending index of patch(face)
                      
  integer(i_kind) :: Nx       ! Element numbers in the x-direction
  integer(i_kind) :: Ny       ! Element numbers in the y-direction
  
  integer(i_kind), parameter :: Nf = 6           ! Number of cube faces
  
  real(r_kind)    :: x_min = -45.   !  start location of x-direction
  real(r_kind)    :: x_max =  45.   !  end location of x-direction
  real(r_kind)    :: y_min = -45.   !  start location of y-direction
  real(r_kind)    :: y_max =  45.   !  end location of y-direction
  
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
                           xhalo               ,&
                           yhalo               ,&
                           nPointsOnEdge
  namelist /dynamic_opt/   reconstruct_scheme,&
                           stencil_width
  
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
    xhalo         = 1
    yhalo         = 1
    nPointsOnEdge = 3
    
    run_days         = 1
    run_hours        = 0
    run_minutes      = 0
    run_seconds      = 0
    history_interval = 360
    integral_scheme  = 'RK4'
    
    ! Read namelist
    call readNamelist
    
    ! Check if dx and dy are avaliable to achieve the integral element number
    if( (x_max - x_min)/dx - int((x_max - x_min)/dx) /= 0. )then
      stop '90 divide dx must be integer(i_kind), choose another dx'
    end if
    
    if( (y_max - y_min)/dy - int((y_max - y_min)/dy) /= 0. )then
        stop '90 divide dy must be integer(i_kind), choose another dy'
    end if
    
    ! Calculate element numbers on x/y direction
    Nx = int((x_max - x_min)/dx)
    Ny = int((y_max - y_min)/dy)
    
    ! Calculate starting and ending index for physical domain
    ids  = 1
    ide  = nx
    jds  = 1
    jde  = ny
    
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
    
    nQuadOrder           = nPointsOnEdge
    nEdgesOnCell         = 4
    nQuadPointsOnCell    = nQuadOrder * nEdgesOnCell
    maxGhostPointsOnCell = 2 * nQuadPointsOnCell
    ! nPointsOnCell = center point + corner points + boundary points + triangle quadrature points + max ghost points
    nPointsOnCell        = 1 + nEdgesOnCell + nEdgesOnCell * nPointsOnEdge + nQuadPointsOnCell + maxGhostPointsOnCell
    
    recBdy = ( stencil_width - 1 ) / 2
    
  end subroutine init_Parameters
  
end module parameters_mod
    