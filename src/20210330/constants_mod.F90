MODULE constants_mod
    
  implicit none
  
  integer  ,parameter :: i2  = 2
  integer  ,parameter :: i4  = 4
  integer  ,parameter :: i8  = 8
  integer  ,parameter :: i16 = 16
  integer  ,parameter :: r2  = 2
  integer  ,parameter :: r4  = 4
  integer  ,parameter :: r8  = 8
  integer  ,parameter :: r16 = 16
  
  integer  ,parameter :: i_kind = i4
  integer  ,parameter :: r_kind = r8
  
  real(r_kind),parameter    :: gravity   = 9.80616
  real(r_kind),parameter    :: pi        = 2.*asin(1.)
  real(r_kind),parameter    :: radius    = 6371220.
  real(r_kind),parameter    :: D2R       = PI/180.    ! convert degree into radian
  real(r_kind),parameter    :: R2D       = 180./PI    ! convert radian into degree
  real(r_kind),parameter    :: Omega     = 7.292E-5
  real(r_kind),parameter    :: FillValue = -9999999999999999.
  real(r_kind),parameter    :: Inf       = huge(Inf)
END MODULE constants_mod
