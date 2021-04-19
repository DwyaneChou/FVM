module temporal_mod
  use constants_mod
  use parameters_mod
  use mesh_mod
  use stat_mod
  use tend_mod
  use spatial_operators_mod
  use io_mod
  implicit none
  
  private
  
  public RK4, RK3_TVD, PC2, SSPRK
  
  integer, parameter :: k1 = -1
  integer, parameter :: k2 = -2
  integer, parameter :: k3 = -3
  integer, parameter :: k4 = -4
  integer, parameter :: new = 1
  integer, parameter :: old = 0
  
    contains

    subroutine RK4(stat_new,stat_old)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(inout) :: stat_old
      
      real(r_kind),dimension(4),parameter :: RK4_weight = [1./6.,1./3.,1./3.,1./6.]
      
      !call copyStat(stat(k1),stat_old)
      
      call spatial_operator (stat_old, tend(k1))
      call update_stat      (stat(k2), stat_old, tend(k1), 0.5_r_kind * dt)
      
      call spatial_operator (stat(k2), tend(k2))
      call update_stat      (stat(k3), stat_old, tend(k2), 0.5_r_kind * dt)
      
      call spatial_operator (stat(k3), tend(k3))
      call update_stat      (stat(k4), stat_old, tend(k3), dt)
      
      call spatial_operator(stat(k4), tend(k4))
            
      tend(new)%q = RK4_weight(1) * tend(k1)%q + RK4_weight(2) * tend(k2)%q + RK4_weight(3) * tend(k3)%q + RK4_weight(4) * tend(k4)%q
      
      call update_stat (stat_new, stat_old, tend(new), dt)
    end subroutine RK4
    
    subroutine RK3_TVD(stat_new,stat_old)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(inout) :: stat_old
      
      !call copyStat(stat(k1),stat_old)
      
      call spatial_operator (stat_old, tend(k1))
      call update_stat      (stat(k2), stat_old, tend(k1), dt)
      
      !call history_write_stat(stat(k2),2)
      stop 'RK3_TVD'
      
      call spatial_operator (stat(k2), tend(k2))
      call update_stat_RK3_TVD_1(stat(k3), stat_old, stat(k2), tend(k2))
      
      !call history_write_stat(stat(k3),2)
      !stop 'RK3_TVD'
      
      call spatial_operator (stat(k3), tend(k3))
      call update_stat_RK3_TVD_2(stat_new, stat_old, stat(k3), tend(k3))
      
      !call history_write_stat(stat_new,2)
      !stop 'RK3_TVD'

    end subroutine RK3_TVD
    
    subroutine SSPRK(stat_new,stat_old)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(inout) :: stat_old
    
      call spatial_operator (stat_old, tend(k1))
      call update_stat      (stat(k1), stat_old, tend(k1), 0.5_r_kind * dt)
      
      call spatial_operator (stat(k1), tend(k2))
      call update_stat      (stat(k2), stat_old, tend(k2), 0.5_r_kind * dt)
      
      call spatial_operator (stat(k2), tend(k3))
      call update_stat_SSPRK(stat(k3), stat_old, stat(k2), tend(k3))
      
      call spatial_operator (stat(k3), tend(k4))
      call update_stat      (stat_new, stat(k3), tend(k4), 0.5_r_kind * dt)
      
    end subroutine SSPRK
    
    subroutine PC2(stat_new,stat_old)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(inout) :: stat_old
      
      integer :: iPoint
      
      !call copyStat(stat(k1),stat_old)
      
      call spatial_operator (stat_old, tend(k1))
      call update_stat      (stat(k2), stat_old, tend(k1), 0.5_r_kind * dt)
      
      call spatial_operator (stat(k2), tend(k2))
      call update_stat      (stat(k3), stat_old, tend(k2), 0.5_r_kind * dt)
      
      call spatial_operator(stat(k3), tend(k3))
      
      call update_stat (stat_new, stat_old, tend(k3), dt)
    end subroutine PC2
    
    subroutine update_stat(stat_new, stat_old, tend, inc_t)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old
      type(tend_field), intent(in   ) :: tend
      real(r_kind)    , intent(in   ) :: inc_t
      
      integer(i_kind) :: iVar,i,j,iPatch

      !$OMP PARALLEL DO PRIVATE(i,j,iVar) COLLAPSE(4)
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            do iVar = 1,nVar
              stat_new%q(iVar,i,j,iPatch) = stat_old%q(iVar,i,j,iPatch) + inc_t * tend%q(iVar,i,j,iPatch)
            enddo
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO

    end subroutine update_stat
    
    subroutine update_stat_RK3_TVD_1(stat_new, stat_old,stat1, tend)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old
      type(stat_field), intent(in   ) :: stat1
      type(tend_field), intent(in   ) :: tend

      integer(i_kind) :: iVar,i,j,iPatch
      
      !$OMP PARALLEL DO PRIVATE(i,j,iVar) COLLAPSE(4)
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            do iVar = 1,nVar
              stat_new%q(iVar,i,j,iPatch) = 0.75 * stat_old%q(iVar,i,j,iPatch) + 0.25 * stat1%q(iVar,i,j,iPatch) + 0.25 * dt * tend%q(iVar,i,j,iPatch)
            enddo
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine update_stat_RK3_TVD_1
    
    subroutine update_stat_RK3_TVD_2(stat_new, stat_old,stat2, tend)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old
      type(stat_field), intent(in   ) :: stat2
      type(tend_field), intent(in   ) :: tend
      
      real(r_kind),dimension(3),parameter :: weight = [1./3., 2./3., 2./3.]
      
      integer(i_kind) :: iVar,i,j,iPatch
      
      !$OMP PARALLEL DO PRIVATE(i,j,iVar) COLLAPSE(4)
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            do iVar = 1,nVar
              stat_new%q(iVar,i,j,iPatch) = weight(1) * stat_old%q(iVar,i,j,iPatch) + weight(2) * stat2%q(iVar,i,j,iPatch) + weight(3) * dt * tend%q(iVar,i,j,iPatch)
            enddo
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine update_stat_RK3_TVD_2
    
    subroutine update_stat_SSPRK(stat_new, stat_old, stat2, tend3)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old
      type(stat_field), intent(in   ) :: stat2
      type(tend_field), intent(in   ) :: tend3
      
      real(r_kind),dimension(3),parameter :: coef = [2./3.,1./3.,1./6.]
      
      stat_new%q = coef(1) * stat_old%q + coef(2) * stat2%q + coef(3) * dt * tend3%q
      
    end subroutine update_stat_SSPRK
    
    subroutine calc_residual(residual, tend_old, tend_new)
      real(r_kind)    , intent(out) :: residual
      type(tend_field), intent(in ) :: tend_old
      type(tend_field), intent(in ) :: tend_new
    
      residual = abs(maxval(tend_new%q - tend_old%q))
    end subroutine calc_residual
    
end module temporal_mod
