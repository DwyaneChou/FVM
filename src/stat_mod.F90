module stat_mod
  use constants_mod
  use parameters_mod
  implicit none
  type stat_field
    real(r_kind), dimension(:,:,:,:), allocatable :: q
  end type stat_field
  
  type(stat_field), dimension(:), allocatable :: stat  ! allocated by n time points, which is used by temporal integration schemes
  
    contains
    subroutine init_stat
      integer :: iT
      
      allocate(stat(-nIntegralSubSteps:1))
      
      do iT = -nIntegralSubSteps, 1
        allocate(stat(iT)%q (nVar,ims:ime,jms:jme,ifs:ife))
        stat(iT)%q = FillValue
      enddo
    end subroutine init_stat

    subroutine copyStat(stat_out,stat_in)
      type(stat_field),intent(inout) :: stat_out
      type(stat_field),intent(in   ) :: stat_in
    
      integer :: iVar,i,j,iPatch

      !$OMP PARALLEL DO PRIVATE(i,j,iVar) COLLAPSE(4)
      do iPatch = ifs,ife
        do j = jms,jme
          do i = ims,ime
            do iVar = 1,nVar
              stat_out%q(iVar,i,j,iPatch) = stat_in%q(iVar,i,j,iPatch)
            enddo
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
    end subroutine copyStat
end module stat_mod
