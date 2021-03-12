module tend_mod
  use constants_mod
  use parameters_mod
  implicit none
  
  type tend_field
    real(r_kind), dimension(:,:,:,:), allocatable :: q
  end type tend_field
  
  type(tend_field), dimension(:), allocatable :: tend ! allocated by n time points, which is used by temporal integration schemes
  
  contains
  
  subroutine init_tend
    integer :: iT
    
    allocate( tend(-nIntegralSubSteps:1) )
    
    do iT = -nIntegralSubSteps, 1
      allocate(tend(iT)%q(nVar,ims:ime,jms:jme,ifs:ife))
      tend(iT)%q = 0.
    enddo
    
  end subroutine init_tend

  subroutine copyTend(tend_out,tend_in)
    type(tend_field),intent(inout) :: tend_out
    type(tend_field),intent(in   ) :: tend_in
    
    integer :: iVar,i,j,iPatch
  
    !tend_out%q = tend_in%q

    !$OMP PARALLEL DO PRIVATE(i,j,iVar) COLLAPSE(4)
    do iPatch = ifs,ife
      do j = jds,jde
        do i = ids,ide
          do iVar = 1,nVar
            tend_out%q(iVar,i,j,iPatch) = tend_in%q(iVar,i,j,iPatch)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine copyTend
end module tend_mod
