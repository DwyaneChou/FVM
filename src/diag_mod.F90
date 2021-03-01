MODULE diag_mod
  use constants_mod
  use mesh_mod
  use parameters_mod
  use stat_mod
  
  implicit none
  
    contains
    subroutine calc_total_mass(total_mass,stat)
      type(stat_field), intent(in ) :: stat
      real            , intent(out) :: total_mass
      
      real    massOnCell(ids:ide,jds:jde,ifs:ife)
      integer i,j,iPatch
      
      do iPatch = ifs, ife
        do j = jds,jde
          do i = ids,ide
            massOnCell(i,j,iPatch) = stat%q(1,i,j,iPatch) * dx * dy
          enddo
        enddo
      enddo
      
      total_mass = sum(massOnCell)
    
    end subroutine calc_total_mass
    
    subroutine calc_total_energy(total_energy,stat)
      type(stat_field), intent(in ) :: stat
      real            , intent(out) :: total_energy
      
      real u (ids:ide,jds:jde,ifs:ife)
      real v (ids:ide,jds:jde,ifs:ife)
      real uc(ids:ide,jds:jde,ifs:ife)
      real vc(ids:ide,jds:jde,ifs:ife)
      real energyOnCell(ids:ide,jds:jde,ifs:ife)
      real KE,PE
      integer i,j,iPatch
      
      do iPatch = ifs, ife
        do j = jds,jde
          do i = ids,ide
            uc(i,j,iPatch) = stat%q(2,i,j,iPatch) / stat%q(1,i,j,iPatch)
            vc(i,j,iPatch) = stat%q(3,i,j,iPatch) / stat%q(1,i,j,iPatch)
            call contrav2cov(u(i,j,iPatch),v(i,j,iPatch),uc(i,j,iPatch),vc(i,j,iPatch),matrixG(:,:,cc,i,j,iPatch))
          enddo
        enddo
      enddo
      
      energyOnCell = 0.
      total_energy = 0.
      do iPatch = ifs, ife
        do j = jds,jde
          do i = ids,ide
            KE = 0.5 * stat%q(1,i,j,iPatch) * ( u(i,j,iPatch) * uc(i,j,iPatch) + v(i,j,iPatch) * vc(i,j,iPatch) )
            PE = 0.5 * ( stat%q(1,i,j,iPatch) + sqrtGC(i,j,iPatch) * zsc(i,j,iPatch) * gravity )**2
            
            energyOnCell(i,j,iPatch) = ( KE + PE ) * dx * dy
          enddo
        enddo
      enddo
      
      total_energy = sum(energyOnCell)
    end subroutine calc_total_energy
END MODULE diag_mod

