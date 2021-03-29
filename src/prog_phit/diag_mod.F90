MODULE diag_mod
  use constants_mod
  use mesh_mod
  use parameters_mod
  use stat_mod
  
  implicit none
  
    contains
    subroutine calc_total_mass(total_mass,stat)
      type(stat_field), intent(in ) :: stat
      real(r_kind)    , intent(out) :: total_mass
      
      real(r_kind)    massOnCell(ids:ide,jds:jde,ifs:ife)
      integer i,j,iPatch
      
      do iPatch = ifs, ife
        do j = jds,jde
          do i = ids,ide
            massOnCell(i,j,iPatch) = ( stat%q(1,i,j,iPatch) - sqrtGC(i,j,iPatch) * zsC(i,j,iPatch) * gravity ) * dx * dy
          enddo
        enddo
      enddo
      
      total_mass = sum(massOnCell)
    
    end subroutine calc_total_mass
    
    subroutine calc_total_energy(total_energy,stat)
      type(stat_field), intent(in ) :: stat
      real(r_kind)    , intent(out) :: total_energy
      
      real(r_kind) u (ids:ide,jds:jde,ifs:ife)
      real(r_kind) v (ids:ide,jds:jde,ifs:ife)
      real(r_kind) uc(ids:ide,jds:jde,ifs:ife)
      real(r_kind) vc(ids:ide,jds:jde,ifs:ife)
      real(r_kind) energyOnCell(ids:ide,jds:jde,ifs:ife)
      real(r_kind) KE,PE
      integer i,j,iPatch
      
      do iPatch = ifs, ife
        do j = jds,jde
          do i = ids,ide
            uc(i,j,iPatch) = stat%q(2,i,j,iPatch) / ( stat%q(1,i,j,iPatch) - sqrtGC(i,j,iPatch) * zsC(i,j,iPatch) * gravity )
            vc(i,j,iPatch) = stat%q(3,i,j,iPatch) / ( stat%q(1,i,j,iPatch) - sqrtGC(i,j,iPatch) * zsC(i,j,iPatch) * gravity )
            call contrav2cov(u(i,j,iPatch),v(i,j,iPatch),uc(i,j,iPatch),vc(i,j,iPatch),matrixG(:,:,cc,i,j,iPatch))
          enddo
        enddo
      enddo
      
      energyOnCell = 0.
      total_energy = 0.
      do iPatch = ifs, ife
        do j = jds,jde
          do i = ids,ide
            KE = 0.5 * ( stat%q(1,i,j,iPatch) - sqrtGC(i,j,iPatch) * zsc(i,j,iPatch) * gravity ) * ( u(i,j,iPatch) * uc(i,j,iPatch) + v(i,j,iPatch) * vc(i,j,iPatch) )
            PE = 0.5 * stat%q(1,i,j,iPatch)**2
            
            energyOnCell(i,j,iPatch) = ( KE + PE ) * dx * dy
          enddo
        enddo
      enddo
      
      total_energy = sum(energyOnCell)
    end subroutine calc_total_energy
END MODULE diag_mod

