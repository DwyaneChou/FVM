MODULE diag_mod
  use constants_mod
  use mesh_mod
  use parameters_mod
  use stat_mod
  use spatial_operators_mod
  
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
            massOnCell(i,j,iPatch) = stat%q(1,i,j,iPatch) * dx * dy
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
            PE = 0.5 * ( stat%q(1,i,j,iPatch) + sqrtGC(i,j,iPatch) * ghsC(i,j,iPatch) )**2
            
            energyOnCell(i,j,iPatch) = ( KE + PE ) * dx * dy
          enddo
        enddo
      enddo
      
      total_energy = sum(energyOnCell)
    end subroutine calc_total_energy
    
    subroutine calc_relative_vorticity(vorticity,stat)
      real(r_kind), dimension(ids:ide,jds:jde,ifs:ife), intent(out) :: vorticity
      type(stat_field), intent(in ) :: stat
      
      integer i,j,iPatch
      
      real(r_kind), dimension(:,:,:  ), allocatable :: u 
      real(r_kind), dimension(:,:,:  ), allocatable :: v 
      real(r_kind), dimension(:,:,:  ), allocatable :: uc
      real(r_kind), dimension(:,:,:  ), allocatable :: vc
      real(r_kind), dimension(:,:,:,:), allocatable :: dudy
      real(r_kind), dimension(:,:,:,:), allocatable :: dvdx
      
      allocate( u   (                  ims:ime,jms:jme,ifs:ife) )
      allocate( v   (                  ims:ime,jms:jme,ifs:ife) )
      allocate( uc  (                  ims:ime,jms:jme,ifs:ife) )
      allocate( vc  (                  ims:ime,jms:jme,ifs:ife) )
      allocate( dudy(nQuadPointsOnCell,ims:ime,jms:jme,ifs:ife) )
      allocate( dvdx(nQuadPointsOnCell,ims:ime,jms:jme,ifs:ife) )
      
      qC = stat%q
      
      call fill_halo(qC)
      
      do iPatch = ifs, ife
        do j = jms,jme
          do i = ims,ime
            uc(i,j,iPatch) = qC(2,i,j,iPatch) / qC(1,i,j,iPatch)
            vc(i,j,iPatch) = qC(3,i,j,iPatch) / qC(1,i,j,iPatch)
            call contrav2cov(u(i,j,iPatch),v(i,j,iPatch),uc(i,j,iPatch),vc(i,j,iPatch),matrixG(:,:,cc,i,j,iPatch))
          enddo
        enddo
      enddo
      
      call reconstruction(v,dqdx=dvdx)
      call reconstruction(u,dqdy=dudy)
      
      do iPatch = ifs, ife
        do j = jds,jde
          do i = ids,ide
            vorticity(i,j,iPatch) = cell_quadrature( ( dvdx(:,i,j,iPatch) - dudy(:,i,j,iPatch) ) / sqrtG(cqs:cqe,i,j,iPatch) )
          enddo
        enddo
      enddo
      
    end subroutine calc_relative_vorticity
    
END MODULE diag_mod

