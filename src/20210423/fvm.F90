    program fvm
      use parameters_mod
      use stat_mod
      use tend_mod
      use reconstruction_mod
      use mesh_mod
      use test_case_mod
      use io_mod
      use temporal_mod
      use spatial_operators_mod
      use diag_mod
      implicit none
      
      integer(i_kind) it
      
      integer(i_kind) :: old = 0
      integer(i_kind) :: new = 1
      
      real(r_kind)    :: total_mass0,total_energy0
      real(r_kind)    :: total_mass ,total_energy
      
      integer(i_kind) :: output_idx, total_output_num
      integer(i_kind) :: output_step
      
      integer(i_kind) :: timeStart,timeEnd
      
      !real(r4) :: timer_coef = 1.e6
      real(r4) :: timer_coef = 1.e4
      
      !Timing start
      call SYSTEM_CLOCK(timeStart)
      
      print*,'Start the Finite Volume Shallow Water model'
      call init_parameters
      call init_reconstruction
      call init_mesh
      call init_stat
      call init_tend
      call init_testcase
      call init_spatial_operator
      
      print*,''
      print*,'Temporal integration scheme is '//trim(adjustl(integral_scheme))
      print*,''
      
      output_step = nint(history_interval/dt)
      
      output_idx       = 1
      total_output_num = nsteps * dt / history_interval
      call history_init      (stat(old))
      call history_write_stat(stat(old),1)
      call calc_total_mass   (total_mass0  ,stat(old))
      call calc_total_energy (total_energy0,stat(old))
      print*,'Initial total mass, total energy',total_mass0,total_energy0
      print*,'output index/total, MCR, ECR :',output_idx-1,'/',total_output_num,' ',0., 0.
      
      ! time integration
      do it = 1,nsteps
        if(trim(adjustl(integral_scheme))=='RK3_TVD')call RK3_TVD(stat(new),stat(old))
        if(trim(adjustl(integral_scheme))=='RK4'    )call RK4    (stat(new),stat(old))
        
        ! Write output
        if( mod(it,output_step)==0 .and. (it>=output_step) )then
          output_idx = output_idx + 1
          
          call history_write_stat(stat(new),output_idx)
          call calc_total_mass  (total_mass  ,stat(new))
          call calc_total_energy(total_energy,stat(old))
          print*,'output index/total, MCR, ECR :',output_idx-1,'/',total_output_num,' ',(total_mass-total_mass0)/total_mass0, (total_energy-total_energy0)/total_energy0
        endif
      
        if(any(isnan(stat(new)%q(1:nVar,ids:ide,jds:jde,ifs:ife))))then
          print*,'Nan appears at iVar i j iPatch: ',maxloc(stat(new)%q(1:nVar,ids:ide,jds:jde,ifs:ife),isnan(stat(new)%q(1:nVar,ids:ide,jds:jde,ifs:ife)))
          print*,'after ',it,' steps'
          stop 'Model blow up, maybe smaller dt would help'
        endif
        
        call copyStat(stat(old), stat(new))
      enddo
      
      ! Timing end
      call SYSTEM_CLOCK(timeEnd)
      
      print*,'It took ',real(timeEnd-timeStart,r4)/timer_coef,' seconds to run this program'
    end program fvm
