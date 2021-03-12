    program fvm
      use parameters_mod
      use spatial_operators_mod
      implicit none
      
      print*,'Start the Finite Volume Shallow Water model'
      call init_parameters
      call init_spatial_operator
      
      call fill_halo
    end program fvm
