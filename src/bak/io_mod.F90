module io_mod
  use netcdf
  use constants_mod
  use parameters_mod
  use mesh_mod
  use stat_mod
  use projection_mod
  use diag_mod
  implicit none
    
    character(13) :: ncFile = 'output.nc'
    
    contains
    subroutine history_init(stat)
      type(stat_field), intent(in) :: stat
      
      integer status
      integer ncid
      integer x_dim_id,y_dim_id
      integer patch_dim_id
      integer time_dim_id
      integer lon_id,lat_id
      integer x_id,y_id
      integer areaCell_id
      integer time_id
      integer phi_id
      integer phis_id
      integer phit_id
      integer zonal_wind_id,meridional_wind_id
      integer vorticity_id
      
      status = nf90_create(ncFile, NF90_CLOBBER + NF90_64BIT_OFFSET , ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_def_dim(ncid,'lon'   ,Nx            ,x_dim_id    )
      status = nf90_def_dim(ncid,'lat'   ,Ny            ,y_dim_id    )
      status = nf90_def_dim(ncid,'nPatch',Nf            ,patch_dim_id)
      status = nf90_def_dim(ncid,'time'  ,NF90_UNLIMITED,time_dim_id )
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_def_var(ncid,'time'           ,NF90_INT   ,(/                               time_dim_id/),time_id           )
      status = nf90_def_var(ncid,'lon'            ,NF90_DOUBLE,(/x_dim_id,y_dim_id,patch_dim_id            /),lon_id            )
      status = nf90_def_var(ncid,'lat'            ,NF90_DOUBLE,(/x_dim_id,y_dim_id,patch_dim_id            /),lat_id            )
      status = nf90_def_var(ncid,'areaCell'       ,NF90_DOUBLE,(/x_dim_id,y_dim_id,patch_dim_id            /),areaCell_id       )
      status = nf90_def_var(ncid,'phis'           ,NF90_DOUBLE,(/x_dim_id,y_dim_id,patch_dim_id            /),phis_id           )
      status = nf90_def_var(ncid,'phi'            ,NF90_DOUBLE,(/x_dim_id,y_dim_id,patch_dim_id,time_dim_id/),phi_id            )
      status = nf90_def_var(ncid,'phit'           ,NF90_DOUBLE,(/x_dim_id,y_dim_id,patch_dim_id,time_dim_id/),phit_id           )
      status = nf90_def_var(ncid,'zonal_wind'     ,NF90_DOUBLE,(/x_dim_id,y_dim_id,patch_dim_id,time_dim_id/),zonal_wind_id     )
      status = nf90_def_var(ncid,'meridional_wind',NF90_DOUBLE,(/x_dim_id,y_dim_id,patch_dim_id,time_dim_id/),meridional_wind_id)
      status = nf90_def_var(ncid,'vorticity      ',NF90_DOUBLE,(/x_dim_id,y_dim_id,patch_dim_id,time_dim_id/),vorticity_id      )
      if(status/=nf90_noerr) call handle_err(status)
      
      !print*,'nf90_put_att'
      status = nf90_put_att(ncid,nf90_global       ,'dx'                 ,real(dx*R2D,r8))
      status = nf90_put_att(ncid,nf90_global       ,'dy'                 ,real(dy*R2D,r8))
      status = nf90_put_att(ncid,nf90_global       ,'dt'                 ,real(dt,r8))
      status = nf90_put_att(ncid,nf90_global       ,'xhalo'              ,xhalo)
      status = nf90_put_att(ncid,nf90_global       ,'yhalo'              ,yhalo)
      status = nf90_put_att(ncid,nf90_global       ,'case_num'           ,case_num)
      status = nf90_put_att(ncid,nf90_global       ,'reconstruct_scheme' ,reconstruct_scheme)
      status = nf90_put_att(ncid,nf90_global       ,'stencil_width'      ,stencil_width)
      status = nf90_put_att(ncid,nf90_global       ,'nPointsOnEdge'      ,nPointsOnEdge)
      status = nf90_put_att(ncid,nf90_global       ,'ims'                ,ims)
      status = nf90_put_att(ncid,nf90_global       ,'ime'                ,ime)
      status = nf90_put_att(ncid,nf90_global       ,'jms'                ,jms)
      status = nf90_put_att(ncid,nf90_global       ,'jme'                ,jme)
      status = nf90_put_att(ncid,nf90_global       ,'ids'                ,ids)
      status = nf90_put_att(ncid,nf90_global       ,'ide'                ,ide)
      status = nf90_put_att(ncid,nf90_global       ,'jds'                ,jds)
      status = nf90_put_att(ncid,nf90_global       ,'jde'                ,jde)
      status = nf90_put_att(ncid,nf90_global       ,'ifs'                ,ifs)
      status = nf90_put_att(ncid,nf90_global       ,'ife'                ,ife)
      
      status = nf90_put_att(ncid,lon_id            ,'units'    ,'degree_east' )
      status = nf90_put_att(ncid,lat_id            ,'units'    ,'degree_north')
      status = nf90_put_att(ncid,areaCell_id       ,'units'    ,'m^2'         )
      status = nf90_put_att(ncid,time_id           ,'units'    ,'seconds'     )
      status = nf90_put_att(ncid,phi_id            ,'units'    ,'m^2/s^2'     )
      status = nf90_put_att(ncid,phis_id           ,'units'    ,'m^2/s^2'     )
      status = nf90_put_att(ncid,phit_id           ,'units'    ,'m^2/s^2'     )
      status = nf90_put_att(ncid,zonal_wind_id     ,'units'    ,'m/s'         )
      status = nf90_put_att(ncid,meridional_wind_id,'units'    ,'m/s'         )
      status = nf90_put_att(ncid,vorticity_id      ,'units'    ,'m/s^2'       )
      
      status = nf90_put_att(ncid,lon_id            ,'long_name','longitude on sphere coordinate for Cells' )
      status = nf90_put_att(ncid,lat_id            ,'long_name','latitude on sphere coordinate for Cells'  )
      status = nf90_put_att(ncid,areaCell_id       ,'long_name','area of cells'                            )
      status = nf90_put_att(ncid,time_id           ,'long_name','time'                                     )
      status = nf90_put_att(ncid,phi_id            ,'long_name','geopotential height on points'            )
      status = nf90_put_att(ncid,phis_id           ,'long_name','surface height on points'                 )
      status = nf90_put_att(ncid,phit_id           ,'long_name','total geopotential height on points'      )
      status = nf90_put_att(ncid,zonal_wind_id     ,'long_name','zonal wind'                               )
      status = nf90_put_att(ncid,meridional_wind_id,'long_name','meridional wind'                          )
      status = nf90_put_att(ncid,vorticity_id      ,'long_name','relative vorticity'                       )
      
      ! Define coordinates
      status = nf90_put_att(ncid, x_id              ,'_CoordinateAxisTypes','lon lat nPatch')
      status = nf90_put_att(ncid, y_id              ,'_CoordinateAxisTypes','lon lat nPatch')
      status = nf90_put_att(ncid, lon_id            ,'_CoordinateAxisTypes','lon lat nPatch')
      status = nf90_put_att(ncid, lat_id            ,'_CoordinateAxisTypes','lon lat nPatch')
      status = nf90_put_att(ncid, phis_id           ,'_CoordinateAxisTypes','lon lat nPatch')
      status = nf90_put_att(ncid, phi_id            ,'_CoordinateAxisTypes','lon lat nPatch time')
      status = nf90_put_att(ncid, phit_id           ,'_CoordinateAxisTypes','lon lat nPatch time')
      status = nf90_put_att(ncid, areaCell_id       ,'_CoordinateAxisTypes','lon lat nPatch time')
      status = nf90_put_att(ncid, zonal_wind_id     ,'_CoordinateAxisTypes','lon lat nPatch time')
      status = nf90_put_att(ncid, meridional_wind_id,'_CoordinateAxisTypes','lon lat nPatch time')
      status = nf90_put_att(ncid, vorticity_id      ,'_CoordinateAxisTypes','lon lat nPatch time')
      
      status = nf90_put_att(ncid,phi_id            ,'_FillValue',real(FillValue,r8))
      status = nf90_put_att(ncid,phis_id           ,'_FillValue',real(FillValue,r8))
      status = nf90_put_att(ncid,phit_id           ,'_FillValue',real(FillValue,r8))
      status = nf90_put_att(ncid,zonal_wind_id     ,'_FillValue',real(FillValue,r8))
      status = nf90_put_att(ncid,meridional_wind_id,'_FillValue',real(FillValue,r8))
      status = nf90_put_att(ncid,vorticity_id      ,'_FillValue',real(FillValue,r8))
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_enddef(ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_put_var(ncid,lon_id     , real(lon     (cc,ids:ide,jds:jde,ifs:ife)*R2D, r8) )
      status = nf90_put_var(ncid,lat_id     , real(lat     (cc,ids:ide,jds:jde,ifs:ife)*R2D, r8) )
      status = nf90_put_var(ncid,areaCell_id, real(areaCell(   ids:ide,jds:jde,ifs:ife), r8))
      status = nf90_put_var(ncid,phis_id    , real(ghsC    (   ids:ide,jds:jde,ifs:ife) / gravity, r8) )
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_close(ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
    end subroutine history_init
    
    subroutine history_write_stat(stat,time_slot_num)
      type(stat_field), intent(in) :: stat
      integer         , intent(in) :: time_slot_num
      
      integer status
      integer ncid
      integer time_id
      integer u_id,v_id
      integer uc_id,vc_id
      integer phi_id
      integer phit_id
      integer zonal_wind_id,meridional_wind_id
      integer vorticity_id
      
      real(r_kind), dimension(:,:,:), allocatable :: phi
      real(r_kind), dimension(:,:,:), allocatable :: u
      real(r_kind), dimension(:,:,:), allocatable :: v
      real(r_kind), dimension(:,:,:), allocatable :: vor
      
      integer :: varid
      real(r8), dimension(ids:ide,jds:jde,ifs:ife) :: varout
      
      integer :: time(1)
      
      integer i,j,iPatch
    
      allocate(phi(ids:ide,jds:jde,ifs:ife))
      allocate(u  (ids:ide,jds:jde,ifs:ife))
      allocate(v  (ids:ide,jds:jde,ifs:ife))
      allocate(vor(ids:ide,jds:jde,ifs:ife))
      
      do iPatch = ifs,ife
        do j = jds,jde
          do i = ids,ide
            phi(i,j,iPatch) = stat%q(1,i,j,iPatch) / sqrtGC(i,j,iPatch)
            u  (i,j,iPatch) = stat%q(2,i,j,iPatch) / stat%q(1,i,j,iPatch)
            v  (i,j,iPatch) = stat%q(3,i,j,iPatch) / stat%q(1,i,j,iPatch)
            call contravProjPlane2Sphere(u(i,j,iPatch), v(i,j,iPatch), u(i,j,iPatch), v(i,j,iPatch), matrixA(:,:,cc,i,j,iPatch))
          enddo
        enddo
      enddo
      
      !call calc_relative_vorticity(vor,stat)
      
      time(1) = time_slot_num
      !print*,'nf90_open'
      status = nf90_open(ncFile,NF90_WRITE,ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
      !print*,'nf90_inq_varid'
      status = nf90_inq_varid(ncid,'time'           , time_id           )
      status = nf90_inq_varid(ncid,'phi'            , phi_id            )
      status = nf90_inq_varid(ncid,'phit'           , phit_id           )
      status = nf90_inq_varid(ncid,'zonal_wind'     , zonal_wind_id     )
      status = nf90_inq_varid(ncid,'meridional_wind', meridional_wind_id)
      status = nf90_inq_varid(ncid,'vorticity'      , vorticity_id      )
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_put_var(ncid, time_id, time, start=(/time_slot_num/),count=(/1/))
      
      ! phi
      varid  = phi_id
      varout = phi
      status = nf90_put_var(ncid, varid, varout, start=(/1,1,1,time_slot_num/),count=(/Nx,Ny,Nf,1/))
      
      !phit
      varid  = phit_id
      varout = phi + ghsC(ids:ide,jds:jde,ifs:ife)
      status = nf90_put_var(ncid, varid, varout, start=(/1,1,1,time_slot_num/),count=(/Nx,Ny,Nf,1/))
      
      ! zonal_wind
      varid  = zonal_wind_id
      varout = u
      status = nf90_put_var(ncid, varid, varout, start=(/1,1,1,time_slot_num/),count=(/Nx,Ny,Nf,1/))
      
      ! meridional_wind
      varid  = meridional_wind_id
      varout = v
      status = nf90_put_var(ncid, varid, varout, start=(/1,1,1,time_slot_num/),count=(/Nx,Ny,Nf,1/))
      if(status/=nf90_noerr) call handle_err(status)
      
      ! vorticity
      varid  = vorticity_id
      varout = vor
      status = nf90_put_var(ncid, varid, varout, start=(/1,1,1,time_slot_num/),count=(/Nx,Ny,Nf,1/))
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_close(ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
    end subroutine history_write_stat
    
    subroutine handle_err(status)
      implicit none
      integer,intent(in)::status
            
      if(status/=nf90_noerr)then
          print*, trim(nf90_strerror(status))
          stop "Stopped by netCDF"
      endif  
    endsubroutine handle_err
end module io_mod
    