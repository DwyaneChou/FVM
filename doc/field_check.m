% This program is the test of Saturated Vapor Pressure
% clc
% clear

var_name = 'phit';
% var_name = 'zonal_wind';
% var_name = 'meridional_wind';
it       = 337;

nc_file = '..\run\output.nc';

dx         = ncreadatt(nc_file,'/','dx');
ids        = ncreadatt(nc_file,'/','ids');
ide        = ncreadatt(nc_file,'/','ide');
jds        = ncreadatt(nc_file,'/','jds');
jde        = ncreadatt(nc_file,'/','jde');
ims        = ncreadatt(nc_file,'/','ims');
ime        = ncreadatt(nc_file,'/','ime');
jms        = ncreadatt(nc_file,'/','jms');
jme        = ncreadatt(nc_file,'/','jme');
ifs        = ncreadatt(nc_file,'/','ifs');
ife        = ncreadatt(nc_file,'/','ife');
Fill_Value = ncreadatt(nc_file,var_name,'_FillValue');

nHalo  = ids-ims;
Nx     = ide-ids+1;
Ny     = jde-jds+1;
Npatch = ife - ifs + 1;

lon = ncread(nc_file,'lon'   );
lat = ncread(nc_file,'lat'   );
var = ncread(nc_file,var_name,[ids,jds,1,it],[Nx,Ny,Npatch,1]);

lon(lon<0) = 360 + lon((lon<0));

lon1d = reshape(lon,[],1);
lat1d = reshape(lat,[],1);
var1d = reshape(var,[],1);

res = dx;
x   = 0:res:360;
y   = -90:res:90;

[lon2d,lat2d] = meshgrid(x,y);

var_plot = griddata(lon1d,lat1d,var1d,lon2d,lat2d,'linear');

figure
pcolor(lon2d,lat2d,var_plot)
shading interp
% set(gca,'CLim',[-16,38])
colormap(jet)