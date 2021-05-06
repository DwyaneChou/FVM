% This program is the test of Saturated Vapor Pressure
% clc
% clear

% var_name = 'phi';
var_name = 'phit';
% var_name = 'zonal_wind';
% var_name = 'meridional_wind';
% var_name = 'vorticity';
it       = 361;

gravity = 9.80616;

nc_file = '..\run\output.nc';
% nc_file = '..\run\output_poly.nc';
% nc_file = '..\run\output_WENO.nc';
% nc_file = '..\run\output_WENO2D.nc';

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
case_num   = ncreadatt(nc_file,'/','case_num');
Fill_Value = ncreadatt(nc_file,var_name,'_FillValue');

nHalo  = ids-ims;
Nx     = ide-ids+1;
Ny     = jde-jds+1;
Npatch = ife - ifs + 1;

lon = ncread(nc_file,'lon'   );
lat = ncread(nc_file,'lat'   );
var = ncread(nc_file,var_name,[ids,jds,1,it],[Nx,Ny,Npatch,1]);
var0= ncread(nc_file,var_name,[ids,jds,1, 1],[Nx,Ny,Npatch,1]);
% var = (var-var0)./var0;
% var = ncread(nc_file,var_name,[ids,jds,1,it],[Nx,Ny,Npatch,1])/gravity;
% var0= ncread(nc_file,var_name,[ids,jds,1,1],[Nx,Ny,Npatch,1])/gravity;

phis = ncread(nc_file,'phis');

coef = 1;
if(case_num==5)
    lon(lon>180) = lon(lon>180) - 360;
    xs = -180;
    xe = 180;
    if(strcmp(var_name,'phit'))
        coef = 1 / gravity;
    end
elseif(case_num==8)
    lon(lon<-270) = 360 - lon(lon<-270);
    lon(lon>90) = lon(lon>90) - 360;
    xs = -270;
    xe = 90;
else
    lon(lon<0) = 360 + lon((lon<0));
    xs = 0;
    xe = 360;
end

lon1d = reshape(lon,[],1);
lat1d = reshape(lat,[],1);
var1d = reshape(var,[],1);
% var1d = reshape(var-var0,[],1);

res = dx;
x   = xs:res:xe;
y   = -90:res:90;

[lon2d,lat2d] = meshgrid(x,y);

var_plot = griddata(lon1d,lat1d,var1d,lon2d,lat2d,'linear');

figure
plt = pcolor(lon2d,lat2d,var_plot*coef);
% shading interp
set(plt,'EdgeColor','None')
% set(gca,'CLim',[ 7.8494e+04,1.0350e+05])
% set(gca,'CLim',[ 4.9e4,5.9e4])
% set(gca,'CLim',[ -7e-9,7e-9])
colormap(jet)

if(case_num==5)
    if(strcmp(var_name,'phit'))
        hold on
        LevelList = 5050:50:5950;
        plt=contour(lon2d,lat2d,var_plot*coef,'LevelList',LevelList,'LineColor','k');
        %     plt=contour(lon2d,lat2d,var_plot,'LevelList',LevelList,'LineColor','k','ShowText','on');
    elseif(strcmp(var_name,'zonal_wind'))
        hold on
        LevelList = -30:10:60;
        plt=contour(lon2d,lat2d,var_plot,'LevelList',LevelList,'LineColor','k');
    elseif(strcmp(var_name,'meridional_wind'))
        hold on
        LevelList = -30:10:30;
        plt=contour(lon2d,lat2d,var_plot,'LevelList',LevelList,'LineColor','k');
    end
end