% clc
% clear

data = importdata('..\run\check_tend.txt');

nx = size(data,2);
ny = nx;
nVar = 3;
nPatch = 6;

var = zeros(nVar,nx,ny,nPatch);

count = 0;
for iVar = 1:nVar
    for iPatch = 1:nPatch
        for j = 1:ny
            count = count + 1;
            var(iVar,:,j,iPatch) = data(count,:);
        end
    end
end

figure
iVar = 3;
iPatch = 5;
ids = 3;
ide = 28;
jds = 3;
jde = 28;
% var_plt = squeeze(var(iVar,ids:ide,jds:jde,iPatch))';
var_plt = squeeze(var(iVar,:,:,iPatch))';
plt = pcolor(var_plt);
colormap(jet)
% set(plt,'edgeColor','none')
% set(gca,'Clim',[-100000,100000])
shading interp