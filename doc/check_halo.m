clc
clear

data = importdata('..\run\check_halo.txt');

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
iVar = 2;
iPatch = 5;
var_plt = squeeze(var(iVar,:,:,iPatch))';
plt = pcolor(var_plt);
colormap(jet)
% set(plt,'edgeColor','none')
shading interp