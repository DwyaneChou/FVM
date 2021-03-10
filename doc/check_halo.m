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

iVar = 3;
iPatch = 6;
var_plt = squeeze(var(iVar,:,:,iPatch))';
plt = pcolor(var_plt);
set(plt,'edgeColor','none')