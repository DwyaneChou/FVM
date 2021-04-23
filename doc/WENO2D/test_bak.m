clc
clear

ngpts = 2;
[pG,wG] = Gaussian_Quadrature(ngpts);
pG = pG / 2;
xG = 0.5;
yG = pG(2);

lack_pos = [];

high_stencil_width = 5;
nxh = high_stencil_width;
nyh = high_stencil_width;
high_maxTerms = high_stencil_width^2;
high_recBdy = ( high_stencil_width - 1 ) / 2;

% Calculate 3th order polynomial
nStencil = 9;
low_stencil_width = 3;
nx            = low_stencil_width;
ny            = low_stencil_width;
low_maxTerms  = low_stencil_width^2;
recBdy        = ( low_stencil_width - 1 ) / 2;

% middle middle
ic(1) = 0;
jc(1) = 0;
% middle left
ic(2) = -recBdy;
jc(2) = 0;
% low left
ic(3) = -recBdy;
jc(3) = -recBdy;
% low middle
ic(4) = 0;
jc(4) = -recBdy;
% low right
ic(5) = recBdy;
jc(5) = -recBdy;
% middle right
ic(6) = recBdy;
jc(6) = 0;
% up right
ic(7) = recBdy;
jc(7) = recBdy;
% up middle
ic(8) = 0;
jc(8) = recBdy;
% up left
ic(9) = -recBdy;
jc(9) = recBdy;

iCOS = 0;
POS  = zeros(high_stencil_width,high_stencil_width);
for j = 1:high_stencil_width
    for i = 1:high_stencil_width
        iCOS = iCOS + 1;
        POS(i,j) = iCOS;
    end
end

A   = zeros(nStencil,low_maxTerms,low_maxTerms);
iA  = zeros(nStencil,low_maxTerms,low_maxTerms);
COS = zeros(nStencil,low_maxTerms); % Cells on Stencil
for iStencil = 1:nStencil
    iPOS = 0;
    for j = -recBdy:recBdy
        for i = -recBdy:recBdy
            iPOS = iPOS + 1;
            iP = i+ic(iStencil);
            jP = j+jc(iStencil);
            x_min = iP-0.5;
            x_max = iP+0.5;
            y_min = jP-0.5;
            y_max = jP+0.5;
            A(iStencil,iPOS,:) = polynomial_integration(nx,ny,x_min,x_max,y_min,y_max);
            COS(iStencil,iPOS) = POS(iP+high_recBdy+1,jP+high_recBdy+1);
        end
    end
%     iA(iStencil,:,:) = inv(squeeze(A(iStencil,:,:)));
end

poly = polynomial(nx,ny,xG,yG);

R = zeros(nStencil,low_maxTerms);
for iStencil = 1:nStencil
%     R(iStencil,:) = poly * squeeze(iA(iStencil,:,:));
    R(iStencil,:) = poly / squeeze(A(iStencil,:,:));
end

% Calculate high order matrix
iPOS = 0;
Ah = zeros(high_maxTerms,high_maxTerms);
for j = -high_recBdy:high_recBdy
    for i = -high_recBdy:high_recBdy
        iPOS = iPOS + 1;
        x_min = i-0.5;
        x_max = i+0.5;
        y_min = j-0.5;
        y_max = j+0.5;
        Ah(iPOS,:) = polynomial_integration(nxh,nyh,x_min,x_max,y_min,y_max);
    end
end
% iAh = inv(Ah);

polyh = polynomial(nxh,nyh,xG,yG);

% Rh = polyh * iAh;
Rh = polyh / Ah;

Rl_on_Rh = zeros(nStencil,high_stencil_width);
for iStencil = 1:nStencil
    iPOS = 0;
    for j = 1:low_stencil_width
        for i = 1:low_stencil_width
            iPOS = iPOS + 1;
            iCOS = COS(iStencil,iPOS);
            Rl_on_Rh(iStencil,iCOS) = R(iStencil,iPOS);
        end
    end
end

C = Rl_on_Rh' \ Rh' ;

for iStencil = 1:nStencil
    post(iStencil,:) = C(iStencil).*Rl_on_Rh(iStencil,:);
end
post_check = sum(post,1);
diff = sum(post_check-Rh)

function poly = polynomial_integration(nx,ny,x_min,x_max,y_min,y_max)
poly = zeros(nx*ny,1);
iCOS = 0;
for j = 0:ny-1
    for i = 0:nx-1
        iCOS = iCOS + 1;
        poly(iCOS) = ( x_max^(i+1) - x_min^(i+1) ) * ( y_max^(j+1) - y_min^(j+1) ) / ( ( i + 1 ) * ( j + 1 ) );
    end
end
end

function poly = polynomial(nx,ny,x,y)
n = length(x);
poly = zeros(n,nx*ny);
iCOS = 0;
for j = 0:ny-1
    for i = 0:nx-1
        iCOS = iCOS + 1;
        poly(:,iCOS) = x.^i .* y.^j;
    end
end
end

function [p,w] = Gaussian_Quadrature(ngpts)
% Calculate Gauss-Laguerre Quadrature Evaluation Points and Weights
syms x
syms t
gpts_c = vpasolve(legendreP(ngpts,x) == 0); % Gaussian points
p      = 1:ngpts;
for k=1:ngpts
    j        = setdiff(p,k);
    wts_c(k) = int(prod((t - gpts_c(j)) ./ (gpts_c(k) - gpts_c(j))), t, -1, 1);
    
%     disp('Gaussian quadrature points')
%     disp(gpts_c(k))
%     
%     disp('Gaussian quadrature weights')
%     disp(wts_c(k))
end

p = double(gpts_c);
w = double(wts_c);
end