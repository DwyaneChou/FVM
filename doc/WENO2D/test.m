clc
clear

ngpts = 2;
[pG,wG] = Gaussian_Quadrature(ngpts);
pG = pG / 2;
xG = 0.5;
yG = pG(1);

% lack_pos = [];
lack_pos = [19,20,24,25];
% lack_pos = [20,25];
% lack_pos = [24,25];
% lack_pos = [25];

high_stencil_width = 5;
nxh = high_stencil_width;
nyh = high_stencil_width;
high_maxTerms = high_stencil_width^2;
high_recBdy = ( high_stencil_width - 1 ) / 2;

% Calculate 3th order polynomial
low_stencil_width = 3;
nx            = low_stencil_width;
ny            = low_stencil_width;
low_maxTerms  = low_stencil_width^2;
recBdy        = ( low_stencil_width - 1 ) / 2;

if(any(lack_pos==19))
    nStencil = 8;
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
%     % up right
%     ic(7) = recBdy;
%     jc(7) = recBdy;
    % up middle
    ic(7) = 0;
    jc(7) = recBdy;
    % up left
    ic(8) = -recBdy;
    jc(8) = recBdy;
else
    nStencil = 9;
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
end

iCOS = 0;
POS  = zeros(high_stencil_width,high_stencil_width);
for j = 1:high_stencil_width
    for i = 1:high_stencil_width
        iCOS = iCOS + 1;
        POS(i,j) = iCOS;
    end
end

A         = zeros(nStencil,low_maxTerms,low_maxTerms);
iA        = zeros(nStencil,low_maxTerms,low_maxTerms);
COS       = zeros(nStencil,low_maxTerms); % Cells on Stencil
existTerm = zeros(nStencil,low_maxTerms);
poly      = zeros(nStencil,low_maxTerms);
for iStencil = 1:nStencil
    iPOS  = 0;
    iRCOS = 0;
    for j = -recBdy:recBdy
        for i = -recBdy:recBdy
            iPOS = iPOS + 1;
            iP = i+ic(iStencil);
            jP = j+jc(iStencil);
            iCOH = POS(iP+high_recBdy+1,jP+high_recBdy+1);
            if ~any(lack_pos==iCOH)
                iRCOS = iRCOS + 1;
                COS(iStencil,iRCOS) = POS(iP+high_recBdy+1,jP+high_recBdy+1);
                existTerm(iStencil,iPOS) = 1;
            end
        end
    end
    nCOS(iStencil) = iRCOS;
    nRC = nCOS(iStencil);
    iRCOS = 0;
    for j = -recBdy:recBdy
        for i = -recBdy:recBdy
            iPOS = iPOS + 1;
            iP = i+ic(iStencil);
            jP = j+jc(iStencil);
            x_min = iP-0.5;
            x_max = iP+0.5;
            y_min = jP-0.5;
            y_max = jP+0.5;
            iCOH  = POS(iP+high_recBdy+1,jP+high_recBdy+1);
            if ~any(lack_pos==iCOH)
                iRCOS = iRCOS + 1;
                A(iStencil,iRCOS,1:nRC) = polynomial_integration(nx,ny,x_min,x_max,y_min,y_max,existTerm(iStencil,:));
            end
        end
    end
    % iA(iStencil,:,:) = inv(squeeze(A(iStencil,:,:)));

    poly(iStencil,1:nRC) = polynomial(nx,ny,xG,yG,existTerm(iStencil,:));
end

R = zeros(nStencil,low_maxTerms);
for iStencil = 1:nStencil
    nRC = nCOS(iStencil);
    R(iStencil,1:nRC) = poly(iStencil,1:nRC) / squeeze(A(iStencil,1:nRC,1:nRC));
end

% Calculate high order matrix
iPOS = 0;
existTermh = zeros(1,high_maxTerms);
nCOSh = 0;
for j = -high_recBdy:high_recBdy
    for i = -high_recBdy:high_recBdy
        iPOS = iPOS + 1;
        iP = i + high_recBdy + 1;
        jP = j + high_recBdy + 1;
        if ~any(lack_pos==POS(iP,jP))
            nCOSh = nCOSh + 1;
            existTermh(iPOS) = 1;
        end
    end
end

Ah = zeros(nCOSh,nCOSh);
iPOS = 0;
for j = -high_recBdy:high_recBdy
    for i = -high_recBdy:high_recBdy
        iPOS = iPOS + 1;
        if existTermh(iPOS)
            x_min = i-0.5;
            x_max = i+0.5;
            y_min = j-0.5;
            y_max = j+0.5;
            
            Ah(iPOS,1:nCOSh) = polynomial_integration(nxh,nyh,x_min,x_max,y_min,y_max,existTermh);
        end
    end
end

polyh = polynomial(nxh,nyh,xG,yG,existTermh);

Rh = polyh / Ah;

Rl_on_Rh = zeros(nStencil,high_stencil_width);
for iStencil = 1:nStencil
    iPOS = 0;
    iRCOS = 0;
    for j = 1:low_stencil_width
        for i = 1:low_stencil_width
            iPOS = iPOS + 1;
            if existTerm(iStencil,iPOS)
                iRCOS = iRCOS + 1;
                iCOS = COS(iStencil,iRCOS);
                Rl_on_Rh(iStencil,iCOS) = R(iStencil,iRCOS);
            end
        end
    end
end

C = Rl_on_Rh' \ Rh' ;

for iStencil = 1:nStencil
    post(iStencil,:) = C(iStencil).*Rl_on_Rh(iStencil,:);
end
post_check = sum(post,1);
diff = sum(post_check-Rh)

function poly = polynomial_integration(nx,ny,x_min,x_max,y_min,y_max,existTerm)
n = sum(existTerm);
poly = zeros(n,1);
iCOS = 0;
iRCOS = 0;
for j = 0:ny-1
    for i = 0:nx-1
        iCOS = iCOS + 1;
        if existTerm(iCOS)
            iRCOS = iRCOS + 1;
            poly(iRCOS) = ( x_max^(i+1) - x_min^(i+1) ) * ( y_max^(j+1) - y_min^(j+1) ) / ( ( i + 1 ) * ( j + 1 ) );
        end
    end
end
end

function poly = polynomial(nx,ny,x,y,existTerm)
np = length(x);
n = sum(existTerm);
poly = zeros(np,n);
iCOS = 0;
iRCOS = 0;
for j = 0:ny-1
    for i = 0:nx-1
        iCOS = iCOS + 1;
        if existTerm(iCOS)
            iRCOS = iRCOS + 1;
            poly(:,iRCOS) = x.^i .* y.^j;
        end
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