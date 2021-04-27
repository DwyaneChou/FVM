clc
clear

digits(36);

ngpts = vpa(2);
[pG,wG] = Gaussian_Quadrature(ngpts);
pG = pG / 2;

nWenoPoints = 12;
% left
x(1) = vpa(-0.5);
y(1) = pG(1);

x(2) = vpa(-0.5);
y(2) = pG(2);

% right
x(3) = vpa(0.5);
y(3) = pG(1);

x(4) = vpa(0.5);
y(4) = pG(2);

% bottom
x(5) = pG(1);
y(5) = vpa(-0.5);

x(6) = pG(2);
y(6) = vpa(-0.5);

% top
x(7) = pG(1);
y(7) = vpa(0.5);

x(8) = pG(2);
y(8) = vpa(0.5);

% quadrature points
x(9) = pG(1);
y(9) = pG(1);

x(10) = pG(2);
y(10) = pG(1);

x(11) = pG(2);
y(11) = pG(2);

x(12) = pG(1);
y(12) = pG(2);

lack = zeros(2,1);
lack(2,1) = 9;

nType = length(lack);

high_stencil_width = 3;
nxh = high_stencil_width;
nyh = high_stencil_width;
high_maxTerms = high_stencil_width^2;
high_recBdy = ( high_stencil_width - 1 ) / 2;

% Calculate 3th order polynomial
low_stencil_width = 2;
nx            = low_stencil_width;
ny            = low_stencil_width;
low_maxTerms  = low_stencil_width^2;
recBdy        = ( low_stencil_width - 1 ) / 2;

fid=fopen('coef3.txt','w');

for iType = 1:2
    lack_pos = lack(iType,:);
    nStencil = 4;
    ic(1:nStencil) = 0;
    jc(1:nStencil) = 0;
    
    iCOS = 0;
    POS  = zeros(high_stencil_width,high_stencil_width);
    for j = 1:high_stencil_width
        for i = 1:high_stencil_width
            iCOS = iCOS + 1;
            POS(i,j) = iCOS;
        end
    end
    
    Coef = vpa( zeros(nWenoPoints,nStencil) );
    for iPoint = 1:nWenoPoints
        xG = x(iPoint);
        yG = y(iPoint);
        
        A         = vpa(zeros(nStencil,low_maxTerms,low_maxTerms));
        iA        = zeros(nStencil,low_maxTerms,low_maxTerms);
        COS       = zeros(nStencil,low_maxTerms); % Cells on Stencil
        existTerm = ones(nStencil,low_maxTerms);
        poly      = vpa( zeros(nStencil,low_maxTerms) );
        
        COS(1,1) = 1;
        COS(1,2) = 2;
        COS(1,3) = 4;
        COS(1,4) = 5;
        COS(2,1) = 2;
        COS(2,2) = 3;
        COS(2,3) = 5;
        COS(2,4) = 6;
        COS(3,1) = 5;
        COS(3,2) = 6;
        COS(3,3) = 8;
        COS(3,4) = 9;
        COS(4,1) = 4;
        COS(4,2) = 5;
        COS(4,3) = 7;
        COS(4,4) = 8;
        xc(1) = -1; yc(1) = -1;
        xc(2) =  0; yc(2) = -1;
        xc(3) =  1; yc(3) = -1;
        xc(4) = -1; yc(4) =  0;
        xc(5) =  0; yc(5) =  0;
        xc(6) =  1; yc(6) =  0;
        xc(7) = -1; yc(7) =  1;
        xc(8) =  0; yc(8) =  1;
        xc(9) =  1; yc(9) =  1;
        if any(lack_pos==COS(3,4))
            existTerm(3,4) = 0;
        end
        
        for iStencil = 1:nStencil
            nCOS(iStencil) = sum(existTerm(iStencil,:));
            nRC = nCOS(iStencil);
            for iRCOS = 1:low_maxTerms
                x_min = xc(COS(iStencil,iRCOS)) - 0.5;
                x_max = xc(COS(iStencil,iRCOS)) + 0.5;
                y_min = yc(COS(iStencil,iRCOS)) - 0.5;
                y_max = yc(COS(iStencil,iRCOS)) + 0.5;
                A(iStencil,iRCOS,1:nRC) = polynomial_integration(nx,ny,x_min,x_max,y_min,y_max,existTerm(iStencil,:));
            end
            % iA(iStencil,:,:) = inv(squeeze(A(iStencil,:,:)));
            
            poly(iStencil,1:nRC) = polynomial(nx,ny,xG,yG,existTerm(iStencil,:));
        end
        
        R = vpa( zeros(nStencil,low_maxTerms) );
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
                    COSh(iPOS) = nCOSh;
                    existTermh(iPOS) = 1;
                end
            end
        end
        
        Ah = vpa( zeros(nCOSh,nCOSh) );
        iPOS = 0;
        iCOS = 0;
        for j = -high_recBdy:high_recBdy
            for i = -high_recBdy:high_recBdy
                iPOS = iPOS + 1;
                if existTermh(iPOS)
                    iCOS = iCOS + 1;
                    x_min = i-0.5;
                    x_max = i+0.5;
                    y_min = j-0.5;
                    y_max = j+0.5;
                    
                    Ah(iCOS,1:nCOSh) = polynomial_integration(nxh,nyh,x_min,x_max,y_min,y_max,existTermh);
                end
            end
        end
        
        polyh = polynomial(nxh,nyh,xG,yG,existTermh);
        
        Rh = polyh / Ah;
        
        Rl_on_Rh = vpa( zeros(nStencil,nCOSh) );
        for iStencil = 1:nStencil
            iPOS = 0;
            iRCOS = 0;
            for j = 1:low_stencil_width
                for i = 1:low_stencil_width
                    iPOS = iPOS + 1;
                    if existTerm(iStencil,iPOS)
                        iRCOS = iRCOS + 1;
                        iCOS = COSh(COS(iStencil,iRCOS));
                        Rl_on_Rh(iStencil,iCOS) = R(iStencil,iRCOS);
                    end
                end
            end
        end
        
        %     C = vpa( Rl_on_Rh' \ Rh' );
        C = inv( Rl_on_Rh * Rl_on_Rh' ) * Rl_on_Rh * Rh';
        
        post = zeros(nStencil,size(Rh,2));
        for iStencil = 1:nStencil
            post(iStencil,:) = C(iStencil).*Rl_on_Rh(iStencil,:);
        end
        post_check = sum(post,1);
        diff = sum(post_check-Rh);
        
        disp(['Type ',num2str(iType),', ','Point ',num2str(iPoint),', diff=',num2str(double(diff))])
        
        if any(lack_pos==19)
            C(9) = C(8);
            C(8) = C(7);
            C(7) = 0;
        end
        
        Coef(iPoint,:) = C;
    end
    
    Coef(abs(Coef)<1.e-12) = 0;
    
    for j = 1:nWenoPoints
        fprintf( fid,['%',num2str(low_maxTerms),'s\n'],char(Coef(j,:)) );
    end
end

fclose(fid);

function poly = polynomial_integration(nx,ny,x_min,x_max,y_min,y_max,existTerm)
digits(36);
n = sum(existTerm);
poly = vpa( zeros(n,1) );
iCOS = 0;
iRCOS = 0;
for j = 0:ny-1
    for i = 0:nx-1
        iCOS = iCOS + 1;
        if existTerm(iCOS)
            iRCOS = iRCOS + 1;
            poly(iRCOS) = vpa( ( x_max^(i+1) - x_min^(i+1) ) * ( y_max^(j+1) - y_min^(j+1) ) / ( ( i + 1 ) * ( j + 1 ) ) );
        end
    end
end
end

function poly = polynomial(nx,ny,x,y,existTerm)
digits(36);
np = length(x);
n = sum(existTerm);
poly = vpa( zeros(np,n) );
iCOS = 0;
iRCOS = 0;
for j = 0:ny-1
    for i = 0:nx-1
        iCOS = iCOS + 1;
        if existTerm(iCOS)
            iRCOS = iRCOS + 1;
            poly(:,iRCOS) = vpa( x.^i .* y.^j );
        end
    end
end
end

function [p,w] = Gaussian_Quadrature(ngpts)
digits(36);
% Calculate Gauss-Laguerre Quadrature Evaluation Points and Weights
syms x
syms t
gpts_c = vpa( vpasolve(legendreP(ngpts,x) == 0) ); % Gaussian points
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

p = vpa(gpts_c);
w = vpa(wts_c);
end