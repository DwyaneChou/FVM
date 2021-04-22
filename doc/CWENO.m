clc
clear

run_seconds      = 2;
dx               = 1/50;
dt               = 0.004;
history_interval = 0.2;
output_path      = 'pictures\';

% Prepare for CWENO by Zhu Jun, 2018
nStencil = 5;
maxTerms = 2 * nStencil - 1;
% Calculate polynomial
poly = zeros(nStencil,maxTerms,maxTerms);
for i = 1:nStencil-1
    n = i*2 + 1;
    for j = -i:i
        poly(i,i+1+j,1:n) = poly_integration(n,j-0.5,j+0.5);
    end
end

invA = zeros(nStencil,maxTerms,maxTerms);
invA(1,1) = 1;
for iStencil = 2:nStencil
    nterm = 2 * iStencil - 1;
    invA(iStencil,1:nterm,1:nterm) = inv( squeeze( poly(iStencil-1,1:nterm,1:nterm) ) );
end

r = zeros(nStencil,nStencil);
for j = 1:nStencil
    for i = 1:j
        r(i,j) = 10^(i-1);
    end
    r(:,j) = r(:,j) / sum(r(:,j));
end

x_min = -1;
x_max = 1;
x_res = dx;

x = x_min:x_res:x_max;

u = ones(size(x));
f = zeros(size(x));

% Square wave
f(x>-0.4&x<0.4) = 2;

% Sine wave
% f = sin( x/(x_max-x_min)*2*pi )+1;

n   = length(f);

nt  = run_seconds/dt;
ht  = history_interval/dt;

stat.u        = u;
stat.f        = f;
stat.x        = x;
stat.dx       = x_res;
stat.n        = n;
stat.ht       = ht;
stat.nt       = nt;
stat.dt       = dt;
stat.r        = r;
stat.invA     = invA;
stat.nStencil = nStencil;
stat.maxTerms = maxTerms;

plot_result(stat,0,output_path)

total_mass0 = sum(stat.f(1:end-1));

iht = 1;
for it = 1:nt
    stat_new = RK4(stat,dt);
    
    if mod(it,ht)==0
        total_mass = sum(stat.f(1:end-1));
        MCR = (total_mass - total_mass0)/total_mass0;
        disp(['Output result at ',num2str(iht*stat.ht*stat.dt),' second(s)',' MCR = ',num2str(MCR,'%e'),...
              ' Max/Min = ',num2str(max(stat.f),'%e'),' ',num2str(min(stat.f),'%e')])
        plot_result(stat,iht,output_path)
        iht = iht + 1;
    end
    
    stat = cp_stat(stat_new);
end

L2 = sqrt(sum((f(1:n-1)-stat.f(1:n-1)).^2)/sum(f(1:n-1).^2));

function tend = spatial_discrete(stat)
u        = stat.u;
f        = stat.f;
dx       = stat.dx;
n        = stat.n;
r        = stat.r;
invA     = stat.invA;
nStencil = stat.nStencil;
maxTerms = stat.maxTerms;

eps = 1.E-16;
nm1 = n-1;
kappa = 4;

stencil_width = 1:4;
max_stencil_width = max(stencil_width);
% Extend halo
nfc = length(f) + 2 * max_stencil_width;
fc = zeros(1,nfc);
fc(1+max_stencil_width:nfc-max_stencil_width) = f(1:n);
fc(1:max_stencil_width                      ) = f(n-max_stencil_width:n-1);
fc(nfc-max_stencil_width+1:nfc              ) = f(2:max_stencil_width+1);

fR = zeros(1,n+1);
fL = zeros(1,n+1);
qc = zeros(5,9);
q  = zeros(5,9);
p  = zeros(5,9);
w = zeros(nStencil,1);
p_final = zeros(1,maxTerms);
for i = 1:n
    idx = i + max_stencil_width;
    qc(1,1) = fc(idx);
    for iStencil = 2:nStencil
        is  = idx - iStencil + 1;
        ie  = idx + iStencil - 1;
        qc(iStencil,1:2*iStencil-1) = fc(is:ie);
    end
    
    for iStencil = 1:nStencil
        nterms = 2 * iStencil - 1;
        q(iStencil,1:nterms) = squeeze( invA(iStencil,1:nterms,1:nterms) ) * qc(iStencil,1:nterms)';
    end
    
    p(1,1) = q(1,1);
    for iStencil = 2:nStencil
        for iterm = 1:2*iStencil-1
            p(iStencil,iterm) = ( q(iStencil,iterm) - dot( r(1:iStencil-1,iStencil), p(1:iStencil-1,iterm) ) ) / r(iStencil,iStencil);
        end
    end
    IS(2) = SI3(p(2,1:3));
    IS(3) = SI5(p(3,1:5));
    IS(4) = SI7(p(4,1:7));
    IS(5) = SI9(p(5,1:9));
    
    xi0 = ( fc(idx) - fc(idx-1) )^2;
    xi1 = ( fc(idx+1) - fc(idx) )^2;
    if xi0>=xi1
        r01 = 1;
    else
        r01 = 10;
    end
    r11 = 11 - r01;
    r01 = r01 / ( r01 + r11 );
    r11 = 1 - r01;
    sigma0 = r01 * ( 1 + abs(xi0-xi1)^kappa / ( xi0 + eps ) );
    sigma1 = r11 * ( 1 + abs(xi0-xi1)^kappa / ( xi1 + eps ) );
    sigma  = sigma0 + sigma1;
    IS(1) = ( sigma0 * ( fc(idx) - fc(idx-1) + sigma1 * ( fc(idx+1) - fc(idx) ) ) )^2 / sigma;
    
    tau = sum( abs( IS(nStencil) - IS(1:nStencil-1) ) / ( nStencil - 1 ) )^( nStencil - 1 );
    
    for iStencil = 1:nStencil
        w(iStencil) = r(iStencil,nStencil) .* ( 1 + tau / ( IS(iStencil) + eps ) );
    end
    w = w ./ sum(w);
    
    nterms = maxTerms;
    for iterm = 1:nterms
        p_final(iterm) = dot( w, p(:,iterm) );
    end
    fR(i+1) = polynomial(p_final,nterms,0.5);
    fL(i  ) = polynomial(p_final,nterms,-0.5);
end
fR(1) = fR(n);
% fL(n+1) = fL(2);

% Calculate tend
dfdx = ( fR(2:n+1) - fR(1:n) ) / dx;
% END of WENO

tend.PV    = -dfdx;
tend.PV    = u .* tend.PV;

% % True solution
% dfdx=cos(stat.x/100.*pi) * pi / 100.;

end

function stat = RK4(stat,dt)
weights = [1/6,1/3,1/3,1/6];

tend1 = spatial_discrete(stat);
stat2 = update_stat(stat,tend1,0.5*dt);

tend2 = spatial_discrete(stat2);
stat3 = update_stat(stat,tend2,0.5*dt);

tend3 = spatial_discrete(stat3);
stat4 = update_stat(stat,tend3,dt);

tend4 = spatial_discrete(stat4);

tend.PV  = weights(1).*tend1.PV  + weights(2).*tend2.PV  + weights(3).*tend3.PV  + weights(4).*tend4.PV;

stat = update_stat(stat,tend,dt);

end

function stat = update_stat(stat,tend,dt)
stat.f    = stat.f    + tend.PV  * dt;
end

function stat_out = cp_stat(stat_in)
stat_out.u        = stat_in.u;
stat_out.f        = stat_in.f;
stat_out.x        = stat_in.x;
stat_out.dx       = stat_in.dx;
stat_out.n        = stat_in.n;
stat_out.ht       = stat_in.ht;
stat_out.nt       = stat_in.nt;
stat_out.dt       = stat_in.dt;
stat_out.r        = stat_in.r;
stat_out.invA     = stat_in.invA;
stat_out.nStencil = stat_in.nStencil;
stat_out.maxTerms = stat_in.maxTerms;
end

function plot_result(stat,time_idx,output_path)
x_min = min(stat.x);
x_max = max(stat.x);

figure('Visible','off')
plot(stat.x,stat.f)
xlim([x_min,x_max])
ylim([-1,2.5])
title(['CWENO at ',num2str(time_idx*stat.ht*stat.dt),' second(s)'])
print(gcf,'-r400','-dpng',[output_path,'CWENO_',num2str(time_idx,'%.4d'),'.png']);
end

function poly = poly_integration(n,x_min,x_max)
% n is term number of the polynomial
% n = degree_of_polynomial + 1
% x_min, x_max are quadrature low/up bounds
poly = zeros(n,1);
for i = 1:n
    poly(i) = ( x_max^i - x_min^i ) / i;
end
end

function poly = polynomial(a,n,x)
% n is term number of the polynomial
% n = degree_of_polynomial + 1
% x is position of points
m = length(x);
poly = zeros(m,n);
for i = 1:n
    poly(:,i) = a(i) * x.^(i-1);
end
poly = sum(poly,2);
end

function SI = SI3(a)
SI = a(2)^2+(13*a(3)^2)/3;
end

function SI = SI5(a)
SI = a(2)^2 + (13*a(3)^2)/3 + (1/2)*a(2)*a(4) + (3129*a(4)^2)/80 + (21/5)*a(3)*a(5) + ...
  (87617*a(5)^2)/140;
end

function SI = SI7(a)
SI = a(2)^2 + (13*a(3)^2)/3 + (3129*a(4)^2)/80 + (21/5)*a(3)*a(5) + (87617*a(5)^2)/140 + ...
  (14127/224)*a(4)*a(6) + (252337135*a(6)^2)/16128 + (1/8)*a(2)*(4*a(4) + a(6)) + ...
  (87/56)*a(3)*a(7) + (508579/336)*a(5)*a(7) + (11102834003*a(7)^2)/19712;
end

function SI = SI9(a)
SI = a(2)^2 + (13*a(3)^2)/3 + (3129*a(4)^2)/80 + (87617*a(5)^2)/140 + (14127/224)*a(4)*a(6) + ...
  (252337135*a(6)^2)/16128 + (508579/336)*a(5)*a(7) + (11102834003*a(7)^2)/19712 + ...
  (12535/384)*a(4)*a(8) + (895099145*a(6)*a(8))/16896 + (16165726308907*a(8)^2)/585728 + ...
  (1/32)*a(2)*(16*a(4) + 4*a(6) + a(8)) + ...
  a(3)*((21*a(5))/5 + (87*a(7))/56 + (37*a(9))/72) + (551543/528)*a(5)*a(9) + ...
  (46545155573*a(7)*a(9))/18304 + (969943578534563*a(9)^2)/549120;
end