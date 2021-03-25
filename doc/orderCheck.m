clc
clear

time_start = 1;
time_end   = 25;

history_path = '..\run';

res_nc{1} = [history_path,'\','output.nc'];

% res_nc{2} = [history_path,'\','output_2p0.nc'];

% res_nc{1} = [history_path,'\','output_4p5.nc'];
% res_nc{2} = [history_path,'\','output_2p25.nc'];
% res_nc{3} = [history_path,'\','output_1p125.nc'];

% res_nc{2} = [history_path,'\','output_1p0.nc'];
% res_nc{2} = [history_path,'\','output_2p0.nc'];

% res_nc{1} = [history_path,'\','output_9p0.nc'];
% res_nc{2} = [history_path,'\','output_4p5.nc'];
% res_nc{3} = [history_path,'\','output_2p25.nc'];
% res_nc{4} = [history_path,'\','output_1p125.nc'];
% res_nc{5} = [history_path,'\','output_0p5625.nc'];

res_num = size(res_nc,2);

for ires = 1:res_num
    ids   = ncreadatt(res_nc{ires},'/','ids');
    ide   = ncreadatt(res_nc{ires},'/','ide');
    jds   = ncreadatt(res_nc{ires},'/','jds');
    jde   = ncreadatt(res_nc{ires},'/','jde');
    ifs   = ncreadatt(res_nc{ires},'/','ifs');
    ife   = ncreadatt(res_nc{ires},'/','ife');
    xhalo = ncreadatt(res_nc{ires},'/','xhalo');
    yhalo = ncreadatt(res_nc{ires},'/','yhalo');
    
    dx (ires) = ncreadatt(res_nc{ires},'/','dx' );
    dy (ires) = ncreadatt(res_nc{ires},'/','dy' );
    
    areaCell = ncread(res_nc{ires},'areaCell');
%     u_end    = ncread(res_nc{ires},'u'      ,[1,1,1,time_end  ],[Inf,Inf,Inf,1]);
%     v_end    = ncread(res_nc{ires},'v'      ,[1,1,1,time_end  ],[Inf,Inf,Inf,1]);
    phi_end  = ncread(res_nc{ires},'phi'    ,[1,1,1,time_end  ],[Inf,Inf,Inf,1]);
%     u_start  = ncread(res_nc{ires},'u'      ,[1,1,1,time_start],[Inf,Inf,Inf,1]);
%     v_start  = ncread(res_nc{ires},'v'      ,[1,1,1,time_start],[Inf,Inf,Inf,1]);
    phi_start= ncread(res_nc{ires},'phi'    ,[1,1,1,time_start],[Inf,Inf,Inf,1]);
    
    L1_phi  (ires) = L1  (phi_end(ids:ide,jds:jde,:),phi_start(ids:ide,jds:jde,:),areaCell(ids:ide,jds:jde,:));
    L2_phi  (ires) = L2  (phi_end(ids:ide,jds:jde,:),phi_start(ids:ide,jds:jde,:),areaCell(ids:ide,jds:jde,:));
    LInf_phi(ires) = LInf(phi_end(ids:ide,jds:jde,:),phi_start(ids:ide,jds:jde,:),areaCell(ids:ide,jds:jde,:));
end

for ires = 2:res_num
    L1order_phi  (ires) = log(L1_phi  (ires)/L1_phi  (ires-1))/log(dx(ires)/dx(ires-1));
    L2order_phi  (ires) = log(L2_phi  (ires)/L2_phi  (ires-1))/log(dx(ires)/dx(ires-1));
    LInforder_phi(ires) = log(LInf_phi(ires)/LInf_phi(ires-1))/log(dx(ires)/dx(ires-1));
end

function reslut = L1(field_model,field_ref,areaCell)

reslut = sum(sum(sum(abs(field_model - field_ref) .* areaCell)))...
       / sum(sum(sum(abs(field_ref) .* areaCell)));

end

function reslut = L2(field_model,field_ref,areaCell)

reslut = sqrt(sum(sum(sum((field_model - field_ref).^2 .* areaCell)))...
             /sum(sum(sum(field_ref.^2 .* areaCell))));

end

function reslut = LInf(field_model,field_ref,areaCell)

reslut = max(max(max(abs(field_model - field_ref) .* areaCell)))...
       / max(max(max(abs(field_ref) .* areaCell)));

end