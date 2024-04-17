
% To get "geom_tab", must run the first part of Z_component.m or R_component.m or
% T_component.m before using this code


file_vel_drop = '50.796_44.746';
dt_syn = 5e-4;
dt_obs = 0.004;
T = 8.0;
nt = 2001;

%
sxpos = unique(geom_tab(:,6));
rxpos = unique(geom_tab(:,7));

% [shot#, trace#]
syn_dict = zeros(size(geom_tab,1),2);
syn_vx = zeros(nt,size(geom_tab,1));
syn_vz = zeros(nt,size(geom_tab,1));
for i=1:size(geom_tab,1)
    if mod(i,100)==0
        fprintf('finish trace # %d\n',i);
    end

    ishot = find(sxpos==geom_tab(i,6));
    itrace = find(rxpos==geom_tab(i,7));
    
    syn_dict(i,1) = ishot;
    syn_dict(i,2) = itrace;
    file_shot_vx = strcat('/home0/cxd170430/codes/FWI_CLDUAN_7_3_2023/yellowstone/new/gardner_rho/',file_vel_drop,'/syn_vx.su.shot',num2str(ishot));
    file_shot_vz = strcat('/home0/cxd170430/codes/FWI_CLDUAN_7_3_2023/yellowstone/new/gardner_rho/',file_vel_drop,'/syn_vy.su.shot',num2str(ishot));
    [vx_here,~,~] = ReadSu(file_shot_vx);
    vx_here_re = interp1(0:dt_syn:T-dt_syn, vx_here, 0:dt_obs:T, 'linear', 0);
    [vz_here,~,~] = ReadSu(file_shot_vz);
    vz_here_re = interp1(0:dt_syn:T-dt_syn, vz_here, 0:dt_obs:T, 'linear', 0);
    syn_vx(:,i) = vx_here_re(:,itrace);
    syn_vz(:,i) = vz_here_re(:,itrace);
end
figure;imagesc(syn_vx);clim([-1e-14,1e-14]);
figure;imagesc(syn_vz);clim([-1e-14,1e-14]);

% OPTIONAL...
output_vx = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vx_',file_vel_drop,'_32bit.bin');
output_vz = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vz_',file_vel_drop,'_32bit.bin');
fp=fopen(output_vx,'wb'); fwrite(fp,syn_vx,'float32'); fclose(fp);
fp=fopen(output_vz,'wb'); fwrite(fp,syn_vz,'float32'); fclose(fp);

% output_vx = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vx_',file_vel_drop,'_64bit.bin');
% output_vz = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vz_',file_vel_drop,'_64bit.bin');
% fp=fopen(output_vx,'wb'); fwrite(fp,syn_vx,'float64'); fclose(fp);
% fp=fopen(output_vz,'wb'); fwrite(fp,syn_vz,'float64'); fclose(fp);

% vx_here_re = interp1(0:dt_syn:T-dt_syn, vx_here, 0:dt_obs:T, 'linear', 0);
% figure;plot(0:dt_syn:T-dt_syn,vx_here(:,50),'b','linewidth',1.5);hold on;plot(0:dt_obs:T,vx_here_re(:,50),'r','linewidth',1.5);
% hold on;plot(0:dt_obs:T,vx_here_re1(:,50),'k','linewidth',1.5);