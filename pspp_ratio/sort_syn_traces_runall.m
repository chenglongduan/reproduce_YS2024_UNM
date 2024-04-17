

%%
clear;
addpath('/home0/cxd170430/codes/matlab/segymat/');
addpath('/home0/cxd170430/codes/matlab/');
addpath('/home0/cxd170430/codes/matlab/util');

[data_Z_raw,header,~] = ReadSu('/data/Chenglong/Yellowstone_vibroseis/su_data_sta_lta/combine.Z.su');

%----remove bad stations which record only one trace----
header([8354,8011,8006,8181,8178]) = [];
data_Z_raw(:,[8354,8011,8006,8181,8178]) = [];
%-------------------------------------------------------

ntr = length(header);
nt = header(1).ns;
dt = header(1).dt * 10^(-6);

% build 2D SR geometry table
% (1) sx, (2) sy, (3) rx, (4) ry, (5) offset, (6) src_new_x, (7) recv_new_x, (8) midpt_new_x_pp
tmp_tab = zeros(ntr,8); %debug only: tmp_tab = zeros(ntr,9);count=0;
for i=1:ntr
    tmp_tab(i,1) = header(i).SourceX;
    tmp_tab(i,2) = header(i).SourceY;
    tmp_tab(i,3) = header(i).GroupX;
    tmp_tab(i,4) = header(i).GroupY;
    sz=header(i).SourceSurfaceElevation;
    rz=header(i).ReceiverGroupElevation;
    tmp_tab(i,5) = sqrt((tmp_tab(i,1)-tmp_tab(i,3))^2+(tmp_tab(i,2)-tmp_tab(i,4))^2+(sz-rz)^2);
    %debug only: count = count+1; tmp_tab(i,9) = count;
end

% sort by offset
[~,order_offset] = sort(tmp_tab(:,5));
geom_tab = tmp_tab(order_offset,:);
data_Z = data_Z_raw(:,order_offset);


% 2-D line projection
x_sr = [geom_tab(:,1);geom_tab(:,3)]; % ntr*2, upper src, lower recv
y_sr = [geom_tab(:,2);geom_tab(:,4)]; % ntr*2, upper src, lower recv

Y_LR = y_sr;
X_LR = [x_sr,ones(length(x_sr),1)];
KB = X_LR\Y_LR;

x_sr_proj = (x_sr+y_sr*KB(1)-KB(1)*KB(2))/(1+KB(1)*KB(1));
y_sr_proj = KB(1)*x_sr_proj+KB(2);

[~,iref] = min(x_sr_proj);
line_dist = zeros(2*ntr,1);
for i=1:2*ntr
    line_dist(i) = sqrt((x_sr_proj(i)-x_sr_proj(iref))^2 + (y_sr_proj(i)-y_sr_proj(iref))^2);
end
  
geom_tab(:,6) = line_dist(1:ntr);
geom_tab(:,7) = line_dist(ntr+1:end);
geom_tab(:,8) = 0.5*(geom_tab(:,6)+geom_tab(:,7));

file_vel_drop = '1.627_1.703';
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

output_vx = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vx_',file_vel_drop,'_32bit.bin');
output_vz = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vz_',file_vel_drop,'_32bit.bin');
fp=fopen(output_vx,'wb'); fwrite(fp,syn_vx,'float32'); fclose(fp);
fp=fopen(output_vz,'wb'); fwrite(fp,syn_vz,'float32'); fclose(fp);


%%
clear;
addpath('/home0/cxd170430/codes/matlab/segymat/');
addpath('/home0/cxd170430/codes/matlab/');
addpath('/home0/cxd170430/codes/matlab/util');

[data_Z_raw,header,~] = ReadSu('/data/Chenglong/Yellowstone_vibroseis/su_data_sta_lta/combine.Z.su');

%----remove bad stations which record only one trace----
header([8354,8011,8006,8181,8178]) = [];
data_Z_raw(:,[8354,8011,8006,8181,8178]) = [];
%-------------------------------------------------------

ntr = length(header);
nt = header(1).ns;
dt = header(1).dt * 10^(-6);

% build 2D SR geometry table
% (1) sx, (2) sy, (3) rx, (4) ry, (5) offset, (6) src_new_x, (7) recv_new_x, (8) midpt_new_x_pp
tmp_tab = zeros(ntr,8); %debug only: tmp_tab = zeros(ntr,9);count=0;
for i=1:ntr
    tmp_tab(i,1) = header(i).SourceX;
    tmp_tab(i,2) = header(i).SourceY;
    tmp_tab(i,3) = header(i).GroupX;
    tmp_tab(i,4) = header(i).GroupY;
    sz=header(i).SourceSurfaceElevation;
    rz=header(i).ReceiverGroupElevation;
    tmp_tab(i,5) = sqrt((tmp_tab(i,1)-tmp_tab(i,3))^2+(tmp_tab(i,2)-tmp_tab(i,4))^2+(sz-rz)^2);
    %debug only: count = count+1; tmp_tab(i,9) = count;
end

% sort by offset
[~,order_offset] = sort(tmp_tab(:,5));
geom_tab = tmp_tab(order_offset,:);
data_Z = data_Z_raw(:,order_offset);


% 2-D line projection
x_sr = [geom_tab(:,1);geom_tab(:,3)]; % ntr*2, upper src, lower recv
y_sr = [geom_tab(:,2);geom_tab(:,4)]; % ntr*2, upper src, lower recv

Y_LR = y_sr;
X_LR = [x_sr,ones(length(x_sr),1)];
KB = X_LR\Y_LR;

x_sr_proj = (x_sr+y_sr*KB(1)-KB(1)*KB(2))/(1+KB(1)*KB(1));
y_sr_proj = KB(1)*x_sr_proj+KB(2);

[~,iref] = min(x_sr_proj);
line_dist = zeros(2*ntr,1);
for i=1:2*ntr
    line_dist(i) = sqrt((x_sr_proj(i)-x_sr_proj(iref))^2 + (y_sr_proj(i)-y_sr_proj(iref))^2);
end
  
geom_tab(:,6) = line_dist(1:ntr);
geom_tab(:,7) = line_dist(ntr+1:end);
geom_tab(:,8) = 0.5*(geom_tab(:,6)+geom_tab(:,7));

file_vel_drop = '4.849_5.164';
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

output_vx = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vx_',file_vel_drop,'_32bit.bin');
output_vz = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vz_',file_vel_drop,'_32bit.bin');
fp=fopen(output_vx,'wb'); fwrite(fp,syn_vx,'float32'); fclose(fp);
fp=fopen(output_vz,'wb'); fwrite(fp,syn_vz,'float32'); fclose(fp);



%%
clear;
addpath('/home0/cxd170430/codes/matlab/segymat/');
addpath('/home0/cxd170430/codes/matlab/');
addpath('/home0/cxd170430/codes/matlab/util');

[data_Z_raw,header,~] = ReadSu('/data/Chenglong/Yellowstone_vibroseis/su_data_sta_lta/combine.Z.su');

%----remove bad stations which record only one trace----
header([8354,8011,8006,8181,8178]) = [];
data_Z_raw(:,[8354,8011,8006,8181,8178]) = [];
%-------------------------------------------------------

ntr = length(header);
nt = header(1).ns;
dt = header(1).dt * 10^(-6);

% build 2D SR geometry table
% (1) sx, (2) sy, (3) rx, (4) ry, (5) offset, (6) src_new_x, (7) recv_new_x, (8) midpt_new_x_pp
tmp_tab = zeros(ntr,8); %debug only: tmp_tab = zeros(ntr,9);count=0;
for i=1:ntr
    tmp_tab(i,1) = header(i).SourceX;
    tmp_tab(i,2) = header(i).SourceY;
    tmp_tab(i,3) = header(i).GroupX;
    tmp_tab(i,4) = header(i).GroupY;
    sz=header(i).SourceSurfaceElevation;
    rz=header(i).ReceiverGroupElevation;
    tmp_tab(i,5) = sqrt((tmp_tab(i,1)-tmp_tab(i,3))^2+(tmp_tab(i,2)-tmp_tab(i,4))^2+(sz-rz)^2);
    %debug only: count = count+1; tmp_tab(i,9) = count;
end

% sort by offset
[~,order_offset] = sort(tmp_tab(:,5));
geom_tab = tmp_tab(order_offset,:);
data_Z = data_Z_raw(:,order_offset);


% 2-D line projection
x_sr = [geom_tab(:,1);geom_tab(:,3)]; % ntr*2, upper src, lower recv
y_sr = [geom_tab(:,2);geom_tab(:,4)]; % ntr*2, upper src, lower recv

Y_LR = y_sr;
X_LR = [x_sr,ones(length(x_sr),1)];
KB = X_LR\Y_LR;

x_sr_proj = (x_sr+y_sr*KB(1)-KB(1)*KB(2))/(1+KB(1)*KB(1));
y_sr_proj = KB(1)*x_sr_proj+KB(2);

[~,iref] = min(x_sr_proj);
line_dist = zeros(2*ntr,1);
for i=1:2*ntr
    line_dist(i) = sqrt((x_sr_proj(i)-x_sr_proj(iref))^2 + (y_sr_proj(i)-y_sr_proj(iref))^2);
end
  
geom_tab(:,6) = line_dist(1:ntr);
geom_tab(:,7) = line_dist(ntr+1:end);
geom_tab(:,8) = 0.5*(geom_tab(:,6)+geom_tab(:,7));

file_vel_drop = '12.771_14.157';
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

output_vx = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vx_',file_vel_drop,'_32bit.bin');
output_vz = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vz_',file_vel_drop,'_32bit.bin');
fp=fopen(output_vx,'wb'); fwrite(fp,syn_vx,'float32'); fclose(fp);
fp=fopen(output_vz,'wb'); fwrite(fp,syn_vz,'float32'); fclose(fp);


%%
clear;
addpath('/home0/cxd170430/codes/matlab/segymat/');
addpath('/home0/cxd170430/codes/matlab/');
addpath('/home0/cxd170430/codes/matlab/util');

[data_Z_raw,header,~] = ReadSu('/data/Chenglong/Yellowstone_vibroseis/su_data_sta_lta/combine.Z.su');

%----remove bad stations which record only one trace----
header([8354,8011,8006,8181,8178]) = [];
data_Z_raw(:,[8354,8011,8006,8181,8178]) = [];
%-------------------------------------------------------

ntr = length(header);
nt = header(1).ns;
dt = header(1).dt * 10^(-6);

% build 2D SR geometry table
% (1) sx, (2) sy, (3) rx, (4) ry, (5) offset, (6) src_new_x, (7) recv_new_x, (8) midpt_new_x_pp
tmp_tab = zeros(ntr,8); %debug only: tmp_tab = zeros(ntr,9);count=0;
for i=1:ntr
    tmp_tab(i,1) = header(i).SourceX;
    tmp_tab(i,2) = header(i).SourceY;
    tmp_tab(i,3) = header(i).GroupX;
    tmp_tab(i,4) = header(i).GroupY;
    sz=header(i).SourceSurfaceElevation;
    rz=header(i).ReceiverGroupElevation;
    tmp_tab(i,5) = sqrt((tmp_tab(i,1)-tmp_tab(i,3))^2+(tmp_tab(i,2)-tmp_tab(i,4))^2+(sz-rz)^2);
    %debug only: count = count+1; tmp_tab(i,9) = count;
end

% sort by offset
[~,order_offset] = sort(tmp_tab(:,5));
geom_tab = tmp_tab(order_offset,:);
data_Z = data_Z_raw(:,order_offset);


% 2-D line projection
x_sr = [geom_tab(:,1);geom_tab(:,3)]; % ntr*2, upper src, lower recv
y_sr = [geom_tab(:,2);geom_tab(:,4)]; % ntr*2, upper src, lower recv

Y_LR = y_sr;
X_LR = [x_sr,ones(length(x_sr),1)];
KB = X_LR\Y_LR;

x_sr_proj = (x_sr+y_sr*KB(1)-KB(1)*KB(2))/(1+KB(1)*KB(1));
y_sr_proj = KB(1)*x_sr_proj+KB(2);

[~,iref] = min(x_sr_proj);
line_dist = zeros(2*ntr,1);
for i=1:2*ntr
    line_dist(i) = sqrt((x_sr_proj(i)-x_sr_proj(iref))^2 + (y_sr_proj(i)-y_sr_proj(iref))^2);
end
  
geom_tab(:,6) = line_dist(1:ntr);
geom_tab(:,7) = line_dist(ntr+1:end);
geom_tab(:,8) = 0.5*(geom_tab(:,6)+geom_tab(:,7));

file_vel_drop = '15.91_17.909';
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

output_vx = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vx_',file_vel_drop,'_32bit.bin');
output_vz = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vz_',file_vel_drop,'_32bit.bin');
fp=fopen(output_vx,'wb'); fwrite(fp,syn_vx,'float32'); fclose(fp);
fp=fopen(output_vz,'wb'); fwrite(fp,syn_vz,'float32'); fclose(fp);


%%
clear;
addpath('/home0/cxd170430/codes/matlab/segymat/');
addpath('/home0/cxd170430/codes/matlab/');
addpath('/home0/cxd170430/codes/matlab/util');

[data_Z_raw,header,~] = ReadSu('/data/Chenglong/Yellowstone_vibroseis/su_data_sta_lta/combine.Z.su');

%----remove bad stations which record only one trace----
header([8354,8011,8006,8181,8178]) = [];
data_Z_raw(:,[8354,8011,8006,8181,8178]) = [];
%-------------------------------------------------------

ntr = length(header);
nt = header(1).ns;
dt = header(1).dt * 10^(-6);

% build 2D SR geometry table
% (1) sx, (2) sy, (3) rx, (4) ry, (5) offset, (6) src_new_x, (7) recv_new_x, (8) midpt_new_x_pp
tmp_tab = zeros(ntr,8); %debug only: tmp_tab = zeros(ntr,9);count=0;
for i=1:ntr
    tmp_tab(i,1) = header(i).SourceX;
    tmp_tab(i,2) = header(i).SourceY;
    tmp_tab(i,3) = header(i).GroupX;
    tmp_tab(i,4) = header(i).GroupY;
    sz=header(i).SourceSurfaceElevation;
    rz=header(i).ReceiverGroupElevation;
    tmp_tab(i,5) = sqrt((tmp_tab(i,1)-tmp_tab(i,3))^2+(tmp_tab(i,2)-tmp_tab(i,4))^2+(sz-rz)^2);
    %debug only: count = count+1; tmp_tab(i,9) = count;
end

% sort by offset
[~,order_offset] = sort(tmp_tab(:,5));
geom_tab = tmp_tab(order_offset,:);
data_Z = data_Z_raw(:,order_offset);


% 2-D line projection
x_sr = [geom_tab(:,1);geom_tab(:,3)]; % ntr*2, upper src, lower recv
y_sr = [geom_tab(:,2);geom_tab(:,4)]; % ntr*2, upper src, lower recv

Y_LR = y_sr;
X_LR = [x_sr,ones(length(x_sr),1)];
KB = X_LR\Y_LR;

x_sr_proj = (x_sr+y_sr*KB(1)-KB(1)*KB(2))/(1+KB(1)*KB(1));
y_sr_proj = KB(1)*x_sr_proj+KB(2);

[~,iref] = min(x_sr_proj);
line_dist = zeros(2*ntr,1);
for i=1:2*ntr
    line_dist(i) = sqrt((x_sr_proj(i)-x_sr_proj(iref))^2 + (y_sr_proj(i)-y_sr_proj(iref))^2);
end
  
geom_tab(:,6) = line_dist(1:ntr);
geom_tab(:,7) = line_dist(ntr+1:end);
geom_tab(:,8) = 0.5*(geom_tab(:,6)+geom_tab(:,7));

file_vel_drop = '19.043_21.767';
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

output_vx = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vx_',file_vel_drop,'_32bit.bin');
output_vz = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vz_',file_vel_drop,'_32bit.bin');
fp=fopen(output_vx,'wb'); fwrite(fp,syn_vx,'float32'); fclose(fp);
fp=fopen(output_vz,'wb'); fwrite(fp,syn_vz,'float32'); fclose(fp);


%%
clear;
addpath('/home0/cxd170430/codes/matlab/segymat/');
addpath('/home0/cxd170430/codes/matlab/');
addpath('/home0/cxd170430/codes/matlab/util');

[data_Z_raw,header,~] = ReadSu('/data/Chenglong/Yellowstone_vibroseis/su_data_sta_lta/combine.Z.su');

%----remove bad stations which record only one trace----
header([8354,8011,8006,8181,8178]) = [];
data_Z_raw(:,[8354,8011,8006,8181,8178]) = [];
%-------------------------------------------------------

ntr = length(header);
nt = header(1).ns;
dt = header(1).dt * 10^(-6);

% build 2D SR geometry table
% (1) sx, (2) sy, (3) rx, (4) ry, (5) offset, (6) src_new_x, (7) recv_new_x, (8) midpt_new_x_pp
tmp_tab = zeros(ntr,8); %debug only: tmp_tab = zeros(ntr,9);count=0;
for i=1:ntr
    tmp_tab(i,1) = header(i).SourceX;
    tmp_tab(i,2) = header(i).SourceY;
    tmp_tab(i,3) = header(i).GroupX;
    tmp_tab(i,4) = header(i).GroupY;
    sz=header(i).SourceSurfaceElevation;
    rz=header(i).ReceiverGroupElevation;
    tmp_tab(i,5) = sqrt((tmp_tab(i,1)-tmp_tab(i,3))^2+(tmp_tab(i,2)-tmp_tab(i,4))^2+(sz-rz)^2);
    %debug only: count = count+1; tmp_tab(i,9) = count;
end

% sort by offset
[~,order_offset] = sort(tmp_tab(:,5));
geom_tab = tmp_tab(order_offset,:);
data_Z = data_Z_raw(:,order_offset);


% 2-D line projection
x_sr = [geom_tab(:,1);geom_tab(:,3)]; % ntr*2, upper src, lower recv
y_sr = [geom_tab(:,2);geom_tab(:,4)]; % ntr*2, upper src, lower recv

Y_LR = y_sr;
X_LR = [x_sr,ones(length(x_sr),1)];
KB = X_LR\Y_LR;

x_sr_proj = (x_sr+y_sr*KB(1)-KB(1)*KB(2))/(1+KB(1)*KB(1));
y_sr_proj = KB(1)*x_sr_proj+KB(2);

[~,iref] = min(x_sr_proj);
line_dist = zeros(2*ntr,1);
for i=1:2*ntr
    line_dist(i) = sqrt((x_sr_proj(i)-x_sr_proj(iref))^2 + (y_sr_proj(i)-y_sr_proj(iref))^2);
end
  
geom_tab(:,6) = line_dist(1:ntr);
geom_tab(:,7) = line_dist(ntr+1:end);
geom_tab(:,8) = 0.5*(geom_tab(:,6)+geom_tab(:,7));

file_vel_drop = '22.179_25.745';
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

output_vx = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vx_',file_vel_drop,'_32bit.bin');
output_vz = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vz_',file_vel_drop,'_32bit.bin');
fp=fopen(output_vx,'wb'); fwrite(fp,syn_vx,'float32'); fclose(fp);
fp=fopen(output_vz,'wb'); fwrite(fp,syn_vz,'float32'); fclose(fp);


%%
clear;
addpath('/home0/cxd170430/codes/matlab/segymat/');
addpath('/home0/cxd170430/codes/matlab/');
addpath('/home0/cxd170430/codes/matlab/util');

[data_Z_raw,header,~] = ReadSu('/data/Chenglong/Yellowstone_vibroseis/su_data_sta_lta/combine.Z.su');

%----remove bad stations which record only one trace----
header([8354,8011,8006,8181,8178]) = [];
data_Z_raw(:,[8354,8011,8006,8181,8178]) = [];
%-------------------------------------------------------

ntr = length(header);
nt = header(1).ns;
dt = header(1).dt * 10^(-6);

% build 2D SR geometry table
% (1) sx, (2) sy, (3) rx, (4) ry, (5) offset, (6) src_new_x, (7) recv_new_x, (8) midpt_new_x_pp
tmp_tab = zeros(ntr,8); %debug only: tmp_tab = zeros(ntr,9);count=0;
for i=1:ntr
    tmp_tab(i,1) = header(i).SourceX;
    tmp_tab(i,2) = header(i).SourceY;
    tmp_tab(i,3) = header(i).GroupX;
    tmp_tab(i,4) = header(i).GroupY;
    sz=header(i).SourceSurfaceElevation;
    rz=header(i).ReceiverGroupElevation;
    tmp_tab(i,5) = sqrt((tmp_tab(i,1)-tmp_tab(i,3))^2+(tmp_tab(i,2)-tmp_tab(i,4))^2+(sz-rz)^2);
    %debug only: count = count+1; tmp_tab(i,9) = count;
end

% sort by offset
[~,order_offset] = sort(tmp_tab(:,5));
geom_tab = tmp_tab(order_offset,:);
data_Z = data_Z_raw(:,order_offset);


% 2-D line projection
x_sr = [geom_tab(:,1);geom_tab(:,3)]; % ntr*2, upper src, lower recv
y_sr = [geom_tab(:,2);geom_tab(:,4)]; % ntr*2, upper src, lower recv

Y_LR = y_sr;
X_LR = [x_sr,ones(length(x_sr),1)];
KB = X_LR\Y_LR;

x_sr_proj = (x_sr+y_sr*KB(1)-KB(1)*KB(2))/(1+KB(1)*KB(1));
y_sr_proj = KB(1)*x_sr_proj+KB(2);

[~,iref] = min(x_sr_proj);
line_dist = zeros(2*ntr,1);
for i=1:2*ntr
    line_dist(i) = sqrt((x_sr_proj(i)-x_sr_proj(iref))^2 + (y_sr_proj(i)-y_sr_proj(iref))^2);
end
  
geom_tab(:,6) = line_dist(1:ntr);
geom_tab(:,7) = line_dist(ntr+1:end);
geom_tab(:,8) = 0.5*(geom_tab(:,6)+geom_tab(:,7));

file_vel_drop = '25.327_29.864';
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

output_vx = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vx_',file_vel_drop,'_32bit.bin');
output_vz = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vz_',file_vel_drop,'_32bit.bin');
fp=fopen(output_vx,'wb'); fwrite(fp,syn_vx,'float32'); fclose(fp);
fp=fopen(output_vz,'wb'); fwrite(fp,syn_vz,'float32'); fclose(fp);


%%
clear;
addpath('/home0/cxd170430/codes/matlab/segymat/');
addpath('/home0/cxd170430/codes/matlab/');
addpath('/home0/cxd170430/codes/matlab/util');

[data_Z_raw,header,~] = ReadSu('/data/Chenglong/Yellowstone_vibroseis/su_data_sta_lta/combine.Z.su');

%----remove bad stations which record only one trace----
header([8354,8011,8006,8181,8178]) = [];
data_Z_raw(:,[8354,8011,8006,8181,8178]) = [];
%-------------------------------------------------------

ntr = length(header);
nt = header(1).ns;
dt = header(1).dt * 10^(-6);

% build 2D SR geometry table
% (1) sx, (2) sy, (3) rx, (4) ry, (5) offset, (6) src_new_x, (7) recv_new_x, (8) midpt_new_x_pp
tmp_tab = zeros(ntr,8); %debug only: tmp_tab = zeros(ntr,9);count=0;
for i=1:ntr
    tmp_tab(i,1) = header(i).SourceX;
    tmp_tab(i,2) = header(i).SourceY;
    tmp_tab(i,3) = header(i).GroupX;
    tmp_tab(i,4) = header(i).GroupY;
    sz=header(i).SourceSurfaceElevation;
    rz=header(i).ReceiverGroupElevation;
    tmp_tab(i,5) = sqrt((tmp_tab(i,1)-tmp_tab(i,3))^2+(tmp_tab(i,2)-tmp_tab(i,4))^2+(sz-rz)^2);
    %debug only: count = count+1; tmp_tab(i,9) = count;
end

% sort by offset
[~,order_offset] = sort(tmp_tab(:,5));
geom_tab = tmp_tab(order_offset,:);
data_Z = data_Z_raw(:,order_offset);


% 2-D line projection
x_sr = [geom_tab(:,1);geom_tab(:,3)]; % ntr*2, upper src, lower recv
y_sr = [geom_tab(:,2);geom_tab(:,4)]; % ntr*2, upper src, lower recv

Y_LR = y_sr;
X_LR = [x_sr,ones(length(x_sr),1)];
KB = X_LR\Y_LR;

x_sr_proj = (x_sr+y_sr*KB(1)-KB(1)*KB(2))/(1+KB(1)*KB(1));
y_sr_proj = KB(1)*x_sr_proj+KB(2);

[~,iref] = min(x_sr_proj);
line_dist = zeros(2*ntr,1);
for i=1:2*ntr
    line_dist(i) = sqrt((x_sr_proj(i)-x_sr_proj(iref))^2 + (y_sr_proj(i)-y_sr_proj(iref))^2);
end
  
geom_tab(:,6) = line_dist(1:ntr);
geom_tab(:,7) = line_dist(ntr+1:end);
geom_tab(:,8) = 0.5*(geom_tab(:,6)+geom_tab(:,7));

file_vel_drop = '31.693_38.644';
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

output_vx = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vx_',file_vel_drop,'_32bit.bin');
output_vz = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vz_',file_vel_drop,'_32bit.bin');
fp=fopen(output_vx,'wb'); fwrite(fp,syn_vx,'float32'); fclose(fp);
fp=fopen(output_vz,'wb'); fwrite(fp,syn_vz,'float32'); fclose(fp);


%%
clear;
addpath('/home0/cxd170430/codes/matlab/segymat/');
addpath('/home0/cxd170430/codes/matlab/');
addpath('/home0/cxd170430/codes/matlab/util');

[data_Z_raw,header,~] = ReadSu('/data/Chenglong/Yellowstone_vibroseis/su_data_sta_lta/combine.Z.su');

%----remove bad stations which record only one trace----
header([8354,8011,8006,8181,8178]) = [];
data_Z_raw(:,[8354,8011,8006,8181,8178]) = [];
%-------------------------------------------------------

ntr = length(header);
nt = header(1).ns;
dt = header(1).dt * 10^(-6);

% build 2D SR geometry table
% (1) sx, (2) sy, (3) rx, (4) ry, (5) offset, (6) src_new_x, (7) recv_new_x, (8) midpt_new_x_pp
tmp_tab = zeros(ntr,8); %debug only: tmp_tab = zeros(ntr,9);count=0;
for i=1:ntr
    tmp_tab(i,1) = header(i).SourceX;
    tmp_tab(i,2) = header(i).SourceY;
    tmp_tab(i,3) = header(i).GroupX;
    tmp_tab(i,4) = header(i).GroupY;
    sz=header(i).SourceSurfaceElevation;
    rz=header(i).ReceiverGroupElevation;
    tmp_tab(i,5) = sqrt((tmp_tab(i,1)-tmp_tab(i,3))^2+(tmp_tab(i,2)-tmp_tab(i,4))^2+(sz-rz)^2);
    %debug only: count = count+1; tmp_tab(i,9) = count;
end

% sort by offset
[~,order_offset] = sort(tmp_tab(:,5));
geom_tab = tmp_tab(order_offset,:);
data_Z = data_Z_raw(:,order_offset);


% 2-D line projection
x_sr = [geom_tab(:,1);geom_tab(:,3)]; % ntr*2, upper src, lower recv
y_sr = [geom_tab(:,2);geom_tab(:,4)]; % ntr*2, upper src, lower recv

Y_LR = y_sr;
X_LR = [x_sr,ones(length(x_sr),1)];
KB = X_LR\Y_LR;

x_sr_proj = (x_sr+y_sr*KB(1)-KB(1)*KB(2))/(1+KB(1)*KB(1));
y_sr_proj = KB(1)*x_sr_proj+KB(2);

[~,iref] = min(x_sr_proj);
line_dist = zeros(2*ntr,1);
for i=1:2*ntr
    line_dist(i) = sqrt((x_sr_proj(i)-x_sr_proj(iref))^2 + (y_sr_proj(i)-y_sr_proj(iref))^2);
end
  
geom_tab(:,6) = line_dist(1:ntr);
geom_tab(:,7) = line_dist(ntr+1:end);
geom_tab(:,8) = 0.5*(geom_tab(:,6)+geom_tab(:,7));

file_vel_drop = '39.897_51.16';
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

output_vx = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vx_',file_vel_drop,'_32bit.bin');
output_vz = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vz_',file_vel_drop,'_32bit.bin');
fp=fopen(output_vx,'wb'); fwrite(fp,syn_vx,'float32'); fclose(fp);
fp=fopen(output_vz,'wb'); fwrite(fp,syn_vz,'float32'); fclose(fp);


%%
clear;
addpath('/home0/cxd170430/codes/matlab/segymat/');
addpath('/home0/cxd170430/codes/matlab/');
addpath('/home0/cxd170430/codes/matlab/util');

[data_Z_raw,header,~] = ReadSu('/data/Chenglong/Yellowstone_vibroseis/su_data_sta_lta/combine.Z.su');

%----remove bad stations which record only one trace----
header([8354,8011,8006,8181,8178]) = [];
data_Z_raw(:,[8354,8011,8006,8181,8178]) = [];
%-------------------------------------------------------

ntr = length(header);
nt = header(1).ns;
dt = header(1).dt * 10^(-6);

% build 2D SR geometry table
% (1) sx, (2) sy, (3) rx, (4) ry, (5) offset, (6) src_new_x, (7) recv_new_x, (8) midpt_new_x_pp
tmp_tab = zeros(ntr,8); %debug only: tmp_tab = zeros(ntr,9);count=0;
for i=1:ntr
    tmp_tab(i,1) = header(i).SourceX;
    tmp_tab(i,2) = header(i).SourceY;
    tmp_tab(i,3) = header(i).GroupX;
    tmp_tab(i,4) = header(i).GroupY;
    sz=header(i).SourceSurfaceElevation;
    rz=header(i).ReceiverGroupElevation;
    tmp_tab(i,5) = sqrt((tmp_tab(i,1)-tmp_tab(i,3))^2+(tmp_tab(i,2)-tmp_tab(i,4))^2+(sz-rz)^2);
    %debug only: count = count+1; tmp_tab(i,9) = count;
end

% sort by offset
[~,order_offset] = sort(tmp_tab(:,5));
geom_tab = tmp_tab(order_offset,:);
data_Z = data_Z_raw(:,order_offset);


% 2-D line projection
x_sr = [geom_tab(:,1);geom_tab(:,3)]; % ntr*2, upper src, lower recv
y_sr = [geom_tab(:,2);geom_tab(:,4)]; % ntr*2, upper src, lower recv

Y_LR = y_sr;
X_LR = [x_sr,ones(length(x_sr),1)];
KB = X_LR\Y_LR;

x_sr_proj = (x_sr+y_sr*KB(1)-KB(1)*KB(2))/(1+KB(1)*KB(1));
y_sr_proj = KB(1)*x_sr_proj+KB(2);

[~,iref] = min(x_sr_proj);
line_dist = zeros(2*ntr,1);
for i=1:2*ntr
    line_dist(i) = sqrt((x_sr_proj(i)-x_sr_proj(iref))^2 + (y_sr_proj(i)-y_sr_proj(iref))^2);
end
  
geom_tab(:,6) = line_dist(1:ntr);
geom_tab(:,7) = line_dist(ntr+1:end);
geom_tab(:,8) = 0.5*(geom_tab(:,6)+geom_tab(:,7));

file_vel_drop = '48.601_67.232';
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

output_vx = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vx_',file_vel_drop,'_32bit.bin');
output_vz = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vz_',file_vel_drop,'_32bit.bin');
fp=fopen(output_vx,'wb'); fwrite(fp,syn_vx,'float32'); fclose(fp);
fp=fopen(output_vz,'wb'); fwrite(fp,syn_vz,'float32'); fclose(fp);






