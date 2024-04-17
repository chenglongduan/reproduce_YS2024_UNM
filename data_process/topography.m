clear;

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
tmp_tab = zeros(ntr,8);
for i=1:ntr
    tmp_tab(i,1) = header(i).SourceX;
    tmp_tab(i,2) = header(i).SourceY;
    tmp_tab(i,3) = header(i).GroupX;
    tmp_tab(i,4) = header(i).GroupY;
    sz=header(i).SourceSurfaceElevation;
    rz=header(i).ReceiverGroupElevation;
    tmp_tab(i,5) = sz/1000; %sqrt((tmp_tab(i,1)-tmp_tab(i,3))^2+(tmp_tab(i,2)-tmp_tab(i,4))^2+(sz-rz)^2)
    tmp_tab(i,8) = rz/1000; %header(i).ReceiverGroupElevation/1000
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
  
geom_tab(:,6) = line_dist(1:ntr)/1000;
geom_tab(:,7) = line_dist(ntr+1:end)/1000;

figure;plot(geom_tab(:,6),zeros(ntr,1),'r*');hold on;plot(geom_tab(:,7),ones(ntr,1)*0.2,'bo');
ylim([-1,1]);set(gca,'Ydir','reverse');xlabel('Distance (m)');legend('source','receiver');
title('2D-projected sources, receivers');

% extract
tmp_rx_rz = geom_tab(:,7:8);
[~,order_rx] = sort(tmp_rx_rz(:,1));
geom_rx_rz = tmp_rx_rz(order_rx,:);

% unique
[~,order_uniq,~] = unique(geom_rx_rz(:,1),'stable');
geom_rx_rz_uniq = geom_rx_rz(order_uniq,:);
geom_rx_rz_uniq(266,:)=[30,geom_rx_rz_uniq(265,2)];

figure;plot(geom_rx_rz_uniq(:,1),geom_rx_rz_uniq(:,2));

% 1D interp
rx_grid = 0:0.01:30;
rz_grid = interp1(geom_rx_rz_uniq(:,1),geom_rx_rz_uniq(:,2),rx_grid,'spline');

figure;plot(geom_rx_rz_uniq(:,1),geom_rx_rz_uniq(:,2));hold on;plot(rx_grid,rz_grid);
figure;plot(rx_grid,rz_grid);ylim([0,10]);

topo = [rx_grid',rz_grid'];
fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/topo.bin','wb');
fwrite(fp,topo,'float64');
fclose(fp);

% 1D smooth
rz_grid_smooth = movmean(rz_grid, 20); % mean
figure;plot(rx_grid,rz_grid_smooth);ylim([0,10]);
figure;plot(geom_rx_rz_uniq(:,1),geom_rx_rz_uniq(:,2));hold on;plot(rx_grid,rz_grid_smooth);

%% To extract discrete source and receiver topography
% Run Lines 1 - 56 first, and then run below:
tmp_rx_rz = geom_tab(:,7:8);
[~,order_rx] = sort(tmp_rx_rz(:,1));
geom_rx_rz = tmp_rx_rz(order_rx,:);
[~,order_uniq,~] = unique(geom_rx_rz(:,1),'stable');
geom_rx_rz_uniq = geom_rx_rz(order_uniq,:);

tmp_sx_sz = [geom_tab(:,6),geom_tab(:,5)];
[~,order_sx] = sort(tmp_sx_sz(:,1));
geom_sx_sz = tmp_sx_sz(order_sx,:);
[~,order_uniq_s,~] = unique(geom_sx_sz(:,1),'stable');
geom_sx_sz_uniq = geom_sx_sz(order_uniq_s,:);

%
writematrix(geom_sx_sz_uniq,'/home0/cxd170430/codes/matlab/yellowstone_imaging/output/topo_src.txt','Delimiter','space');
writematrix(geom_rx_rz_uniq,'/home0/cxd170430/codes/matlab/yellowstone_imaging/output/topo_recv.txt','Delimiter','space');


