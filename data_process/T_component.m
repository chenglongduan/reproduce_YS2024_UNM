%% T-component data reflection imaging

% Run '01_segy2su.sh' and '02_combine2one.sh' before running this one

clear;
addpath('/home0/cxd170430/codes/matlab/segymat/');
addpath('/home0/cxd170430/codes/matlab/');
addpath('/home0/cxd170430/codes/matlab/util');


%% T component

[data_T_raw,header,~] = ReadSu('/data/Chenglong/Yellowstone_vibroseis/su_data_sta_lta/combine.T.su');

%----remove bad stations which record only one trace----
header([8354,8011,8006,8181,8178]) = [];
data_T_raw(:,[8354,8011,8006,8181,8178]) = [];
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
    tmp_tab(i,5) = sqrt((tmp_tab(i,1)-tmp_tab(i,3))^2+(tmp_tab(i,2)-tmp_tab(i,4))^2+(sz-rz)^2);
end

% sort by offset
[~,order_offset] = sort(tmp_tab(:,5));
geom_tab = tmp_tab(order_offset,:);
data_T = data_T_raw(:,order_offset);


figure;imagesc(data_T);colormap(colorbar_bwr);colorbar;clim([0.5,2.5]);


% 2-D line projection
x_sr = [geom_tab(:,1);geom_tab(:,3)]; % ntr*2, upper src, lower recv
y_sr = [geom_tab(:,2);geom_tab(:,4)]; % ntr*2, upper src, lower recv

Y_LR = y_sr;
X_LR = [x_sr,ones(length(x_sr),1)];
KB = X_LR\Y_LR;

x_sr_proj = (x_sr+y_sr*KB(1)-KB(1)*KB(2))/(1+KB(1)*KB(1));
y_sr_proj = KB(1)*x_sr_proj+KB(2);

figure;plot(geom_tab(:,1),geom_tab(:,2),'r*');hold on;plot(geom_tab(:,3),geom_tab(:,4),'bo');hold on;
plot(x_sr_proj(1:ntr),y_sr_proj(1:ntr),'ro');hold on;plot(x_sr_proj(ntr+1:end),y_sr_proj(ntr+1:end),'ko');
axis equal;

[~,iref] = min(x_sr_proj);
line_dist = zeros(2*ntr,1);
for i=1:2*ntr
    line_dist(i) = sqrt((x_sr_proj(i)-x_sr_proj(iref))^2 + (y_sr_proj(i)-y_sr_proj(iref))^2);
end
  
geom_tab(:,6) = line_dist(1:ntr);
geom_tab(:,7) = line_dist(ntr+1:end);
geom_tab(:,8) = 0.5*(geom_tab(:,6)+geom_tab(:,7));

figure;plot(geom_tab(:,6),zeros(ntr,1),'r*');hold on;plot(geom_tab(:,7),ones(ntr,1)*0.2,'bo');hold on;
plot(geom_tab(:,8),ones(ntr,1)*0.4,'mo');
ylim([-1,1]);set(gca,'Ydir','reverse');xlabel('Distance (m)');legend('source','receiver','PP reflection point');
title('2D-projected sources, receivers, PP reflection points');


%% mute
% (x,y) points of the mute polygons
%muteD_ipt1=[1,110]; muteD_ipt2=[1667,270]; muteD_ipt3=[11218,900]; %780cut direct
muteD_ipt1=[1,110]; muteD_ipt2=[1667,480]; muteD_ipt3=[11218,900]; %cut direct+pp reflection
muteA_ipt1=[1,40]; muteA_ipt2=[3560,2001];
% assign (index) offset range
ioffset_range = []; % 2000:11218
% assign clip value
clip_value = 8.0; % 10.0
% read stack data
T_stack = read_stack('/data/Chenglong/Yellowstone_vibroseis/sta_lta_binstack/raw/EHT.stalta.txt',2001,100);
% value for muted part
value4muted = 0; % 0 or NaN, where 0 is for adjoint, NaN is for kirchhoff

figure;subplot(1,2,1);
imagesc(10:100:10000,1:2001,T_stack);colormap(colorbar_bwr);clim([0.5,1.5]);%2.5
hold on;plot([geom_tab(muteD_ipt1(1),5),geom_tab(muteD_ipt2(1),5)],[muteD_ipt1(2),muteD_ipt2(2)],'k-');
hold on;plot([geom_tab(muteD_ipt2(1),5),geom_tab(muteD_ipt3(1),5)],[muteD_ipt2(2),muteD_ipt3(2)],'k-');
hold on;plot([geom_tab(muteA_ipt1(1),5),geom_tab(muteA_ipt2(1),5)],[muteA_ipt1(2),muteA_ipt2(2)],'k-');title('mute: stacked T');
subplot(1,2,2);
imagesc(data_T);colormap(colorbar_bwr);clim([0.5,clip_value]);
hold on;plot([muteD_ipt1(1),muteD_ipt2(1)],[muteD_ipt1(2),muteD_ipt2(2)],'k-');
hold on;plot([muteD_ipt2(1),muteD_ipt3(1)],[muteD_ipt2(2),muteD_ipt3(2)],'k-');
hold on;plot([muteA_ipt1(1),muteA_ipt2(1)],[muteA_ipt1(2),muteA_ipt2(2)],'k-');title('mute: raw T');
%**************************************************************************


% do not change below!!!
it_muteD = floor(muteD_ipt1(2)+(muteD_ipt2(2)-muteD_ipt1(2))*((muteD_ipt1(1):muteD_ipt2(1))-muteD_ipt1(1))/(muteD_ipt2(1)-muteD_ipt1(1)));
it_muteD = [it_muteD(1:end-1),floor(muteD_ipt2(2)+(muteD_ipt3(2)-muteD_ipt2(2))*((muteD_ipt2(1):muteD_ipt3(1))-muteD_ipt2(1))/(muteD_ipt3(1)-muteD_ipt2(1)))];
data_T_muteD = data_T;
for i=1:length(it_muteD)
    data_T_muteD(1:it_muteD(i),muteD_ipt1(1)-1+i)=value4muted;
end
it_muteA = floor(muteA_ipt1(2)+(muteA_ipt2(2)-muteA_ipt1(2))*((muteA_ipt1(1):muteA_ipt2(1))-muteA_ipt1(1))/(muteA_ipt2(1)-muteA_ipt1(1)));
data_T_muteA = data_T_muteD;
for i=1:length(it_muteA)
    data_T_muteA(it_muteA(i):nt,muteA_ipt1(1)-1+i)=value4muted;
end
if ~isempty(ioffset_range)
    data_T_muteA(:,setdiff(1:ntr,ioffset_range))=value4muted;
end
if ~isempty(clip_value)
    data_T_muteA(data_T_muteA>clip_value)=clip_value;
end

figure;imagesc(data_T_muteA);colormap(colorbar_bwr);colorbar;title('final input data for imaging');


%% (Optional) generate adjoint imaging files
src_filename = '/home0/cxd170430/codes/ACOUSTIC2D_RTM/yellowstone_field_data/source_file.txt';
rec_dir = '/home0/cxd170430/codes/ACOUSTIC2D_RTM/yellowstone_field_data/receiver/';
dat_dir = '/home0/cxd170430/codes/ACOUSTIC2D_RTM/yellowstone_field_data/data/stalta_t_ps/';
tshift = 0.1;

%cell_rec_dat = extract_shot_gathers(data_T_muteA, geom_tab, src_filename, rec_dir, dat_dir, tshift, dt);
cell_rec_dat = extract_shot_gathers(data_T_muteA, geom_tab, '', '', dat_dir, tshift, dt);


%% PP PS imaging

data_mat = data_T_muteA;
nx = 81; %121  81   61
dx = 375; %250  375  500
nz = 201;
dz = 50;
dt1 = dt;
tshift = 0.0;
vp = 4500;
vs = 2500;
snr_qc = 1;
snr_threshold = 0.9;
su_file_env='/data/Chenglong/Yellowstone_vibroseis/su_data_envelope/combine.Z.su';

% (optional)
if snr_qc==0
    [A_all_T, A_each_T, data_mat_sg] = PP_imaging_srcloop(data_mat, geom_tab, vp, nx, dx, nz, dz, dt, tshift);
else
    data_mat_qc = SNR(su_file_env, data_mat, 'Z', snr_threshold);
    [A_all_T, A_each_T, data_mat_sg] = PP_imaging_srcloop(data_mat_qc, geom_tab, vp, nx, dx, nz, dz, dt, tshift);
end

fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/img375_pp_t_nosig.bin','wb');
fwrite(fp,A_all_T,'float64');fclose(fp);


% (required)
if snr_qc==0
    [A1_all_T, A1_each_T, ccp_x_section] = PS_imaging_srcloop(data_mat, geom_tab, vp, vs, nx, dx, nz, dz, dt, tshift);
else
    data_mat_qc = SNR(su_file_env, data_mat, 'Z', snr_threshold);
    [A1_all_T, A1_each_T, ccp_x_section] = PS_imaging_srcloop(data_mat_qc, geom_tab, vp, vs, nx, dx, nz, dz, dt, tshift);
end

fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/img375_ps_t.bin','wb');
fwrite(fp,A1_all_T,'float64');fclose(fp);


%[~,~] = PP_imaging_traceloop(data_mat, 'stalta', geom_tab, vp, nx, dx, nz, dz, dt1, tshift);
%[~,~,~] = PS_imaging_traceloop(data_mat, 'stalta', 0, geom_tab, vp, vs, nx, dx, nz, dz, dt1, tshift);

%figure;imagesc(PPsum_z,'AlphaData',~isnan(PPsum_z));colormap jet;colorbar;


%*******************imaging condition*************************************
axis_x = 0:dx/1000:(nx-1)*dx/1000; axis_z = 0:dz/1000:(nz-1)*dz/1000;
figure;imagesc(axis_x,axis_z,A_all_T+A1_all_T);colormap jet;colorbar;title('PP + PS');xlabel('Distance (km)');ylabel('Depth (km)');
figure;imagesc(axis_x,axis_z,A_all_T.*A1_all_T);colormap jet;colorbar;title('PP * PS');xlabel('Distance (km)');ylabel('Depth (km)');
%*************************************************************************



%========================plot shot gathers=============================
isrc_plot = 10:10:60;
trace_num = [];
for i=1:length(isrc_plot)
    isrc = isrc_plot(i);
    if isrc==1
        itr_beg = 1;
        itr_end = ntr_vector(isrc);
    else
        itr_beg = sum(ntr_vector(1:isrc-1))+1;
        itr_end = sum(ntr_vector(1:isrc));
    end
    trace_num = [trace_num, itr_beg:itr_end];
end
figure;imagesc(data_mat_sg(:,trace_num));colormap(colorbar_bwr);clim([0,4]);title('source gathers');
%=======================================================================





