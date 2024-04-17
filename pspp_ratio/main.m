
clear;
addpath('/home0/cxd170430/codes/matlab/segymat/');
addpath('/home0/cxd170430/codes/matlab/');
addpath('/home0/cxd170430/codes/matlab/util');
addpath('/home0/cxd170430/codes/matlab/yellowstone_imaging/pspp_ratio_matlab');
addpath('/home0/cxd170430/codes/matlab/yellowstone_imaging');

%% 1. load envelope data, do SNR and offset selection (observed)

su_fileZ_env='/data/Chenglong/Yellowstone_vibroseis/su_data_envelope/combine.Z.su';
su_fileR_env='/data/Chenglong/Yellowstone_vibroseis/su_data_envelope/combine.R.su';
su_fileT_env='/data/Chenglong/Yellowstone_vibroseis/su_data_envelope/combine.T.su';
su_fileRT_env='/data/Chenglong/Yellowstone_vibroseis/su_data_envelope/combine.RT.su';
% SNR params
yn_snr_qc = 1;
snr_threshold = 1.0;
% offset range (meter)
offset_min = 2200; % 2000, 4000
offset_max = 5200; % 5200,5500 5250 6200

[Z_env_new, R_env_new, T_env_new, RT_env_new, geom_tab_new] = ...
load_env_so(su_fileZ_env, su_fileR_env, su_fileT_env, su_fileRT_env, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);

ZRT_env_new = Z_env_new + R_env_new + T_env_new;

%% 2. load stacked STA/LTA data, do +0.1s time correction (observed)

Z_stack = read_stack('/data/Chenglong/Yellowstone_vibroseis/sta_lta_binstack/raw/EHZ.stalta.txt',2001,100);
R_stack = read_stack('/data/Chenglong/Yellowstone_vibroseis/sta_lta_binstack/raw/EHR.stalta.txt',2001,100);
T_stack = read_stack('/data/Chenglong/Yellowstone_vibroseis/sta_lta_binstack/raw/EHT.stalta.txt',2001,100);
RT_stack = read_stack('/data/Chenglong/Yellowstone_vibroseis/sta_lta_binstack/raw/EHRT.stalta.txt',2001,100);

tshift = 0.1;
nt = 2001;
dt = 0.004;
ntr_stack_stalta = 100;
if tshift~=0
    itshift = round(abs(tshift)/dt);
    data_tmp1 = Z_stack;
    Z_stack = zeros(nt,ntr_stack_stalta);
    data_tmp2 = R_stack;
    R_stack = zeros(nt,ntr_stack_stalta);
    data_tmp3 = T_stack;
    T_stack = zeros(nt,ntr_stack_stalta);
    data_tmp4 = RT_stack;
    RT_stack = zeros(nt,ntr_stack_stalta);
    if(tshift>0)
        Z_stack(itshift+1:nt,:) = data_tmp1(1:nt-itshift,:);
        R_stack(itshift+1:nt,:) = data_tmp2(1:nt-itshift,:);
        T_stack(itshift+1:nt,:) = data_tmp3(1:nt-itshift,:);
        RT_stack(itshift+1:nt,:) = data_tmp4(1:nt-itshift,:);
    else
        Z_stack(1:nt-itshift,:) = data_tmp1(itshift+1:nt,:);
        R_stack(1:nt-itshift,:) = data_tmp2(itshift+1:nt,:);
        T_stack(1:nt-itshift,:) = data_tmp3(itshift+1:nt,:);
        RT_stack(1:nt-itshift,:) = data_tmp4(itshift+1:nt,:);
    end
end

figure;imagesc(0.1:0.1:10,0:dt:(nt-1)*dt,Z_stack);title('Z1');colormap(colorbar_bwr);clim([0.5,1.5]);
xlabel('offset (km)');ylabel('time (s)');hold on;
plot([offset_min,offset_min]/1000,[0,8],'k--');hold on;plot([offset_max,offset_max]/1000,[0,8],'k--');
figure;imagesc(0.1:0.1:10,0:dt:(nt-1)*dt,R_stack);title('R1');colormap(colorbar_bwr);clim([0.5,1.5]);
xlabel('offset (km)');ylabel('time (s)');hold on;
plot([offset_min,offset_min]/1000,[0,8],'k--');hold on;plot([offset_max,offset_max]/1000,[0,8],'k--');
figure;imagesc(0.1:0.1:10,0:dt:(nt-1)*dt,T_stack);title('T1');colormap(colorbar_bwr);clim([0.5,1.5]);
xlabel('offset (km)');ylabel('time (s)');hold on;
plot([offset_min,offset_min]/1000,[0,8],'k--');hold on;plot([offset_max,offset_max]/1000,[0,8],'k--');


%% 3. compute predicted PP, PS arrival times
% compute traveltimes based on selected envelope traces
% and then project them to stacked stalta traces based on the offset values

vp = 4500;
vs = 2500;
z = 3890;
timewin_p_prior = -0.1;
timewin_p_post = 0.596;
timewin_pp_prior = 0.0; % 0.05
timewin_pp_post = 0.3; % 0.15
timewin_pp_noise = 0.1;
timewin_ps_prior = 0.08; % 0.1
timewin_ps_post = 0.45; % 0.3
timewin_ps_noise = 0.1;

ntr_new = size(Z_env_new,2);
t_p = zeros(ntr_new,1);
t_pp = zeros(ntr_new,1);
t_ps = zeros(ntr_new,1);

% Tpp
for i=1:ntr_new
    x_offset = geom_tab_new(i,5);
    t_p(i) = x_offset/vp;
    t_pp(i) = sqrt(x_offset * x_offset + 4.0 * z * z)/vp;
end

% Tps
gamma = vp/vs;
for i=1:ntr_new
    syms mx;
    x_offset = geom_tab_new(i,5);
    x_ccp = vpasolve(sqrt(gamma^2 + (gamma^2-1)*mx*mx/(z*z)) * x_offset / (1+sqrt(gamma^2 + (gamma^2-1)*mx*mx/(z*z))) +...
                     mx/sqrt(gamma^2 + (gamma^2-1)*mx*mx/(z*z)) == x_offset, mx);
    t_ps(i) = sqrt(x_ccp*x_ccp + z*z)/vp + sqrt((x_offset-x_ccp)*(x_offset-x_ccp) + z*z)/vs;
end

% visualize
figure;plot(t_p);hold on;plot(t_pp);hold on;plot(t_ps);legend('P','PP','PS');ylabel('Time (s)');xlabel('Trace number');
plot_check_timewin(ZRT_env_new,Z_env_new,R_env_new,T_env_new,RT_env_new,Z_stack,R_stack,T_stack,RT_stack,geom_tab_new,...
    nt,dt,ntr_new,t_p,t_pp,t_ps,timewin_p_prior,timewin_p_post,timewin_pp_prior,timewin_pp_post,timewin_ps_prior,timewin_ps_post,'combine');

% (Optional) measure time window uncertainty
UQ_tw_pp_prior = [0.00,0.01,0.02,0.03,0.04,0.05]; % 0.05
UQ_tw_pp_post = [0.2,0.22,0.24,0.26,0.28,0.3]; % 0.15
UQ_tw_ps_prior = [0.08,0.084,0.088,0.092,0.096,0.1]; % 0.1
UQ_tw_ps_post = [0.3,0.33,0.36,0.39,0.42,0.45]; % 0.3

plot_check_timewin(ZRT_env_new,Z_env_new,R_env_new,T_env_new,RT_env_new,Z_stack,R_stack,T_stack,RT_stack,geom_tab_new,...
    nt,dt,ntr_new,t_p,t_pp,t_ps,timewin_p_prior,timewin_p_post,min(UQ_tw_pp_prior),min(UQ_tw_pp_post),...
    min(UQ_tw_ps_prior),min(UQ_tw_ps_post),'combine');

plot_check_timewin(ZRT_env_new,Z_env_new,R_env_new,T_env_new,RT_env_new,Z_stack,R_stack,T_stack,RT_stack,geom_tab_new,...
    nt,dt,ntr_new,t_p,t_pp,t_ps,timewin_p_prior,timewin_p_post,max(UQ_tw_pp_prior),max(UQ_tw_pp_post),...
    max(UQ_tw_ps_prior),max(UQ_tw_ps_post),'combine');


%% 4. compute AVO (observed)
% average offset range=2~6km by bins (bins: [2, 2.1] (2.1, 2.2] ... (5.9, 6] (6, 6.1])
offset_bin_length = 200; % meter  100,200,500

% ====PS/PP====
bin_psppratio = amp_ratio(ZRT_env_new,ZRT_env_new,geom_tab_new,dt,t_ps,t_pp,timewin_ps_prior,timewin_ps_post,...
    timewin_pp_prior,timewin_pp_post,offset_min,offset_bin_length,offset_max,'median',1500);

% std_twin_psppratio = amp_ratio_std_twin(ZRT_env_new,ZRT_env_new,geom_tab_new,dt,t_ps,t_pp,UQ_tw_ps_prior,UQ_tw_ps_post,...
%     UQ_tw_pp_prior,UQ_tw_pp_post,offset_min,offset_bin_length,offset_max,'median');
% figure;errorbar(bin_psppratio(:,1),bin_psppratio(:,2),std_twin_psppratio,'o');
% xlabel('offset (m)');ylabel('PS/PP ratio');

N_bootstrap = 1000;
std_bsp_psppratio = amp_ratio_std_bsp(ZRT_env_new,ZRT_env_new,geom_tab_new,dt,t_ps,t_pp,timewin_ps_prior,timewin_ps_post,...
    timewin_pp_prior,timewin_pp_post,offset_min,offset_bin_length,offset_max,'median',N_bootstrap);

figure;errorbar(bin_psppratio(:,1),bin_psppratio(:,2),std_bsp_psppratio,'o');
xlabel('offset (m)');ylabel('PS/PP ratio');
ylim([0.5,1.1]);xlim([offset_min-100,offset_max]);

% (Optional) PP,PS
bin_dps_dpp = amp_pp_ps(ZRT_env_new,geom_tab_new,dt,t_ps,t_pp,timewin_ps_prior,timewin_ps_post,timewin_ps_noise,...
    timewin_pp_prior,timewin_pp_post,timewin_pp_noise,offset_min,offset_bin_length,offset_max,1500);


%% 5. (Optional) experiment PP amplitude reduction than P arrival (UK=unknown)
% SNR params
%yn_snr_qc_UK = 1;
%snr_threshold_UK = 1.0;
% offset range (meter)
offset_min_UK = 2000;
offset_max_UK = 5200;
offset_bin_UK = 500;
bin_UK = offset_min_UK:offset_bin_UK:offset_max_UK;
nbin_UK = length(bin_UK);
%[Z_env_new_UK, R_env_new_UK, T_env_new_UK, RT_env_new_UK, geom_tab_new_UK] = ...
%load_env_so(su_fileZ_env, su_fileR_env, su_fileT_env, su_fileRT_env, 'Z', yn_snr_qc_UK, snr_threshold_UK, offset_min_UK, offset_max_UK);
%ZRT_env_new_UK = Z_env_new_UK + R_env_new_UK + T_env_new_UK;

for itmp=1:nbin_UK-1
    if itmp==1
        index_UK = find(geom_tab_new(:,5)>=bin_UK(itmp) & geom_tab_new(:,5)<=bin_UK(itmp+1));
    else
        index_UK = find(geom_tab_new(:,5)>bin_UK(itmp) & geom_tab_new(:,5)<=bin_UK(itmp+1));
    end
    data_target1 = Z_env_new(:,index_UK);
    data_target2 = ZRT_env_new(:,index_UK);
    
    tp_med_UK = median(t_p(index_UK)); tp_mean_UK = mean(t_p(index_UK));
    tpp_med_UK = median(t_pp(index_UK)); tpp_mean_UK = mean(t_pp(index_UK));
    tps_med_UK = median(t_ps(index_UK)); tps_mean_UK = mean(t_ps(index_UK));

    exp_trace1_median = zeros(nt,1);
    exp_trace1_mean = zeros(nt,1);
    exp_trace2_median = zeros(nt,1);
    exp_trace2_mean = zeros(nt,1);
    for it=1:nt
        exp_trace1_median(it) = median(data_target1(it,:));
        exp_trace1_mean(it) = mean(data_target1(it,:));
        exp_trace2_median(it) = median(data_target2(it,:));
        exp_trace2_mean(it) = mean(data_target2(it,:));
    end
    
    figure;
    subplot(2,2,1);plot(0:dt:(nt-1)*dt,exp_trace1_mean,'k');hold on;
    %plot([tp_mean_UK,tp_mean_UK],[min(exp_trace1_mean),max(exp_trace1_mean)],'r--');hold on;
    plot([tpp_mean_UK,tpp_mean_UK],[min(exp_trace1_mean),max(exp_trace1_mean)],'r--');hold on;
    plot([tpp_mean_UK-timewin_pp_prior,tpp_mean_UK-timewin_pp_prior],[min(exp_trace1_mean),max(exp_trace1_mean)],'r--');hold on;
    plot([tpp_mean_UK+timewin_pp_post,tpp_mean_UK+timewin_pp_post],[min(exp_trace1_mean),max(exp_trace1_mean)],'r--');
    title(['Z: ',num2str(bin_UK(itmp)/1000),'-',num2str(bin_UK(itmp+1)/1000),' km']);
    
    subplot(2,2,3);plot(0:dt:(nt-1)*dt,exp_trace1_median,'k');hold on;
    %plot([tp_med_UK,tp_med_UK],[min(exp_trace1_median),max(exp_trace1_median)],'r--');hold on;
    plot([tpp_med_UK,tpp_med_UK],[min(exp_trace1_median),max(exp_trace1_median)],'r--');
    plot([tpp_med_UK-timewin_pp_prior,tpp_med_UK-timewin_pp_prior],[min(exp_trace1_median),max(exp_trace1_median)],'r--');
    plot([tpp_med_UK+timewin_pp_post,tpp_med_UK+timewin_pp_post],[min(exp_trace1_median),max(exp_trace1_median)],'r--');
    
    subplot(2,2,2);plot(0:dt:(nt-1)*dt,exp_trace2_mean,'k');hold on;
    plot([tpp_mean_UK,tpp_mean_UK],[min(exp_trace2_mean),max(exp_trace2_mean)],'b-');hold on;
    plot([tpp_mean_UK-timewin_pp_prior,tpp_mean_UK-timewin_pp_prior],[min(exp_trace2_mean),max(exp_trace2_mean)],'b--');hold on;
    plot([tpp_mean_UK+timewin_pp_post,tpp_mean_UK+timewin_pp_post],[min(exp_trace2_mean),max(exp_trace2_mean)],'b--');hold on;
    plot([tps_mean_UK,tps_mean_UK],[min(exp_trace2_mean),max(exp_trace2_mean)],'r-');hold on;
    plot([tps_mean_UK-timewin_ps_prior,tps_mean_UK-timewin_ps_prior],[min(exp_trace2_mean),max(exp_trace2_mean)],'r--');hold on;
    plot([tps_mean_UK+timewin_ps_post,tps_mean_UK+timewin_ps_post],[min(exp_trace2_mean),max(exp_trace2_mean)],'r--');
    title(['Z+R+T: ',num2str(bin_UK(itmp)/1000),'-',num2str(bin_UK(itmp+1)/1000),' km']);
    
    subplot(2,2,4);plot(0:dt:(nt-1)*dt,exp_trace2_median,'k');hold on;
    plot([tpp_med_UK,tpp_med_UK],[min(exp_trace2_median),max(exp_trace2_median)],'b-');hold on;
    plot([tpp_med_UK-timewin_pp_prior,tpp_med_UK-timewin_pp_prior],[min(exp_trace2_median),max(exp_trace2_median)],'b--');hold on;
    plot([tpp_med_UK+timewin_pp_post,tpp_med_UK+timewin_pp_post],[min(exp_trace2_median),max(exp_trace2_median)],'b--');hold on;
    plot([tps_med_UK,tps_med_UK],[min(exp_trace2_median),max(exp_trace2_median)],'r-');hold on;
    plot([tps_med_UK-timewin_ps_prior,tps_med_UK-timewin_ps_prior],[min(exp_trace2_median),max(exp_trace2_median)],'r--');hold on;
    plot([tps_med_UK+timewin_ps_post,tps_med_UK+timewin_ps_post],[min(exp_trace2_median),max(exp_trace2_median)],'r--');
end


%% 6. time window (synthetic)
syn_timewin_p_prior = 0.28;
syn_timewin_p_post = 0.3;
syn_timewin_pp_prior = 0.05;
syn_timewin_pp_post = 0.3;
syn_timewin_ps_prior = 0.08;
syn_timewin_ps_post = 0.45;

% 10-*
su_ez_10_10 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_10_10_32bit.su';
su_ex_10_10 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_10_10_32bit.su';
[~, ~, ez_10_10, ex_10_10, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_10_10, su_ex_10_10, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);

% 15-*
su_ez_15_10 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_15_10_32bit.su';
su_ex_15_10 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_15_10_32bit.su';
su_ez_15_15 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_15_15_32bit.su';
su_ex_15_15 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_15_15_32bit.su';
[~, ~, ez_15_10, ex_15_10, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_15_10, su_ex_15_10, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_15_15, ex_15_15, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_15_15, su_ex_15_15, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);

% 20-*
su_ez_20_10 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_20_10_32bit.su';
su_ex_20_10 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_20_10_32bit.su';
su_ez_20_15 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_20_15_32bit.su';
su_ex_20_15 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_20_15_32bit.su';
su_ez_20_20 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_20_20_32bit.su';
su_ex_20_20 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_20_20_32bit.su';
[~, ~, ez_20_10, ex_20_10, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_20_10, su_ex_20_10, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_20_15, ex_20_15, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_20_15, su_ex_20_15, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_20_20, ex_20_20, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_20_20, su_ex_20_20, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);

% 25-*
su_ez_25_10 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_25_10_32bit.su';
su_ex_25_10 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_25_10_32bit.su';
su_ez_25_15 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_25_15_32bit.su';
su_ex_25_15 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_25_15_32bit.su';
su_ez_25_20 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_25_20_32bit.su';
su_ex_25_20 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_25_20_32bit.su';
su_ez_25_25 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_25_25_32bit.su';
su_ex_25_25 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_25_25_32bit.su';
[~, ~, ez_25_10, ex_25_10, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_25_10, su_ex_25_10, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_25_15, ex_25_15, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_25_15, su_ex_25_15, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_25_20, ex_25_20, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_25_20, su_ex_25_20, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_25_25, ex_25_25, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_25_25, su_ex_25_25, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);

% 30-*
su_ez_30_10 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_30_10_32bit.su';
su_ex_30_10 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_30_10_32bit.su';
su_ez_30_15 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_30_15_32bit.su';
su_ex_30_15 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_30_15_32bit.su';
su_ez_30_20 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_30_20_32bit.su';
su_ex_30_20 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_30_20_32bit.su';
su_ez_30_25 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_30_25_32bit.su';
su_ex_30_25 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_30_25_32bit.su';
su_ez_30_30 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_30_30_32bit.su';
su_ex_30_30 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_30_30_32bit.su';
[~, ~, ez_30_10, ex_30_10, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_30_10, su_ex_30_10, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_30_15, ex_30_15, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_30_15, su_ex_30_15, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_30_20, ex_30_20, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_30_20, su_ex_30_20, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_30_25, ex_30_25, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_30_25, su_ex_30_25, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_30_30, ex_30_30, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_30_30, su_ex_30_30, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);

% 35-*
su_ez_35_10 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_35_10_32bit.su';
su_ex_35_10 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_35_10_32bit.su';
su_ez_35_15 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_35_15_32bit.su';
su_ex_35_15 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_35_15_32bit.su';
su_ez_35_20 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_35_20_32bit.su';
su_ex_35_20 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_35_20_32bit.su';
su_ez_35_25 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_35_25_32bit.su';
su_ex_35_25 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_35_25_32bit.su';
su_ez_35_30 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_35_30_32bit.su';
su_ex_35_30 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_35_30_32bit.su';
su_ez_35_35 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_35_35_32bit.su';
su_ex_35_35 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_35_35_32bit.su';
[~, ~, ez_35_10, ex_35_10, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_35_10, su_ex_35_10, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_35_15, ex_35_15, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_35_15, su_ex_35_15, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_35_20, ex_35_20, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_35_20, su_ex_35_20, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_35_25, ex_35_25, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_35_25, su_ex_35_25, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_35_30, ex_35_30, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_35_30, su_ex_35_30, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_35_35, ex_35_35, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_35_35, su_ex_35_35, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);

% 40-*
su_ez_40_10 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_40_10_32bit.su';
su_ex_40_10 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_40_10_32bit.su';
su_ez_40_15 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_40_15_32bit.su';
su_ex_40_15 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_40_15_32bit.su';
su_ez_40_20 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_40_20_32bit.su';
su_ex_40_20 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_40_20_32bit.su';
su_ez_40_25 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_40_25_32bit.su';
su_ex_40_25 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_40_25_32bit.su';
su_ez_40_30 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_40_30_32bit.su';
su_ex_40_30 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_40_30_32bit.su';
su_ez_40_35 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_40_35_32bit.su';
su_ex_40_35 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_40_35_32bit.su';
su_ez_40_40 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_40_40_32bit.su';
su_ex_40_40 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_40_40_32bit.su';
[~, ~, ez_40_10, ex_40_10, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_40_10, su_ex_40_10, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_40_15, ex_40_15, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_40_15, su_ex_40_15, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_40_20, ex_40_20, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_40_20, su_ex_40_20, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_40_25, ex_40_25, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_40_25, su_ex_40_25, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_40_30, ex_40_30, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_40_30, su_ex_40_30, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_40_35, ex_40_35, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_40_35, su_ex_40_35, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_40_40, ex_40_40, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_40_40, su_ex_40_40, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);

% 45-*
su_ez_45_10 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_45_10_32bit.su';
su_ex_45_10 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_45_10_32bit.su';
su_ez_45_15 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_45_15_32bit.su';
su_ex_45_15 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_45_15_32bit.su';
su_ez_45_20 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_45_20_32bit.su';
su_ex_45_20 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_45_20_32bit.su';
su_ez_45_25 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_45_25_32bit.su';
su_ex_45_25 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_45_25_32bit.su';
su_ez_45_30 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_45_30_32bit.su';
su_ex_45_30 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_45_30_32bit.su';
su_ez_45_35 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_45_35_32bit.su';
su_ex_45_35 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_45_35_32bit.su';
su_ez_45_40 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_45_40_32bit.su';
su_ex_45_40 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_45_40_32bit.su';
su_ez_45_45 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_45_45_32bit.su';
su_ex_45_45 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_45_45_32bit.su';
[~, ~, ez_45_10, ex_45_10, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_45_10, su_ex_45_10, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_45_15, ex_45_15, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_45_15, su_ex_45_15, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_45_20, ex_45_20, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_45_20, su_ex_45_20, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_45_25, ex_45_25, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_45_25, su_ex_45_25, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_45_30, ex_45_30, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_45_30, su_ex_45_30, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_45_35, ex_45_35, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_45_35, su_ex_45_35, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_45_40, ex_45_40, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_45_40, su_ex_45_40, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_45_45, ex_45_45, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_45_45, su_ex_45_45, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);

% 50-*
su_ez_50_15 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_50_15_32bit.su';
su_ex_50_15 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_50_15_32bit.su';
su_ez_50_20 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_50_20_32bit.su';
su_ex_50_20 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_50_20_32bit.su';
su_ez_50_25 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_50_25_32bit.su';
su_ex_50_25 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_50_25_32bit.su';
su_ez_50_30 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_50_30_32bit.su';
su_ex_50_30 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_50_30_32bit.su';
su_ez_50_35 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_50_35_32bit.su';
su_ex_50_35 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_50_35_32bit.su';
su_ez_50_40 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_50_40_32bit.su';
su_ex_50_40 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_50_40_32bit.su';
su_ez_50_45 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_50_45_32bit.su';
su_ex_50_45 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_50_45_32bit.su';
su_ez_50_50 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_50_50_32bit.su';
su_ex_50_50 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_50_50_32bit.su';
[~, ~, ez_50_15, ex_50_15, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_50_15, su_ex_50_15, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_50_20, ex_50_20, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_50_20, su_ex_50_20, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_50_25, ex_50_25, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_50_25, su_ex_50_25, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_50_30, ex_50_30, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_50_30, su_ex_50_30, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_50_35, ex_50_35, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_50_35, su_ex_50_35, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_50_40, ex_50_40, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_50_40, su_ex_50_40, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_50_45, ex_50_45, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_50_45, su_ex_50_45, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_50_50, ex_50_50, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_50_50, su_ex_50_50, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);

% 55-*
su_ez_55_20 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_55_20_32bit.su';
su_ex_55_20 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_55_20_32bit.su';
su_ez_55_25 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_55_25_32bit.su';
su_ex_55_25 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_55_25_32bit.su';
su_ez_55_30 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_55_30_32bit.su';
su_ex_55_30 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_55_30_32bit.su';
su_ez_55_35 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_55_35_32bit.su';
su_ex_55_35 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_55_35_32bit.su';
su_ez_55_40 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_55_40_32bit.su';
su_ex_55_40 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_55_40_32bit.su';
su_ez_55_45 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_55_45_32bit.su';
su_ex_55_45 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_55_45_32bit.su';
su_ez_55_50 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_55_50_32bit.su';
su_ex_55_50 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_55_50_32bit.su';
su_ez_55_55 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_55_55_32bit.su';
su_ex_55_55 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_55_55_32bit.su';
[~, ~, ez_55_20, ex_55_20, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_55_20, su_ex_55_20, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_55_25, ex_55_25, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_55_25, su_ex_55_25, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_55_30, ex_55_30, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_55_30, su_ex_55_30, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_55_35, ex_55_35, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_55_35, su_ex_55_35, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_55_40, ex_55_40, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_55_40, su_ex_55_40, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_55_45, ex_55_45, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_55_45, su_ex_55_45, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_55_50, ex_55_50, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_55_50, su_ex_55_50, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_55_55, ex_55_55, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_55_55, su_ex_55_55, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);

% 60-*
su_ez_60_30 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_60_30_32bit.su';
su_ex_60_30 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_60_30_32bit.su';
su_ez_60_35 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_60_35_32bit.su';
su_ex_60_35 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_60_35_32bit.su';
su_ez_60_40 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_60_40_32bit.su';
su_ex_60_40 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_60_40_32bit.su';
su_ez_60_45 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_60_45_32bit.su';
su_ex_60_45 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_60_45_32bit.su';
su_ez_60_50 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_60_50_32bit.su';
su_ex_60_50 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_60_50_32bit.su';
su_ez_60_55 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_60_55_32bit.su';
su_ex_60_55 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_60_55_32bit.su';
su_ez_60_60 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_60_60_32bit.su';
su_ex_60_60 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_60_60_32bit.su';
[~, ~, ez_60_30, ex_60_30, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_60_30, su_ex_60_30, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_60_35, ex_60_35, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_60_35, su_ex_60_35, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_60_40, ex_60_40, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_60_40, su_ex_60_40, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_60_45, ex_60_45, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_60_45, su_ex_60_45, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_60_50, ex_60_50, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_60_50, su_ex_60_50, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_60_55, ex_60_55, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_60_55, su_ex_60_55, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_60_60, ex_60_60, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_60_60, su_ex_60_60, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);

% 65-*
su_ez_65_40 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_65_40_32bit.su';
su_ex_65_40 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_65_40_32bit.su';
su_ez_65_45 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_65_45_32bit.su';
su_ex_65_45 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_65_45_32bit.su';
su_ez_65_50 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_65_50_32bit.su';
su_ex_65_50 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_65_50_32bit.su';
su_ez_65_55 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_65_55_32bit.su';
su_ex_65_55 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_65_55_32bit.su';
su_ez_65_60 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_65_60_32bit.su';
su_ex_65_60 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_65_60_32bit.su';
su_ez_65_65 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_65_65_32bit.su';
su_ex_65_65 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_65_65_32bit.su';
[~, ~, ez_65_40, ex_65_40, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_65_40, su_ex_65_40, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_65_45, ex_65_45, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_65_45, su_ex_65_45, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_65_50, ex_65_50, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_65_50, su_ex_65_50, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_65_55, ex_65_55, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_65_55, su_ex_65_55, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_65_60, ex_65_60, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_65_60, su_ex_65_60, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_65_65, ex_65_65, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_65_65, su_ex_65_65, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);

% 70-*
su_ez_70_50 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_70_50_32bit.su';
su_ex_70_50 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_70_50_32bit.su';
su_ez_70_55 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_70_55_32bit.su';
su_ex_70_55 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_70_55_32bit.su';
su_ez_70_60 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_70_60_32bit.su';
su_ex_70_60 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_70_60_32bit.su';
su_ez_70_65 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_70_65_32bit.su';
su_ex_70_65 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_70_65_32bit.su';
su_ez_70_70 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_70_70_32bit.su';
su_ex_70_70 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_70_70_32bit.su';
[~, ~, ez_70_50, ex_70_50, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_70_50, su_ex_70_50, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_70_55, ex_70_55, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_70_55, su_ex_70_55, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_70_60, ex_70_60, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_70_60, su_ex_70_60, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_70_65, ex_70_65, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_70_65, su_ex_70_65, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_70_70, ex_70_70, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_70_70, su_ex_70_70, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);

% others
su_ez_10_30 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_10_30_32bit.su';
su_ex_10_30 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_10_30_32bit.su';
su_ez_10_50 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_10_50_32bit.su';
su_ex_10_50 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_10_50_32bit.su';
su_ez_30_50 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_30_50_32bit.su';
su_ex_30_50 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_30_50_32bit.su';
[~, ~, ez_10_30, ex_10_30, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_10_30, su_ex_10_30, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_10_50, ex_10_50, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_10_50, su_ex_10_50, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_30_50, ex_30_50, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_30_50, su_ex_30_50, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);


figure;subplot(2,2,1);
imagesc(1:ntr_new,0:dt:(nt-1)*dt,ez_10_10);clim([0,2e-14]);colormap(colorbar_bwr);
xlabel(['trace# (offset ',num2str(offset_min/1000),'-',num2str(offset_max/1000),' km)']);ylabel('time (s)');
hold on;plot(1:ntr_new,t_p-syn_timewin_p_prior,'k');
hold on;plot(1:ntr_new,t_p+syn_timewin_p_post,'k');
hold on;plot(1:ntr_new,t_pp-syn_timewin_pp_prior,'k');
hold on;plot(1:ntr_new,t_pp+syn_timewin_pp_post,'k');
hold on;plot(1:ntr_new,t_ps-syn_timewin_ps_prior,'k');
hold on;plot(1:ntr_new,t_ps+syn_timewin_ps_post,'k');title('ez (10-10)');
subplot(2,2,2);
imagesc(1:ntr_new,0:dt:(nt-1)*dt,ex_10_10);clim([0,2e-14]);colormap(colorbar_bwr);
xlabel(['trace# (offset ',num2str(offset_min/1000),'-',num2str(offset_max/1000),' km)']);ylabel('time (s)');
hold on;plot(1:ntr_new,t_p-syn_timewin_p_prior,'k');
hold on;plot(1:ntr_new,t_p+syn_timewin_p_post,'k');
hold on;plot(1:ntr_new,t_pp-syn_timewin_pp_prior,'k');
hold on;plot(1:ntr_new,t_pp+syn_timewin_pp_post,'k');
hold on;plot(1:ntr_new,t_ps-syn_timewin_ps_prior,'k');
hold on;plot(1:ntr_new,t_ps+syn_timewin_ps_post,'k');title('ex (10-10)');
subplot(2,2,3);
imagesc(1:ntr_new,0:dt:(nt-1)*dt,ez_10_10+ex_10_10);clim([0,2e-14]);colormap(colorbar_bwr);
xlabel(['trace# (offset ',num2str(offset_min/1000),'-',num2str(offset_max/1000),' km)']);ylabel('time (s)');
hold on;plot(1:ntr_new,t_p-syn_timewin_p_prior,'k');
hold on;plot(1:ntr_new,t_p+syn_timewin_p_post,'k');
hold on;plot(1:ntr_new,t_pp-syn_timewin_pp_prior,'k');
hold on;plot(1:ntr_new,t_pp+syn_timewin_pp_post,'k');
hold on;plot(1:ntr_new,t_ps-syn_timewin_ps_prior,'k');
hold on;plot(1:ntr_new,t_ps+syn_timewin_ps_post,'k');title('ez+ex (10-10)');

itr=4000;
figure;plot(0:dt:(nt-1)*dt,ex_10_10(:,itr),'k');hold on;
plot(t_pp(itr)-syn_timewin_pp_prior,0,'r*');hold on;plot(t_pp(itr)+syn_timewin_pp_post,0,'r*');


%% 7 compute AVO - trace by trace (synthetic)
% average offset range=2~6km by bins (bins: [2, 2.1] (2.1, 2.2] ... (5.9, 6] (6, 6.1])

% ====PS/PP====
% 10-*
ez_ex_10_10 = ez_10_10 + ex_10_10;
bin_psppratio_syn1010 = amp_ratio(ez_ex_10_10,ez_ex_10_10,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

% 15-*
ez_ex_15_10 = ez_15_10 + ex_15_10;
bin_psppratio_syn1510 = amp_ratio(ez_ex_15_10,ez_ex_15_10,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_15_15 = ez_15_15 + ex_15_15;
bin_psppratio_syn1515 = amp_ratio(ez_ex_15_15,ez_ex_15_15,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

% 20-*
ez_ex_20_10 = ez_20_10 + ex_20_10;
bin_psppratio_syn2010 = amp_ratio(ez_ex_20_10,ez_ex_20_10,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_20_15 = ez_20_15 + ex_20_15;
bin_psppratio_syn2015 = amp_ratio(ez_ex_20_15,ez_ex_20_15,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_20_20 = ez_20_20 + ex_20_20;
bin_psppratio_syn2020 = amp_ratio(ez_ex_20_20,ez_ex_20_20,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

% 25-*
ez_ex_25_10 = ez_25_10 + ex_25_10;
bin_psppratio_syn2510 = amp_ratio(ez_ex_25_10,ez_ex_25_10,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_25_15 = ez_25_15 + ex_25_15;
bin_psppratio_syn2515 = amp_ratio(ez_ex_25_15,ez_ex_25_15,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_25_20 = ez_25_20 + ex_25_20;
bin_psppratio_syn2520 = amp_ratio(ez_ex_25_20,ez_ex_25_20,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_25_25 = ez_25_25 + ex_25_25;
bin_psppratio_syn2525 = amp_ratio(ez_ex_25_25,ez_ex_25_25,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

% 30-*
ez_ex_30_10 = ez_30_10 + ex_30_10;
bin_psppratio_syn3010 = amp_ratio(ez_ex_30_10,ez_ex_30_10,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_30_15 = ez_30_15 + ex_30_15;
bin_psppratio_syn3015 = amp_ratio(ez_ex_30_15,ez_ex_30_15,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_30_20 = ez_30_20 + ex_30_20;
bin_psppratio_syn3020 = amp_ratio(ez_ex_30_20,ez_ex_30_20,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_30_25 = ez_30_25 + ex_30_25;
bin_psppratio_syn3025 = amp_ratio(ez_ex_30_25,ez_ex_30_25,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_30_30 = ez_30_30 + ex_30_30;
bin_psppratio_syn3030 = amp_ratio(ez_ex_30_30,ez_ex_30_30,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

% 35-*
ez_ex_35_10 = ez_35_10 + ex_35_10;
bin_psppratio_syn3510 = amp_ratio(ez_ex_35_10,ez_ex_35_10,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_35_15 = ez_35_15 + ex_35_15;
bin_psppratio_syn3515 = amp_ratio(ez_ex_35_15,ez_ex_35_15,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_35_20 = ez_35_20 + ex_35_20;
bin_psppratio_syn3520 = amp_ratio(ez_ex_35_20,ez_ex_35_20,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_35_25 = ez_35_25 + ex_35_25;
bin_psppratio_syn3525 = amp_ratio(ez_ex_35_25,ez_ex_35_25,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_35_30 = ez_35_30 + ex_35_30;
bin_psppratio_syn3530 = amp_ratio(ez_ex_35_30,ez_ex_35_30,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_35_35 = ez_35_35 + ex_35_35;
bin_psppratio_syn3535 = amp_ratio(ez_ex_35_35,ez_ex_35_35,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

% 40-*
ez_ex_40_10 = ez_40_10 + ex_40_10;
bin_psppratio_syn4010 = amp_ratio(ez_ex_40_10,ez_ex_40_10,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_40_15 = ez_40_15 + ex_40_15;
bin_psppratio_syn4015 = amp_ratio(ez_ex_40_15,ez_ex_40_15,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_40_20 = ez_40_20 + ex_40_20;
bin_psppratio_syn4020 = amp_ratio(ez_ex_40_20,ez_ex_40_20,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_40_25 = ez_40_25 + ex_40_25;
bin_psppratio_syn4025 = amp_ratio(ez_ex_40_25,ez_ex_40_25,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_40_30 = ez_40_30 + ex_40_30;
bin_psppratio_syn4030 = amp_ratio(ez_ex_40_30,ez_ex_40_30,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_40_35 = ez_40_35 + ex_40_35;
bin_psppratio_syn4035 = amp_ratio(ez_ex_40_35,ez_ex_40_35,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_40_40 = ez_40_40 + ex_40_40;
bin_psppratio_syn4040 = amp_ratio(ez_ex_40_40,ez_ex_40_40,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

% 45-*
ez_ex_45_10 = ez_45_10 + ex_45_10;
bin_psppratio_syn4510 = amp_ratio(ez_ex_45_10,ez_ex_45_10,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_45_15 = ez_45_15 + ex_45_15;
bin_psppratio_syn4515 = amp_ratio(ez_ex_45_15,ez_ex_45_15,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_45_20 = ez_45_20 + ex_45_20;
bin_psppratio_syn4520 = amp_ratio(ez_ex_45_20,ez_ex_45_20,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_45_25 = ez_45_25 + ex_45_25;
bin_psppratio_syn4525 = amp_ratio(ez_ex_45_25,ez_ex_45_25,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_45_30 = ez_45_30 + ex_45_30;
bin_psppratio_syn4530 = amp_ratio(ez_ex_45_30,ez_ex_45_30,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_45_35 = ez_45_35 + ex_45_35;
bin_psppratio_syn4535 = amp_ratio(ez_ex_45_35,ez_ex_45_35,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_45_40 = ez_45_40 + ex_45_40;
bin_psppratio_syn4540 = amp_ratio(ez_ex_45_40,ez_ex_45_40,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_45_45 = ez_45_45 + ex_45_45;
bin_psppratio_syn4545 = amp_ratio(ez_ex_45_45,ez_ex_45_45,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

% 50-*
ez_ex_50_15 = ez_50_15 + ex_50_15;
bin_psppratio_syn5015 = amp_ratio(ez_ex_50_15,ez_ex_50_15,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_50_20 = ez_50_20 + ex_50_20;
bin_psppratio_syn5020 = amp_ratio(ez_ex_50_20,ez_ex_50_20,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_50_25 = ez_50_25 + ex_50_25;
bin_psppratio_syn5025 = amp_ratio(ez_ex_50_25,ez_ex_50_25,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_50_30 = ez_50_30 + ex_50_30;
bin_psppratio_syn5030 = amp_ratio(ez_ex_50_30,ez_ex_50_30,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_50_35 = ez_50_35 + ex_50_35;
bin_psppratio_syn5035 = amp_ratio(ez_ex_50_35,ez_ex_50_35,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_50_40 = ez_50_40 + ex_50_40;
bin_psppratio_syn5040 = amp_ratio(ez_ex_50_40,ez_ex_50_40,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_50_45 = ez_50_45 + ex_50_45;
bin_psppratio_syn5045 = amp_ratio(ez_ex_50_45,ez_ex_50_45,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_50_50 = ez_50_50 + ex_50_50;
bin_psppratio_syn5050 = amp_ratio(ez_ex_50_50,ez_ex_50_50,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

% 55-*
ez_ex_55_20 = ez_55_20 + ex_55_20;
bin_psppratio_syn5520 = amp_ratio(ez_ex_55_20,ez_ex_55_20,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_55_25 = ez_55_25 + ex_55_25;
bin_psppratio_syn5525 = amp_ratio(ez_ex_55_25,ez_ex_55_25,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_55_30 = ez_55_30 + ex_55_30;
bin_psppratio_syn5530 = amp_ratio(ez_ex_55_30,ez_ex_55_30,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_55_35 = ez_55_35 + ex_55_35;
bin_psppratio_syn5535 = amp_ratio(ez_ex_55_35,ez_ex_55_35,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_55_40 = ez_55_40 + ex_55_40;
bin_psppratio_syn5540 = amp_ratio(ez_ex_55_40,ez_ex_55_40,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_55_45 = ez_55_45 + ex_55_45;
bin_psppratio_syn5545 = amp_ratio(ez_ex_55_45,ez_ex_55_45,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_55_50 = ez_55_50 + ex_55_50;
bin_psppratio_syn5550 = amp_ratio(ez_ex_55_50,ez_ex_55_50,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_55_55 = ez_55_55 + ex_55_55;
bin_psppratio_syn5555 = amp_ratio(ez_ex_55_55,ez_ex_55_55,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

% 60-*
ez_ex_60_30 = ez_60_30 + ex_60_30;
bin_psppratio_syn6030 = amp_ratio(ez_ex_60_30,ez_ex_60_30,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_60_35 = ez_60_35 + ex_60_35;
bin_psppratio_syn6035 = amp_ratio(ez_ex_60_35,ez_ex_60_35,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_60_40 = ez_60_40 + ex_60_40;
bin_psppratio_syn6040 = amp_ratio(ez_ex_60_40,ez_ex_60_40,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_60_45 = ez_60_45 + ex_60_45;
bin_psppratio_syn6045 = amp_ratio(ez_ex_60_45,ez_ex_60_45,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_60_50 = ez_60_50 + ex_60_50;
bin_psppratio_syn6050 = amp_ratio(ez_ex_60_50,ez_ex_60_50,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_60_55 = ez_60_55 + ex_60_55;
bin_psppratio_syn6055 = amp_ratio(ez_ex_60_55,ez_ex_60_55,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_60_60 = ez_60_60 + ex_60_60;
bin_psppratio_syn6060 = amp_ratio(ez_ex_60_60,ez_ex_60_60,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

% 65-*
ez_ex_65_40 = ez_65_40 + ex_65_40;
bin_psppratio_syn6540 = amp_ratio(ez_ex_65_40,ez_ex_65_40,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_65_45 = ez_65_45 + ex_65_45;
bin_psppratio_syn6545 = amp_ratio(ez_ex_65_45,ez_ex_65_45,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_65_50 = ez_65_50 + ex_65_50;
bin_psppratio_syn6550 = amp_ratio(ez_ex_65_50,ez_ex_65_50,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_65_55 = ez_65_55 + ex_65_55;
bin_psppratio_syn6555 = amp_ratio(ez_ex_65_55,ez_ex_65_55,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_65_60 = ez_65_60 + ex_65_60;
bin_psppratio_syn6560 = amp_ratio(ez_ex_65_60,ez_ex_65_60,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_65_65 = ez_65_65 + ex_65_65;
bin_psppratio_syn6565 = amp_ratio(ez_ex_65_65,ez_ex_65_65,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

% 70-*
ez_ex_70_50 = ez_70_50 + ex_70_50;
bin_psppratio_syn7050 = amp_ratio(ez_ex_70_50,ez_ex_70_50,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_70_55 = ez_70_55 + ex_70_55;
bin_psppratio_syn7055 = amp_ratio(ez_ex_70_55,ez_ex_70_55,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_70_60 = ez_70_60 + ex_70_60;
bin_psppratio_syn7060 = amp_ratio(ez_ex_70_60,ez_ex_70_60,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_70_65 = ez_70_65 + ex_70_65;
bin_psppratio_syn7065 = amp_ratio(ez_ex_70_65,ez_ex_70_65,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_70_70 = ez_70_70 + ex_70_70;
bin_psppratio_syn7070 = amp_ratio(ez_ex_70_70,ez_ex_70_70,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

% others
ez_ex_10_30 = ez_10_30 + ex_10_30; % bad
bin_psppratio_syn1030 = amp_ratio(ez_ex_10_30,ez_ex_10_30,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_10_50 = ez_10_50 + ex_10_50; % bad
bin_psppratio_syn1050 = amp_ratio(ez_ex_10_50,ez_ex_10_50,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_30_50 = ez_30_50 + ex_30_50;
bin_psppratio_syn3050 = amp_ratio(ez_ex_30_50,ez_ex_30_50,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);


%% 8. plots (fix dvp)
figure;
plot(bin_psppratio_syn1010(:,1),bin_psppratio_syn1010(:,2),'o-');hold on;
plot(bin_psppratio_syn1510(:,1),bin_psppratio_syn1510(:,2),'o-');hold on;
plot(bin_psppratio_syn1515(:,1),bin_psppratio_syn1515(:,2),'o-');hold on;
errorbar(bin_psppratio(:,1),bin_psppratio(:,2),std_bsp_psppratio,'ko');legend('10-10','15-10','15-15');
xlabel('offset (m)');ylabel('PS/PP ratio');xlim([offset_min-100,offset_max-100]);
% figure;subplot(1,2,1);
% plot(bin_psppratio_syn1010(:,1),bin_psppratio_syn1010(:,5),'o-');hold on;
% plot(bin_psppratio_syn1510(:,1),bin_psppratio_syn1510(:,5),'o-');hold on;
% plot(bin_psppratio_syn1515(:,1),bin_psppratio_syn1515(:,5),'o-');hold on;
% legend('10-10','15-10','15-15');
% xlabel('offset (m)');ylabel('dPS');xlim([offset_min-100,offset_max-100]);
% subplot(1,2,2);
% plot(bin_psppratio_syn1010(:,1),bin_psppratio_syn1010(:,6),'o-');hold on;
% plot(bin_psppratio_syn1510(:,1),bin_psppratio_syn1510(:,6),'o-');hold on;
% plot(bin_psppratio_syn1515(:,1),bin_psppratio_syn1515(:,6),'o-');hold on;
% legend('10-10','15-10','15-15');
% xlabel('offset (m)');ylabel('dPP');xlim([offset_min-100,offset_max-100]);

figure;
plot(bin_psppratio_syn2010(:,1),bin_psppratio_syn2010(:,2),'o-');hold on;
plot(bin_psppratio_syn2015(:,1),bin_psppratio_syn2015(:,2),'o-');hold on;
plot(bin_psppratio_syn2020(:,1),bin_psppratio_syn2020(:,2),'o-');hold on;
errorbar(bin_psppratio(:,1),bin_psppratio(:,2),std_bsp_psppratio,'ko');legend('20-10','20-15','20-20');
xlabel('offset (m)');ylabel('PS/PP ratio');xlim([offset_min-100,offset_max-100]);
% figure;subplot(1,2,1);
% plot(bin_psppratio_syn2010(:,1),bin_psppratio_syn2010(:,5),'o-');hold on;
% plot(bin_psppratio_syn2015(:,1),bin_psppratio_syn2015(:,5),'o-');hold on;
% plot(bin_psppratio_syn2020(:,1),bin_psppratio_syn2020(:,5),'o-');hold on;
% legend('20-10','20-15','20-20');
% xlabel('offset (m)');ylabel('dPS');xlim([offset_min-100,offset_max-100]);
% subplot(1,2,2);
% plot(bin_psppratio_syn2010(:,1),bin_psppratio_syn2010(:,6),'o-');hold on;
% plot(bin_psppratio_syn2015(:,1),bin_psppratio_syn2015(:,6),'o-');hold on;
% plot(bin_psppratio_syn2020(:,1),bin_psppratio_syn2020(:,6),'o-');hold on;
% legend('20-10','20-15','20-20');
% xlabel('offset (m)');ylabel('dPP');xlim([offset_min-100,offset_max-100]);

figure;
plot(bin_psppratio_syn2510(:,1),bin_psppratio_syn2510(:,2),'o-');hold on;
plot(bin_psppratio_syn2515(:,1),bin_psppratio_syn2515(:,2),'o-');hold on;
plot(bin_psppratio_syn2520(:,1),bin_psppratio_syn2520(:,2),'o-');hold on;
plot(bin_psppratio_syn2525(:,1),bin_psppratio_syn2525(:,2),'o-');hold on;
errorbar(bin_psppratio(:,1),bin_psppratio(:,2),std_bsp_psppratio,'ko');legend('25-10','25-15','25-20','25-25');
xlabel('offset (m)');ylabel('PS/PP ratio');xlim([offset_min-100,offset_max-100]);
% figure;subplot(1,2,1);
% plot(bin_psppratio_syn2510(:,1),bin_psppratio_syn2510(:,5),'o-');hold on;
% plot(bin_psppratio_syn2515(:,1),bin_psppratio_syn2515(:,5),'o-');hold on;
% plot(bin_psppratio_syn2520(:,1),bin_psppratio_syn2520(:,5),'o-');hold on;
% plot(bin_psppratio_syn2525(:,1),bin_psppratio_syn2525(:,5),'o-');hold on;
% legend('25-10','25-15','25-20','25-25');
% xlabel('offset (m)');ylabel('dPS');xlim([offset_min-100,offset_max-100]);
% subplot(1,2,2);
% plot(bin_psppratio_syn2510(:,1),bin_psppratio_syn2510(:,6),'o-');hold on;
% plot(bin_psppratio_syn2515(:,1),bin_psppratio_syn2515(:,6),'o-');hold on;
% plot(bin_psppratio_syn2520(:,1),bin_psppratio_syn2520(:,6),'o-');hold on;
% plot(bin_psppratio_syn2525(:,1),bin_psppratio_syn2525(:,6),'o-');hold on;
% legend('25-10','25-15','25-20','25-25');
% xlabel('offset (m)');ylabel('dPP');xlim([offset_min-100,offset_max-100]);

figure;
plot(bin_psppratio_syn3010(:,1),bin_psppratio_syn3010(:,2),'o-');hold on;
plot(bin_psppratio_syn3015(:,1),bin_psppratio_syn3015(:,2),'o-');hold on;
plot(bin_psppratio_syn3020(:,1),bin_psppratio_syn3020(:,2),'o-');hold on;
plot(bin_psppratio_syn3025(:,1),bin_psppratio_syn3025(:,2),'o-');hold on;
plot(bin_psppratio_syn3030(:,1),bin_psppratio_syn3030(:,2),'o-');hold on;
errorbar(bin_psppratio(:,1),bin_psppratio(:,2),std_bsp_psppratio,'ko');legend('30-10','30-15','30-20','30-25','30-30');
xlabel('offset (m)');ylabel('PS/PP ratio');xlim([offset_min-100,offset_max-100]);
% figure;subplot(1,2,1);
% plot(bin_psppratio_syn3010(:,1),bin_psppratio_syn3010(:,5),'o-');hold on;
% plot(bin_psppratio_syn3015(:,1),bin_psppratio_syn3015(:,5),'o-');hold on;
% plot(bin_psppratio_syn3020(:,1),bin_psppratio_syn3020(:,5),'o-');hold on;
% plot(bin_psppratio_syn3025(:,1),bin_psppratio_syn3025(:,5),'o-');hold on;
% plot(bin_psppratio_syn3030(:,1),bin_psppratio_syn3030(:,5),'o-');hold on;
% legend('30-10','30-15','30-20','30-25','30-30');
% xlabel('offset (m)');ylabel('dPS');xlim([offset_min-100,offset_max-100]);
% subplot(1,2,2);
% plot(bin_psppratio_syn3010(:,1),bin_psppratio_syn3010(:,6),'o-');hold on;
% plot(bin_psppratio_syn3015(:,1),bin_psppratio_syn3015(:,6),'o-');hold on;
% plot(bin_psppratio_syn3020(:,1),bin_psppratio_syn3020(:,6),'o-');hold on;
% plot(bin_psppratio_syn3025(:,1),bin_psppratio_syn3025(:,6),'o-');hold on;
% plot(bin_psppratio_syn3030(:,1),bin_psppratio_syn3030(:,6),'o-');hold on;
% legend('30-10','30-15','30-20','30-25','30-30');
% xlabel('offset (m)');ylabel('dPP');xlim([offset_min-100,offset_max-100]);

figure;
plot(bin_psppratio_syn3510(:,1),bin_psppratio_syn3510(:,2),'o-');hold on;
plot(bin_psppratio_syn3515(:,1),bin_psppratio_syn3515(:,2),'o-');hold on;
plot(bin_psppratio_syn3520(:,1),bin_psppratio_syn3520(:,2),'o-');hold on;
plot(bin_psppratio_syn3525(:,1),bin_psppratio_syn3525(:,2),'o-');hold on;
plot(bin_psppratio_syn3530(:,1),bin_psppratio_syn3530(:,2),'o-');hold on;
plot(bin_psppratio_syn3535(:,1),bin_psppratio_syn3535(:,2),'o-');hold on;
errorbar(bin_psppratio(:,1),bin_psppratio(:,2),std_bsp_psppratio,'ko');legend('35-10','35-15','35-20','35-25','35-30','35-35');
xlabel('offset (m)');ylabel('PS/PP ratio');xlim([offset_min-100,offset_max-100]);
% figure;subplot(1,2,1);
% plot(bin_psppratio_syn3510(:,1),bin_psppratio_syn3510(:,5),'o-');hold on;
% plot(bin_psppratio_syn3515(:,1),bin_psppratio_syn3515(:,5),'o-');hold on;
% plot(bin_psppratio_syn3520(:,1),bin_psppratio_syn3520(:,5),'o-');hold on;
% plot(bin_psppratio_syn3525(:,1),bin_psppratio_syn3525(:,5),'o-');hold on;
% plot(bin_psppratio_syn3530(:,1),bin_psppratio_syn3530(:,5),'o-');hold on;
% plot(bin_psppratio_syn3535(:,1),bin_psppratio_syn3535(:,5),'o-');hold on;
% legend('35-10','35-15','35-20','35-25','35-30','35-35');
% xlabel('offset (m)');ylabel('dPS');xlim([offset_min-100,offset_max-100]);
% subplot(1,2,2);
% plot(bin_psppratio_syn3510(:,1),bin_psppratio_syn3510(:,6),'o-');hold on;
% plot(bin_psppratio_syn3515(:,1),bin_psppratio_syn3515(:,6),'o-');hold on;
% plot(bin_psppratio_syn3520(:,1),bin_psppratio_syn3520(:,6),'o-');hold on;
% plot(bin_psppratio_syn3525(:,1),bin_psppratio_syn3525(:,6),'o-');hold on;
% plot(bin_psppratio_syn3530(:,1),bin_psppratio_syn3530(:,6),'o-');hold on;
% plot(bin_psppratio_syn3535(:,1),bin_psppratio_syn3535(:,6),'o-');hold on;
% legend('35-10','35-15','35-20','35-25','35-30','35-35');
% xlabel('offset (m)');ylabel('dPP');xlim([offset_min-100,offset_max-100]);

figure;
plot(bin_psppratio_syn4010(:,1),bin_psppratio_syn4010(:,2),'o-');hold on;
plot(bin_psppratio_syn4015(:,1),bin_psppratio_syn4015(:,2),'o-');hold on;
plot(bin_psppratio_syn4020(:,1),bin_psppratio_syn4020(:,2),'o-');hold on;
plot(bin_psppratio_syn4025(:,1),bin_psppratio_syn4025(:,2),'o-');hold on;
plot(bin_psppratio_syn4030(:,1),bin_psppratio_syn4030(:,2),'o-');hold on;
plot(bin_psppratio_syn4035(:,1),bin_psppratio_syn4035(:,2),'o-');hold on;
plot(bin_psppratio_syn4040(:,1),bin_psppratio_syn4040(:,2),'o-');hold on;
errorbar(bin_psppratio(:,1),bin_psppratio(:,2),std_bsp_psppratio,'ko');legend('40-10','40-15','40-20','40-25','40-30','40-35','40-40');
xlabel('offset (m)');ylabel('PS/PP ratio');xlim([offset_min-100,offset_max-100]);
% figure;subplot(1,2,1);
% plot(bin_psppratio_syn4010(:,1),bin_psppratio_syn4010(:,5),'o-');hold on;
% plot(bin_psppratio_syn4015(:,1),bin_psppratio_syn4015(:,5),'o-');hold on;
% plot(bin_psppratio_syn4020(:,1),bin_psppratio_syn4020(:,5),'o-');hold on;
% plot(bin_psppratio_syn4025(:,1),bin_psppratio_syn4025(:,5),'o-');hold on;
% plot(bin_psppratio_syn4030(:,1),bin_psppratio_syn4030(:,5),'o-');hold on;
% plot(bin_psppratio_syn4035(:,1),bin_psppratio_syn4035(:,5),'o-');hold on;
% plot(bin_psppratio_syn4040(:,1),bin_psppratio_syn4040(:,5),'o-');hold on;
% legend('40-10','40-15','40-20','40-25','40-30','40-35','40-40');
% xlabel('offset (m)');ylabel('dPS');xlim([offset_min-100,offset_max-100]);
% subplot(1,2,2);
% plot(bin_psppratio_syn4010(:,1),bin_psppratio_syn4010(:,6),'o-');hold on;
% plot(bin_psppratio_syn4015(:,1),bin_psppratio_syn4015(:,6),'o-');hold on;
% plot(bin_psppratio_syn4020(:,1),bin_psppratio_syn4020(:,6),'o-');hold on;
% plot(bin_psppratio_syn4025(:,1),bin_psppratio_syn4025(:,6),'o-');hold on;
% plot(bin_psppratio_syn4030(:,1),bin_psppratio_syn4030(:,6),'o-');hold on;
% plot(bin_psppratio_syn4035(:,1),bin_psppratio_syn4035(:,6),'o-');hold on;
% plot(bin_psppratio_syn4040(:,1),bin_psppratio_syn4040(:,6),'o-');hold on;
% legend('40-10','40-15','40-20','40-25','40-30','40-35','40-40');
% xlabel('offset (m)');ylabel('dPP');xlim([offset_min-100,offset_max-100]);

figure;
plot(bin_psppratio_syn4510(:,1),bin_psppratio_syn4510(:,2),'o-');hold on;
plot(bin_psppratio_syn4515(:,1),bin_psppratio_syn4515(:,2),'o-');hold on;
plot(bin_psppratio_syn4520(:,1),bin_psppratio_syn4520(:,2),'o-');hold on;
plot(bin_psppratio_syn4525(:,1),bin_psppratio_syn4525(:,2),'o-');hold on;
plot(bin_psppratio_syn4530(:,1),bin_psppratio_syn4530(:,2),'o-');hold on;
plot(bin_psppratio_syn4535(:,1),bin_psppratio_syn4535(:,2),'o-');hold on;
plot(bin_psppratio_syn4540(:,1),bin_psppratio_syn4540(:,2),'o-');hold on;
plot(bin_psppratio_syn4545(:,1),bin_psppratio_syn4545(:,2),'o-');hold on;
errorbar(bin_psppratio(:,1),bin_psppratio(:,2),std_bsp_psppratio,'ko');legend('45-10','45-15','45-20','45-25','45-30','45-35','45-40','45-45');
xlabel('offset (m)');ylabel('PS/PP ratio');xlim([offset_min-100,offset_max-100]);
% figure;subplot(1,2,1);
% plot(bin_psppratio_syn4510(:,1),bin_psppratio_syn4510(:,5),'o-');hold on;
% plot(bin_psppratio_syn4515(:,1),bin_psppratio_syn4515(:,5),'o-');hold on;
% plot(bin_psppratio_syn4520(:,1),bin_psppratio_syn4520(:,5),'o-');hold on;
% plot(bin_psppratio_syn4525(:,1),bin_psppratio_syn4525(:,5),'o-');hold on;
% plot(bin_psppratio_syn4530(:,1),bin_psppratio_syn4530(:,5),'o-');hold on;
% plot(bin_psppratio_syn4535(:,1),bin_psppratio_syn4535(:,5),'o-');hold on;
% plot(bin_psppratio_syn4540(:,1),bin_psppratio_syn4540(:,5),'o-');hold on;
% plot(bin_psppratio_syn4545(:,1),bin_psppratio_syn4545(:,5),'o-');hold on;
% legend('45-10','45-15','45-20','45-25','45-30','45-35','45-40','45-45');
% xlabel('offset (m)');ylabel('dPS');xlim([offset_min-100,offset_max-100]);
% subplot(1,2,2);
% plot(bin_psppratio_syn4510(:,1),bin_psppratio_syn4510(:,6),'o-');hold on;
% plot(bin_psppratio_syn4515(:,1),bin_psppratio_syn4515(:,6),'o-');hold on;
% plot(bin_psppratio_syn4520(:,1),bin_psppratio_syn4520(:,6),'o-');hold on;
% plot(bin_psppratio_syn4525(:,1),bin_psppratio_syn4525(:,6),'o-');hold on;
% plot(bin_psppratio_syn4530(:,1),bin_psppratio_syn4530(:,6),'o-');hold on;
% plot(bin_psppratio_syn4535(:,1),bin_psppratio_syn4535(:,6),'o-');hold on;
% plot(bin_psppratio_syn4540(:,1),bin_psppratio_syn4540(:,6),'o-');hold on;
% plot(bin_psppratio_syn4545(:,1),bin_psppratio_syn4545(:,6),'o-');hold on;
% legend('45-10','45-15','45-20','45-25','45-30','45-35','45-40','45-45');
% xlabel('offset (m)');ylabel('dPP');xlim([offset_min-100,offset_max-100]);

figure;
plot(bin_psppratio_syn5015(:,1),bin_psppratio_syn5015(:,2),'o-');hold on;
plot(bin_psppratio_syn5020(:,1),bin_psppratio_syn5020(:,2),'o-');hold on;
plot(bin_psppratio_syn5025(:,1),bin_psppratio_syn5025(:,2),'o-');hold on;
plot(bin_psppratio_syn5030(:,1),bin_psppratio_syn5030(:,2),'o-');hold on;
plot(bin_psppratio_syn5035(:,1),bin_psppratio_syn5035(:,2),'o-');hold on;
plot(bin_psppratio_syn5040(:,1),bin_psppratio_syn5040(:,2),'o-');hold on;
plot(bin_psppratio_syn5045(:,1),bin_psppratio_syn5045(:,2),'o-');hold on;
plot(bin_psppratio_syn5050(:,1),bin_psppratio_syn5050(:,2),'o-');hold on;
errorbar(bin_psppratio(:,1),bin_psppratio(:,2),std_bsp_psppratio,'ko');legend('50-15','50-20','50-25','50-30','50-35','50-40','50-45','50-50');
xlabel('offset (m)');ylabel('PS/PP ratio');xlim([offset_min-100,offset_max-100]);
% figure;subplot(1,2,1);
% plot(bin_psppratio_syn5015(:,1),bin_psppratio_syn5015(:,5),'o-');hold on;
% plot(bin_psppratio_syn5020(:,1),bin_psppratio_syn5020(:,5),'o-');hold on;
% plot(bin_psppratio_syn5025(:,1),bin_psppratio_syn5025(:,5),'o-');hold on;
% plot(bin_psppratio_syn5030(:,1),bin_psppratio_syn5030(:,5),'o-');hold on;
% plot(bin_psppratio_syn5035(:,1),bin_psppratio_syn5035(:,5),'o-');hold on;
% plot(bin_psppratio_syn5040(:,1),bin_psppratio_syn5040(:,5),'o-');hold on;
% plot(bin_psppratio_syn5045(:,1),bin_psppratio_syn5045(:,5),'o-');hold on;
% plot(bin_psppratio_syn5050(:,1),bin_psppratio_syn5050(:,5),'o-');hold on;
% legend('50-15','50-20','50-25','50-30','50-35','50-40','50-45','50-50');
% xlabel('offset (m)');ylabel('dPS');xlim([offset_min-100,offset_max-100]);
% subplot(1,2,2);
% plot(bin_psppratio_syn5015(:,1),bin_psppratio_syn5015(:,6),'o-');hold on;
% plot(bin_psppratio_syn5020(:,1),bin_psppratio_syn5020(:,6),'o-');hold on;
% plot(bin_psppratio_syn5025(:,1),bin_psppratio_syn5025(:,6),'o-');hold on;
% plot(bin_psppratio_syn5030(:,1),bin_psppratio_syn5030(:,6),'o-');hold on;
% plot(bin_psppratio_syn5035(:,1),bin_psppratio_syn5035(:,6),'o-');hold on;
% plot(bin_psppratio_syn5040(:,1),bin_psppratio_syn5040(:,6),'o-');hold on;
% plot(bin_psppratio_syn5045(:,1),bin_psppratio_syn5045(:,6),'o-');hold on;
% plot(bin_psppratio_syn5050(:,1),bin_psppratio_syn5050(:,6),'o-');hold on;
% legend('50-15','50-20','50-25','50-30','50-35','50-40','50-45','50-50');
% xlabel('offset (m)');ylabel('dPP');xlim([offset_min-100,offset_max-100]);

figure;
plot(bin_psppratio_syn5520(:,1),bin_psppratio_syn5520(:,2),'o-');hold on;
plot(bin_psppratio_syn5525(:,1),bin_psppratio_syn5525(:,2),'o-');hold on;
plot(bin_psppratio_syn5530(:,1),bin_psppratio_syn5530(:,2),'o-');hold on;
plot(bin_psppratio_syn5535(:,1),bin_psppratio_syn5535(:,2),'o-');hold on;
plot(bin_psppratio_syn5540(:,1),bin_psppratio_syn5540(:,2),'o-');hold on;
plot(bin_psppratio_syn5545(:,1),bin_psppratio_syn5545(:,2),'o-');hold on;
plot(bin_psppratio_syn5550(:,1),bin_psppratio_syn5550(:,2),'o-');hold on;
plot(bin_psppratio_syn5555(:,1),bin_psppratio_syn5555(:,2),'o-');hold on;
errorbar(bin_psppratio(:,1),bin_psppratio(:,2),std_bsp_psppratio,'ko');legend('55-20','55-25','55-30','55-35','55-40','55-45','55-50','55-55');
xlabel('offset (m)');ylabel('PS/PP ratio');xlim([offset_min-100,offset_max-100]);
% figure;subplot(1,2,1);
% plot(bin_psppratio_syn5520(:,1),bin_psppratio_syn5520(:,5),'o-');hold on;
% plot(bin_psppratio_syn5525(:,1),bin_psppratio_syn5525(:,5),'o-');hold on;
% plot(bin_psppratio_syn5530(:,1),bin_psppratio_syn5530(:,5),'o-');hold on;
% plot(bin_psppratio_syn5535(:,1),bin_psppratio_syn5535(:,5),'o-');hold on;
% plot(bin_psppratio_syn5540(:,1),bin_psppratio_syn5540(:,5),'o-');hold on;
% plot(bin_psppratio_syn5545(:,1),bin_psppratio_syn5545(:,5),'o-');hold on;
% plot(bin_psppratio_syn5550(:,1),bin_psppratio_syn5550(:,5),'o-');hold on;
% plot(bin_psppratio_syn5555(:,1),bin_psppratio_syn5555(:,5),'o-');hold on;
% legend('55-20','55-25','55-30','55-35','55-40','55-45','55-50','55-55');
% xlabel('offset (m)');ylabel('dPS');xlim([offset_min-100,offset_max-100]);
% subplot(1,2,2);
% plot(bin_psppratio_syn5520(:,1),bin_psppratio_syn5520(:,6),'o-');hold on;
% plot(bin_psppratio_syn5525(:,1),bin_psppratio_syn5525(:,6),'o-');hold on;
% plot(bin_psppratio_syn5530(:,1),bin_psppratio_syn5530(:,6),'o-');hold on;
% plot(bin_psppratio_syn5535(:,1),bin_psppratio_syn5535(:,6),'o-');hold on;
% plot(bin_psppratio_syn5540(:,1),bin_psppratio_syn5540(:,6),'o-');hold on;
% plot(bin_psppratio_syn5545(:,1),bin_psppratio_syn5545(:,6),'o-');hold on;
% plot(bin_psppratio_syn5550(:,1),bin_psppratio_syn5550(:,6),'o-');hold on;
% plot(bin_psppratio_syn5555(:,1),bin_psppratio_syn5555(:,6),'o-');hold on;
% legend('55-20','55-25','55-30','55-35','55-40','55-45','55-50','55-55');
% xlabel('offset (m)');ylabel('dPP');xlim([offset_min-100,offset_max-100]);

figure;
plot(bin_psppratio_syn6030(:,1),bin_psppratio_syn6030(:,2),'o-');hold on;
plot(bin_psppratio_syn6035(:,1),bin_psppratio_syn6035(:,2),'o-');hold on;
plot(bin_psppratio_syn6040(:,1),bin_psppratio_syn6040(:,2),'o-');hold on;
plot(bin_psppratio_syn6045(:,1),bin_psppratio_syn6045(:,2),'o-');hold on;
plot(bin_psppratio_syn6050(:,1),bin_psppratio_syn6050(:,2),'o-');hold on;
plot(bin_psppratio_syn6055(:,1),bin_psppratio_syn6055(:,2),'o-');hold on;
plot(bin_psppratio_syn6060(:,1),bin_psppratio_syn6060(:,2),'o-');hold on;
errorbar(bin_psppratio(:,1),bin_psppratio(:,2),std_bsp_psppratio,'ko');legend('60-30','60-35','60-40','60-45','60-50','60-55','60-60');
xlabel('offset (m)');ylabel('PS/PP ratio');xlim([offset_min-100,offset_max-100]);

figure;
plot(bin_psppratio_syn6540(:,1),bin_psppratio_syn6540(:,2),'o-');hold on;
plot(bin_psppratio_syn6545(:,1),bin_psppratio_syn6545(:,2),'o-');hold on;
plot(bin_psppratio_syn6550(:,1),bin_psppratio_syn6550(:,2),'o-');hold on;
plot(bin_psppratio_syn6555(:,1),bin_psppratio_syn6555(:,2),'o-');hold on;
plot(bin_psppratio_syn6560(:,1),bin_psppratio_syn6560(:,2),'o-');hold on;
plot(bin_psppratio_syn6565(:,1),bin_psppratio_syn6565(:,2),'o-');hold on;
errorbar(bin_psppratio(:,1),bin_psppratio(:,2),std_bsp_psppratio,'ko');legend('65-40','65-45','65-50','65-55','65-60','65-65');
xlabel('offset (m)');ylabel('PS/PP ratio');xlim([offset_min-100,offset_max-100]);

figure;
plot(bin_psppratio_syn7050(:,1),bin_psppratio_syn7050(:,2),'o-');hold on;
plot(bin_psppratio_syn7055(:,1),bin_psppratio_syn7055(:,2),'o-');hold on;
plot(bin_psppratio_syn7060(:,1),bin_psppratio_syn7060(:,2),'o-');hold on;
plot(bin_psppratio_syn7065(:,1),bin_psppratio_syn7065(:,2),'o-');hold on;
plot(bin_psppratio_syn7070(:,1),bin_psppratio_syn7070(:,2),'o-');hold on;
errorbar(bin_psppratio(:,1),bin_psppratio(:,2),std_bsp_psppratio,'ko');legend('70-50','70-55','70-60','70-65','70-70');
xlabel('offset (m)');ylabel('PS/PP ratio');xlim([offset_min-100,offset_max-100]);

% good collection
figure;
plot(bin_psppratio_syn2015(:,1),bin_psppratio_syn2015(:,2),'o-');hold on;
plot(bin_psppratio_syn3525(:,1),bin_psppratio_syn3525(:,2),'o-');hold on;
plot(bin_psppratio_syn4030(:,1),bin_psppratio_syn4030(:,2),'o-');hold on;
plot(bin_psppratio_syn5035(:,1),bin_psppratio_syn5035(:,2),'o-');hold on;
errorbar(bin_psppratio(:,1),bin_psppratio(:,2),std_bsp_psppratio,'ko');legend('20-15','35-25','40-30','50-35');
xlabel('offset (m)');ylabel('PS/PP ratio');xlim([offset_min-100,offset_max-100]);
% figure;subplot(1,2,1);
% plot(bin_psppratio_syn2015(:,1),bin_psppratio_syn2015(:,5),'o-');hold on;
% plot(bin_psppratio_syn3525(:,1),bin_psppratio_syn3525(:,5),'o-');hold on;
% plot(bin_psppratio_syn4030(:,1),bin_psppratio_syn4030(:,5),'o-');hold on;
% plot(bin_psppratio_syn5035(:,1),bin_psppratio_syn5035(:,5),'o-');hold on;
% legend('20-15','35-25','40-30','50-35');
% xlabel('offset (m)');ylabel('dPS');xlim([offset_min-100,offset_max-100]);
% subplot(1,2,2);
% plot(bin_psppratio_syn2015(:,1),bin_psppratio_syn2015(:,6),'o-');hold on;
% plot(bin_psppratio_syn3525(:,1),bin_psppratio_syn3525(:,6),'o-');hold on;
% plot(bin_psppratio_syn4030(:,1),bin_psppratio_syn4030(:,6),'o-');hold on;
% plot(bin_psppratio_syn5035(:,1),bin_psppratio_syn5035(:,6),'o-');hold on;
% legend('20-15','35-25','40-30','50-35');
% xlabel('offset (m)');ylabel('dPP');xlim([offset_min-100,offset_max-100]);

% bad results
figure;
plot(bin_psppratio_syn3050(:,1),bin_psppratio_syn3050(:,2),'o-');hold on;
plot(bin_psppratio_syn1030(:,1),bin_psppratio_syn1030(:,2),'o-');hold on;
plot(bin_psppratio_syn1050(:,1),bin_psppratio_syn1050(:,2),'o-');hold on;
errorbar(bin_psppratio(:,1),bin_psppratio(:,2),std_bsp_psppratio,'ko');
xlabel('offset (m)');ylabel('PS/PP ratio');xlim([offset_min-100,offset_max-100]);

%% compute RMS and Chi-square misfit
dvp=[10,15,15,20,20,20,25,25,25,25,30,30,30,30,30,35,35,35,35,35,35,40,40,40,40,40,40,40,...
    45,45,45,45,45,45,45,45,50,50,50,50,50,50,50,50,55,55,55,55,55,55,55,55,60,60,60,60,60,60,60,...
    65,65,65,65,65,65,70,70,70,70,70];
dvs=[10,10,15,10,15,20,10,15,20,25,10,15,20,25,30,10,15,20,25,30,35,10,15,20,25,30,35,40,...
    10,15,20,25,30,35,40,45,15,20,25,30,35,40,45,50,20,25,30,35,40,45,50,55,30,35,40,45,50,55,60,...
    40,45,50,55,60,65,50,55,60,65,70];
figure;plot(dvp,dvs,'bo');hold on;
plot([10,15,20,25,30,35,40,45,50,55,60,65,70],[10,15,20,25,30,35,40,45,50,55,60,65,70],'k--');axis equal;
xlim([10,72]);ylim([8,70]);xlabel('dvp/vp (%)');ylabel('dvs/vs (%)');grid on;

N_model = length(dvp);
misfit_tab = zeros(N_model,5);

im = 1;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn1010(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn1010(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 2;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn1510(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn1510(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 3;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn1515(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn1515(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 4;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn2010(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn2010(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 5;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn2015(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn2015(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 6;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn2020(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn2020(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 7;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn2510(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn2510(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 8;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn2515(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn2515(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 9;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn2520(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn2520(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 10;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn2525(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn2525(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 11;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn3010(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn3010(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 12;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn3015(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn3015(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 13;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn3020(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn3020(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 14;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn3025(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn3025(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=15;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn3030(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn3030(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=16;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn3510(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn3510(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=17;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn3515(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn3515(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=18;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn3520(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn3520(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=19;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn3525(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn3525(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=20;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn3530(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn3530(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=21;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn3535(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn3535(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=22;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn4010(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn4010(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=23;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn4015(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn4015(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=24;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn4020(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn4020(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=25;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn4025(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn4025(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=26;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn4030(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn4030(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=27;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn4035(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn4035(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=28;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn4040(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn4040(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=29;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn4510(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn4510(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=30;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn4515(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn4515(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=31;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn4520(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn4520(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=32;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn4525(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn4525(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=33;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn4530(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn4530(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=34;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn4535(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn4535(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=35;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn4540(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn4540(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=36;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn4545(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn4545(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=37;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn5015(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn5015(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=38;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn5020(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn5020(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=39;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn5025(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn5025(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=40;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn5030(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn5030(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=41;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn5035(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn5035(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=42;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn5040(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn5040(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=43;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn5045(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn5045(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=44;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn5050(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn5050(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=45;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn5520(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn5520(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=46;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn5525(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn5525(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=47;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn5530(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn5530(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=48;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn5535(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn5535(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=49;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn5540(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn5540(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=50;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn5545(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn5545(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=51;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn5550(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn5550(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=52;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn5555(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn5555(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=53;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn6030(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn6030(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=54;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn6035(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn6035(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=55;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn6040(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn6040(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=56;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn6045(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn6045(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=57;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn6050(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn6050(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=58;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn6055(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn6055(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=59;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn6060(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn6060(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=60;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn6540(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn6540(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=61;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn6545(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn6545(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=62;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn6550(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn6550(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=63;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn6555(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn6555(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=64;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn6560(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn6560(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=65;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn6565(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn6565(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=66;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn7050(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn7050(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=67;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn7055(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn7055(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=68;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn7060(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn7060(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=69;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn7065(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn7065(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im=70;
fprintf('[dvp,dvs]=[%d, %d]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_syn7070(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_syn7070(:,2), bin_psppratio(:,2), std_bsp_psppratio);

% plot misfits
figure;
semilogx(misfit_tab(:,3),misfit_tab(:,4),'bo');
xlim([1.01,4.55]);
xlabel('\deltaVp/\deltaVs');ylabel('RMS misfit');

figure;
plot(misfit_tab(:,3),misfit_tab(:,5),'bo');
xlim([1.01,3.5]);
xlabel('\deltaVp/\deltaVs');ylabel('Chi-square misfit');

%% output for Python
% 1. 1D misfit curve------
% sort by dvp/dvs, and exclude dvp/dvs=1
misfit_tab_1 = sortrows(misfit_tab,3);
misfit_tab_2 = misfit_tab_1(misfit_tab_1(:,3)>1,:);

% (optional)
fp = fopen('/home0/cxd170430/NatureGeo/avo/misfit_rms_chi.bin','wb');
fwrite(fp,misfit_tab_2,'float32');
fclose(fp);

% 2. 2D misfit matrix------
ox = 10; dx = 5; ex = 70;
oy = 10; dy = 5; ey = 70;
nx = length(ox:dx:ex);
ny = length(oy:dy:ey);
rms_mat = NaN(ny,nx);
chi_mat = NaN(ny,nx);
for i=1:N_model
    ix = (misfit_tab(i,1)-ox)/dx+1;
    iy = ny - (misfit_tab(i,2)-oy)/dy;
    rms_mat(iy,ix) = misfit_tab(i,4);
    chi_mat(iy,ix) = misfit_tab(i,5);
end

% (optional)
fp = fopen('/home0/cxd170430/NatureGeo/avo/rms_mat.bin','wb');
fwrite(fp,rms_mat,'float32');
fclose(fp);

fp = fopen('/home0/cxd170430/NatureGeo/avo/chi_mat.bin','wb');
fwrite(fp,chi_mat,'float32');
fclose(fp);







