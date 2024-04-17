
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




%% 6. time window (synthetic)
syn_timewin_p_prior = 0.28;
syn_timewin_p_post = 0.3;
syn_timewin_pp_prior = 0.05;
syn_timewin_pp_post = 0.3;
syn_timewin_ps_prior = 0.08;
syn_timewin_ps_post = 0.45;

% fluid
su_ez_f1 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_2.081_1.273_32bit.su';
su_ex_f1 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_2.081_1.273_32bit.su';

su_ez_f2 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_6.106_3.9_32bit.su';
su_ex_f2 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_6.106_3.9_32bit.su';

su_ez_f3 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_15.687_10.998_32bit.su';
su_ex_f3 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_15.687_10.998_32bit.su';

su_ez_f4 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_19.44_14.078_32bit.su';
su_ex_f4 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_19.44_14.078_32bit.su';

su_ez_f5 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_23.204_17.32_32bit.su';
su_ex_f5 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_23.204_17.32_32bit.su';

su_ez_f6 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_27.014_20.746_32bit.su';
su_ex_f6 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_27.014_20.746_32bit.su';

su_ez_f7 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_30.907_24.385_32bit.su';
su_ex_f7 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_30.907_24.385_32bit.su';

su_ez_f8 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_39.113_32.46_32bit.su';
su_ex_f8 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_39.113_32.46_32bit.su';

su_ez_f9 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_50.796_44.746_32bit.su';
su_ex_f9 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_50.796_44.746_32bit.su';

su_ez_f10 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_65.994_61.838_32bit.su';
su_ex_f10 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_65.994_61.838_32bit.su';

[~, ~, ez_f1, ex_f1, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_f1, su_ex_f1, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_f2, ex_f2, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_f2, su_ex_f2, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_f3, ex_f3, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_f3, su_ex_f3, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_f4, ex_f4, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_f4, su_ex_f4, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_f5, ex_f5, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_f5, su_ex_f5, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_f6, ex_f6, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_f6, su_ex_f6, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_f7, ex_f7, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_f7, su_ex_f7, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_f8, ex_f8, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_f8, su_ex_f8, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_f9, ex_f9, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_f9, su_ex_f9, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_f10, ex_f10, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_f10, su_ex_f10, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);

% melt
su_ez_m1 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_1.627_1.703_32bit.su';
su_ex_m1 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_1.627_1.703_32bit.su';

su_ez_m2 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_4.849_5.164_32bit.su';
su_ex_m2 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_4.849_5.164_32bit.su';

su_ez_m3 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_12.771_14.157_32bit.su';
su_ex_m3 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_12.771_14.157_32bit.su';

su_ez_m4 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_15.91_17.909_32bit.su';
su_ex_m4 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_15.91_17.909_32bit.su';

su_ez_m5 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_19.043_21.767_32bit.su';
su_ex_m5 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_19.043_21.767_32bit.su';

su_ez_m6 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_22.179_25.745_32bit.su';
su_ex_m6 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_22.179_25.745_32bit.su';

su_ez_m7 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_25.327_29.864_32bit.su';
su_ex_m7 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_25.327_29.864_32bit.su';

su_ez_m8 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_31.693_38.644_32bit.su';
su_ex_m8 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_31.693_38.644_32bit.su';

su_ez_m9 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_39.897_51.16_32bit.su';
su_ex_m9 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_39.897_51.16_32bit.su';

su_ez_m10 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_48.601_67.232_32bit.su';
su_ex_m10 = '/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_48.601_67.232_32bit.su';

[~, ~, ez_m1, ex_m1, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_m1, su_ex_m1, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_m2, ex_m2, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_m2, su_ex_m2, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_m3, ex_m3, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_m3, su_ex_m3, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_m4, ex_m4, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_m4, su_ex_m4, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_m5, ex_m5, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_m5, su_ex_m5, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_m6, ex_m6, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_m6, su_ex_m6, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_m7, ex_m7, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_m7, su_ex_m7, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_m8, ex_m8, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_m8, su_ex_m8, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_m9, ex_m9, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_m9, su_ex_m9, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);
[~, ~, ez_m10, ex_m10, ~] = ...
load_synenv_so(su_fileZ_env, su_fileR_env, su_ez_m10, su_ex_m10, 'Z', yn_snr_qc, snr_threshold, offset_min, offset_max);



%% 7 compute AVO - trace by trace (synthetic)
% average offset range=2~6km by bins (bins: [2, 2.1] (2.1, 2.2] ... (5.9, 6] (6, 6.1])

% ====PS/PP====
% fluid
ez_ex_f1 = ez_f1 + ex_f1;
bin_psppratio_synf1 = amp_ratio(ez_ex_f1,ez_ex_f1,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_f2 = ez_f2 + ex_f2;
bin_psppratio_synf2 = amp_ratio(ez_ex_f2,ez_ex_f2,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_f3 = ez_f3 + ex_f3;
bin_psppratio_synf3 = amp_ratio(ez_ex_f3,ez_ex_f3,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_f4 = ez_f4 + ex_f4;
bin_psppratio_synf4 = amp_ratio(ez_ex_f4,ez_ex_f4,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_f5 = ez_f5 + ex_f5;
bin_psppratio_synf5 = amp_ratio(ez_ex_f5,ez_ex_f5,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_f6 = ez_f6 + ex_f6;
bin_psppratio_synf6 = amp_ratio(ez_ex_f6,ez_ex_f6,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_f7 = ez_f7 + ex_f7;
bin_psppratio_synf7 = amp_ratio(ez_ex_f7,ez_ex_f7,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_f8 = ez_f8 + ex_f8;
bin_psppratio_synf8 = amp_ratio(ez_ex_f8,ez_ex_f8,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_f9 = ez_f9 + ex_f9;
bin_psppratio_synf9 = amp_ratio(ez_ex_f9,ez_ex_f9,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_f10 = ez_f10 + ex_f10;
bin_psppratio_synf10 = amp_ratio(ez_ex_f10,ez_ex_f10,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

% melt
ez_ex_m1 = ez_m1 + ex_m1;
bin_psppratio_synm1 = amp_ratio(ez_ex_m1,ez_ex_m1,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_m2 = ez_m2 + ex_m2;
bin_psppratio_synm2 = amp_ratio(ez_ex_m2,ez_ex_m2,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_m3 = ez_m3 + ex_m3;
bin_psppratio_synm3 = amp_ratio(ez_ex_m3,ez_ex_m3,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_m4 = ez_m4 + ex_m4;
bin_psppratio_synm4 = amp_ratio(ez_ex_m4,ez_ex_m4,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_m5 = ez_m5 + ex_m5;
bin_psppratio_synm5 = amp_ratio(ez_ex_m5,ez_ex_m5,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_m6 = ez_m6 + ex_m6;
bin_psppratio_synm6 = amp_ratio(ez_ex_m6,ez_ex_m6,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_m7 = ez_m7 + ex_m7;
bin_psppratio_synm7 = amp_ratio(ez_ex_m7,ez_ex_m7,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_m8 = ez_m8 + ex_m8;
bin_psppratio_synm8 = amp_ratio(ez_ex_m8,ez_ex_m8,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_m9 = ez_m9 + ex_m9;
bin_psppratio_synm9 = amp_ratio(ez_ex_m9,ez_ex_m9,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_m10 = ez_m10 + ex_m10;
bin_psppratio_synm10 = amp_ratio(ez_ex_m10,ez_ex_m10,geom_tab_new,dt,t_ps,t_pp,syn_timewin_ps_prior,syn_timewin_ps_post,...
    syn_timewin_pp_prior,syn_timewin_pp_post,offset_min,offset_bin_length,offset_max,'max',2e-14);


%% 8. final plots
% fluid (good)
figure;
plot(bin_psppratio_synf1(:,1),bin_psppratio_synf1(:,2),'o-');hold on;
plot(bin_psppratio_synf2(:,1),bin_psppratio_synf2(:,2),'o-');hold on;
plot(bin_psppratio_synf3(:,1),bin_psppratio_synf3(:,2),'o-');hold on;
plot(bin_psppratio_synf4(:,1),bin_psppratio_synf4(:,2),'o-');hold on;
plot(bin_psppratio_synf5(:,1),bin_psppratio_synf5(:,2),'o-');hold on;
plot(bin_psppratio_synf6(:,1),bin_psppratio_synf6(:,2),'o-');hold on;
plot(bin_psppratio_synf7(:,1),bin_psppratio_synf7(:,2),'o-');hold on;
plot(bin_psppratio_synf8(:,1),bin_psppratio_synf8(:,2),'o-');hold on;
plot(bin_psppratio_synf9(:,1),bin_psppratio_synf9(:,2),'o-');hold on;
plot(bin_psppratio_synf10(:,1),bin_psppratio_synf10(:,2),'o-');hold on;
errorbar(bin_psppratio(:,1),bin_psppratio(:,2),std_bsp_psppratio,'ko');legend('0.01','0.03','0.08','0.1','0.12','0.14','0.16','0.2','0.25','0.3');
xlabel('offset (m)');ylabel('PS/PP ratio');xlim([offset_min-100,offset_max-100]);

% melt (bad)
figure;
plot(bin_psppratio_synm1(:,1),bin_psppratio_synm1(:,2),'o-');hold on;
plot(bin_psppratio_synm2(:,1),bin_psppratio_synm2(:,2),'o-');hold on;
plot(bin_psppratio_synm3(:,1),bin_psppratio_synm3(:,2),'o-');hold on;
plot(bin_psppratio_synm4(:,1),bin_psppratio_synm4(:,2),'o-');hold on;
plot(bin_psppratio_synm5(:,1),bin_psppratio_synm5(:,2),'o-');hold on;
plot(bin_psppratio_synm6(:,1),bin_psppratio_synm6(:,2),'o-');hold on;
plot(bin_psppratio_synm7(:,1),bin_psppratio_synm7(:,2),'o-');hold on;
plot(bin_psppratio_synm8(:,1),bin_psppratio_synm8(:,2),'o-');hold on;
plot(bin_psppratio_synm9(:,1),bin_psppratio_synm9(:,2),'o-');hold on;
plot(bin_psppratio_synm10(:,1),bin_psppratio_synm10(:,2),'o-');hold on;
errorbar(bin_psppratio(:,1),bin_psppratio(:,2),std_bsp_psppratio,'ko');legend('0.01','0.03','0.08','0.1','0.12','0.14','0.16','0.2','0.25','0.3');
xlabel('offset (m)');ylabel('PS/PP ratio');xlim([offset_min-100,offset_max-100]);


%% compute RMS and Chi-square misfit
% ---fluid---
dvp=[-2.081,-6.106,-15.687,-19.44,-23.204,-27.014,-30.907,-39.113,-50.796,-65.994];
dvs=[-1.273,-3.9,-10.998,-14.078,-17.32,-20.746,-24.385,-32.46,-44.746,-61.838];

N_model = length(dvp);
misfit_tab = zeros(N_model,5);

im = 1;
fprintf('[dvp,dvs]=[%.3f, %.3f]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_synf1(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_synf1(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 2;
fprintf('[dvp,dvs]=[%.3f, %.3f]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_synf2(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_synf2(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 3;
fprintf('[dvp,dvs]=[%.3f, %.3f]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_synf3(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_synf3(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 4;
fprintf('[dvp,dvs]=[%.3f, %.3f]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_synf4(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_synf4(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 5;
fprintf('[dvp,dvs]=[%.3f, %.3f]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_synf5(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_synf5(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 6;
fprintf('[dvp,dvs]=[%.3f, %.3f]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_synf6(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_synf6(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 7;
fprintf('[dvp,dvs]=[%.3f, %.3f]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_synf7(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_synf7(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 8;
fprintf('[dvp,dvs]=[%.3f, %.3f]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_synf8(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_synf8(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 9;
fprintf('[dvp,dvs]=[%.3f, %.3f]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_synf9(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_synf9(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 10;
fprintf('[dvp,dvs]=[%.3f, %.3f]\n',dvp(im),dvs(im));
misfit_tab(im,1) = dvp(im); misfit_tab(im,2) = dvs(im); misfit_tab(im,3) = dvp(im)/dvs(im);
misfit_tab(im,4) = RMSE(bin_psppratio_synf10(:,2), bin_psppratio(:,2));
misfit_tab(im,5) = Chi_square(bin_psppratio_synf10(:,2), bin_psppratio(:,2), std_bsp_psppratio);

% plot fluid fraction vs. misfit
fluid_frac = [0.01, 0.03, 0.08, 0.1, 0.12, 0.14, 0.16, 0.2, 0.25, 0.3];
figure;
plot(fluid_frac,misfit_tab(:,5),'bo-');
xlabel('Fluid fraction');ylabel('Chi-square misfit');

% ---melt---
dvp=[-1.627,-4.849,-12.771,-15.91,-19.043,-22.179,-25.327,-31.693,-39.897,-48.601];
dvs=[-1.703,-5.164,-14.157,-17.909,-21.767,-25.745,-29.864,-38.644,-51.16,-67.232];

N_model = length(dvp);
misfit_tab1 = zeros(N_model,5);

im = 1;
fprintf('[dvp,dvs]=[%.3f, %.3f]\n',dvp(im),dvs(im));
misfit_tab1(im,1) = dvp(im); misfit_tab1(im,2) = dvs(im); misfit_tab1(im,3) = dvp(im)/dvs(im);
misfit_tab1(im,4) = RMSE(bin_psppratio_synm1(:,2), bin_psppratio(:,2));
misfit_tab1(im,5) = Chi_square(bin_psppratio_synm1(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 2;
fprintf('[dvp,dvs]=[%.3f, %.3f]\n',dvp(im),dvs(im));
misfit_tab1(im,1) = dvp(im); misfit_tab1(im,2) = dvs(im); misfit_tab1(im,3) = dvp(im)/dvs(im);
misfit_tab1(im,4) = RMSE(bin_psppratio_synm2(:,2), bin_psppratio(:,2));
misfit_tab1(im,5) = Chi_square(bin_psppratio_synm2(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 3;
fprintf('[dvp,dvs]=[%.3f, %.3f]\n',dvp(im),dvs(im));
misfit_tab1(im,1) = dvp(im); misfit_tab1(im,2) = dvs(im); misfit_tab1(im,3) = dvp(im)/dvs(im);
misfit_tab1(im,4) = RMSE(bin_psppratio_synm3(:,2), bin_psppratio(:,2));
misfit_tab1(im,5) = Chi_square(bin_psppratio_synm3(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 4;
fprintf('[dvp,dvs]=[%.3f, %.3f]\n',dvp(im),dvs(im));
misfit_tab1(im,1) = dvp(im); misfit_tab1(im,2) = dvs(im); misfit_tab1(im,3) = dvp(im)/dvs(im);
misfit_tab1(im,4) = RMSE(bin_psppratio_synm4(:,2), bin_psppratio(:,2));
misfit_tab1(im,5) = Chi_square(bin_psppratio_synm4(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 5;
fprintf('[dvp,dvs]=[%.3f, %.3f]\n',dvp(im),dvs(im));
misfit_tab1(im,1) = dvp(im); misfit_tab1(im,2) = dvs(im); misfit_tab1(im,3) = dvp(im)/dvs(im);
misfit_tab1(im,4) = RMSE(bin_psppratio_synm5(:,2), bin_psppratio(:,2));
misfit_tab1(im,5) = Chi_square(bin_psppratio_synm5(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 6;
fprintf('[dvp,dvs]=[%.3f, %.3f]\n',dvp(im),dvs(im));
misfit_tab1(im,1) = dvp(im); misfit_tab1(im,2) = dvs(im); misfit_tab1(im,3) = dvp(im)/dvs(im);
misfit_tab1(im,4) = RMSE(bin_psppratio_synm6(:,2), bin_psppratio(:,2));
misfit_tab1(im,5) = Chi_square(bin_psppratio_synm6(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 7;
fprintf('[dvp,dvs]=[%.3f, %.3f]\n',dvp(im),dvs(im));
misfit_tab1(im,1) = dvp(im); misfit_tab1(im,2) = dvs(im); misfit_tab1(im,3) = dvp(im)/dvs(im);
misfit_tab1(im,4) = RMSE(bin_psppratio_synm7(:,2), bin_psppratio(:,2));
misfit_tab1(im,5) = Chi_square(bin_psppratio_synm7(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 8;
fprintf('[dvp,dvs]=[%.3f, %.3f]\n',dvp(im),dvs(im));
misfit_tab1(im,1) = dvp(im); misfit_tab1(im,2) = dvs(im); misfit_tab1(im,3) = dvp(im)/dvs(im);
misfit_tab1(im,4) = RMSE(bin_psppratio_synm8(:,2), bin_psppratio(:,2));
misfit_tab1(im,5) = Chi_square(bin_psppratio_synm8(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 9;
fprintf('[dvp,dvs]=[%.3f, %.3f]\n',dvp(im),dvs(im));
misfit_tab1(im,1) = dvp(im); misfit_tab1(im,2) = dvs(im); misfit_tab1(im,3) = dvp(im)/dvs(im);
misfit_tab1(im,4) = RMSE(bin_psppratio_synm9(:,2), bin_psppratio(:,2));
misfit_tab1(im,5) = Chi_square(bin_psppratio_synm9(:,2), bin_psppratio(:,2), std_bsp_psppratio);

im = 10;
fprintf('[dvp,dvs]=[%.3f, %.3f]\n',dvp(im),dvs(im));
misfit_tab1(im,1) = dvp(im); misfit_tab1(im,2) = dvs(im); misfit_tab1(im,3) = dvp(im)/dvs(im);
misfit_tab1(im,4) = RMSE(bin_psppratio_synm10(:,2), bin_psppratio(:,2));
misfit_tab1(im,5) = Chi_square(bin_psppratio_synm10(:,2), bin_psppratio(:,2), std_bsp_psppratio);

% plot melt fraction vs. misfit
melt_frac = [0.01, 0.03, 0.08, 0.1, 0.12, 0.14, 0.16, 0.2, 0.25, 0.3];
figure;
plot(melt_frac,misfit_tab1(:,5),'bo-');
xlabel('Melt fraction');ylabel('Chi-square misfit');

% plot together
figure;
plot(fluid_frac,misfit_tab(:,5),'bo-');hold on;
plot(melt_frac,misfit_tab1(:,5),'ro-');
xlabel('Porosity');ylabel('Chi-square misfit');


%% output for Python
% 1. AVO curve------
writematrix([bin_psppratio,std_bsp_psppratio],'/home0/cxd170430/NatureGeo/interpretation/plot_data/obs.txt','Delimiter',' ');

writematrix(bin_psppratio_synf1,'/home0/cxd170430/NatureGeo/interpretation/plot_data/f0.01.txt','Delimiter',' ');
writematrix(bin_psppratio_synf2,'/home0/cxd170430/NatureGeo/interpretation/plot_data/f0.03.txt','Delimiter',' ');
writematrix(bin_psppratio_synf3,'/home0/cxd170430/NatureGeo/interpretation/plot_data/f0.08.txt','Delimiter',' ');
writematrix(bin_psppratio_synf4,'/home0/cxd170430/NatureGeo/interpretation/plot_data/f0.1.txt','Delimiter',' ');
writematrix(bin_psppratio_synf5,'/home0/cxd170430/NatureGeo/interpretation/plot_data/f0.12.txt','Delimiter',' ');
writematrix(bin_psppratio_synf6,'/home0/cxd170430/NatureGeo/interpretation/plot_data/f0.14.txt','Delimiter',' ');
writematrix(bin_psppratio_synf7,'/home0/cxd170430/NatureGeo/interpretation/plot_data/f0.16.txt','Delimiter',' ');
writematrix(bin_psppratio_synf8,'/home0/cxd170430/NatureGeo/interpretation/plot_data/f0.2.txt','Delimiter',' ');
writematrix(bin_psppratio_synf9,'/home0/cxd170430/NatureGeo/interpretation/plot_data/f0.25.txt','Delimiter',' ');
writematrix(bin_psppratio_synf10,'/home0/cxd170430/NatureGeo/interpretation/plot_data/f0.3.txt','Delimiter',' ');

writematrix(bin_psppratio_synm1,'/home0/cxd170430/NatureGeo/interpretation/plot_data/m0.01.txt','Delimiter',' ');
writematrix(bin_psppratio_synm2,'/home0/cxd170430/NatureGeo/interpretation/plot_data/m0.03.txt','Delimiter',' ');
writematrix(bin_psppratio_synm3,'/home0/cxd170430/NatureGeo/interpretation/plot_data/m0.08.txt','Delimiter',' ');
writematrix(bin_psppratio_synm4,'/home0/cxd170430/NatureGeo/interpretation/plot_data/m0.1.txt','Delimiter',' ');
writematrix(bin_psppratio_synm5,'/home0/cxd170430/NatureGeo/interpretation/plot_data/m0.12.txt','Delimiter',' ');
writematrix(bin_psppratio_synm6,'/home0/cxd170430/NatureGeo/interpretation/plot_data/m0.14.txt','Delimiter',' ');
writematrix(bin_psppratio_synm7,'/home0/cxd170430/NatureGeo/interpretation/plot_data/m0.16.txt','Delimiter',' ');
writematrix(bin_psppratio_synm8,'/home0/cxd170430/NatureGeo/interpretation/plot_data/m0.2.txt','Delimiter',' ');
writematrix(bin_psppratio_synm9,'/home0/cxd170430/NatureGeo/interpretation/plot_data/m0.25.txt','Delimiter',' ');
writematrix(bin_psppratio_synm10,'/home0/cxd170430/NatureGeo/interpretation/plot_data/m0.3.txt','Delimiter',' ');

% 2. misfit curve------
writematrix([fluid_frac',misfit_tab],'/home0/cxd170430/NatureGeo/interpretation/plot_data/misfit_fluid.txt','Delimiter',' ');
writematrix([melt_frac',misfit_tab1],'/home0/cxd170430/NatureGeo/interpretation/plot_data/misfit_melt.txt','Delimiter',' ');








