

gamma=vp/vs;

sqrt(gamma^2 + (gamma^2-1)*mx*mx/(z*z)) * x_offset / (1+sqrt(gamma^2 + (gamma^2-1)*mx*mx/(z*z)))

mx/sqrt(gamma^2 + (gamma^2-1)*mx*mx/(z*z))

xp+xs=x_offset;



t_ps_seg = zeros(100,1);
sx=0;
rx=100:100:10000;
gamma=vp/vs;
% Tps
for i=1:100
    syms mx;
    x_offset = rx(i);
    x_ccp = vpasolve(sqrt(gamma^2 + (gamma^2-1)*mx*mx/(z*z)) * x_offset / (1+sqrt(gamma^2 + (gamma^2-1)*mx*mx/(z*z)))+...
                     mx/sqrt(gamma^2 + (gamma^2-1)*mx*mx/(z*z)) == x_offset, mx);
    t_ps_seg(i) = sqrt(x_ccp*x_ccp + z*z) / vp + sqrt((x_offset-x_ccp)*(x_offset-x_ccp) + z*z) / vs;
end

figure;plot(t_ps_duan,'b');hold on;plot(t_ps_seg,'r--');
legend('Chenglong','SEG');



t_ps_new = zeros(100,1);
sx=0;
rx=100:100:10000;
% Tps
gamma = vp/vs;
for i=1:100
    syms mx;
    x_offset = rx(i);
    x_ccp = vpasolve(sqrt(gamma^2 + (gamma^2-1)*mx*mx/(z*z)) * x_offset / (1+sqrt(gamma^2 + (gamma^2-1)*mx*mx/(z*z))) +...
                     mx/sqrt(gamma^2 + (gamma^2-1)*mx*mx/(z*z)) == x_offset, mx);
    t_ps_new(i) = sqrt(x_ccp*x_ccp + z*z)/vp + sqrt((x_offset-x_ccp)*(x_offset-x_ccp) + z*z)/vs;
end

figure;plot(t_ps_duan,'b');hold on;plot(t_ps_new,'r--');legend('raw_duan','new_duan');


vp=4500;vs=2500;z=4000;
incident_angle = 0:0.01:90;
N = length(incident_angle);
pp_angle = zeros(N,1);
ps_angle = zeros(N,1);
sp_angle = zeros(N,1);
ss_angle = zeros(N,1);
offset_pp = zeros(N,1);offset_ps = zeros(N,1);offset_ss = zeros(N,1);offset_sp = zeros(N,1);
tt_pp = zeros(N,1);tt_ps = zeros(N,1);tt_ss = zeros(N,1);tt_sp = zeros(N,1);
for i=1:N
    pp_angle(i) = incident_angle(i);
    offset_pp(i) = z*(tand(incident_angle(i))+tand(pp_angle(i)));
    tt_pp(i) = z/cosd(incident_angle(i))/vp+z/cosd(pp_angle(i))/vp;
    
    ss_angle(i) = incident_angle(i);
    offset_ss(i) = z*(tand(incident_angle(i))+tand(ss_angle(i)));
    tt_ss(i) = z/cosd(incident_angle(i))/vs+z/cosd(ss_angle(i))/vs;
    
    crit1 = (vs/vp)*sind(incident_angle(i));
    if crit1<=1
        ps_angle(i) = asind(crit1);
        offset_ps(i) = z*(tand(incident_angle(i))+tand(ps_angle(i)));
        tt_ps(i) = z/cosd(incident_angle(i))/vp+z/cosd(ps_angle(i))/vs;
    else
        ps_angle(i) = NaN;
        offset_ps(i) = NaN;
        tt_ps(i) = NaN;
    end
    
    crit2 = (vp/vs)*sind(incident_angle(i));
    if crit2<=1
        sp_angle(i) = asind(crit2);
        offset_sp(i) = z*(tand(incident_angle(i))+tand(sp_angle(i)));
        tt_sp(i) = z/cosd(incident_angle(i))/vs+z/cosd(sp_angle(i))/vp;
    else
        sp_angle(i) = NaN;
        offset_sp(i) = NaN;
        tt_sp(i) = NaN;
    end
end
figure;plot(incident_angle,pp_angle,'b','linewidth',1.5);
hold on;plot(incident_angle,ss_angle,'r--','linewidth',1.5);
hold on;plot(incident_angle,ps_angle,'k','linewidth',1.5);
hold on;plot(incident_angle,sp_angle,'m','linewidth',1.5);
legend('PP','SS','PS','SP');
xlabel('incident angle (deg)');ylabel('reflection angle (deg)');
figure;plot(offset_pp(1:7508)/1000,tt_pp(1:7508),'linewidth',1.5);
hold on;plot(offset_ss(1:7508)/1000,tt_ss(1:7508),'linewidth',1.5);
hold on;plot(offset_ps(1:8170)/1000,tt_ps(1:8170),'linewidth',1.5);
hold on;plot(offset_sp(1:3336)/1000,tt_sp(1:3336),'--','linewidth',1.5);
legend('PP','SS','PS','SP');xlim([0,31]);ylim([0,13]);
xlabel('offset (km)');ylabel('traveltime (s)');


%% old 'velocity_drop.m' codes

% measure envelope amplitudes ratio due to velocity drops
% normalized Epp vs. offset
% normalized Eps vs. offset

% use Z to measure PP, use R to measure PS

% load Z
fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/Z_envelope.bin','rb');
z_env = fread(fp,[2001 11218],'float64');
fclose(fp);

fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/Z_wiggle.bin','rb');
z_wiggle = fread(fp,[2001 11218],'float64');
fclose(fp);

fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/Z_stack_stalta.bin','rb');
z_stack = fread(fp,[2001 100],'float64');
fclose(fp);
% load R
fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/R_envelope.bin','rb');
r_env = fread(fp,[2001 11218],'float64');
fclose(fp);

fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/R_wiggle.bin','rb');
r_wiggle = fread(fp,[2001 11218],'float64');
fclose(fp);

fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/R_stack_stalta.bin','rb');
r_stack = fread(fp,[2001 100],'float64');
fclose(fp);
% load T
fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/T_envelope.bin','rb');
t_env = fread(fp,[2001 11218],'float64');
fclose(fp);

fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/T_wiggle.bin','rb');
t_wiggle = fread(fp,[2001 11218],'float64');
fclose(fp);

fp =fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/T_stack_stalta.bin','rb');
t_stack = fread(fp,[2001 100],'float64');
fclose(fp);

figure;imagesc(z_env);clim([0,1000]);colorbar;
figure;imagesc(z_wiggle);clim([-500,500]);colorbar;
figure;imagesc(z_stack);colormap(colorbar_bwr);clim([0.5,1.5]);
figure;imagesc(r_env);clim([0,1000]);colorbar;
figure;imagesc(r_wiggle);clim([-500,500]);colorbar;
figure;imagesc(r_stack);colormap(colorbar_bwr);clim([0.5,1.5]);
figure;imagesc(t_env);clim([0,1000]);colorbar;
figure;imagesc(t_wiggle);clim([-500,500]);colorbar;
figure;imagesc(t_stack);colormap(colorbar_bwr);clim([0.5,1.5]);

% (obs) extract windowed offset traces

offset_min = 5000; % meter
offset_max = 7000; % meter

geom_tab_new = geom_tab(geom_tab(:,5)>=offset_min & geom_tab(:,5)<=offset_max, 5:7);
trace_index = find(geom_tab(:,5)>=offset_min & geom_tab(:,5)<=offset_max);

z_env_offsetwin = z_env(:,trace_index);
r_env_offsetwin = r_env(:,trace_index);
t_env_offsetwin = t_env(:,trace_index);

figure;imagesc(z_env_offsetwin);clim([0,1000]);colorbar;

% OPTIONAL...
% fp = fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/obs_wiggle_5_7km','wb');
% fwrite(,'float32');
% fclose(fp);
% 
% fp = fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/obs_env_5_7km','wb');
% fwrite(,'float32');
% fclose(fp);

% load velocity traces
% 10/10
fp = fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vx_10_10_64bit.bin','rb');
vx_10_10 = fread(fp,[2001 11218],'float64');
fclose(fp);
fp = fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vz_10_10_64bit.bin','rb');
vz_10_10 = fread(fp,[2001 11218],'float64');
fclose(fp);
% 10/30
fp = fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vx_10_30_64bit.bin','rb');
vx_10_30 = fread(fp,[2001 11218],'float64');
fclose(fp);
fp = fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vz_10_30_64bit.bin','rb');
vz_10_30 = fread(fp,[2001 11218],'float64');
fclose(fp);
% 10/50
fp = fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vx_10_50_64bit.bin','rb');
vx_10_50 = fread(fp,[2001 11218],'float64');
fclose(fp);
fp = fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vz_10_50_64bit.bin','rb');
vz_10_50 = fread(fp,[2001 11218],'float64');
fclose(fp);
% 30/10
fp = fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vx_30_10_64bit.bin','rb');
vx_30_10 = fread(fp,[2001 11218],'float64');
fclose(fp);
fp = fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vz_30_10_64bit.bin','rb');
vz_30_10 = fread(fp,[2001 11218],'float64');
fclose(fp);
% 30/30
fp = fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vx_30_30_64bit.bin','rb');
vx_30_30 = fread(fp,[2001 11218],'float64');
fclose(fp);
fp = fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vz_30_30_64bit.bin','rb');
vz_30_30 = fread(fp,[2001 11218],'float64');
fclose(fp);
% 30/50
fp = fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vx_30_50_64bit.bin','rb');
vx_30_50 = fread(fp,[2001 11218],'float64');
fclose(fp);
fp = fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vz_30_50_64bit.bin','rb');
vz_30_50 = fread(fp,[2001 11218],'float64');
fclose(fp);
% 50/30
fp = fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vx_50_30_64bit.bin','rb');
vx_50_30 = fread(fp,[2001 11218],'float64');
fclose(fp);
fp = fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vz_50_30_64bit.bin','rb');
vz_50_30 = fread(fp,[2001 11218],'float64');
fclose(fp);
% 50/50
fp = fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vx_50_50_64bit.bin','rb');
vx_50_50 = fread(fp,[2001 11218],'float64');
fclose(fp);
fp = fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vz_50_50_64bit.bin','rb');
vz_50_50 = fread(fp,[2001 11218],'float64');
fclose(fp);

itr1=2000;
figure;
subplot(2,1,1);plot(vx_10_10(:,itr1));hold on;plot(vx_10_30(:,itr1));hold on;plot(vx_10_50(:,itr1));hold on;plot(vx_30_10(:,itr1));hold on;
plot(vx_30_30(:,itr1));hold on;plot(vx_30_50(:,itr1));hold on;plot(vx_50_30(:,itr1));hold on;plot(vx_50_50(:,itr1));xlim([1,2001]);title(strcat('R',num2str(itr1)));
subplot(2,1,2);plot(vz_10_10(:,itr1));hold on;plot(vz_10_30(:,itr1));hold on;plot(vz_10_50(:,itr1));hold on;plot(vz_30_10(:,itr1));hold on;
plot(vz_30_30(:,itr1));hold on;plot(vz_30_50(:,itr1));hold on;plot(vz_50_30(:,itr1));hold on;plot(vz_50_50(:,itr1));xlim([1,2001]);title(strcat('Z',num2str(itr1)));

itr2=4000;
figure;
subplot(2,1,1);plot(vx_10_10(:,itr2));hold on;plot(vx_10_30(:,itr2));hold on;plot(vx_10_50(:,itr2));hold on;plot(vx_30_10(:,itr2));hold on;
plot(vx_30_30(:,itr2));hold on;plot(vx_30_50(:,itr2));hold on;plot(vx_50_30(:,itr2));hold on;plot(vx_50_50(:,itr2));xlim([1,2001]);title(strcat('R',num2str(itr2)));
subplot(2,1,2);plot(vz_10_10(:,itr2));hold on;plot(vz_10_30(:,itr2));hold on;plot(vz_10_50(:,itr2));hold on;plot(vz_30_10(:,itr2));hold on;
plot(vz_30_30(:,itr2));hold on;plot(vz_30_50(:,itr2));hold on;plot(vz_50_30(:,itr2));hold on;plot(vz_50_50(:,itr2));xlim([1,2001]);title(strcat('Z',num2str(itr2)));

itr3=6000;
figure;
subplot(2,1,1);plot(vx_10_10(:,itr3));hold on;plot(vx_10_30(:,itr3));hold on;plot(vx_10_50(:,itr3));hold on;plot(vx_30_10(:,itr3));hold on;
plot(vx_30_30(:,itr3));hold on;plot(vx_30_50(:,itr3));hold on;plot(vx_50_30(:,itr3));hold on;plot(vx_50_50(:,itr3));xlim([1,2001]);title(strcat('R',num2str(itr3)));
subplot(2,1,2);plot(vz_10_10(:,itr3));hold on;plot(vz_10_30(:,itr3));hold on;plot(vz_10_50(:,itr3));hold on;plot(vz_30_10(:,itr3));hold on;
plot(vz_30_30(:,itr3));hold on;plot(vz_30_50(:,itr3));hold on;plot(vz_50_30(:,itr3));hold on;plot(vz_50_50(:,itr3));xlim([1,2001]);title(strcat('Z',num2str(itr3)));

itr4=8000;
figure;
subplot(2,1,1);plot(vx_10_10(:,itr4));hold on;plot(vx_10_30(:,itr4));hold on;plot(vx_10_50(:,itr4));hold on;plot(vx_30_10(:,itr4));hold on;
plot(vx_30_30(:,itr4));hold on;plot(vx_30_50(:,itr4));hold on;plot(vx_50_30(:,itr4));hold on;plot(vx_50_50(:,itr4));xlim([1,2001]);title(strcat('R',num2str(itr4)));
subplot(2,1,2);plot(vz_10_10(:,itr4));hold on;plot(vz_10_30(:,itr4));hold on;plot(vz_10_50(:,itr4));hold on;plot(vz_30_10(:,itr4));hold on;
plot(vz_30_30(:,itr4));hold on;plot(vz_30_50(:,itr4));hold on;plot(vz_50_30(:,itr4));hold on;plot(vz_50_50(:,itr4));xlim([1,2001]);title(strcat('Z',num2str(itr4)));

itr5=10000;
figure;
subplot(2,1,1);plot(vx_10_10(:,itr5));hold on;plot(vx_10_30(:,itr5));hold on;plot(vx_10_50(:,itr5));hold on;plot(vx_30_10(:,itr5));hold on;
plot(vx_30_30(:,itr5));hold on;plot(vx_30_50(:,itr5));hold on;plot(vx_50_30(:,itr5));hold on;plot(vx_50_50(:,itr5));xlim([1,2001]);title(strcat('R',num2str(itr5)));
subplot(2,1,2);plot(vz_10_10(:,itr5));hold on;plot(vz_10_30(:,itr5));hold on;plot(vz_10_50(:,itr5));hold on;plot(vz_30_10(:,itr5));hold on;
plot(vz_30_30(:,itr5));hold on;plot(vz_30_50(:,itr5));hold on;plot(vz_50_30(:,itr5));hold on;plot(vz_50_50(:,itr5));xlim([1,2001]);title(strcat('Z',num2str(itr5)));


% load envelope traces
% 10/10
fp = fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_10_10_64bit.bin','rb');
ex_10_10 = fread(fp,[2001 11218],'float64');
fclose(fp);
fp = fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_10_10_64bit.bin','rb');
ez_10_10 = fread(fp,[2001 11218],'float64');
fclose(fp);

% 50/50
fp = fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ex_50_50_64bit.bin','rb');
ex_50_50 = fread(fp,[2001 11218],'float64');
fclose(fp);
fp = fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/envelope/ez_50_50_64bit.bin','rb');
ez_50_50 = fread(fp,[2001 11218],'float64');
fclose(fp);

%fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/test64bit_vx.bin','rb');test2_vx=fread(fp,[2001 11218],'float64');
%fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/test64bit_vz.bin','rb');test2_vz=fread(fp,[2001 11218],'float64');

itr=7228;
figure;subplot(2,1,1);plot(data_env_offset(:,itr)/max(data_env_offset(:,itr)));hold on;plot(Z_stack(:,60)/max(Z_stack(:,60)));
subplot(2,1,2);plot(ez_50_50(:,itr),'r');%plot(vz_50_50(:,itr));hold on;
itr=6000;
figure;subplot(2,1,1);plot(data_env_offset(:,itr));hold on;plot(Z_stack(:,48));
subplot(2,1,2);plot(ez_50_50(:,itr),'r');
itr=3091;
figure;subplot(2,1,1);plot(data_env_offset(:,itr));
subplot(2,1,2);plot(ez_50_50(:,itr),'r');
itr=6087;
figure;subplot(2,1,1);plot(data_env_offset(:,itr));
subplot(2,1,2);plot(ez_50_50(:,itr),'r');

%-----


%
figure;imagesc(data_env_offset);clim([0,1000])
figure;imagesc(10:100:10000,1:2001,Z_stack);colormap(colorbar_bwr);clim([0.5,1.4]);
figure;imagesc(ex_50_50);clim([0,1e-15])


%% for combining
file_vel_drop = '50_50';
dt_syn = 5e-4;
dt_obs = 0.004;
T = 8.0;
% [shot#, trace#]
syn_dict = zeros(size(geom_tab,1),2);
syn_vx = zeros(2001,size(geom_tab,1));
syn_vz = zeros(2001,size(geom_tab,1));
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
% OPTIONAL...
output_vx = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vx_',file_vel_drop,'_32bit.bin');
output_vz = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vz_',file_vel_drop,'_32bit.bin');
fp=fopen(output_vx,'wb'); fwrite(fp,syn_vx,'float32'); fclose(fp);
fp=fopen(output_vz,'wb'); fwrite(fp,syn_vz,'float32'); fclose(fp);

output_vx = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vx_',file_vel_drop,'_64bit.bin');
output_vz = strcat('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/velocity_drop/gardner/vz_',file_vel_drop,'_64bit.bin');
fp=fopen(output_vx,'wb'); fwrite(fp,syn_vx,'float64'); fclose(fp);
fp=fopen(output_vz,'wb'); fwrite(fp,syn_vz,'float64'); fclose(fp);




