index = snr < snr_threshold;
fprintf('\n SNR < %f removes %d (%.1f percent) traces\n\n',snr_threshold,sum(index),sum(index)/ntr*100);

data_env_snr = data_env;
geom_tab_snr = geom_tab;


itr_remove = zeros(sum(index),1);
    %offset_threshold = zeros(sum(index),1);
    count=1;
    for i=1:ntr_snr
        if snr(i) < snr_threshold
            %data_env_snr(:,itr_snr(i)) = NaN;
            %geom_tab_snr(itr_snr(i),:) = NaN;
            %offset_threshold(count) = geom_tab(itr_snr(i),5);
            itr_remove(count) = itr_snr(i);
            count=count+1;
        end
    end

    data_env_snr(:,itr_remove) = [];
    geom_tab_snr(itr_remove,:) = [];





%% 7.1 compute PS/PP ratio - trace by trace (synthetic)
% average offset range=2~6km by bins (bins: [2, 2.1] (2.1, 2.2] ... (5.9, 6] (6, 6.1])

pptmp_syn = zeros(ntr_new,1);
pstmp_syn = zeros(ntr_new,1);
ps_pp_ratio_syn = zeros(ntr_new,1);
for i=1:ntr_new
    ipp1 = floor((t_pp(i)-timewin_pp_prior)/dt);
    ipp2 = ceil((t_pp(i)+timewin_pp_post)/dt);
    pptmp_syn(i) = max(ez_10_10(ipp1:ipp2,i));

    ips1 = floor((t_ps(i)-timewin_ps_prior)/dt);
    ips2 = ceil((t_ps(i)+timewin_ps_post)/dt);
    pstmp_syn(i) = max(ex_10_10(ips1:ips2,i));

    ps_pp_ratio_syn(i) = pstmp_syn(i)/pptmp_syn(i);
end

offset_bin_length = 100; % meter
bins = offset_min:offset_bin_length:offset_max;
nbins = length(bins);
ratio_bin_syn = zeros(nbins-1,1);
for j=1:nbins-1
    if j==1
        index = find(geom_tab_new(:,5)>=bins(j) & geom_tab_new(:,5)<=bins(j+1));
    else
        index = find(geom_tab_new(:,5)>bins(j) & geom_tab_new(:,5)<=bins(j+1));
    end

    ratio_bin_syn(j) = median(ps_pp_ratio_syn(index));
end

figure;plot(bins(1:nbins-1),ratio_bin_syn,'o-');xlabel('offset (m)');

%
bin_psppratio_syn = amp_ratio(ex_10_10,ez_10_10,geom_tab_new,dt,t_ps,t_pp,timewin_ps_prior,timewin_ps_post,...
    timewin_pp_prior,timewin_pp_post,offset_min,offset_bin_length,offset_max,'max');

%% backup (trace by trace synthetic)
ez_ex_10_10 = ez_10_10+ex_10_10;
pptmp_syn = zeros(ntr_new,1);
pstmp_syn = zeros(ntr_new,1);
ps_pp_ratio_syn = zeros(ntr_new,1);
for i=1:ntr_new
    ipp1 = floor((t_pp(i)-timewin_pp_prior)/dt);
    ipp2 = ceil((t_pp(i)+timewin_pp_post)/dt);
    pptmp_syn(i) = max(ez_ex_10_10(ipp1:ipp2,i));

    ips1 = floor((t_ps(i)-timewin_ps_prior)/dt);
    ips2 = ceil((t_ps(i)+timewin_ps_post)/dt);
    pstmp_syn(i) = max(ez_ex_10_10(ips1:ips2,i));

    ps_pp_ratio_syn(i) = pstmp_syn(i)/pptmp_syn(i);
end

offset_bin_length = 100; % meter
bins = offset_min:offset_bin_length:offset_max;
nbins = length(bins);
ratio_bin_syn = zeros(nbins-1,1);
for j=1:nbins-1
    if j==1
        index = find(geom_tab_new(:,5)>=bins(j) & geom_tab_new(:,5)<=bins(j+1));
    else
        index = find(geom_tab_new(:,5)>bins(j) & geom_tab_new(:,5)<=bins(j+1));
    end

    ratio_bin_syn(j) = median(ps_pp_ratio_syn(index));
end

figure;plot(bins(1:nbins-1),ratio_bin_syn,'o-');xlabel('offset (m)');

%
bin_psppratio_syn = amp_ratio(ez_ex_10_10,ez_ex_10_10,geom_tab_new,dt,t_ps,t_pp,timewin_ps_prior,timewin_ps_post,...
    timewin_pp_prior,timewin_pp_post,offset_min,offset_bin_length,offset_max,'max');



%% 7.2 compute PS/PP ratio - bin by bin (synthetic)
% offset-bin traces stack and then compute PS/PP ratio
offset_bin_length = 100; % meter
bins = offset_min:offset_bin_length:offset_max;
nbins = length(bins);
ez_stack_syn = zeros(nt,nbins-1);
ex_stack_syn = zeros(nt,nbins-1);
tpp_env_stack = zeros(nbins-1,1);
tps_env_stack = zeros(nbins-1,1);
ratio_bin_syn = zeros(nbins-1,1);
for j=1:nbins-1
    if j==1
        index = find(geom_tab_new(:,5)>=bins(j) & geom_tab_new(:,5)<=bins(j+1));
    else
        index = find(geom_tab_new(:,5)>bins(j) & geom_tab_new(:,5)<=bins(j+1));
    end

    ez_stack_syn(:,j) = sum(ez_10_10(:,index),2);
    ex_stack_syn(:,j) = sum(ex_10_10(:,index),2);
    tpp_env_stack(j) = median(t_pp(index));
    tps_env_stack(j) = median(t_ps(index));
end

figure;subplot(1,2,1);
imagesc(1:nbins-1,0:dt:(nt-1)*dt,ez_stack_syn);clim([0,1e-12]);colormap(colorbar_bwr);
hold on;plot(1:nbins-1,tpp_env_stack-timewin_pp_prior,'k');
hold on;plot(1:nbins-1,tpp_env_stack+timewin_pp_post,'k');
hold on;plot(1:nbins-1,tps_env_stack-timewin_ps_prior,'k');
hold on;plot(1:nbins-1,tps_env_stack+timewin_ps_post,'k');
subplot(1,2,2);
imagesc(1:nbins-1,0:dt:(nt-1)*dt,ex_stack_syn);clim([0,1e-12]);colormap(colorbar_bwr);
hold on;plot(1:nbins-1,tpp_env_stack-timewin_pp_prior,'k');
hold on;plot(1:nbins-1,tpp_env_stack+timewin_pp_post,'k');
hold on;plot(1:nbins-1,tps_env_stack-timewin_ps_prior,'k');
hold on;plot(1:nbins-1,tps_env_stack+timewin_ps_post,'k');

for k=1:nbins-1
    ipp1 = floor((tpp_env_stack(k)-timewin_pp_prior)/dt);
    ipp2 = ceil((tpp_env_stack(k)+timewin_pp_post)/dt);
    ips1 = floor((tps_env_stack(k)-timewin_ps_prior)/dt);
    ips2 = ceil((tps_env_stack(k)+timewin_ps_post)/dt);
    ratio_bin_syn(k) = max(ex_stack_syn(ips1:ips2,k))/max(ez_stack_syn(ipp1:ipp2,k));
end

figure;plot(bins(1:nbins-1),ratio_bin_syn,'o-');xlabel('offset (m)');

%% backup (bin by bin synthetic)
ez_ex_stack_syn = ez_stack_syn+ex_stack_syn;
ratio_bin_syn = zeros(nbins-1,1);
for k=1:nbins-1
    ipp1 = floor((tpp_env_stack(k)-timewin_pp_prior)/dt);
    ipp2 = ceil((tpp_env_stack(k)+timewin_pp_post)/dt);
    ips1 = floor((tps_env_stack(k)-timewin_ps_prior)/dt);
    ips2 = ceil((tps_env_stack(k)+timewin_ps_post)/dt);
    ratio_bin_syn(k) = max(ez_ex_stack_syn(ips1:ips2,k))/max(ez_ex_stack_syn(ipp1:ipp2,k));
end

figure;plot(bins(1:nbins-1),ratio_bin_syn,'o-');xlabel('offset (m)');



%% PS/PP
pptmp = zeros(ntr_new,1);
pstmp = zeros(ntr_new,1);
ps_pp_ratio = zeros(ntr_new,1);
for i=1:ntr_new
    ipp1 = floor((t_pp(i)-timewin_pp_prior)/dt);
    ipp2 = ceil((t_pp(i)+timewin_pp_post)/dt);
    pptmp(i) = median(ZRT_env_new(ipp1:ipp2,i));
    %pptmp(i) = mean(Z_env_new(ipp1:ipp2,i));

    ips1 = floor((t_ps(i)-timewin_ps_prior)/dt);
    ips2 = ceil((t_ps(i)+timewin_ps_post)/dt);
    pstmp(i) = median(ZRT_env_new(ips1:ips2,i));
    %pstmp(i) = median(R_env_new(ips1:ips2,i));
    %pstmp(i) = median(T_env_new(ips1:ips2,i));
    %pstmp(i) = median(RT_env_new(ips1:ips2,i));

    ps_pp_ratio(i) = pstmp(i)/pptmp(i);
end

offset_bin_length = 100; % meter
bins = offset_min:offset_bin_length:offset_max;
nbins = length(bins);
ratio_bin = zeros(nbins-1,1);
for j=1:nbins-1
    if j==1
        index = find(geom_tab_new(:,5)>=bins(j) & geom_tab_new(:,5)<=bins(j+1));
    else
        index = find(geom_tab_new(:,5)>bins(j) & geom_tab_new(:,5)<=bins(j+1));
    end

    ratio_bin(j) = median(ps_pp_ratio(index));
end

figure;plot(bins(1:nbins-1),ratio_bin,'o-');xlabel('offset (m)');ylabel('PS/PP');


%% backup PP/P and PS/P
% OBSERVED
% ====PP/P====(bad)
bin_pppratio = amp_ratio(ZRT_env_new,ZRT_env_new,geom_tab_new,dt,t_pp,t_p,timewin_pp_prior,timewin_pp_post,...
    timewin_p_prior,timewin_p_post,offset_min,offset_bin_length,offset_max,'median',1500);

% ====PS/P====(bad)
bin_pspratio = amp_ratio(ZRT_env_new,ZRT_env_new,geom_tab_new,dt,t_ps,t_p,timewin_ps_prior,timewin_ps_post,...
    timewin_p_prior,timewin_p_post,offset_min,offset_bin_length,offset_max,'median',1500);

% 5. (Optional) spline curve fitting
% open curve fitting toolbox
curveFitter(bin_psppratio(:,1),bin_psppratio(:,2));

[fitresult, gof] = createFitpspp_tmp(bin_psppratio(:,1), bin_psppratio(:,2));
spline_curve = feval(fitresult,bin_psppratio(:,1));
figure;plot(bin_psppratio(:,1),spline_curve,'r');xlabel('offset (m)');ylabel('PS/PP');

[fitresult1, gof1] = createFitpspp_tmp1(bin_psppratio(:,1), bin_psppratio(:,2));
spline_curve1 = feval(fitresult1,bin_psppratio(:,1));
figure;plot(bin_psppratio(:,1),spline_curve1,'r');xlabel('offset (m)');ylabel('PS/PP');
%curveFitter(bin_pppratio(:,1),bin_pppratio(:,2));
%curveFitter(bin_pspratio(:,1),bin_pspratio(:,2));


% SYNTHETIC
% ====PP/P====(bad)
ez_ex_10_10 = ez_10_10 + ex_10_10;
bin_pppratio_syn1010 = amp_ratio(ez_ex_10_10,ez_ex_10_10,geom_tab_new,dt,t_pp,t_p,timewin_pp_prior,timewin_pp_post,...
    timewin_p_prior,timewin_p_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_30_10 = ez_30_10 + ex_30_10;
bin_pppratio_syn3010 = amp_ratio(ez_ex_30_10,ez_ex_30_10,geom_tab_new,dt,t_pp,t_p,timewin_pp_prior,timewin_pp_post,...
    timewin_p_prior,timewin_p_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_30_30 = ez_30_30 + ex_30_30;
bin_pppratio_syn3030 = amp_ratio(ez_ex_30_30,ez_ex_30_30,geom_tab_new,dt,t_pp,t_p,timewin_pp_prior,timewin_pp_post,...
    timewin_p_prior,timewin_p_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_30_50 = ez_30_50 + ex_30_50;
bin_pppratio_syn3050 = amp_ratio(ez_ex_30_50,ez_ex_30_50,geom_tab_new,dt,t_pp,t_p,timewin_pp_prior,timewin_pp_post,...
    timewin_p_prior,timewin_p_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_50_30 = ez_50_30 + ex_50_30;
bin_pppratio_syn5030 = amp_ratio(ez_ex_50_30,ez_ex_50_30,geom_tab_new,dt,t_pp,t_p,timewin_pp_prior,timewin_pp_post,...
    timewin_p_prior,timewin_p_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_50_50 = ez_50_50 + ex_50_50;
bin_pppratio_syn5050 = amp_ratio(ez_ex_50_50,ez_ex_50_50,geom_tab_new,dt,t_pp,t_p,timewin_pp_prior,timewin_pp_post,...
    timewin_p_prior,timewin_p_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

figure;plot(bin_pppratio_syn1010(:,1),bin_pppratio_syn1010(:,2),'o-');
hold on;plot(bin_pppratio_syn3010(:,1),bin_pppratio_syn3010(:,2),'o-');
hold on;plot(bin_pppratio_syn3030(:,1),bin_pppratio_syn3030(:,2),'o-');
hold on;plot(bin_pppratio_syn3050(:,1),bin_pppratio_syn3050(:,2),'o-');
hold on;plot(bin_pppratio_syn5030(:,1),bin_pppratio_syn5030(:,2),'o-');
hold on;plot(bin_pppratio_syn5050(:,1),bin_pppratio_syn5050(:,2),'o-');
legend('10-10','30-10','30-30','30-50','50-30','50-50');
xlim([1950,6150]);
grid on;

% ====PS/P====(bad)
ez_ex_10_10 = ez_10_10 + ex_10_10;
bin_pspratio_syn1010 = amp_ratio(ez_ex_10_10,ez_ex_10_10,geom_tab_new,dt,t_ps,t_p,timewin_ps_prior,timewin_ps_post,...
    timewin_p_prior,timewin_p_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_30_10 = ez_30_10 + ex_30_10;
bin_pspratio_syn3010 = amp_ratio(ez_ex_30_10,ez_ex_30_10,geom_tab_new,dt,t_ps,t_p,timewin_ps_prior,timewin_ps_post,...
    timewin_p_prior,timewin_p_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_30_30 = ez_30_30 + ex_30_30;
bin_pspratio_syn3030 = amp_ratio(ez_ex_30_30,ez_ex_30_30,geom_tab_new,dt,t_ps,t_p,timewin_ps_prior,timewin_ps_post,...
    timewin_p_prior,timewin_p_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_30_50 = ez_30_50 + ex_30_50;
bin_pspratio_syn3050 = amp_ratio(ez_ex_30_50,ez_ex_30_50,geom_tab_new,dt,t_ps,t_p,timewin_ps_prior,timewin_ps_post,...
    timewin_p_prior,timewin_p_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_50_30 = ez_50_30 + ex_50_30;
bin_pspratio_syn5030 = amp_ratio(ez_ex_50_30,ez_ex_50_30,geom_tab_new,dt,t_ps,t_p,timewin_ps_prior,timewin_ps_post,...
    timewin_p_prior,timewin_p_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

ez_ex_50_50 = ez_50_50 + ex_50_50;
bin_pspratio_syn5050 = amp_ratio(ez_ex_50_50,ez_ex_50_50,geom_tab_new,dt,t_ps,t_p,timewin_ps_prior,timewin_ps_post,...
    timewin_p_prior,timewin_p_post,offset_min,offset_bin_length,offset_max,'max',2e-14);

figure;plot(bin_pspratio_syn1010(:,1),bin_pspratio_syn1010(:,2),'o-');
hold on;plot(bin_pspratio_syn3010(:,1),bin_pspratio_syn3010(:,2),'o-');
hold on;plot(bin_pspratio_syn3030(:,1),bin_pspratio_syn3030(:,2),'o-');
hold on;plot(bin_pspratio_syn3050(:,1),bin_pspratio_syn3050(:,2),'o-');
hold on;plot(bin_pspratio_syn5030(:,1),bin_pspratio_syn5030(:,2),'o-');
hold on;plot(bin_pspratio_syn5050(:,1),bin_pspratio_syn5050(:,2),'o-');
legend('10-10','30-10','30-30','30-50','50-30','50-50');
xlim([1950,6150]);
grid on;





