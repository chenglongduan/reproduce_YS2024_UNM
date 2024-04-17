function bin_ratio_5col = amp_pp_ps(obs_datamat,geom_tab_new,dt,t_ps,t_pp,twin_ps_prior,twin_ps_post,twin_ps_noise,...
    twin_pp_prior,twin_pp_post,twin_pp_noise,bin_min,bin_incre,bin_max,plot_maxamp)

% Real data: measure PS amplitudes, PP amplitudes, dPS, dPP evolving with offsets

ntr_new = length(t_ps);
if ntr_new ~= size(obs_datamat,2)
    error('wrong trace number!');
end

% check time windows on data
nt = size(obs_datamat,1);
figure;
imagesc(1:ntr_new,0:dt:(nt-1)*dt,obs_datamat);clim([0,plot_maxamp]);colormap(colorbar_bwr);xlabel('trace# by offset');ylabel('time (s)');
hold on;plot(1:ntr_new,t_pp-twin_pp_prior-twin_pp_noise,'y--');
hold on;plot(1:ntr_new,t_pp-twin_pp_prior,'y');
hold on;plot(1:ntr_new,t_pp+twin_pp_post,'y');
hold on;plot(1:ntr_new,t_ps-twin_ps_prior-twin_ps_noise,'k--');
hold on;plot(1:ntr_new,t_ps-twin_ps_prior,'k');
hold on;plot(1:ntr_new,t_ps+twin_ps_post,'k');
title('real env data');


% main body
amp_pp = zeros(ntr_new,1);
amp_ps = zeros(ntr_new,1);

for i=1:ntr_new
    it0_pp = floor((t_pp(i)-twin_pp_prior-twin_pp_noise)/dt);
    it1_pp = floor((t_pp(i)-twin_pp_prior)/dt);
    it2_pp = ceil((t_pp(i)+twin_pp_post)/dt);

    it0_ps = floor((t_ps(i)-twin_ps_prior-twin_ps_noise)/dt);
    it1_ps = floor((t_ps(i)-twin_ps_prior)/dt);
    it2_ps = ceil((t_ps(i)+twin_ps_post)/dt);

    %
    %amp_pp(i) = median(obs_datamat(it1_pp:it2_pp,i)) - median(obs_datamat(it0_pp:it1_pp,i));
    %amp_ps(i) = median(obs_datamat(it1_ps:it2_ps,i)) - median(obs_datamat(it0_ps:it1_ps,i));

    %
    %amp_pp(i) = median(obs_datamat(it1_pp:it2_pp,i));
    %amp_ps(i) = median(obs_datamat(it1_ps:it2_ps,i));

    % 1
    amp_pp(i) = max(obs_datamat(it1_pp:it2_pp,i)) - min(obs_datamat(it0_pp:it1_pp,i));
    amp_ps(i) = max(obs_datamat(it1_ps:it2_ps,i)) - min(obs_datamat(it0_ps:it1_ps,i));

    
%     if strcmp(method,'max')
%         amp_numer_tmp(i) = max(numer_datamat(it1_numer:it2_numer,i));
%         amp_denom_tmp(i) = max(denom_datamat(it1_denom:it2_denom,i));
%     end
%     if strcmp(method,'median')
%         amp_numer_tmp(i) = median(numer_datamat(it1_numer:it2_numer,i));
%         amp_denom_tmp(i) = median(denom_datamat(it1_denom:it2_denom,i));
%     end
%     if strcmp(method,'mean')
%         amp_numer_tmp(i) = mean(numer_datamat(it1_numer:it2_numer,i));
%         amp_denom_tmp(i) = mean(denom_datamat(it1_denom:it2_denom,i));
%     end

end


bins = bin_min:bin_incre:bin_max;
nbins = length(bins);
bin_ratio_5col = zeros(nbins-1,5);
for j=1:nbins-1
    if j==1
        index = find(geom_tab_new(:,5)>=bins(j) & geom_tab_new(:,5)<=bins(j+1));
    else
        index = find(geom_tab_new(:,5)>bins(j) & geom_tab_new(:,5)<=bins(j+1));
    end

    bin_ratio_5col(j,1) = 0.5*(bins(j)+bins(j+1)); % bins(j)
    bin_ratio_5col(j,2) = median(amp_ps(index)); % PS
    bin_ratio_5col(j,3) = median(amp_pp(index)); % PP
    bin_ratio_5col(j,4) = bin_ratio_5col(j,2)/bin_ratio_5col(1,2) -1; % dPS
    bin_ratio_5col(j,5) = bin_ratio_5col(j,3)/bin_ratio_5col(1,3) -1; % dPP
end

figure;
subplot(1,2,1);plot(bin_ratio_5col(:,1),bin_ratio_5col(:,4),'o-');
xlabel('offset (m)');ylabel('dPS');
subplot(1,2,2);plot(bin_ratio_5col(:,1),bin_ratio_5col(:,5),'o-');
xlabel('offset (m)');ylabel('dPP');


end