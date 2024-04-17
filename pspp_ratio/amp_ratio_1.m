function bin_ratio_6col = amp_ratio(numer_datamat,denom_datamat,geom_tab_new,dt,t_numer,t_denom,twin_numer_prior,twin_numer_post,...
    twin_denom_prior,twin_denom_post,bin_min,bin_incre,bin_max,method,plot_maxamp)

% Real data: measure PS/PP ratios evolving with offsets
% Synthetic data: measure PS/PP ratios, PS amplitudes, PP amplitudes, dPS, dPP evolving with offsets

ntr_new = length(t_numer);
if ntr_new ~= size(numer_datamat,2)
    error('wrong trace number!');
end

% check time windows on data
nt = size(numer_datamat,1);
figure;subplot(1,2,1);
imagesc(1:ntr_new,0:dt:(nt-1)*dt,numer_datamat);clim([0,plot_maxamp]);colormap(colorbar_bwr);xlabel('trace# by offset');ylabel('time (s)');
hold on;plot(1:ntr_new,t_numer-twin_numer_prior,'y');
hold on;plot(1:ntr_new,t_numer+twin_numer_post,'y');
hold on;plot(1:ntr_new,t_denom-twin_denom_prior,'k');
hold on;plot(1:ntr_new,t_denom+twin_denom_post,'k');title('numerator data');
subplot(1,2,2);
imagesc(1:ntr_new,0:dt:(nt-1)*dt,denom_datamat);clim([0,plot_maxamp]);colormap(colorbar_bwr);xlabel('trace# by offset');ylabel('time (s)');
hold on;plot(1:ntr_new,t_numer-twin_numer_prior,'y');
hold on;plot(1:ntr_new,t_numer+twin_numer_post,'y');
hold on;plot(1:ntr_new,t_denom-twin_denom_prior,'k');
hold on;plot(1:ntr_new,t_denom+twin_denom_post,'k');title('denominator data');


% main body
amp_numer_tmp = zeros(ntr_new,1);
amp_denom_tmp = zeros(ntr_new,1);
ratio_tmp = zeros(ntr_new,1);

for i=1:ntr_new
    it1_numer = floor((t_numer(i)-twin_numer_prior)/dt);
    it2_numer = ceil((t_numer(i)+twin_numer_post)/dt);

    it1_denom = floor((t_denom(i)-twin_denom_prior)/dt);
    it2_denom = ceil((t_denom(i)+twin_denom_post)/dt);
    
    if strcmp(method,'max')
        amp_numer_tmp(i) = max(numer_datamat(it1_numer:it2_numer,i));
        amp_denom_tmp(i) = max(denom_datamat(it1_denom:it2_denom,i));
    end
    if strcmp(method,'median')
        amp_numer_tmp(i) = median(numer_datamat(it1_numer:it2_numer,i));
        amp_denom_tmp(i) = median(denom_datamat(it1_denom:it2_denom,i));
    end
    if strcmp(method,'mean')
        amp_numer_tmp(i) = mean(numer_datamat(it1_numer:it2_numer,i));
        amp_denom_tmp(i) = mean(denom_datamat(it1_denom:it2_denom,i));
    end

    ratio_tmp(i) = amp_numer_tmp(i)/amp_denom_tmp(i);
end


bins = bin_min:bin_incre:bin_max;
nbins = length(bins);
bin_ratio_6col = zeros(nbins-1,6);
for j=1:nbins-1
    if j==1
        index = find(geom_tab_new(:,5)>=bins(j) & geom_tab_new(:,5)<=bins(j+1));
    else
        index = find(geom_tab_new(:,5)>bins(j) & geom_tab_new(:,5)<=bins(j+1));
    end

    bin_ratio_6col(j,1) = 0.5*(bins(j)+bins(j+1)); % bins(j)
    bin_ratio_6col(j,2) = median(ratio_tmp(index)); % mean
    bin_ratio_6col(j,3) = median(amp_numer_tmp(index)); % PS
    bin_ratio_6col(j,4) = median(amp_denom_tmp(index)); % PP
    bin_ratio_6col(j,5) = bin_ratio_6col(j,3)/bin_ratio_6col(1,3) -1; % dPS
    bin_ratio_6col(j,6) = bin_ratio_6col(j,4)/bin_ratio_6col(1,4) -1; % dPP
end

figure;
subplot(1,3,1);plot(bin_ratio_6col(:,1),bin_ratio_6col(:,2),'o-');
xlabel('offset (m)');ylabel('ratio');
subplot(1,3,2);plot(bin_ratio_6col(:,1),bin_ratio_6col(:,5),'o-');
xlabel('offset (m)');ylabel('dPS');
subplot(1,3,3);plot(bin_ratio_6col(:,1),bin_ratio_6col(:,6),'o-');
xlabel('offset (m)');ylabel('dPP');


end