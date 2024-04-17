function data_mat_qc = SNR(su_file_env, data_mat, component, snr_threshold)
% Plot SNR histogram
% Data quality control based on SNR


[data_env,header,~] = ReadSu(su_file_env);

%----remove bad stations which record only one trace----
header([8354,8011,8006,8181,8178]) = [];
data_env(:,[8354,8011,8006,8181,8178]) = [];
%-------------------------------------------------------

ntr = length(header);
nt = header(1).ns;
%dt = header(1).dt * 10^(-6);

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
data_env_offset = data_env(:,order_offset);
geom_tab_offset = tmp_tab(order_offset,:);


% (x,y) points of the signal time window
if strcmp(component,'Z')
    muteD_ipt2=[1667,229]; muteD_ipt3=[11218,732];
    muteD_ipt4=[1667,105]; muteD_ipt5=[11218,608];
end
if strcmp(component,'R')
    muteD_ipt2=[1667,270]; muteD_ipt3=[11218,900];
    muteD_ipt4=[1667,146]; muteD_ipt5=[11218,776];
end
if strcmp(component,'T')
    muteD_ipt2=[1667,270]; muteD_ipt3=[11218,900];
    muteD_ipt4=[1667,146]; muteD_ipt5=[11218,776];
end

itr_snr = muteD_ipt2(1):muteD_ipt3(1);
ntr_snr = length(itr_snr);

it_mute_big = floor(muteD_ipt2(2)+(muteD_ipt3(2)-muteD_ipt2(2))*((muteD_ipt2(1):muteD_ipt3(1))-muteD_ipt2(1))/(muteD_ipt3(1)-muteD_ipt2(1)));
it_mute_small = floor(muteD_ipt4(2)+(muteD_ipt5(2)-muteD_ipt4(2))*((muteD_ipt4(1):muteD_ipt5(1))-muteD_ipt4(1))/(muteD_ipt5(1)-muteD_ipt4(1)));


figure;
imagesc(data_env_offset);colormap(colorbar_bwr);clim([0,1500]);
hold on;plot([muteD_ipt2(1),muteD_ipt3(1)],[muteD_ipt2(2),muteD_ipt3(2)],'k-');
hold on;plot([muteD_ipt4(1),muteD_ipt5(1)],[muteD_ipt4(2),muteD_ipt5(2)],'k-');

% figure;
% imagesc(data_offset);colormap(colorbar_bwr);clim([0,1500]);
% hold on;plot(itr_snr,it_mute_big,'k-','LineWidth',3);
% hold on;plot(itr_snr,it_mute_small,'k-','LineWidth',3);



snr = zeros(ntr_snr,1);
for i=1:ntr_snr
    signal = data_env_offset(it_mute_small(i):it_mute_big(i),itr_snr(i));
    noise = data_env_offset(setdiff(1:nt,it_mute_small(i):it_mute_big(i)),itr_snr(i));
    %snr(i) = max(signal)/median(noise);
    snr(i) = median(signal)/median(noise);
    %snr(i) = max(signal)/max(noise);
end


figure;
histogram(snr);
xlim([0,30]);
xlabel('Signal-to-noise ratio (SNR)');ylabel('Count');


%% quality control
index = snr < snr_threshold;
fprintf('\n SNR < %f removes %d (%.1f percent) traces\n\n',snr_threshold,sum(index),sum(index)/ntr*100);

data_mat_qc = data_mat;
offset_threshold = zeros(sum(index),1); count=1;
for i=1:ntr_snr
    if snr(i) < snr_threshold
        data_mat_qc(:,itr_snr(i)) = NaN;
        offset_threshold(count) = geom_tab_offset(itr_snr(i),5);count=count+1;
    end
end

figure;
histogram(offset_threshold);
xlabel('offset (m)');ylabel('Count');


end