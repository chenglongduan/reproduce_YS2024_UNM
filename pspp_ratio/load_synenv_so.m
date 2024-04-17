function [Zenv_final,Renv_final,Tenv_final,RTenv_final,geom_tab_final] = ...
    load_synenv_so(su_file1, su_file2, su_file3, su_file4, component, yn_snr_qc, snr_threshold, offset_min, offset_max)
% load envelope (env) data
% select traces by SNR (s) and offset (o)

% major
[Zenv_raw,header,~] = ReadSu(su_file1);
header([8354,8011,8006,8181,8178]) = [];
Zenv_raw(:,[8354,8011,8006,8181,8178]) = [];
ntr = length(header);
nt = header(1).ns;
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
[~,order_offset] = sort(tmp_tab(:,5));
geom_tab = tmp_tab(order_offset,:);
Zenv = Zenv_raw(:,order_offset);
%------2-D line projection
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
%------
% auxil
[Renv_raw,~,~] = ReadSu(su_file2);
Renv_raw(:,[8354,8011,8006,8181,8178]) = [];
Renv = Renv_raw(:,order_offset);
[Tenv,~,~] = ReadSu(su_file3);
[RTenv,~,~] = ReadSu(su_file4);


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
itr_snr = itr_snr';
ntr_snr = length(itr_snr);

it_mute_big = floor(muteD_ipt2(2)+(muteD_ipt3(2)-muteD_ipt2(2))*((muteD_ipt2(1):muteD_ipt3(1))-muteD_ipt2(1))/(muteD_ipt3(1)-muteD_ipt2(1)));
it_mute_small = floor(muteD_ipt4(2)+(muteD_ipt5(2)-muteD_ipt4(2))*((muteD_ipt4(1):muteD_ipt5(1))-muteD_ipt4(1))/(muteD_ipt5(1)-muteD_ipt4(1)));

% figure;
% imagesc(Zenv);colormap(colorbar_bwr);clim([0,1500]);
% hold on;plot([muteD_ipt2(1),muteD_ipt3(1)],[muteD_ipt2(2),muteD_ipt3(2)],'k-');
% hold on;plot([muteD_ipt4(1),muteD_ipt5(1)],[muteD_ipt4(2),muteD_ipt5(2)],'k-');


%% compute SNR of each trace
snr = zeros(ntr_snr,1);
for i=1:ntr_snr
    signal = Zenv(it_mute_small(i):it_mute_big(i),itr_snr(i));
    noise = Zenv(setdiff(1:nt,it_mute_small(i):it_mute_big(i)),itr_snr(i));
    snr(i) = median(signal)/median(noise);
end

% figure;
% histogram(snr);
% xlim([0,30]);
% xlabel('Signal-to-noise ratio (SNR)');ylabel('Count');


%% select SNR
% major
Zenv_snr = Zenv;
geom_tab_snr = geom_tab;
% auxil
Renv_snr = Renv;
Tenv_snr = Tenv;
RTenv_snr = RTenv;

if yn_snr_qc
    index = snr < snr_threshold;
    fprintf('\n SNR < %f removes %d (%.1f percent) traces\n\n',snr_threshold,sum(index),sum(index)/ntr*100);

    itr_remove = itr_snr(index);
    
    % major
    Zenv_snr(:,itr_remove) = [];
    geom_tab_snr(itr_remove,:) = [];
    % auxil
    Renv_snr(:,itr_remove) = [];
    Tenv_snr(:,itr_remove) = [];
    RTenv_snr(:,itr_remove) = [];
end


%% select offset
index_o = geom_tab_snr(:,5) >= offset_min & geom_tab_snr(:,5) <= offset_max;
% major
Zenv_final = Zenv_snr(:,index_o);
geom_tab_final = geom_tab_snr(index_o,:);
% auxil
Renv_final = Renv_snr(:,index_o);
Tenv_final = Tenv_snr(:,index_o);
RTenv_final = RTenv_snr(:,index_o);



end