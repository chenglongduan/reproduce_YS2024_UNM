%% Z-component

% 1. save STA/LTA (all traces)
fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/Z_stalta.bin','wb');
fwrite(fp,data_Z,'float64');
fclose(fp);


% 2. save stacked STA/LTA
fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/Z_stack_stalta.bin','wb');
fwrite(fp,Z_stack,'float64');
fclose(fp);


% 3. save envelope (all traces)
su_file_env='/data/Chenglong/Yellowstone_vibroseis/su_data_envelope/combine.Z.su';
[data_env,header,~] = ReadSu(su_file_env);
header([8354,8011,8006,8181,8178]) = [];
data_env(:,[8354,8011,8006,8181,8178]) = [];
ntr = length(header);
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
data_env_offset = data_env(:,order_offset);
fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/Z_envelope.bin','wb');
fwrite(fp,data_env_offset,'float64');
fclose(fp);

% bin stack
bins = 100:100:10000; % meters
data_env_stack = zeros(2001,length(bins));
for j=1:length(bins)
    if j==1
        index = find(tmp_tab_offset<=bins(1));
    else
        index = find(tmp_tab_offset>bins(j-1) & tmp_tab_offset<=bins(j));
    end
    data_env_stack(:,j) = sum(data_env_offset(:,index),2);
end
figure;imagesc(data_env_stack);clim([1e3,1e5]);
fp=fopen('/home0/cxd170430/NatureGeo/seismic_data/data/B/Z_env_stack.bin','wb');
fwrite(fp,data_env_stack,'float32');
fclose(fp);


% 4. save velocity wiggles (all traces)
su_file_wiggle='/data/Chenglong/Yellowstone_vibroseis/su_data_vel_traces/combine.Z.su';
[data_wiggle,header,~] = ReadSu(su_file_wiggle);
header([8354,8011,8006,8181,8178]) = [];
data_wiggle(:,[8354,8011,8006,8181,8178]) = [];
ntr = length(header);
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
[tmp_tab_offset,order_offset] = sort(tmp_tab(:,5));
data_wiggle_offset = data_wiggle(:,order_offset);
figure;imagesc(data_wiggle_offset);clim([-1e2,1e2]);

fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/Z_wiggle.bin','wb');
fwrite(fp,data_wiggle_offset,'float64');
fclose(fp);
% tr=4128;
% figure;plot(data_wiggle_offset(:,tr)/max(abs(data_wiggle_offset(:,tr))));
% hold on;plot(data_env_offset(:,tr)/max(data_env_offset(:,tr)));
% hold on;plot(data_Z(:,tr)/max(data_Z(:,tr)));

% bin stack
bins = 100:50:10000; % meters
data_wiggle_stack = zeros(2001,length(bins));
for j=1:length(bins)
    if j==1
        index = find(tmp_tab_offset<=bins(1));
    else
        index = find(tmp_tab_offset>bins(j-1) & tmp_tab_offset<=bins(j));
    end
    data_wiggle_stack(:,j) = sum(data_wiggle_offset(:,index),2);
end
figure;imagesc(data_wiggle_stack);clim([-0.6e4,0.6e4]);


%% R-component

% 1. save STA/LTA (all traces)
fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/R_stalta.bin','wb');
fwrite(fp,data_R,'float64');
fclose(fp);

% 2. save stacked STA/LTA
fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/R_stack_stalta.bin','wb');
fwrite(fp,R_stack,'float64');
fclose(fp);

% 3. save envelope (all traces)
su_file_env='/data/Chenglong/Yellowstone_vibroseis/su_data_envelope/combine.R.su';
[data_env,header,~] = ReadSu(su_file_env);
header([8354,8011,8006,8181,8178]) = [];
data_env(:,[8354,8011,8006,8181,8178]) = [];
ntr = length(header);
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
data_env_offset = data_env(:,order_offset);
fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/R_envelope.bin','wb');
fwrite(fp,data_env_offset,'float64');
fclose(fp);


% 4. save velocity wiggles (all traces)
su_file_wiggle='/data/Chenglong/Yellowstone_vibroseis/su_data_vel_traces/combine.R.su';
[data_wiggle,header,~] = ReadSu(su_file_wiggle);
header([8354,8011,8006,8181,8178]) = [];
data_wiggle(:,[8354,8011,8006,8181,8178]) = [];
ntr = length(header);
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
[tmp_tab_offset,order_offset] = sort(tmp_tab(:,5));
data_wiggle_offset = data_wiggle(:,order_offset);
fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/R_wiggle.bin','wb');
fwrite(fp,data_wiggle_offset,'float64');
fclose(fp);
% tr=4128;
% figure;plot(data_wiggle_offset(:,tr)/max(abs(data_wiggle_offset(:,tr))));
% hold on;plot(data_env_offset(:,tr)/max(data_env_offset(:,tr)));
% hold on;plot(data_R(:,tr)/max(data_R(:,tr)));

% bin stack
bins = 100:50:10000; % meters
data_wiggle_stack = zeros(2001,length(bins));
for j=1:length(bins)
    if j==1
        index = find(tmp_tab_offset<=bins(1));
    else
        index = find(tmp_tab_offset>bins(j-1) & tmp_tab_offset<=bins(j));
    end
    data_wiggle_stack(:,j) = sum(data_wiggle_offset(:,index),2);
end
figure;imagesc(data_wiggle_stack);clim([-0.6e4,0.6e4]);


%% T-component

% 1. save STA/LTA (all traces)
fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/T_stalta.bin','wb');
fwrite(fp,data_T,'float64');
fclose(fp);

% 2. save stacked STA/LTA
fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/T_stack_stalta.bin','wb');
fwrite(fp,T_stack,'float64');
fclose(fp);

% 3. save envelope (all traces)
su_file_env='/data/Chenglong/Yellowstone_vibroseis/su_data_envelope/combine.T.su';
[data_env,header,~] = ReadSu(su_file_env);
header([8354,8011,8006,8181,8178]) = [];
data_env(:,[8354,8011,8006,8181,8178]) = [];
ntr = length(header);
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
data_env_offset = data_env(:,order_offset);
fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/T_envelope.bin','wb');
fwrite(fp,data_env_offset,'float64');
fclose(fp);


% 4. save velocity wiggles (all traces)
su_file_wiggle='/data/Chenglong/Yellowstone_vibroseis/su_data_vel_traces/combine.T.su';
[data_wiggle,header,~] = ReadSu(su_file_wiggle);
header([8354,8011,8006,8181,8178]) = [];
data_wiggle(:,[8354,8011,8006,8181,8178]) = [];
ntr = length(header);
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
[tmp_tab_offset,order_offset] = sort(tmp_tab(:,5));
data_wiggle_offset = data_wiggle(:,order_offset);
fp=fopen('/home0/cxd170430/codes/matlab/yellowstone_imaging/output/T_wiggle.bin','wb');
fwrite(fp,data_wiggle_offset,'float64');
fclose(fp);
% tr=4128;
% figure;plot(data_wiggle_offset(:,tr)/max(abs(data_wiggle_offset(:,tr))));
% hold on;plot(data_env_offset(:,tr)/max(data_env_offset(:,tr)));
% hold on;plot(data_T(:,tr)/max(data_T(:,tr)));

% bin stack
bins = 100:50:10000; % meters
data_wiggle_stack = zeros(2001,length(bins));
for j=1:length(bins)
    if j==1
        index = find(tmp_tab_offset<=bins(1));
    else
        index = find(tmp_tab_offset>bins(j-1) & tmp_tab_offset<=bins(j));
    end
    data_wiggle_stack(:,j) = sum(data_wiggle_offset(:,index),2);
end
figure;imagesc(data_wiggle_stack);clim([-0.6e4,0.6e4]);

