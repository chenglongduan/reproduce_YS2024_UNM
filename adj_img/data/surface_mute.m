
%% STA/LTA PP image
clear;

fp=fopen('./vz_pp_dt0_final.rsf@','rb');
pp_dt0=fread(fp,'float32');fclose(fp);pp_dt0=reshape(pp_dt0,251,751);

pp_top = medfilt2(pp_dt0(1:18,:),[15 1],"symmetric");
pp_top(pp_top<-6) = 0;
pp_top1 = medfilt2(pp_top,[3 9],"symmetric");
pp_final = [pp_top1;pp_dt0(19:251,:)];

fp=fopen('./vz_pp_dt0_final.rsf@','wb');
fwrite(fp,pp_final,'float32');fclose(fp);

% double check
fp=fopen('./vz_pp_dt0_final.rsf@','rb');
pp_dt0_new=fread(fp,'float32');fclose(fp);pp_dt0_new=reshape(pp_dt0_new,251,751);
diff1 = pp_dt0(1:18,:) - pp_dt0_new(1:18,:);
diff2 = pp_dt0(19:251,:) - pp_dt0_new(19:251,:);
[min(diff1(:)),max(diff1(:))]
[min(diff2(:)),max(diff2(:))]

%% STA/LTA PS_R image
clear;

fp=fopen('./vr_ps_dt0_final.rsf@','rb');
ps_dt0=fread(fp,'float32');fclose(fp);ps_dt0=reshape(ps_dt0,251,751);

ps_top = medfilt2(ps_dt0(1:8,:),[15 1],"symmetric");
ps_top(ps_top<-3) = 0;
ps_top1 = medfilt2(ps_top,[51 7],"symmetric");
ps_final = [ps_top1;ps_dt0(9:251,:)];

fp=fopen('./vr_ps_dt0_final.rsf@','wb');
fwrite(fp,ps_final,'float32');fclose(fp);

% double check
fp=fopen('./vr_ps_dt0_final.rsf@','rb');
ps_dt0_new=fread(fp,'float32');fclose(fp);ps_dt0_new=reshape(ps_dt0_new,251,751);
diff1 = ps_dt0(1:8,:) - ps_dt0_new(1:8,:);
diff2 = ps_dt0(9:251,:) - ps_dt0_new(9:251,:);
[min(diff1(:)),max(diff1(:))]
[min(diff2(:)),max(diff2(:))]

%% STA/LTA PS_T image
clear;

fp=fopen('./vt_ps_dt0_final.rsf@','rb');
ps_dt0=fread(fp,'float32');fclose(fp);ps_dt0=reshape(ps_dt0,251,751);

ps_top = medfilt2(ps_dt0(1:8,:),[15 1],"symmetric");
ps_top(ps_top<-3) = 0;
ps_top1 = medfilt2(ps_top,[51 7],"symmetric");
ps_final = [ps_top1;ps_dt0(9:251,:)];

fp=fopen('./vt_ps_dt0_final.rsf@','wb');
fwrite(fp,ps_final,'float32');fclose(fp);

% double check
fp=fopen('./vt_ps_dt0_final.rsf@','rb');
ps_dt0_new=fread(fp,'float32');fclose(fp);ps_dt0_new=reshape(ps_dt0_new,251,751);
diff1 = ps_dt0(1:8,:) - ps_dt0_new(1:8,:);
diff2 = ps_dt0(9:251,:) - ps_dt0_new(9:251,:);
[min(diff1(:)),max(diff1(:))]
[min(diff2(:)),max(diff2(:))]

%% STA/LTA PS_R+T image
clear;

fp=fopen('./R_plus_T_final.rsf@','rb');
ps_dt0=fread(fp,'float32');fclose(fp);ps_dt0=reshape(ps_dt0,251,751);

ps_top = medfilt2(ps_dt0(1:8,:),[15 1],"symmetric");
ps_top(ps_top<-3) = 0;
ps_top1 = medfilt2(ps_top,[51 7],"symmetric");
ps_final = [ps_top1;ps_dt0(9:251,:)];

fp=fopen('./R_plus_T_final.rsf@','wb');
fwrite(fp,ps_final,'float32');fclose(fp);

% double check
fp=fopen('./R_plus_T_final.rsf@','rb');
ps_dt0_new=fread(fp,'float32');fclose(fp);ps_dt0_new=reshape(ps_dt0_new,251,751);
diff1 = ps_dt0(1:8,:) - ps_dt0_new(1:8,:);
diff2 = ps_dt0(9:251,:) - ps_dt0_new(9:251,:);
[min(diff1(:)),max(diff1(:))]
[min(diff2(:)),max(diff2(:))]



