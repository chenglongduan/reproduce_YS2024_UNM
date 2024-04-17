

%% select source wavelet
% obs sta/lta data
geom_tab(itr,5)   % raw trace offset
ceil(geom_tab(itr,5)/100)   % stacked trace offset number
itr=5000;
figure;subplot(2,1,1);plot(data_Z(:,itr));hold on;plot(Z_stack(:,41));
subplot(2,1,2);plot(ez_50_50(:,itr),'r');
cp1=304;cp2=403;test_wlt = data_Z(:,itr);test_wlt(1:cp1-1)=test_wlt(cp1);test_wlt(cp2+1:end)=test_wlt(cp2);
figure;plot(test_wlt);
test_wlt=test_wlt-min(test_wlt);
figure;plot(test_wlt);
%fp=fopen('/home0/cxd170430/codes/ACOUSTIC2D_RTM/spectra/test_wlt.bin','wb');fwrite(fp,test_wlt,'float32');fclose(fp);

itr=6000;
figure;subplot(2,1,1);plot(data_Z(:,itr));hold on;plot(Z_stack(:,49));
subplot(2,1,2);plot(ez_50_50(:,itr),'r');
itr=7000;
figure;subplot(2,1,1);plot(data_Z(:,itr));hold on;plot(Z_stack(:,58));
subplot(2,1,2);plot(ez_50_50(:,itr),'r');


% match gaussian wavelet
amp=4.8e3;td=1.08;f0=5;
for it=1:nt
    tau = pi*(it*dt-1.5/f0-td)*f0;
    wlt(it) = exp(-tau*tau)/(2.0*pi*pi*f0*f0) * amp; % + 0.866
end
figure;plot(wlt)
hold on;plot(test_wlt,'r');

%fp=fopen('/home0/cxd170430/codes/ACOUSTIC2D_RTM/spectra/syn_Z_envelope.bin','wb');fwrite(fp,ez_50_50,'float32');fclose(fp);


% obs raw wiggle data
% match ricker wavelet
amp=1;td=1;f0=10;
for it=1:nt
    tau = pi*(it*dt-1.5/f0-td)*f0;
    wlt1(it) = ((1.0-2.0*tau*tau)*exp(-tau*tau)) * amp;
end
wlt1_env = envelope(wlt1);
figure;plot(wlt1);hold on;plot(wlt1_env);
%fp=fopen('/home0/cxd170430/codes/ACOUSTIC2D_RTM/spectra/test_wlt.bin','wb');fwrite(fp,wlt1_env,'float32');fclose(fp);


%% 



