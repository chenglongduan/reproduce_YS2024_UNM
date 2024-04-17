% Make 1D numerical models
% Author: Chenglong Duan, 2023
% This script should work together with "velocity_drop.m" to generate geometry files for modeling

clear;
addpath('/home0/cxd170430/codes/matlab/util');

%% constant models

xmax = 30000; % unit: m
zmax = 10000; % unit: m
dh = 20; % unit: m
vp_bg = 4500; vs_bg = 2500; % background velocity: m/s
rho_model = 'gardner'; % 'bg' or 'gardner' or 'barton'

if strcmp(rho_model,'gardner')
    rho_bg = gardner_rho(vp_bg); % unit: kg/m3
end
if strcmp(rho_model,'barton')
    rho_bg = barton_rho(vp_bg,'mean');
end

[nx,nz,~,~] = dist2grid(xmax, zmax, dh, 0, 0);
vp_homo = ones(nz,nx)*vp_bg;
vs_homo = ones(nz,nx)*vs_bg;
rho_homo = ones(nz,nx)*rho_bg;

figure;subplot(2,2,1);imagesc(vp_homo);colorbar;
subplot(2,2,2);imagesc(vs_homo);colorbar;
subplot(2,2,3);imagesc(rho_homo);colorbar;

% OPTIONAL...
fp=fopen('/home0/cxd170430/codes/ACOUSTIC2D_RTM/yellowstone_field_data/vphomo_ys.bin','wb');
fwrite(fp,vp_homo,'float32');fclose(fp);
fp=fopen('/home0/cxd170430/codes/ACOUSTIC2D_RTM/yellowstone_field_data/vshomo_ys.bin','wb');
fwrite(fp,vs_homo,'float32');fclose(fp);
fp=fopen('/home0/cxd170430/codes/ACOUSTIC2D_RTM/yellowstone_field_data/rhohomo_ys.bin','wb');
fwrite(fp,rho_homo,'float32');fclose(fp);



%% low velocity fracture model: 2 layers

xmax = 30000; % unit: m
zmax = 10000; % unit: m
dh = 5; % unit: m
vp_bg = 4500; vs_bg = 2500; % background velocity: m/s
dln_vp = 0.48601;  % vp perturbation
dln_vs = 0.67232;  % vs perturbation
rho_model = 'gardner'; % 'bg' or 'gardner' or 'barton'
interface = 4000; % unit: m
%=========================================================================

if strcmp(rho_model,'gardner')
    rho_bg = gardner_rho(vp_bg); % unit: kg/m3
end
if strcmp(rho_model,'barton')
    rho_bg = barton_rho(vp_bg,'mean');
end

vp_drop = vp_bg*(1-dln_vp);
vs_drop = vs_bg*(1-dln_vs);
if vs_drop >= vp_drop
    error('vs_drop >= vp_drop, not physically hold!');
end
if strcmp(rho_model,'gardner')
    rho_drop = gardner_rho(vp_drop);
end
if strcmp(rho_model,'barton')
    rho_drop = barton_rho(vp_drop,'mean');
end
if strcmp(rho_model,'bg')
    rho_drop = rho_bg;
end
dln_rho = (1-rho_drop/rho_bg)*100;
fprintf('drho/rho = %f %%\n',dln_rho);

nx = floor(xmax/dh); nz = floor(zmax/dh); iz = floor(interface/dh);
vp_low1=zeros(nz,nx); vs_low1=zeros(nz,nx); rho_low1=zeros(nz,nx);

vp_low1(1:iz,:) = vp_bg; vp_low1(iz+1:nz,:) = vp_drop;
vs_low1(1:iz,:) = vs_bg; vs_low1(iz+1:nz,:) = vs_drop;
rho_low1(1:iz,:) = rho_bg; rho_low1(iz+1:nz,:) = rho_drop;

figure;subplot(2,2,1);imagesc(vp_low1);colorbar;
subplot(2,2,2);imagesc(vs_low1);colorbar;
subplot(2,2,3);imagesc(rho_low1);colorbar;
[min(vp_low1(:)),max(vp_low1(:));min(vs_low1(:)),max(vs_low1(:));min(rho_low1(:)),max(rho_low1(:))]

% OPTIONAL...
NPOINTS_PML = 40;
freesurface = 0;
fout1 = strcat('/home0/cxd170430/codes/FWI_CLDUAN_7_3_2023/ys_',num2str(dln_vp*100),'_',num2str(dln_vs*100),'.vp');
fout2 = strcat('/home0/cxd170430/codes/FWI_CLDUAN_7_3_2023/ys_',num2str(dln_vp*100),'_',num2str(dln_vs*100),'.vs');
fout3 = strcat('/home0/cxd170430/codes/FWI_CLDUAN_7_3_2023/ys_',num2str(dln_vp*100),'_',num2str(dln_vs*100),'.rho');
vp_pad = pad_bin(vp_low1, NPOINTS_PML, freesurface, fout1);
vs_pad = pad_bin(vs_low1, NPOINTS_PML, freesurface, fout2);
rho_pad = pad_bin(rho_low1, NPOINTS_PML, freesurface, fout3);


%% low velocity fracture model: 1 thin layer

xmax = 30000; % unit: m
zmax = 10000; % unit: m
dh = 5; % unit: m
vp_bg = 4500; vs_bg = 2500; % background velocity: m/s
dln_vp = 0.23204;  % vp perturbation
dln_vs = 0.1732;  % vs perturbation
interface = 3800; % unit: m
thickness = 1000; % unit: m
gap = 1000; % unit: m
%=========================================================================

%rho_bg = 310*vp_bg^0.25; % unit: kg/m3
rho_bg = gardner_rho(vp_bg); % unit: kg/m3

vp_drop = vp_bg*(1-dln_vp);
vs_drop = vs_bg*(1-dln_vs);
if vs_drop >= vp_drop
    error('vs_drop >= vp_drop, not physically hold!');
end
%rho_drop = 310*vp_drop^0.25;
rho_drop = gardner_rho(vp_drop);
%rho_drop = rho_bg;
dln_rho = (1-rho_drop/rho_bg)*100;
fprintf('drho/rho = %f %%\n',dln_rho);

nx = floor(xmax/dh); nz = floor(zmax/dh);
iz1 = floor(interface/dh); iz2 = iz1+floor(thickness/dh);
vp_low2=zeros(nz,nx); vs_low2=zeros(nz,nx); rho_low2=zeros(nz,nx);

vp_low2(:,:)=vp_bg; vs_low2(:,:)=vs_bg; rho_low2(:,:)=rho_bg;
if gap>0
    ixx = floor((1000:gap:xmax)/dh);
    for i=1:floor(length(ixx)/2)
        ix1 = ixx((i-1)*2+1); ix2 = ixx(i*2);
        if ix2==nx
            break;
        end
        vp_low2(iz1:iz2,ix1:ix2)=vp_drop;
        vs_low2(iz1:iz2,ix1:ix2)=vs_drop;
        rho_low2(iz1:iz2,ix1:ix2)=rho_drop;
    end
else
    vp_low2(iz1:iz2,:)=vp_drop;
    vs_low2(iz1:iz2,:)=vs_drop;
    rho_low2(iz1:iz2,:)=rho_drop;
end

figure;subplot(2,2,1);imagesc(vp_low2);colorbar;
subplot(2,2,2);imagesc(vs_low2);colorbar;
subplot(2,2,3);imagesc(rho_low2);colorbar;
[min(vp_low2(:)),max(vp_low2(:));min(vs_low2(:)),max(vs_low2(:));min(rho_low2(:)),max(rho_low2(:))]

% OPTIONAL...
NPOINTS_PML = 40;
freesurface = 0;
fout1 = '/data/Chenglong/backup_hopper/FWI_CLDUAN_updated/project_twolayer/model/pad/synthetic/test_1kmgap.vp';
fout2 = '/data/Chenglong/backup_hopper/FWI_CLDUAN_updated/project_twolayer/model/pad/synthetic/test_1kmgap.vs';
fout3 = '/data/Chenglong/backup_hopper/FWI_CLDUAN_updated/project_twolayer/model/pad/synthetic/test_1kmgap.rho';
vp_pad = pad_bin(vp_low2, NPOINTS_PML, freesurface, fout1);
vs_pad = pad_bin(vs_low2, NPOINTS_PML, freesurface, fout2);
rho_pad = pad_bin(rho_low2, NPOINTS_PML, freesurface, fout3);


%% Point Spread Function (PSF)
xmax = 30000; % unit: m
zmax = 10000; % unit: m
dh = 5; % unit: m
vp_bg = 4500; vs_bg = 2500; % background velocity: m/s
dln_vp = 0.23204;  % vp perturbation
dln_vs = 0.1732;  % vs perturbation
interface = 3800; % unit: m
x_points = 29000; % unit: m
radius = 20; % unit: m
%=========================================================================

%rho_bg = 310*vp_bg^0.25; % unit: kg/m3
rho_bg = gardner_rho(vp_bg); % unit: kg/m3

vp_drop = vp_bg*(1-dln_vp);
vs_drop = vs_bg*(1-dln_vs);
if vs_drop >= vp_drop
    error('vs_drop >= vp_drop, not physically hold!');
end
%rho_drop = 310*vp_drop^0.25;
rho_drop = gardner_rho(vp_drop);
%rho_drop = rho_bg;
dln_rho = (1-rho_drop/rho_bg)*100;
fprintf('drho/rho = %f %%\n',dln_rho);

nx = floor(xmax/dh); nz = floor(zmax/dh);
iz_pts = floor(interface/dh)+1; ix_pts = floor(x_points/dh)+1;
vp_low2=zeros(nz,nx); vs_low2=zeros(nz,nx); rho_low2=zeros(nz,nx);

vp_low2(:,:)=vp_bg; vs_low2(:,:)=vs_bg; rho_low2(:,:)=rho_bg;
for i=1:length(ix_pts)
    ix1=ix_pts(i)-floor(radius/dh); ix2=ix_pts(i)+floor(radius/dh);
    iz1=iz_pts-floor(radius/dh); iz2=iz_pts+floor(radius/dh);
    vp_low2(iz1:iz2,ix1:ix2)=vp_drop;
    vs_low2(iz1:iz2,ix1:ix2)=vs_drop;
    rho_low2(iz1:iz2,ix1:ix2)=rho_drop;
end

figure;subplot(2,2,1);imagesc(vp_low2);colorbar;
subplot(2,2,2);imagesc(vs_low2);colorbar;
subplot(2,2,3);imagesc(rho_low2);colorbar;
[min(vp_low2(:)),max(vp_low2(:));min(vs_low2(:)),max(vs_low2(:));min(rho_low2(:)),max(rho_low2(:))]

% OPTIONAL...
NPOINTS_PML = 20;
freesurface = 0;
fout1 = '/home0/cxd170430/NatureGeo/synthetic/model/syn_29km.vp';
fout2 = '/home0/cxd170430/NatureGeo/synthetic/model/syn_29km.vs';
fout3 = '/home0/cxd170430/NatureGeo/synthetic/model/syn_29km.rho';
vp_pad = pad_bin(vp_low2, NPOINTS_PML, freesurface, fout1);
vs_pad = pad_bin(vs_low2, NPOINTS_PML, freesurface, fout2);
rho_pad = pad_bin(rho_low2, NPOINTS_PML, freesurface, fout3);



%% 1 dipping layer model

xmax = 15001; % unit: m
zmax = 8001;  % unit: m
dh = 10;  % unit: m
vp1 = 2000; % 1st layer Vp
vp2 = 2500; % 2nd layer Vp
slope = (4500-3500)/xmax;
%=========================================================================

nx = floor((xmax-1)/dh); nz = floor((zmax-1)/dh);

vs1 = vp1/sqrt(3); rho1 = 310*vp1^0.25;
vs2 = vp2/sqrt(3); rho2 = 310*vp2^0.25;

% true model
z = @(x) 3500+slope*x;

vp_dip=zeros(nz,nx); vs_dip=zeros(nz,nx); rho_dip=zeros(nz,nx);
for ix=1:nx
    value_x = (ix-1)*dh;
    for iz=1:nz
        value_z = (iz-1)*dh;
        if value_z < z(value_x)
            vp_dip(iz,ix) = vp1;
            vs_dip(iz,ix) = vs1;
            rho_dip(iz,ix) = rho1;
        else
            vp_dip(iz,ix) = vp2;
            vs_dip(iz,ix) = vs2;
            rho_dip(iz,ix) = rho2;
        end
    end
end

figure;subplot(2,2,1);imagesc(vp_dip);colorbar;
subplot(2,2,2);imagesc(vs_dip);colorbar;
subplot(2,2,3);imagesc(rho_dip);colorbar;
[min(vp_dip(:)),max(vp_dip(:));min(vs_dip(:)),max(vs_dip(:));min(rho_dip(:)),max(rho_dip(:))]

% OPTIONAL...
NPOINTS_PML = 40;
freesurface = 0;
vp_pad = pad_bin(vp_dip, NPOINTS_PML, freesurface, '/home0/cxd170430/codes/FWI_CLDUAN_7_3_2023/ys_dip.vp');
vs_pad = pad_bin(vs_dip, NPOINTS_PML, freesurface, '/home0/cxd170430/codes/FWI_CLDUAN_7_3_2023/ys_dip.vs');
rho_pad = pad_bin(rho_dip, NPOINTS_PML, freesurface, '/home0/cxd170430/codes/FWI_CLDUAN_7_3_2023/ys_dip.rho');


%% fracture model (yellowstone)


fp=fopen('/data/Chenglong/fwi_test/model/vpinit_homo.bin','wb');fwrite(fp,vp,'float32');
fp=fopen('/data/Chenglong/fwi_test/model/vsinit_homo.bin','wb');fwrite(fp,vs,'float32');
fp=fopen('/data/Chenglong/fwi_test/model/rhoinit_homo.bin','wb');fwrite(fp,rho,'float32');

%2. true (trucated)
vp1=vp; vp1(395:405,400:1100)=2500; %2000;
vs1=vs; vs1(395:405,400:1100)=1400; %1000;
rho1=rho; rho1(395:405,400:1100)=310*2500^0.25; %2073;

figure;imagesc(vp1)
figure;imagesc(vs1)
figure;imagesc(rho1)

fp=fopen('/data/Chenglong/fwi_test/model/vp.bin','wb');fwrite(fp,vp1,'float32');
fp=fopen('/data/Chenglong/fwi_test/model/vs.bin','wb');fwrite(fp,vs1,'float32');
fp=fopen('/data/Chenglong/fwi_test/model/rho.bin','wb');fwrite(fp,rho1,'float32');

%3. true (continuous)
vp2=vp; vp2(395:405,:)=2500;
vs2=vs; vs2(395:405,:)=1400;
rho2=rho; rho2(395:405,:)=310*2500^0.25;

figure;imagesc(vp2)
figure;imagesc(vs2)
figure;imagesc(rho2)

fp=fopen('/data/Chenglong/fwi_test/model/vp_c.bin','wb');fwrite(fp,vp2,'float32');
fp=fopen('/data/Chenglong/fwi_test/model/vs_c.bin','wb');fwrite(fp,vs2,'float32');
fp=fopen('/data/Chenglong/fwi_test/model/rho_c.bin','wb');fwrite(fp,rho2,'float32');

%4. true (point)
vp3=vp; vp3(395:405,740:760)=2500;
vs3=vs; vs3(395:405,740:760)=1400;
rho3=rho; rho3(395:405,740:760)=310*2500^0.25;

figure;imagesc(vp3)
figure;imagesc(vs3)
figure;imagesc(rho3)

fp=fopen('/data/Chenglong/fwi_test/model/vp_p.bin','wb');fwrite(fp,vp3,'float32');
fp=fopen('/data/Chenglong/fwi_test/model/vs_p.bin','wb');fwrite(fp,vs3,'float32');
fp=fopen('/data/Chenglong/fwi_test/model/rho_p.bin','wb');fwrite(fp,rho3,'float32');

% smooth true
%vp2 = imgaussfilt(vp1,20); figure;imagesc(vp2);



