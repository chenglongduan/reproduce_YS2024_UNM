% source, receiver UTM X,Y to Latitude,Longitude

%% true S-R locations
ns = size(geom_tab,1);
src_latlon = zeros(ns,2);

for i=1:ns

    [src_latlon(i,1),src_latlon(i,2)] = utm2deg(geom_tab(i,1),geom_tab(i,2),'12 T');

end


nr = ns;
rec_latlon = zeros(nr,2);

for i=1:nr

    [rec_latlon(i,1),rec_latlon(i,2)] = utm2deg(geom_tab(i,3),geom_tab(i,4),'12 T');

end

writematrix(src_latlon,'/home0/cxd170430/codes/matlab/yellowstone_imaging/figs/src_latlon.txt','Delimiter','space');
writematrix(rec_latlon,'/home0/cxd170430/codes/matlab/yellowstone_imaging/figs/rec_latlon.txt','Delimiter','space');


%% 2D projected S-R locations
% start and end points
[xmin,ixmin] = min(x_sr_proj);
[ymax,iymax] = max(y_sr_proj);
if ixmin ~= iymax
    error('ixmin not equal iymax!');
end

[xmax,ixmax] = max(x_sr_proj);
[ymin,iymin] = min(y_sr_proj);
if ixmax ~= iymin
    error('ixmax not equal iymin!');
end

[Lat_start,Lon_start] = utm2deg(xmin,ymax,'12 T');

[Lat_end,Lon_end] = utm2deg(xmax,ymin,'12 T');

fprintf('start point:\n %f  %f\n',Lat_start,Lon_start);
fprintf('end point:\n %f  %f\n',Lat_end,Lon_end);

% plot
figure;plot(geom_tab(:,1),geom_tab(:,2),'r*');hold on;plot(geom_tab(:,3),geom_tab(:,4),'bo');hold on;
plot(x_sr_proj(1:ntr),y_sr_proj(1:ntr),'ro');hold on;plot(x_sr_proj(ntr+1:end),y_sr_proj(ntr+1:end),'ko');
axis equal;

% all points
sr_line_latlon = zeros(length(x_sr_proj),2);
for i=1:length(x_sr_proj)
    [sr_line_latlon(i,1),sr_line_latlon(i,2)] = utm2deg(x_sr_proj(i),y_sr_proj(i),'12 T');
end

writematrix(sr_line_latlon,'/home0/cxd170430/codes/matlab/yellowstone_imaging/figs/sr_line_latlon.txt','Delimiter','space');


