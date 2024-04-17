% convert from line distance (in km) to UTM XY (in m)

utm_rx_ry_dist_r = [geom_tab(:,3),geom_tab(:,4),geom_tab(:,7)/1000];


dist1 = 5;
index = find(abs(utm_rx_ry_dist_r(:,3)-dist1)<0.05);
fprintf('The nearest to %d is %f\n',dist1,utm_rx_ry_dist_r(index(1),3));
[utm_rx_ry_dist_r(index(1),1),utm_rx_ry_dist_r(index(1),2)]
[Lat,Lon] = utm2deg(utm_rx_ry_dist_r(index(1),1),utm_rx_ry_dist_r(index(1),2),'12 T')


dist1 = 8;
index = find(abs(utm_rx_ry_dist_r(:,3)-dist1)<0.02);
fprintf('The nearest to %d is %f\n',dist1,utm_rx_ry_dist_r(index(1),3));
[utm_rx_ry_dist_r(index(1),1),utm_rx_ry_dist_r(index(1),2)]
[Lat,Lon] = utm2deg(utm_rx_ry_dist_r(index(1),1),utm_rx_ry_dist_r(index(1),2),'12 T')


dist1 = 10;
index = find(abs(utm_rx_ry_dist_r(:,3)-dist1)<0.06);
fprintf('The nearest to %d is %f\n',dist1,utm_rx_ry_dist_r(index(1),3));
[utm_rx_ry_dist_r(index(1),1),utm_rx_ry_dist_r(index(1),2)]
[Lat,Lon] = utm2deg(utm_rx_ry_dist_r(index(1),1),utm_rx_ry_dist_r(index(1),2),'12 T')


dist1 = 12;
index = find(abs(utm_rx_ry_dist_r(:,3)-dist1)<0.02);
fprintf('The nearest to %d is %f\n',dist1,utm_rx_ry_dist_r(index(1),3));
[utm_rx_ry_dist_r(index(1),1),utm_rx_ry_dist_r(index(1),2)]
[Lat,Lon] = utm2deg(utm_rx_ry_dist_r(index(1),1),utm_rx_ry_dist_r(index(1),2),'12 T')


dist1 = 13.5;
index = find(abs(utm_rx_ry_dist_r(:,3)-dist1)<0.005);
fprintf('The nearest to %d is %f\n',dist1,utm_rx_ry_dist_r(index(1),3));
[utm_rx_ry_dist_r(index(1),1),utm_rx_ry_dist_r(index(1),2)]
[Lat,Lon] = utm2deg(utm_rx_ry_dist_r(index(1),1),utm_rx_ry_dist_r(index(1),2),'12 T')


dist1 = 15;
index = find(abs(utm_rx_ry_dist_r(:,3)-dist1)<0.002);
fprintf('The nearest to %d is %f\n',dist1,utm_rx_ry_dist_r(index(1),3));
[utm_rx_ry_dist_r(index(1),1),utm_rx_ry_dist_r(index(1),2)]
[Lat,Lon] = utm2deg(utm_rx_ry_dist_r(index(1),1),utm_rx_ry_dist_r(index(1),2),'12 T')


dist1 = 17;
index = find(abs(utm_rx_ry_dist_r(:,3)-dist1)<0.05);
fprintf('The nearest to %d is %f\n',dist1,utm_rx_ry_dist_r(index(1),3));
[utm_rx_ry_dist_r(index(1),1),utm_rx_ry_dist_r(index(1),2)]
[Lat,Lon] = utm2deg(utm_rx_ry_dist_r(index(1),1),utm_rx_ry_dist_r(index(1),2),'12 T')


dist1 = 18;
index = find(abs(utm_rx_ry_dist_r(:,3)-dist1)<0.03);
fprintf('The nearest to %d is %f\n',dist1,utm_rx_ry_dist_r(index(1),3));
[utm_rx_ry_dist_r(index(1),1),utm_rx_ry_dist_r(index(1),2)]
[Lat,Lon] = utm2deg(utm_rx_ry_dist_r(index(1),1),utm_rx_ry_dist_r(index(1),2),'12 T')


dist1 = 21;
index = find(abs(utm_rx_ry_dist_r(:,3)-dist1)<0.15);
fprintf('The nearest to %d is %f\n',dist1,utm_rx_ry_dist_r(index(1),3));
[utm_rx_ry_dist_r(index(1),1),utm_rx_ry_dist_r(index(1),2)]
[Lat,Lon] = utm2deg(utm_rx_ry_dist_r(index(1),1),utm_rx_ry_dist_r(index(1),2),'12 T')