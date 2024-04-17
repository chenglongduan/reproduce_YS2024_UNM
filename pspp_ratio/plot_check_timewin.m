function plot_check_timewin(ZRT_env_new,Z_env_new,R_env_new,T_env_new,RT_env_new,Z_stack,R_stack,T_stack,RT_stack,geom_tab_new,...
    nt,dt,ntr_new,t_p,t_pp,t_ps,timewin_p_prior,timewin_p_post,timewin_pp_prior,timewin_pp_post,timewin_ps_prior,timewin_ps_post,method)


% combine
if strcmp(method,'combine')
    figure;subplot(3,2,1);
    imagesc(100:100:10000,0:dt:(nt-1)*dt,Z_stack);colormap(colorbar_bwr);clim([0.5,1.5]);xlabel('offset (m)');ylabel('time (s)');
    hold on;plot(geom_tab_new(:,5),t_p-timewin_p_prior,'k');
    hold on;plot(geom_tab_new(:,5),t_p+timewin_p_post,'k');
    hold on;plot(geom_tab_new(:,5),t_pp-timewin_pp_prior,'k');
    hold on;plot(geom_tab_new(:,5),t_pp+timewin_pp_post,'k');title('Z stacked STA/LTA');
    subplot(3,2,2);
    imagesc(100:100:10000,0:dt:(nt-1)*dt,R_stack);colormap(colorbar_bwr);clim([0.5,1.5]);xlabel('offset (m)');ylabel('time (s)');
    hold on;plot(geom_tab_new(:,5),t_ps-timewin_ps_prior,'k');
    hold on;plot(geom_tab_new(:,5),t_ps+timewin_ps_post,'k');title('R stacked STA/LTA');
    subplot(3,2,3);
    imagesc(100:100:10000,0:dt:(nt-1)*dt,T_stack);colormap(colorbar_bwr);clim([0.5,1.5]);xlabel('offset (m)');ylabel('time (s)');
    hold on;plot(geom_tab_new(:,5),t_ps-timewin_ps_prior,'k');
    hold on;plot(geom_tab_new(:,5),t_ps+timewin_ps_post,'k');title('T stacked STA/LTA');
    subplot(3,2,4);
    imagesc(100:100:10000,0:dt:(nt-1)*dt,Z_stack+R_stack+T_stack);colormap(colorbar_bwr);clim([2.0,4.0]);xlabel('offset (m)');ylabel('time (s)');
    hold on;plot(geom_tab_new(:,5),t_pp-timewin_pp_prior,'k');
    hold on;plot(geom_tab_new(:,5),t_pp+timewin_pp_post,'k');
    hold on;plot(geom_tab_new(:,5),t_ps-timewin_ps_prior,'k');
    hold on;plot(geom_tab_new(:,5),t_ps+timewin_ps_post,'k');title('Z+R+T stacked STA/LTA');
    subplot(3,2,5);
    imagesc(1:ntr_new,0:dt:(nt-1)*dt,ZRT_env_new);colormap(colorbar_bwr);clim([0,1500]);xlabel('trace# (offset 2-6km)');ylabel('time (s)');
    hold on;plot(1:ntr_new,t_p-timewin_p_prior,'k');
    hold on;plot(1:ntr_new,t_p+timewin_p_post,'k');
    hold on;plot(1:ntr_new,t_pp-timewin_pp_prior,'k');
    hold on;plot(1:ntr_new,t_pp+timewin_pp_post,'k');
    hold on;plot(1:ntr_new,t_ps-timewin_ps_prior,'k');
    hold on;plot(1:ntr_new,t_ps+timewin_ps_post,'k');title('ZRT envelope');
end


% individual
if strcmp(method,'individual')
    figure;subplot(1,2,1);
    imagesc(100:100:10000,0:dt:(nt-1)*dt,Z_stack);colormap(colorbar_bwr);clim([0.5,1.5]);xlabel('offset (m)');ylabel('time (s)');
    hold on;plot(geom_tab_new(:,5),t_pp-timewin_pp_prior,'k');
    hold on;plot(geom_tab_new(:,5),t_pp+timewin_pp_post,'k');title('Z stacked STA/LTA');
    subplot(1,2,2);
    imagesc(1:ntr_new,0:dt:(nt-1)*dt,Z_env_new);colormap(colorbar_bwr);clim([0,1500]);xlabel('trace# (offset 2-6km)');ylabel('time (s)');
    hold on;plot(1:ntr_new,t_pp-timewin_pp_prior,'k');
    hold on;plot(1:ntr_new,t_pp+timewin_pp_post,'k');title('Z envelope');
    
    figure;subplot(1,2,1);
    imagesc(100:100:10000,0:dt:(nt-1)*dt,R_stack);colormap(colorbar_bwr);clim([0.5,1.5]);xlabel('offset (m)');ylabel('time (s)');
    hold on;plot(geom_tab_new(:,5),t_ps-timewin_ps_prior,'k');
    hold on;plot(geom_tab_new(:,5),t_ps+timewin_ps_post,'k');title('R stacked STA/LTA');
    subplot(1,2,2);
    imagesc(1:ntr_new,0:dt:(nt-1)*dt,R_env_new);colormap(colorbar_bwr);clim([0,1500]);xlabel('trace# (offset 2-6km)');ylabel('time (s)');
    hold on;plot(1:ntr_new,t_ps-timewin_ps_prior,'k');
    hold on;plot(1:ntr_new,t_ps+timewin_ps_post,'k');title('R envelope');
    
    figure;subplot(1,2,1);
    imagesc(100:100:10000,0:dt:(nt-1)*dt,T_stack);colormap(colorbar_bwr);clim([0.5,1.5]);xlabel('offset (m)');ylabel('time (s)');
    hold on;plot(geom_tab_new(:,5),t_ps-timewin_ps_prior,'k');
    hold on;plot(geom_tab_new(:,5),t_ps+timewin_ps_post,'k');title('T stacked STA/LTA');
    subplot(1,2,2);
    imagesc(1:ntr_new,0:dt:(nt-1)*dt,T_env_new);colormap(colorbar_bwr);clim([0,1500]);xlabel('trace# (offset 2-6km)');ylabel('time (s)');
    hold on;plot(1:ntr_new,t_ps-timewin_ps_prior,'k');
    hold on;plot(1:ntr_new,t_ps+timewin_ps_post,'k');title('T envelope');
    
    figure;subplot(1,2,1);
    imagesc(100:100:10000,0:dt:(nt-1)*dt,RT_stack);colormap(colorbar_bwr);clim([0.5,1.5]);xlabel('offset (m)');ylabel('time (s)');
    hold on;plot(geom_tab_new(:,5),t_ps-timewin_ps_prior,'k');
    hold on;plot(geom_tab_new(:,5),t_ps+timewin_ps_post,'k');title('RT stacked STA/LTA');
    subplot(1,2,2);
    imagesc(1:ntr_new,0:dt:(nt-1)*dt,RT_env_new);colormap(colorbar_bwr);clim([0,1500]);xlabel('trace# (offset 2-6km)');ylabel('time (s)');
    hold on;plot(1:ntr_new,t_ps-timewin_ps_prior,'k');
    hold on;plot(1:ntr_new,t_ps+timewin_ps_post,'k');title('RT envelope');
end


end