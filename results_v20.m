  %% RESULTS


for wIdx = 1:number_of_robots    
    robots(wIdx).EMC_x = [];
    robots(wIdx).EMC_y = [];    
    robots(wIdx).EMC_z = [];
    
    robots(wIdx).EMC_v = [];
    
    robots(wIdx).EMC_x_meas = [];
    robots(wIdx).EMC_y_meas = [];
    robots(wIdx).EMC_z_meas = [];
    
    robots(wIdx).EMC_x_double = [];
    robots(wIdx).EMC_y_double = [];    
    robots(wIdx).EMC_z_double = [];   
    
        
    robots(wIdx).EMC_x_double_2 = [];
    robots(wIdx).EMC_y_double_2 = [];    
    robots(wIdx).EMC_z_double_2 = [];
end


EMC_X = [];
EMC_Y = [];
EMC_Z = [];

EV_T = [];


for wIdx = 1:number_of_robots
    figure(510+wIdx)    
    hold on;
    
    for uIdx = 1:number_of_robots
            if uIdx == 1
                plot3(robots(uIdx).pos_hist(1,end-1), robots(uIdx).pos_hist(2,end-1), robots(uIdx).pos_hist(3,end-1),...
                    'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',color_1(uIdx,:));
            else
                plot3(robots(uIdx).pos_hist(1,end), robots(uIdx).pos_hist(2,end), robots(uIdx).pos_hist(3,end),...
                    'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',color_1(uIdx,:));
            end
    end
    
%         % Plot pos_est_double positions
    for uIdx = 1: 6: length(robots(wIdx).pose_est_double(1,:))
        x_data = uIdx;
        y_data = uIdx+1;
        z_data = uIdx+2;
        
        x_data_2 = uIdx+3;
        y_data_2 = uIdx+4;
        z_data_2 = uIdx+5;
        
        % Show Estimated position relative to robot wIdx
        
        for cont = 1 : number_of_robots-1
            if cont >= wIdx
                plot3(robots(wIdx).pose_est_double(cont,x_data), robots(wIdx).pose_est_double(cont,y_data),...
                    robots(wIdx).pose_est_double(cont,z_data), 'o','Markersize', 10,'LineWidth',1.5,'Color',color_2(cont+1,:));
                
                plot3(robots(wIdx).pose_est_double(cont,x_data_2), robots(wIdx).pose_est_double(cont,y_data_2),...
                    robots(wIdx).pose_est_double(cont,z_data_2), 'o','Markersize', 10,'LineWidth',1.5,'Color',color_2(cont+1,:));
            else
                plot3(robots(wIdx).pose_est_double(cont,x_data), robots(wIdx).pose_est_double(cont,y_data),...
                    robots(wIdx).pose_est_double(cont,z_data), 'o','Markersize', 10,'LineWidth',1.5,'Color',color_2(cont,:));
%                 
                plot3(robots(wIdx).pose_est_double(cont,x_data_2), robots(wIdx).pose_est_double(cont,y_data_2),...
                    robots(wIdx).pose_est_double(cont,z_data_2), 'o','Markersize', 10,'LineWidth',1.5,'Color',color_2(cont,:));
            end
            %             uIdx
            %             pause
        end
        

        
    end
    
    % Error absoluto doble soluciones
    robots(wIdx).EMC_x_double = robots(wIdx).pose_real(:,(1:3:end)) - robots(wIdx).pose_est_double(:,(1:6:end));
    robots(wIdx).EMC_y_double = robots(wIdx).pose_real(:,(2:3:end)) - robots(wIdx).pose_est_double(:,(2:6:end));
    robots(wIdx).EMC_z_double = robots(wIdx).pose_real(:,(3:3:end)) - robots(wIdx).pose_est_double(:,(3:6:end));
    
    robots(wIdx).EMC_x_double_2 = robots(wIdx).pose_real(:,(1:3:end)) - robots(wIdx).pose_est_double(:,(4:6:end));
    robots(wIdx).EMC_y_double_2 = robots(wIdx).pose_real(:,(2:3:end)) - robots(wIdx).pose_est_double(:,(5:6:end));
    robots(wIdx).EMC_z_double_2 = robots(wIdx).pose_real(:,(3:3:end)) - robots(wIdx).pose_est_double(:,(6:6:end));
    
%             robots(wIdx).EMC_x_meas = [robots(wIdx).EMC_x_meas, E_X];
%         robots(wIdx).EMC_y_meas = [robots(wIdx).EMC_y_meas, E_Y];
%         robots(wIdx).EMC_z_meas = [robots(wIdx).EMC_z_meas, E_Z];
    
    %
    
    EV = [];
    for uIdx = 1: 3: length(robots(wIdx).pose_real(1,:))
        x_data = uIdx;
        y_data = uIdx+1;
        z_data = uIdx+2;

                        
        % Show Real positions on trayectory
        
        for cont = 1 : number_of_robots-1
            if cont >= wIdx
                plot3(robots(wIdx).pose_real(cont,x_data), robots(wIdx).pose_real(cont,y_data),...
                    robots(wIdx).pose_real(cont,z_data), '*','Color',color_2(cont+1,:),'LineWidth',1.5);
            else
                plot3(robots(wIdx).pose_real(cont,x_data), robots(wIdx).pose_real(cont,y_data),...
                    robots(wIdx).pose_real(cont,z_data), '*','Color',color_2(cont,:),'LineWidth',1.5);
            end
        end

% 
%         Show Estimated position relative to robot wIdx
        
        for cont = 1 : number_of_robots-1
            if cont >= wIdx
                plot3(robots(wIdx).pose_est(cont,x_data), robots(wIdx).pose_est(cont,y_data),...
                    robots(wIdx).pose_est(cont,z_data), 'o','Markersize', 10,'LineWidth',1.5,'Color',color_2(cont+1,:));
            else
                plot3(robots(wIdx).pose_est(cont,x_data), robots(wIdx).pose_est(cont,y_data),...
                    robots(wIdx).pose_est(cont,z_data), 'o','Markersize', 10,'LineWidth',1.5,'Color',color_2(cont,:));
            end
%             pause 
%             uIdx
        end

        set(gca,'FontSize',11); set(gcf,'Color','White');
        
        % Error Real position and Estimation         
                       
        
        
        % Error absoluto
        E_X = robots(wIdx).pose_real(:,x_data) - robots(wIdx).pose_est(:,x_data);
        E_Y = robots(wIdx).pose_real(:,y_data) - robots(wIdx).pose_est(:,y_data);
        E_Z = robots(wIdx).pose_real(:,z_data) - robots(wIdx).pose_est(:,z_data);
        
        
        robots(wIdx).EMC_x = [robots(wIdx).EMC_x, E_X];
        robots(wIdx).EMC_y = [robots(wIdx).EMC_y, E_Y];
        robots(wIdx).EMC_z = [robots(wIdx).EMC_z, E_Z];
        
        % Error en la medida
        E_X = robots(wIdx).pose_real(:,x_data) - robots(wIdx).pose_meas(:,x_data);
        E_Y = robots(wIdx).pose_real(:,y_data) - robots(wIdx).pose_meas(:,y_data);
        E_Z = robots(wIdx).pose_real(:,z_data) - robots(wIdx).pose_meas(:,z_data);
        
        robots(wIdx).EMC_x_meas = [robots(wIdx).EMC_x_meas, E_X];
        robots(wIdx).EMC_y_meas = [robots(wIdx).EMC_y_meas, E_Y];
        robots(wIdx).EMC_z_meas = [robots(wIdx).EMC_z_meas, E_Z];        
        
    end
    
    
   
%     
%     
    % Plot history positions
    for uIdx = 1:number_of_robots
        %         if wIdx ~= uIdx
        if uIdx == 1
            plot3(robots(uIdx).pos_hist(1,1:end-1), robots(uIdx).pos_hist(2,1:end-1), robots(uIdx).pos_hist(3,1:end-1),...
                'Color',color_1(uIdx,:),'LineWidth',1);
        else
            plot3(robots(uIdx).pos_hist(1,:), robots(uIdx).pos_hist(2,:), robots(uIdx).pos_hist(3,:),...
                'Color',color_1(uIdx,:),'LineWidth',1);
        end
    end
    
    
    xlabel('x(m)')
    ylabel('y(m)')
    zlabel('z(m)')
    axis square
    view(45,45)
    grid on
    
    
    legend(['Robot ',num2str(wIdx)],'Real','Estimation');
    
    EV_T = [];
    for i=1 : number_of_robots-1
        E_V = [];
        for j=1 : length(robots(wIdx).EMC_y(1,:))
            %             E_V = [E_V, norm([robots(wIdx).EMC_x(i,j), robots(wIdx).EMC_y(i,j), robots(wIdx).EMC_z(i,j)])];
            
            % Magnitud Vector Error
            E_V = [E_V, norm([robots(wIdx).EMC_x(i,j), robots(wIdx).EMC_y(i,j),robots(wIdx).EMC_z(i,j)])];
        end
        EV_T = [EV_T;E_V];
    end
    robots(wIdx).EMC_v = [robots(wIdx).EMC_v; EV_T];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % SAVE
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    % for wIdx = 1 : number_of_robots
    %     cadena_filename = ['errores',num2str(wIdx),'_v.mat'];    
    %     Ev = robots(wIdx).EMC_v;
    %     save(cadena_filename,'Ev')       
    % end    
    
    % MOSTRAR PANTALLA
    % Save first time
    
%     if test == 1
%         ECM_d_estimacion_v = sqrt( mean( (robots(wIdx).EMC_v).^2,2))     % MSE
%         E_Vector = robots(wIdx).EMC_v;                          % Vector 4 puntos
%         cadena_filename = ['A_ERR_BLE', num2str(wIdx),'_',num2str(sigma_d),'.mat'];
%         save (cadena_filename,'E_Vector')
%     else        
%         % Append data
%         cadena_filename = ['A_ERR_BLE', num2str(wIdx),'_',num2str(sigma_d),'.mat'];
%         errores(wIdx).V = load (cadena_filename);
%         errores(wIdx).V.E_Vector = [errores(wIdx).V.E_Vector ; robots(wIdx).EMC_v];
%         E_Vector = errores(wIdx).V.E_Vector;
%         save(cadena_filename,'E_Vector')
%     end

%     Z _ noise
%     if test == 1
%         ECM_d_estimacion_v = sqrt( mean( (robots(wIdx).EMC_v).^2,2))     % MSE
%         E_Vector = robots(wIdx).EMC_v;                          % Vector 4 puntos
%         cadena_filename = ['A_ERR_BLE_Z', num2str(wIdx),'_',num2str(sigma_d),'.mat'];
%         save (cadena_filename,'E_Vector')
%     else        
% %         Append data
%         cadena_filename = ['A_ERR_BLE_Z', num2str(wIdx),'_',num2str(sigma_d),'.mat'];
%         errores(wIdx).V = load (cadena_filename);
%         errores(wIdx).V.E_Vector = [errores(wIdx).V.E_Vector ; robots(wIdx).EMC_v];
%         E_Vector = errores(wIdx).V.E_Vector;
%         save(cadena_filename,'E_Vector')
%     end
%     
    % X and Y data
    
%     if test == 1
% %         ECM_d_estimacion_v = sqrt( mean( (robots(wIdx).EMC_v).^2,2))     % MSE
%         E_Vector = robots(wIdx).EMC_x;                          % Vector 4 puntos
%         cadena_filename = ['A_ERR_BLE_x', num2str(wIdx),'_',num2str(sigma_d),'.mat'];
%         save (cadena_filename,'E_Vector')
%         
%         E_Vector = robots(wIdx).EMC_y;                          % Vector 4 puntos
%         cadena_filename = ['A_ERR_BLE_y', num2str(wIdx),'_',num2str(sigma_d),'.mat'];
%         save (cadena_filename,'E_Vector')
%     else
%         % Append data X
%         cadena_filename = ['A_ERR_BLE_x', num2str(wIdx),'_',num2str(sigma_d),'.mat'];
%         errores(wIdx).V = load (cadena_filename);
%         errores(wIdx).V.E_Vector = [errores(wIdx).V.E_Vector ; robots(wIdx).EMC_x];
%         E_Vector = errores(wIdx).V.E_Vector;
%         save(cadena_filename,'E_Vector')
%         
%         cadena_filename = ['A_ERR_BLE_y', num2str(wIdx),'_',num2str(sigma_d),'.mat'];
%         errores(wIdx).V = load (cadena_filename);
%         errores(wIdx).V.E_Vector = [errores(wIdx).V.E_Vector ; robots(wIdx).EMC_y];
%         E_Vector = errores(wIdx).V.E_Vector;
%         save(cadena_filename,'E_Vector')
%     end
% 
    % RUIDO ODOMETRIA
%     if test == 1
%         ECM_d_estimacion_v = sqrt( mean( (robots(wIdx).EMC_v).^2,2));     % MSE
%         E_Vector = robots(wIdx).EMC_v;                          % Vector 4 puntos
%         cadena_filename = ['OVE', num2str(wIdx),'_',num2str(odometry_noise),'.mat'];
%         save (cadena_filename,'E_Vector')
%     else        
%         % Append data
%         cadena_filename = ['OVE', num2str(wIdx),'_',num2str(odometry_noise),'.mat'];
%         errores(wIdx).V = load (cadena_filename);
%         errores(wIdx).V.E_Vector = [errores(wIdx).V.E_Vector ; robots(wIdx).EMC_v];
%         E_Vector = errores(wIdx).V.E_Vector;
%         save(cadena_filename,'E_Vector')
%     end
    
    
%     % RUIDO COMBINADo
%     if test == 1
%         ECM_d_estimacion_v = sqrt( mean( (robots(wIdx).EMC_v).^2,2));     % MSE
%         E_Vector = robots(wIdx).EMC_v;                          % Vector 4 puntos
%         cadena_filename = ['CVE', num2str(wIdx),'_',num2str(odometry_noise),'.mat'];
%         save (cadena_filename,'E_Vector')
%     else        
%         % Append data
%         cadena_filename = ['CVE', num2str(wIdx),'_',num2str(odometry_noise),'.mat'];
%         errores(wIdx).V = load (cadena_filename);
%         errores(wIdx).V.E_Vector = [errores(wIdx).V.E_Vector ; robots(wIdx).EMC_v];
%         E_Vector = errores(wIdx).V.E_Vector;
%         save(cadena_filename,'E_Vector')
%     end
    
        % MOSTRAR PANTALLA
%     ECM_d_estimacion_v = sqrt( mean(robots(wIdx).EMC_v,2))
%     cadena_filename = ['ECM_od', num2str(wIdx),num2str(odometry_noise),'.mat'];
%     save (cadena_filename,'ECM_d_estimacion_v')
    
            % MOSTRAR PANTALLA
%     ECM_d_estimacion_v = sqrt( mean(robots(wIdx).EMC_v,2))
%     cadena_filename = ['ECM_od', num2str(wIdx),num2str(odometry_noise),'.mat'];
%     save (cadena_filename,'ECM_d_estimacion_v')
    
%     error_cuadra_x = robots(wIdx).EMC_x.^2;
%     error_cuadra_y = robots(wIdx).EMC_y.^2;
%     error_cuadra_z = robots(wIdx).EMC_z.^2;
    
%     emcx(:,wIdx) = sqrt( mean(error_cuadra_x,2))
%     emcy(:,wIdx) = sqrt( mean(error_cuadra_y,2))
%     emcz(:,wIdx) = sqrt( mean(error_cuadra_z,2))
    
    

    
%     text( 0, 0, sprintf( ' MSE estimacion = %f', ECM_d_estimacion_v),   'BackgroundColor',[1  .6 .6]);
    
end

%%
%%%%%%%%%%%%%%%%%%%%%% 
%     # of estimations vs error
%     Time vs Error

    fig3 = 810;
    fig4 = 910;
    fig5 = 1000;
    fig6 = 1100;
    n_points_estimated = timesteps/(inc*number_of_robots)



t_est_r = [];
for wIdx = 1 : number_of_robots
%     uIdx = 
    for uIdx = 1 : number_of_robots-1

         %%%%%%%%%%%
        %   error vs time
        %%%% ERRROES 2D
%         figure(fig3+wIdx)
%         
%         set(gca,'FontSize',11); set(gcf,'Color','White');
%         subplot(2,1,1)
%         
%         
%         t_est_r(:,wIdx) = (0:n_points_estimated-1)*inc*number_of_robots + inc*(wIdx);   % Tiempos en que estima robot
%         
%         plot(t_est_r(:,wIdx),robots(wIdx).EMC_x(uIdx,:),'-*')
% %         title(['Error de estimación Robot ',num2str(wIdx)]);
%         ylabel('Distancia (m)')
%         hold on
%         title('Error Posicion X')
%         subplot(2,1,2)
%         plot(t_est_r(:,wIdx),robots(wIdx).EMC_y(uIdx,:),'-*')
%         ylabel('Distancia (m)')
%         hold on
%         title('Error Posicion Y')
%         xlabel('Time')
        
        
        %%%% ERRROES 3D
        
%         figure(fig4+wIdx)
%         
%         set(gca,'FontSize',11); set(gcf,'Color','White');
% %         subplot(2,1,1)
%         
%         
%         t_est_r(:,wIdx) = (0:n_points_estimated-1)*inc*number_of_robots + inc*(wIdx);   % Tiempos en que estima robot
%         
%         plot(t_est_r(:,wIdx),robots(wIdx).EMC_x(uIdx,:),'-*',t_est_r(:,wIdx),robots(wIdx).EMC_y(uIdx,:),'-*',t_est_r(:,wIdx),robots(wIdx).EMC_z(uIdx,:),'-*')
% %         title(['Error de estimación Robot ',num2str(wIdx)]);
%         ylabel('Error (m)')
%         xlabel('Time (seg)')
%         hold on
%         grid on
        
        %%%%%%%%%%%%%%%%%%%%%%%%
                 %   error vs time
                %%%% ERRROES 2D double
                figure(fig5+wIdx)

                set(gca,'FontSize',11); set(gcf,'Color','White');
                subplot(2,1,1)


                t_est_r(:,wIdx) = (0:n_points_estimated-1)*inc*number_of_robots + inc*(wIdx);   % Tiempos en que estima robot

                plot(t_est_r(:,wIdx),robots(wIdx).EMC_x_double(uIdx,:),'-*')
        %         title(['Error de estimación Robot ',num2str(wIdx)]);
                ylabel('Distancia (m)')
                hold on
                title('Error Posicion X')
                subplot(2,1,2)
                plot(t_est_r(:,wIdx),robots(wIdx).EMC_y_double(uIdx,:),'-*')
                ylabel('Distancia (m)')
                hold on
                title('Error Posicion Y')
                xlabel('Time')

                         %   error vs time
                %%%% ERRROES 2D double
                figure(fig5+wIdx)

                set(gca,'FontSize',11); set(gcf,'Color','White');
                subplot(2,1,1)

                t_est_r(:,wIdx) = (0:n_points_estimated-1)*inc*number_of_robots + inc*(wIdx);   % Tiempos en que estima robot

                plot(t_est_r(:,wIdx),robots(wIdx).EMC_x_double_2(uIdx,:),'-*')
        %         title(['Error de estimación Robot ',num2str(wIdx)]);
                ylabel('Distancia (m)')
                hold on
                title('Error Posicion X')
                subplot(2,1,2)
                plot(t_est_r(:,wIdx),robots(wIdx).EMC_y_double_2(uIdx,:),'-*')
                ylabel('Distancia (m)')
                hold on
                title('Error Posicion Y')
                xlabel('Time')
        %%%%%%%%%%%%%%%%%%%%%%%%
        pause
    end
end
