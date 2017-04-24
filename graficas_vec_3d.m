% GRAFICAS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT ROBOTS & DISTANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(100+uIdx)
clf
%     subplot(2,number_of_robots,uIdx)
hold on;


% Posiciòn de Robot moviendose W



% if wIdx == 1
    plot3(robots(wIdx).position(1), robots(wIdx).position(2), robots(wIdx).position(3),...
        'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',color(wIdx,:));
    % Posiciòn de Robot Estatico U
    plot3(k_1(1,uIdx),k_1(2,uIdx),k_1(3,uIdx),...
        'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',color(uIdx,:));
    
% elseif wIdx == 2
%     plot3(robots(wIdx).position(1), robots(wIdx).position(2), robots(wIdx).position(3),...
%         'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1]);
    % Posiciòn de Robot Estatico U
%     plot3(k_1(1,uIdx),k_1(2,uIdx),k_1(3,uIdx),...
%         'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0]);
% end

% Vector Distancia K
real_landmark = robots(wIdx).position;
real_landmark2 = robots(uIdx).position;
plot3([real_landmark(1) real_landmark2(1)], [real_landmark(2) real_landmark2(2)], [real_landmark(3) real_landmark2(3)],...
        'LineStyle','--',...
        'Color',[0 0.5 1],...
        'LineWidth',1.2)

% Each robot circle MAGENTA K_1
% Vector Distance K_1
real_landmark = k_1(:,wIdx);
real_landmark2 = robots(uIdx).position;
plot3([real_landmark(1) real_landmark2(1)], [real_landmark(2) real_landmark2(2)], [real_landmark(3) real_landmark2(3)],...
        'LineStyle','--',...
        'Color',[1 0 1],...
        'LineWidth',1.2)


% plot3(robots(wIdx).pos_hist(1,:), robots(wIdx).pos_hist(2,:), robots(wIdx).pos_hist(3,:),'r');

% medida K
% plot3(x_meas, y_meas, z_meas, '*','Color',[0.5,0.8,0.2],'LineWidth',1.5);

% FRAME
MT_1 = robots(wIdx).MT_1;
%              frame_mod(MT,'c',0.5);
MT = robots(wIdx).MT;
%              frame_mod(MT_1,'m',0.5);

% Vector desplazamiento
plot3([MT_1(1,4),MT(1,4)], [MT_1(2,4),MT(2,4)], [MT_1(3,4),MT(3,4)],...
    'k',...
    'LineWidth',1);

% MOVEMENT
% for j =1:number_of_robots
% if wIdx == 1
plot3(robots(wIdx).pos_hist(1,:), robots(wIdx).pos_hist(2,:), robots(wIdx).pos_hist(3,:),...
    'Color',color(wIdx,:),'LineWidth',1);

% plot3(robots(wIdx).pos_hist(1,:), robots(wIdx).pos_hist(2,:), robots(wIdx).pos_hist(3,:),...
%     'Color',color(uIdx,:),'LineWidth',1);

%     plot3(robots(wIdx).pos_hist(1,:), robots(wIdx).pos_hist(2,:), robots(wIdx).pos_hist(3,:),'r',...
%         'LineWidth',1);
%     plot3(robots(uIdx).pos_hist(1,:), robots(uIdx).pos_hist(2,:), robots(uIdx).pos_hist(3,:),'b',...
%         'LineWidth',1);
    
% % elseif wIdx == 2
%     plot3(robots(wIdx).pos_hist(1,:), robots(wIdx).pos_hist(2,:), robots(wIdx).pos_hist(3,:),'b',...
%         'LineWidth',1);
%     plot3(robots(uIdx).pos_hist(1,:), robots(uIdx).pos_hist(2,:), robots(uIdx).pos_hist(3,:),'r',...
%         'LineWidth',1);
% end


% end





% title(['Robot ',num2str(wIdx),' Distance and Travel to Robot ', num2str(uIdx)])
% title(['Robot ',num2str(wIdx),' Distancia y Movimiento con respecto al Robot ', num2str(uIdx)])
set(gca,'FontSize',11); set(gcf,'Color','White');
legend(['Robot ',num2str(wIdx)], ['Robot ',num2str(uIdx)], 'Distancia K','Distancia K-1','Desplazamiento','Movimiento')
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
% view(45,45)
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT ESTIMATED POSITION W NEIGHBOUR ROBOTS EACH REFERENCE FRAME
%%%%%%%%%%%%%%%%%%%%
figure(300+uIdx)
clf
%     subplot(2,number_of_robots,uIdx)
hold on;

% Posiciòn de Robot moviendose W
% if wIdx == 1
%     plot3(0,0,0,...
%         'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0]);
% %     plot3(robots(wIdx).pos_hist(1,:), robots(wIdx).pos_hist(2,:), robots(wIdx).pos_hist(3,:),'r');
%     
% elseif wIdx == 2
%     plot3(0,0,0,...
%         'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1]);
% %     plot3(robots(wIdx).pos_hist(1,:), robots(wIdx).pos_hist(2,:), robots(wIdx).pos_hist(3,:),'b');
% end

    plot3(0,0,0,...
        'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',color(wIdx,:));


% Estimacion 1
plot3(x_estimado, y_estimado, z_estimado, 'o','MarkerSize',10,...
            'MarkerEdgeColor','k','MarkerFaceColor',[255/255 255/255 0/255]);  
        
% Estimación 2
plot3(x_estimado_2, y_estimado_2, z_estimado, 'o','MarkerSize',10,...
            'MarkerEdgeColor','k','MarkerFaceColor',[255/255 198/255 0/255]);     
    
% Real position
% if wIdx == 1
%     plot3(x_real, y_real, z_real, 'b*','MarkerSize',6);   
% elseif wIdx == 2
%     plot3(x_real, y_real, z_real, 'r*','MarkerSize',6);
% end

plot3(x_real, y_real, z_real, '*','MarkerEdgeColor',color(uIdx,:),'MarkerSize',6); 


% Medida
% plot3(x_meas_loc, y_meas_loc, z_meas_loc, '*','Color',[0.5,0.8,0.2],'LineWidth',1.5);

% Resultado
plot3(x_resolve, y_resolve, z_estimado, 'go', 'Markersize', 10,'LineWidth',1.5);  



% Error medio cuadratico
% ECM_d_estimacion_x = [sqrt( mean( ( x_estimado - x_real)^2  ) ) , sqrt( mean( ( x_estimado_2 - x_real)^2  ) )];
% ECM_d_estimacion_y = [sqrt( mean( ( y_estimado - y_real)^2  ) ) , sqrt( mean( ( y_estimado_2 - y_real)^2  ) )];
% ECM_d_estimacion_z = [sqrt( mean( ( z_estimado - z_real)^2  ) ) , sqrt( mean( ( z_estimado - z_real)^2  ) )];

% text( 0, -3, sprintf( 'MSE measurements = %f', ECM_d_estimacion_x),   'BackgroundColor',[1  .6 .6]);
% text( 0, -4, sprintf( 'MSE measurements = %f', ECM_d_estimacion_y),   'BackgroundColor',[1  .6 .6]);
% text( 0, -5, sprintf( 'MSE measurements = %f', ECM_d_estimacion_z),   'BackgroundColor',[1  .6 .6]);

% axis([-10, 10, -10, 10]);
% title(['Robot ',num2str(wIdx),' Estimation to Robot ', num2str(uIdx)])
%              title(['Robot ',num2str(wIdx),' Estimación del Robot ', num2str(uIdx)])
set(gca,'FontSize',11); set(gcf,'Color','White');
xlabel('x')
ylabel('y')
zlabel('z')
legend(['Robot ',num2str(wIdx)],'Estimación 1','Estimación 2',['Robot ',num2str(uIdx)],'Selección')
%              axis([-5, 2, -0, 5, 0, 5]);
% axis([-6, 1, -5, 2, -1,5]);
axis square
axis equal
% view(45,45)
grid on


%%
%% SUBPLOT
% 
% figure(60+wIdx)
% subplot(2,number_of_robots,uIdx)
% hold on;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % PLOT ROBOTS & DISTANCE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %       % Show the sensor measurement as an arrow
% % History and circle W
% plot(robots(wIdx).pos_hist(1,:), robots(wIdx).pos_hist(2,:),'r');
% plot(robots(wIdx).position(1), robots(wIdx).position(2), 'ro');
% 
% % Real Distance K
% line([robots(wIdx).position(1), robots(wIdx).position(1) + cos(robots(wIdx).anglek)*robots(wIdx).distk], ...
%     [robots(wIdx).position(2), robots(wIdx).position(2) + sin(robots(wIdx).anglek)*robots(wIdx).distk],'Color','g');
% 
% % Each robot circle MAGENTA K_1
% %                       plot(robots(wIdx).pk_1(1,uIdx),robots(wIdx).pk_1(2,uIdx),'mo')
% plot(k_1(1,uIdx),k_1(2,uIdx),'mo')
% 
% % Real Distance K_1
% line([k_1(1,wIdx), k_1(1,wIdx) + robots(wIdx).distk_1*cos(robots(wIdx).anglek_1)], ...
%     [k_1(2,wIdx), k_1(2,wIdx) + robots(wIdx).distk_1*sin(robots(wIdx).anglek_1)],'Color','m');
% 
% % uk_wk
% %               line([robot_u.position(1), robot_u.position(1) + 8*cos(theta_uk_wk(1))], ...
% %                   [robot_u.position(2), robot_u.position(2) + 8*sin(theta_uk_wk(1))],'Color','c');
% %               line([robot_u.position(1), robot_u.position(1) + 8*cos(theta_uk_wk(2))], ...
% %                   [robot_u.position(2), robot_u.position(2) + 8*sin(theta_uk_wk(2))],'Color','y');
% 
% %               % wk_uk
% %               line([robots(wIdx).position(1), robots(wIdx).position(1) + robots(wIdx).distk*cos(robots(wIdx).anglek) + 3*cos(phi_wk_uk(1))], ...
% %                   [robots(wIdx).position(2), robots(wIdx).position(2) + robots(wIdx).distk*sin(robots(wIdx).anglek)+ 3*cos(phi_wk_uk(1))],'Color','k');
% %               line([robots(wIdx).position(1), robots(wIdx).position(1) + robots(wIdx).distk*cos(robots(wIdx).anglek) + 3*cos(phi_wk_uk(2))], ...
% %                   [robots(wIdx).position(2), robots(wIdx).position(2) + robots(wIdx).distk*sin(robots(wIdx).anglek)+ 3*cos(phi_wk_uk(2))],'Color','g');
% 
% axis([-10, 10, -10, 10]);
% grid on
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % PLOT ESTIMATED POSITION W NEIGHBOUR ROBOTS EACH REFERENCE FRAME
% %%%%%%%%%%%%%%%%%%%%
% subplot(2,number_of_robots,uIdx+number_of_robots)
% %                       figure(30+uIdx)
% 
% hold on;
% 
% 
% plot(0, 0, 'ro');
% % plot(robots(wIdx).pose_est_k(1,1),robots(wIdx).pose_est_k(1,2),'b*');
% % plot(robots(wIdx).pose_est_k(2,1),robots(wIdx).pose_est_k(2,2),'b*');
% r = 2;
% % Show Vector of Real Initial Orientation of Robots
% M = [r*cos(robots(wIdx).pos_hist(3,1)); r*sin(robots(wIdx).pos_hist(3,1))];
% plotv(M,'-')
% %                       if isempty(robots(wIdx).theta_wk_uk) ~= 1
% M = [r*cos(robots(wIdx).position(3)), robots(wIdx).distk*cos(robots(wIdx).theta_wk_uk(1)+robots(wIdx).epsi), robots(wIdx).distk*cos(robots(wIdx).theta_wk_uk(2)+robots(wIdx).epsi);...
%     r*sin(robots(wIdx).position(3)), robots(wIdx).distk*sin(robots(wIdx).theta_wk_uk(1)+robots(wIdx).epsi), robots(wIdx).distk*sin(robots(wIdx).theta_wk_uk(2)+robots(wIdx).epsi)];
% plotv(M,'-')
% %                       end
% 
% axis([-10, 10, -10, 10]);
% %               title ('Robot', wIdx);
% title(['Robot ',num2str(uIdx),' Vectors'])
% grid on;
%                       