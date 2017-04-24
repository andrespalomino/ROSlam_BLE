function graficas_3d (robot_wIdx, robot_uIdx, MT)

figure(10)
% PLot distances
subplot(2,2,1);
hold on

plot3(robot_wIdx(1), robot_wIdx(2), robot_wIdx(3), 'ro');
% plot3(robots(wIdx).pos_hist(1,:), robots(wIdx).pos_hist(2,:), robots(wIdx).pos_hist(3,:),'r');

plot3([robot_wIdx(1) robot_uIdx(1)], [robot_wIdx(2) robot_uIdx(2)], [robot_wIdx(3) robot_uIdx(3)],'c','LineWidth',1.5);

% Marco referencia
frame_mod(MT,'r',1);
title('Initial positions and distances');
axis([-10, 10, -10, 10, -10, 10]);
xlabel('x')
ylabel('y')
zlabel('z')
view(45,45)
grid on

subplot(2,2,2);
hold on
plot3([robot_wIdx(1) robot_uIdx(1)], [robot_wIdx(2) robot_uIdx(2)], [robot_wIdx(3) robot_uIdx(3)],'c','LineWidth',1.5);
plot3(robot_wIdx(1), robot_wIdx(2), robot_wIdx(3), 'ro');
% Marco referencia
frame_mod(MT,'r',1);
axis([-10, 10, -10, 10, -10, 10]);
xlabel('x')
ylabel('y')
zlabel('z')
view(0,90)
grid on

subplot(2,2,3);
hold on
plot3([robot_wIdx(1) robot_uIdx(1)], [robot_wIdx(2) robot_uIdx(2)], [robot_wIdx(3) robot_uIdx(3)],'c','LineWidth',1.5);
plot3(robot_wIdx(1), robot_wIdx(2), robot_wIdx(3), 'ro');
% Marco referencia
frame_mod(MT,'r',1);
axis([-10, 10, -10, 10, -10, 10]);
xlabel('x')
ylabel('y')
zlabel('z')
view(0,0)
grid on

subplot(2,2,4);
hold on
plot3([robot_wIdx(1) robot_uIdx(1)], [robot_wIdx(2) robot_uIdx(2)], [robot_wIdx(3) robot_uIdx(3)],'c','LineWidth',1.5)
plot3(robot_wIdx(1), robot_wIdx(2), robot_wIdx(3), 'ro');
% Marco referencia
frame_mod(MT,'r',1);
axis([-10, 10, -10, 10, -10, 10]);
xlabel('x')
ylabel('y')
zlabel('z')
view(90,0)
grid on


figure(20)
hold on
plot3([robot_wIdx(1) robot_uIdx(1)], [robot_wIdx(2) robot_uIdx(2)], [robot_wIdx(3) robot_uIdx(3)],'c','LineWidth',1.5);
plot3(robot_wIdx(1), robot_wIdx(2), robot_wIdx(3), 'ro');
plot3(robot_wIdx(1), robot_wIdx(2), robot_wIdx(3),'o','MarkerSize',10,...
            'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1]);    
% Marco referencia
% frame_mod(MT,'r',1);
title('Posiciones, orientaciones y grafico de comunicaciones Inicial ');
% title('Initial Position');
set(gca,'FontSize',11); set(gcf,'Color','White');
% title('Initial positions, orientation and communication graph ');

% axis([-10, 10, -10, 10, -10, 10]);
axis square
% axis equal
xlabel('x')
ylabel('y')
zlabel('z')
view(45,45)
% axis([-4, 2, -2, 4, 1, 5]);
grid on


% figure(30)
% hold on
fx = robot_uIdx(1)-robot_wIdx(1);
fy = robot_uIdx(2)-robot_wIdx(2);
fz = robot_uIdx(3)-robot_wIdx(3);

F = sqrt(fx^2 + fy^2 + fz^2);

ang_xy = (fx/F);
ang_yz = (fy/F);
ang_xz = (fz/F);

% 
% ang_xy = atan2(robot_wIdx(2)-robot_uIdx(2),robot_wIdx(1)-robot_uIdx(1))
% ang_yz = atan2(robot_wIdx(3)-robot_uIdx(3),robot_wIdx(2)-robot_uIdx(2))
% ang_xz = atan2(robot_wIdx(3)-robot_uIdx(3),robot_wIdx(1)-robot_uIdx(1))

% r = F;
% 
% plot3([robot_wIdx(1) r*(ang_xy)], [robot_wIdx(2) r*(ang_yz)], [robot_wIdx(3)  r*(ang_xz)],'g','LineWidth',1.5);
% 
% pause 
% plot3(robot_wIdx(1), robot_wIdx(2), robot_wIdx(3), 'ro');
% scatter3(robot_wIdx(1), robot_wIdx(2), robot_wIdx(3),'MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 .75 .75]);
% % Marco referencia
% frame_mod(MT,'r',1);
% title('Initial positions and distances');
% % axis([-10, 10, -10, 10, -10, 10]);
% axis equal
% xlabel('x')
% ylabel('y')
% zlabel('z')
% view(45,45)
% grid on

end