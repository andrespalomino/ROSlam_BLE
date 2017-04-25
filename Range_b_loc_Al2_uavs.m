
% break

% LOOP

v_sigma = [0 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8]
% for sigma = 5:length(v_sigma)
    

    


% for test = 1:5
    
% test = 1;
% break;
% ------------------------------
% ------------------------------
% TRAYECTORIA 4 ROBOTS FLOR
% ------------------------------
% ------------------------------

% clear everything
% clear 
close all
% clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
error_resolve = 0.4

% Max 1.5 BLE noise
% sigma_d = v_sigma(sigma)
sigma_d = 0.02
sigma_z = 0.02

measurement_variance_z = normrnd(0, sigma_z); 

odometry_noise = 0;   

% v_sigma = [0 0.001 0.002 0.005 0.01 0.015 0.02 0.03 0.05]

% Number of Robots
number_of_robots = 5;

% The number of timesteps for the simulation
timesteps = 80;

% The Gaussian variance of our sensor readings
measurement_variance = [sigma_d;
                        0;
                        0];   % Distance_x

% The maximum distance from which our sensor can sense a landmark or other
% robot
max_read_distance = 80;

color = [1 0 0;0 0 1; 0 1 1; [127 198 78]/255]; 
color_1 = [1 0 0;0 0 1;
    [65 176 135]/255;
    [153 0 153]/255;
    [255 105 0]/255;
    [140 205 96]/255;    
    [131 156 159]/255;
    [204 0 153]/255];
    
% color_2 = [[251 128 105]/255;[78 179 211]/255; [140 205 96]/255; [204 0 153]/255]; 

color_2 = [[251 128 105]/255;[78 179 211]/255; [140 205 96]/255; [204 0 153]/255; [0.3756 0.7007 0.3032]]; 
% color_1 = color_2;
color_1 = color_1
color_2 = color_1

% color = [1 0 0;0 0 1; 0.6135 0.7261 0.2074; 0.9386 0.0713 0.3635]; 

% round to take odometry and plot
inc = 4;
count_round = 4;
count_round_flip = 2; 
% count_round_flip_2 = count_round_flip + 2;

% % round to take odometry and plot
% inc = 5;
% count_round = 5;
% count_round_flip = 2; 

count_round_plot = count_round;

% Take odometry every 3 rounds 
flag_round = 0;

% Rounds to apply algorithm 
j = 1;

% break
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROBOT U
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The initial starting position of the robot
real_position = [0;      % x
                0;     % y
                 0];  % rotation
                                     
% The movement command given to he robot at each timestep                 
movement_command = [0;     % Distance
                    0];    % Rotation      
                
% The Gaussian variance of the movement commands
movement_variance = [0;   % Distance
                     0]; % Rotation
Q = [movement_variance(1), 0.0;
     0.0, movement_variance(2)];
                 
R = [measurement_variance(1), 0.0, 0.0;
     0.0, measurement_variance(2), 0.0;
     0.0, 0.0, measurement_variance(3)];
                
             
  robot_u.position = real_position;
  robot_u.pos_hist = [real_position];  
  robot_u.odometry = [];
  robot_u.odometry_hist = [];
  robot_u.pos_est = [];
  
% Create the robots and initialize them all to be in the same initial
% position. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROBOTS W
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RANDOM
% a = -10;
% b = 10;
% for i = 1:number_of_robots
%     rx = a + (b-a)*rand;
%     ry = a + (b-a)*rand;
%     rt = a + (b-a)*rand;
%     robots(i).position = [rx; ry; rt];
%     robots(i).pos_hist = robots(i).position;
% end


% FIXED
% The initial starting position of the robot
real_position_w1 = [0;      % x
                 0;     % y
                 3];     %z
                   % rotation
                   
% initial euler angles phi,theta,psi
real_angle_w1 = [0,0,0];
T = R_euler(real_angle_w1(1),real_angle_w1(2),real_angle_w1(3),0);
MT = [[T,real_position_w1];[0 0 0 1]];
robots(1).MT = MT;

% The initial starting position of the robot
real_position_w2 = [1;      % x
                0.5;     % y
                2.5];      % z
                  % rotation
                  
% % The initial starting position of the robot
% real_position_w2 = [-0.5;      % x
%                 0;     % y
%                 3];      % z
%                   % rotation
                  
% initial euler angles phi,theta,psi
real_angle_w2 = [0,0,0];
T = R_euler(real_angle_w2(1),real_angle_w2(2),real_angle_w2(3),0);
MT = [[T,real_position_w2];[0 0 0 1]];
robots(2).MT = MT;

% The initial starting position of the robot
real_position_w3 = [2;      % x
                 1;     % y
                 2];     % z
                   % rotation
                   
%                    % The initial starting position of the robot
% real_position_w3 = [0;      % x
%                  0.5;     % y
%                  3];     % z
%                    % rotation
                   
% initial euler angles phi,theta,psi
real_angle_w3 = [0,0,0];  
T = R_euler(real_angle_w3(1),real_angle_w3(2),real_angle_w3(3),0);
MT = [[T,real_position_w3];[0 0 0 1]];
robots(3).MT = MT;
             
% The initial starting position of the robot
real_position_w4 = [3;      % x
                1;     % y
                2];      %z
                  % rotation
                  
%                   % The initial starting position of the robot
% real_position_w4 = [-0.5;      % x
%                 0.5;     % y
%                 3];      %z
%                   % rotation
                  
% initial euler angles phi,theta,psi
real_angle_w4 = [0,0,0];   
T = R_euler(real_angle_w4(1),real_angle_w4(2),real_angle_w4(3),0);
MT = [[T,real_position_w4];[0 0 0 1]];
robots(4).MT = MT;

                 
movement_command2 = [0.1;   % Distance x
                     0.1;   % Distance y
                     0;   % Distance z
                     0.0;   % phi
                     0.0;   % theta
                     0];    % psi        
                 
              
                 
                 
%%%%% INITIAL POSE
% The initial starting position of the robot

for n=1:number_of_robots
    
%     Random pos 0 - 5
%     robots(n).position =  (6*rand(1,3)-3)';
%     Fixed pos
    robots(n).position =  [n,n,n]';
    
    robots(n).pos_hist = robots(n).position;
    robots(n).pos_hist_all = robots(n).position;
    
    real_angle_ = [0,0,0];   
    T = R_euler(real_angle_(1),real_angle_(2),real_angle_(3),0);
    MT = [[T,robots(n).position];[0 0 0 1]];
    robots(n).MT = MT;

%     color(n,:) = [rand rand rand];
%     color_2(n,:) = [rand rand rand];
%     cadena_ = ['real_position_w' , num2str(n)]
%     var = str2num(cadena_)
%     ['A_ERR_BLE', num2str(wIdx),'_',num2str(sigma_d),'.mat']
end
%%%--------------------------------


                
% robots(1).position = real_position_w1;
% robots(2).position = real_position_w2;
% robots(3).position = real_position_w3;
% robots(4).position = real_position_w4;
% robots(1).pos_hist = real_position_w1;
% robots(2).pos_hist = real_position_w2;
% robots(3).pos_hist = real_position_w3;
% robots(4).pos_hist = real_position_w4;
% 
% robots(1).pos_hist_all = real_position_w1;
% robots(2).pos_hist_all = real_position_w2;
% robots(3).pos_hist_all = real_position_w3;
% robots(4).pos_hist_all = real_position_w4;



% Initialize variables of each robot W
  for wIdx=1:number_of_robots
%     robots(wIdx).position = [0;0;0];  
%     robots(wIdx).pos_hist = []; 

    robots(wIdx).read_distance = [];
    robots(wIdx).read_ang_x = [];
    robots(wIdx).read_ang_y = [];
    robots(wIdx).read_ang_z = [];    
    robots(wIdx).read_angle  = [];
%     robots(wIdx).vector_x = [];
%     robots(wIdx).vector_y = [];
    robots(wIdx).odometry  = [];
    robots(wIdx).odometry_hist  = [];
    robots(wIdx).pose_est_double = [];
    robots(wIdx).pose_est = [];
    robots(wIdx).pose_real = [];
    robots(wIdx).pose_meas = [];
%     robots(wIdx).x_est = [];
    robots(wIdx).gamma = [];
    robots(wIdx).beta = [];           
    robots(wIdx).distk_1 = [];
    robots(wIdx).anglek_1 = [];
    robots(wIdx).distk = [];
    robots(wIdx).theta_wk_uk = [];
    robots(wIdx).epsi = [];                    
%     robots(wIdx).EMC_x = [];
%     robots(wIdx).EMC_y = [];
%     robots(wIdx).EMC_z = [];
%     robots(wIdx).EMC_v = [];
    robots(wIdx).pos_flip = zeros(number_of_robots,4);
    robots(wIdx).theta_flip = zeros(number_of_robots,2);
    if wIdx == 1
        robots(wIdx).state = 1;
    else 
        robots(wIdx).state = 0;
    end
  end

  
%   Adjacency Matrix
  AM = ones(number_of_robots);
  AM(logical(eye(size(AM)))) = 0;
  
  
  next = 1;
  
  
%   INITIAL FOR THETA UK_1_WK and PK_1 DISTANCE 

% figure(1)
% clf;
% hold on;
 
 
  for wIdx=1:number_of_robots
        
      real_landmark = robots(wIdx).position;
      MT = robots(wIdx).MT; 
      robots(wIdx).MT_1 = MT;
      for uIdx = 1: number_of_robots
              
              real_landmark2 = robots(uIdx).position;
              % Take a real (noisy) measurement from the robot to the landmark
              [z_real,a_xy,a_z] = get_Distance_3d(real_landmark, real_landmark2, [0;0]);
%               robots(wIdx).theta_uk_1_wk(uIdx) = z_real(2);
              
              
              % If the robot is close enough, then we can spot it
              if(z_real(1) < max_read_distance)                  
                  aux(uIdx) = z_real(1);
%                   aux_ang(uIdx) = z_real(2);
              else
                  aux(uIdx) = 0;
%                   aux_ang(uIdx) = 0;
              end
              
%               Plot graficas varias vistas ejes
%               graficas_3d (real_landmark, real_landmark2, MT);
              
            
      end
      aux_pk_1(:,wIdx) = robots(wIdx).position; 
%       robots(wIdx).distk_1 = aux;      
  end
  
% Previous position  
k_1 = aux_pk_1;

% pause

% break
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SIMULATION
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for round = 0:timesteps
%       if round <= timesteps/4
%           movement_command2 = [0.1;   % Distance
%               -0.04];    % Rotation
%           %             % Different control actions over time
%       elseif round >= timesteps/4 && round < timesteps/2
%           movement_command2 = [0.05;     % Distance
%               0];    % Rotation
%       elseif round >= timesteps/2 && round < (timesteps/4)*3
%           movement_command2 = [0.06;     % Distance
%               0.02];    % Rotation
%       end
      
      
%       round
%       count_round
      
      %   For every neighboring robot w or robots
      % para cada robot w en el alcance de u Get Distance
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %       Get Distance
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      for wIdx = 1: number_of_robots
          
          real_landmark = robots(wIdx).position;
          %           MT = robots(wIdx).MT;          
              for uIdx = 1: number_of_robots
                  
                  real_landmark2 = robots(uIdx).position;
                  % Take a real (noisy) measurement from the robot to the landmark
                  %               [z_real, G] = get_Distance(real_landmark, real_landmark2, [0;0]);
                  [z_real,a_x,a_y,a_z] = get_Distance_3d(real_landmark, real_landmark2, measurement_variance);
                  
                  % If the robot is close enough, then we can spot it
                  %               if(z_real(1) < max_read_distance)
                  aux(uIdx) = z_real;
                  %                   aux_ang(uIdx) = z_real(2);
                  %               else
                  %                   aux(uIdx) = 0;
                  aux_ang_x(uIdx) = a_x;
                  aux_ang_y(uIdx) = a_y;
                  aux_ang_z(uIdx) = a_z;
                  %               end
                  
              end          
          
          robots(wIdx).read_distance = [robots(wIdx).read_distance, aux'];
          robots(wIdx).read_ang_x = [robots(wIdx).read_ang_x, aux_ang_x'];
          robots(wIdx).read_ang_y = [robots(wIdx).read_ang_y, aux_ang_y'];
          robots(wIdx).read_ang_z = [robots(wIdx).read_ang_z, aux_ang_z'];
          dist = robots(wIdx).read_distance;
          %           pause
          %           robots(wIdx).read_angle = [robots(wIdx).read_angle,aux_ang'];
          
      end %for number_robots
%       
      
  
        
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %       Get Odometry WHEN round = count_round
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %         % Get Odometry real position and past position
      %         %     Robots are capable of using odometry to estimate their pose change between
      %         % rounds in their own local coordinate system.
      %         %     It is assumed
      %         % that odometry estimates are reliable over intervals of two or three rounds
      %         % (i.e. i >= j ?3), but suffer from drift over longer time intervals.
      %       %% Odometry of u
      %       robot_u.pos_hist
      %       robots(wIdx).pos_hist
      %       robots(wIdx).read_distance
      % theta_uk_wk_1
      %       round
      %       pause
      
%       if round == 1
%           
%           %           robot_u.odometry = 0;
%           %           robot_u.odometry_hist = [0;0;0];
%           %           %               robot_u.odometry2 = 0;
          %
                    
      if round == count_round-count_round_flip
          
          % Execute program to estimate 2 solutions beforehand to later
          % compare for consistency
          
          flip_estimation;
         

          
      elseif round == count_round

           % Get Odometry real position and past position ROBOT U
           
          
%           robot_u.odometry = get_Odometry(robot_u.position, robot_u.pos_hist(:,round-(inc-1)));
%           robot_u.odometry_hist = [robot_u.odometry_hist,robot_u.odometry];
%           
%           robot_u.pos_est = [robot_u.pos_est,robot_u.position];
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K Estimation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          


          % Get Odometry real position and past position ROBOT W
          for wIdx = 1:number_of_robots
              
              %----------------------------------------
              % 3: broadcast puk?1|uk, ?uk?1|uk, ?uk?1
              %----------------------------------------
              robots(wIdx).odometry =    get_Odometry_pos( ...
                  (robots(wIdx).position - robots(wIdx).pos_hist(:,1)), (k_1(:,wIdx) - robots(wIdx).pos_hist(:,1)), odometry_noise);              
              
              % Get measurement of altitud
              dreal_z_k_1 = k_1(3,wIdx);
              dreal_z_k   = k_1(3,wIdx) + robots(wIdx).odometry(3);  
              
              measurement_variance_z = normrnd(0, sigma_z); 

              % traveled Distance with odometry
              % X-Y-Z 
              luk = sqrt(robots(wIdx).odometry(1)^2 + robots(wIdx).odometry(2)^2 + robots(wIdx).odometry(3)^2);
              % X-Y
              luk_xy = sqrt(robots(wIdx).odometry(1)^2 + robots(wIdx).odometry(2)^2);
              % X-Z
              luk_xz = sqrt(robots(wIdx).odometry(1)^2 + robots(wIdx).odometry(3)^2);
              % Y-Z
              luk_yz = sqrt(robots(wIdx).odometry(2)^2 + robots(wIdx).odometry(3)^2);  
              
              
              robots(wIdx).epsi = atan2(robots(wIdx).odometry(2),robots(wIdx).odometry(1));
              
%               robots(wIdx).odometry_hist  = [robots(wIdx).odometry_hist, [robots(wIdx).odometry;robots(wIdx).epsi]];
              
              % inicialize var aux to calculate error
              EMC_x = 0;
              EMC_y = 0;
              EMC_z = 0;
              theta_flip = [];
              pos_real = [];
              pos_est = [];
              pos_meas = [];
              pos_est_double = [];
              for uIdx = 1: number_of_robots
                  %               Take a real (noisy) measurement (DISTANCE)
                  %               from the robot to the landmark
                  %----------------------------------------
                  % 4: receive pwk?1|wk, ?wk?1|wk, ?uk?1 for w ? Nuk
                  %----------------------------------------
                  % if state = mobile then
                  if (robots(wIdx).state == 1) && (wIdx ~= uIdx)
                      
                      % Get delta_z of altitud + noise                                            
                       
                      delta_z_k_1 = k_1(3,uIdx) - dreal_z_k_1 + measurement_variance_z;
                      
                      delta_z_k   = k_1(3,uIdx) - dreal_z_k + measurement_variance_z;
                      
                      
                      
                      % Measurement Distance and angle K_1
                      real_landmark = k_1(:,uIdx);
%                       [z_real, G] = get_Distance(k_1(:,wIdx), real_landmark, [0;0]);
                      [z_real] = get_Distance_3d(k_1(:,wIdx), real_landmark, measurement_variance);
                      robots(wIdx).distk_1  = z_real;
%                       robots(wIdx).anglek_1 = z_real(2);
                      
                      real_landmark2 = robots(uIdx).position;
                      % Measurement Distance and angle K
%                       [z_real, G] = get_Distance(robots(wIdx).position, real_landmark2, [0;0]);
                      [z_real,a_x,a_y,a_z] = get_Distance_3d(robots(wIdx).position, real_landmark2, [0,0]);
                      sin_varianza = z_real
                      [z_real,a_x,a_y,a_z] = get_Distance_3d(robots(wIdx).position, real_landmark2, measurement_variance);
                      
                      robots(wIdx).distk    = z_real;
                      con_varianza = z_real;
                      
%                       pause
%                       aux(uIdx) = z_real(1);
%                       robots(wIdx).anglek   = z_real(2);
                        
                      % Vector en x_y 
                      v_xy_k   = sqrt(robots(wIdx).distk^2 - delta_z_k^2);
                      v_xy_k_1 = sqrt(robots(wIdx).distk_1^2 - delta_z_k_1^2);
                      
%                       aux(uIdx) = z_real(1);
%                       aux_ang_x(uIdx) = a_x;
%                       aux_ang_y(uIdx) = a_y;
%                       aux_ang_z(uIdx) = a_z;
                      
                      
%                       aux(uIdx) = z_real(1);
%                       aux_ang(uIdx) = a_xy;
%                       aux_ang_z(uIdx) = a_z;

%%%%%%            


                      % 6: ?uk ? ? ˆwk|uk through Eq. (4-5)
                      % 7: ? ˆwk|uk ? use Eq. (6-7) ?w ? Nuk
                      % 8: use previous state resolve flip in ?uk
                      
                      %
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      % % RELATIVE POSITION AND ORIENTATION ESTIMATION 
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      dk = v_xy_k;
%                       robots(wIdx).dk = dk;
                      dk_1 = v_xy_k_1;
                      %               dk_1 = robots(wIdx).distk_1;
                      %
                      luk = luk_xy;
                      
                      robots(wIdx).beta = acos((luk^2 + dk^2 - dk_1^2)/(2*luk * dk));
                      
                      if (luk^2 + dk^2 - dk_1^2)/(2*luk * dk) < -1
                          robots(wIdx).beta = acos(-1);
                      elseif (luk^2 + dk^2 - dk_1^2)/(2*luk * dk) > 1
                          robots(wIdx).beta = acos(1);
                      end
                      
                      
                      robots(wIdx).gamma = acos((dk^2 + dk_1^2 - luk^2)/(2*dk* dk_1));
                      %
                      % angle(Puk_1|uk)
                      %               alfa_uk = robots(wIdx).pos_hist(3,round-(inc-1));
                      %               alfa_uk = -robots(wIdx).odometry(3)+pi;
                      %               alfa_uk = robots(wIdx).epsi+pi;
                      %
                      alfa_uk = 0+pi;
                      % 2 posible solutions
                      theta_wk_uk = [alfa_uk + robots(wIdx).beta, alfa_uk - robots(wIdx).beta];
                      robots(wIdx).theta_wk_uk = theta_wk_uk;
                      
%                       
                      
%                       
%                       theta_flip = [theta_flip, theta_wk_uk'];
                      
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      % ESTIMATED RELATIVE POSITION
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      x_estimado = dk*cos(robots(wIdx).theta_wk_uk(1) + robots(wIdx).epsi);
                      y_estimado = dk*sin(robots(wIdx).theta_wk_uk(1) + robots(wIdx).epsi);
                      z_estimado = delta_z_k;
                      
%                       vxy_estimado = sqrt(x_estimado^2 + y_estimado^2);
                      
%                       pos_s1_x(uIdx) = x_estimado;
%                       
                      x_estimado_2 = dk*cos(robots(wIdx).theta_wk_uk(2) + robots(wIdx).epsi);
                      y_estimado_2 = dk*sin(robots(wIdx).theta_wk_uk(2) + robots(wIdx).epsi);
                      
%                       vxy_estimado_2 = sqrt(x_estimado_2^2 + y_estimado_2^2);
                      

                      %%%%%%%%%%%%%%%%%%%%%%%
                      pos_ahora = [x_estimado,y_estimado,x_estimado_2,y_estimado_2]
                      pos_ant = robots(wIdx).pos_flip(uIdx,:);
                      
                      % pos_dif should be near to previous pos
                      pos_dif = [robots(wIdx).odometry(1) + x_estimado, robots(wIdx).odometry(2) + y_estimado, robots(wIdx).odometry(1) + x_estimado_2,...
                          robots(wIdx).odometry(2) + y_estimado_2];
                      
                      
                      
%                       pos_flip_error =  sqrt ((pos_ant - pos_dif).^2)
                      pos_flip_error =  abs (pos_ant - pos_dif);
                      pos_ant = [pos_ant(3:4),pos_ant(1:2)];
%                       pos_flip_error_2 =  sqrt ((pos_ant - pos_dif).^2)
                      pos_flip_error_2 =  abs (pos_ant - pos_dif);
                      
                      theta_antes = robots(wIdx).theta_flip(uIdx,:);
                      theta_ahora = theta_wk_uk;
                      %%%%%%%%%%%%%%%%%%%%%%%
                      % Resolve flip Ambiguity
                      
                      m_flip_error1 = sqrt(pos_flip_error(1)^2 + pos_flip_error(2)^2) ;                     
                      m_flip_error2 = sqrt(pos_flip_error(3)^2 + pos_flip_error(4)^2);
                      
                      m_flip_error_21 = sqrt(pos_flip_error_2(1)^2 + pos_flip_error_2(2)^2);                  
                      m_flip_error_22 = sqrt(pos_flip_error_2(3)^2 + pos_flip_error_2(4)^2);
%                       
%                       if (pos_flip_error(1) <= error_resolve) || (pos_flip_error(2) <= error_resolve)
%                           pos_resolve = [x_estimado, y_estimado];
%                       elseif(pos_flip_error_2(1) <= error_resolve) && (pos_flip_error_2(2) <= error_resolve)
%                           pos_resolve = [x_estimado, y_estimado];
%                       else
%                           pos_resolve = [x_estimado_2, y_estimado_2];
%                       end

                      
                      if ((m_flip_error1 <= m_flip_error2) && (m_flip_error1 <= error_resolve))
                          
                              pos_resolve = [x_estimado, y_estimado];
                          
%                               pos_resolve = [x_estimado_2, y_estimado_2];
                          
                      elseif((m_flip_error_21 <= m_flip_error_22) && (m_flip_error_21 <= error_resolve))
%                           if 
                              pos_resolve = [x_estimado, y_estimado];
%                           else
%                               pos_resolve = [x_estimado_2, y_estimado_2];
%                           end
                      else
                          pos_resolve = [x_estimado_2, y_estimado_2];
                      end                      
                      
                      
                      % Solución Doble 
%                       pos_est_double = pos_ahora;
                      x_resolve_double = [pos_ahora(1),pos_ahora(3)];
                      y_resolve_double = [pos_ahora(2),pos_ahora(4)];
                      
                      % una Solucion escogida

                      x_resolve = pos_resolve(1);
                      y_resolve = pos_resolve(2);
                      z_resolve = z_estimado;
%                       z_resolve_double = z_estimado;
                      
                      
                      % Pos estimada GLOBAL pos_resolve + robot_pos 
                      pos_resolve
                      pos_est = [pos_est; (robots(wIdx).position' + [x_resolve, y_resolve, z_resolve])]
                      pos_est_double = [pos_est_double; (robots(wIdx).position' + ...
                          [x_resolve_double(1), y_resolve_double(1), z_resolve]),...
                          (robots(wIdx).position' + [x_resolve_double(2), y_resolve_double(2), z_resolve])]
                      
%                       pause

                      
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      % REAL RELATIVE POSITION
                      real_landmark = robots(wIdx).position;
                      real_landmark2 = robots(uIdx).position;
                      x_real = real_landmark2(1) - real_landmark(1);
                      y_real = real_landmark2(2) - real_landmark(2);
                      z_real = real_landmark2(3) - real_landmark(3);
                      %%%%%%%
                      % Distancia Medida Global
                      x_meas = real_landmark(1) + cos(a_x) * robots(wIdx).distk;
                      y_meas = real_landmark(2) + cos(a_y) * robots(wIdx).distk;
                      z_meas = real_landmark(3) + cos(a_z) * robots(wIdx).distk;
                      % Distancia Medida Local
                      x_meas_loc = cos(a_x) * robots(wIdx).distk;
                      y_meas_loc = cos(a_y) * robots(wIdx).distk;
                      z_meas_loc = cos(a_z) * robots(wIdx).distk;
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      
                      pos_real = [pos_real; (robots(wIdx).position' + [x_real, y_real, z_real])];
                      pos_meas = [pos_meas; [x_meas,y_meas,z_meas]];
                      
                      
                      % ERROR
                      error_absoluto_medida_x = abs((x_real - x_meas_loc)/x_real)*100;
                      error_absoluto_medida_y = abs((y_real - y_meas_loc)/y_real)*100;
                      error_absoluto_medida_z = abs((z_real - z_meas_loc)/z_real)*100;
                      
                      error_absoluto_estimacion_x = abs((x_meas_loc - x_resolve)/x_meas_loc)*100;
                      error_absoluto_estimacion_y = abs((y_meas_loc - y_resolve)/y_meas_loc)*100;
                      error_absoluto_estimacion_z = abs((z_meas_loc - z_resolve)/z_meas_loc)*100;
                      
                      arg_coseno = (luk^2 + dk^2 - dk_1^2)/(2*luk * dk);
                      
%                       pos_glo_real = robots(wIdx).position';
                      
%                       ECM_d_estimacion_x = [sqrt( mean( ( x_estimado - x_real)^2  ) ) , sqrt( mean( ( x_estimado_2 - x_real)^2  ) )];
%                       ECM_d_estimacion_y = [sqrt( mean( ( y_estimado - y_real)^2  ) ) , sqrt( mean( ( y_estimado_2 - y_real)^2  ) )];
%                       ECM_d_estimacion_z = [sqrt( mean( ( z_estimado - z_real)^2  ) ) , sqrt( mean( ( z_estimado - z_real)^2  ) )];
                      

    

%                       EMC_x(uIdx) = x_resolve - x_real;                
%                       EMC_y(uIdx) = y_resolve - y_real;
%                       EMC_z(uIdx) = z_estimado - z_real;
                      
                      %
                      %
                      %               % 2 posible solutions
                      %               %               theta_uk_wk = [theta_uk_wk_1 + robots(wIdx).gamma; theta_uk_wk_1 - robots(wIdx).gamma];
                      %
                      %               theta_uk_wk = [robots(wIdx).theta_uk_1_wk + robots(wIdx).gamma, robots(wIdx).theta_uk_1_wk - robots(wIdx).gamma];
                      %               robots(wIdx).theta_uk_wk = theta_uk_wk;
                      %               %               robots(wIdx).theta_uk_wk = theta_uk_wk;
                      %               %               t_uk_1_wk = robots(wIdx).theta_uk_1_wk*180/pi
                      %               %               robots(wIdx).anglek_1*180/pi
                      %               %               pause

                      %               % 8 posible solutions
                      %               phi_wk_uk = [theta_wk_uk(1) - theta_uk_wk(1) + pi; theta_wk_uk(1) - theta_uk_wk(2) + pi];
                      %               robots(wIdx).phi_wk_uk = phi_wk_uk;
                      %               p_uk_wk = phi_wk_uk*180/pi;
                      %
                      %               %               theta_uk_wk_1 = theta_uk_wk;
                      %               %               pos = dk * [cos(theta_wk_uk(2)) ; sin(theta_wk_uk(2))];
                      %               %               pose = [pos ; phi_uk_wk];
                      %               %
                      %               %               robots(wIdx).pose_est = [robots(wIdx).pose_est , pose];

                      %               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      % %               t_uk_wk = theta_uk_wk*180/pi
                      %
                      
                      %                   else
                      %                       theta_wk_uk =
                      
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      % GRAFICAS DISTANCIA Y VECTORES ESTIMACION
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       graficas_vec_3d;
%                       pause
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                     
%                   else
%                       
%                       for wwIdx = 1:number_of_robots
%                           for uIdx = 1:number_of_robots
%                               robots(wwIdx).pk_1(:,uIdx) = robots(uIdx).position;
%                               pk_1 = robots(wwIdx).pk_1(:,uIdx)
%                           end
%                       end
                      
                  end 
              end
              
              
              if (robots(wIdx).state == 1)
%                   robots(wIdx).EMC_x = [robots(wIdx).EMC_x; EMC_x];
%                   EMC_x_acu = robots(wIdx).EMC_x;
%                   robots(wIdx).EMC_y = [robots(wIdx).EMC_y; EMC_y];
%                   robots(wIdx).EMC_z = [robots(wIdx).EMC_z; EMC_z];
                  robots(wIdx).pose_est = [robots(wIdx).pose_est,pos_est];
                  robots(wIdx).pose_real = [robots(wIdx).pose_real,pos_real];
                  robots(wIdx).pose_meas = [robots(wIdx).pose_meas,pos_meas];
                  
                  robots(wIdx).pose_est_double = [robots(wIdx).pose_est_double,pos_est_double];
%                   robots(wIdx).read_distance = [robots(wIdx).read_distance, aux'];
%                   robots(wIdx).read_ang_x = [robots(wIdx).read_ang_x, aux_ang_x'];
%                   robots(wIdx).read_ang_y = [robots(wIdx).read_ang_y, aux_ang_y'];
%                   robots(wIdx).read_ang_z = [robots(wIdx).read_ang_z, aux_ang_z'];
                  
%                   robots(wIdx).read_distance = [robots(wIdx).read_distance, aux'];
%                   robots(wIdx).read_ang_xy = [robots(wIdx).read_ang_xy, aux_ang'];
%                   robots(wIdx).read_ang_z = [robots(wIdx).read_ang_z, aux_ang_z'];
%                   dist = robots(wIdx).read_distance;
              end
              
              % Actualizar Pk_1 y theta_uk_1_wk    %% ESTIMATION 
              k_1(:,wIdx) = robots(wIdx).position;
              robots(wIdx).MT_1 = robots(wIdx).MT;
%               robots(wIdx).distk_1 = aux;
%               robots(1).pk_1(:,1) = robots(wIdx).position;
%               robots(2).pk_1(:,1) = robots(wIdx).position;
%               robots(3).pk_1(:,1) = robots(wIdx).position;
              % robots(wIdx).theta_uk_1_wk = robots(wIdx).anglek;
                            
          end

          % MOSTRAR VALORES PROGRAMA
          %------------------------------------------------------
          %           u_pos = robot_u.position
          %           w_pos = robots(wIdx).position
          %           u_od = robot_u.odometry
          %           w_od = robots_odometry
          %           d_j = c
          %           d_k = d
          %           x_est =[x_est;x]
          %           posss = robots(wIdx).pose_est
          %           angle = robots(wIdx).read_angle(round)
          %                       dk
          %                       dk_1
          %             luk
          %           gamma = robots(wIdx).gamma*180/pi
          %           beta = robots(wIdx).beta*180/pi
          
          %
          %             theta_uk_wk_1
          %             phi_uk_wk
          %           alfa = alfa_uk*180/pi
          %           theta_wk_uk*180/pi
          %           epsi = robots(wIdx).epsi*180/pi
          %           t_uk_wk = theta_uk_wk*180/pi
          %             theta_uk_wk*180/pi
          %             phi_uk_wk
          %             robots(wIdx).pose_est
          %             alfa_uk3
          %------------------------------------------------------
          
          
          
          
          %------------------------------------------------------
          % state ? motion-scheduler
          %------------------------------------------------------
          for wIdx = 1:number_of_robots
              robots(wIdx).state = 0;
          end
          
          next = next +1;
          if next == number_of_robots + 1;
              next = 1;
          end
          robots(next).state = 1;
          

          
          
          % aumentar count_round para comparar
          count_round = count_round + inc
          
%           pause
          
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % PLOTTING
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      figure(50)
      clf;
      hold on;
     
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % PLOT REAL ROBOT
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % History and circle U
%       plot(robot_u.pos_hist(1,:), robot_u.pos_hist(2,:), 'b');
%       plot(robot_u.position(1), robot_u.position(2), 'bo');
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % PLOT ROBOTS & DISTANCE
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % Show the sensor measurement as an arrow
      %       last = size(robot_u.pos_est);
      %       last_ph =size(robot_u.pos_hist(1,:));
               
      for wIdx=1: number_of_robots
         real_landmark = robots(wIdx).position;    
         MT = robots(wIdx).MT;
         
         for uIdx = 1: number_of_robots
             real_landmark2 = robots(uIdx).position;
             % History and circle W
             %              graficas_3d (real_landmark, real_landmark2, MT)
             
             %              plot3(robots(wIdx).position(1), robots(wIdx).position(2), robots(wIdx).position(3), 'ro');
             % BLUE
%              wIdx = 1
%              plot3(robots(wIdx).position(1), robots(wIdx).position(2), robots(wIdx).position(3),'o','MarkerSize',8,...
%                  'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1]);
%              plot3(0,0, robots(wIdx).position(3),'o','MarkerSize',8,...
%                  'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1]);
             
%              if robots(wIdx).state == 1
%                  % RED
%                  plot3(robots(wIdx).position(1), robots(wIdx).position(2), robots(wIdx).position(3),'o','MarkerSize',8,...
%                      'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0]);
% %              end
%             if wIdx == 1
%                 plot3(robots(wIdx).position(1), robots(wIdx).position(2), robots(wIdx).position(3),'o','MarkerSize',8,...
%                     'MarkerEdgeColor','k','MarkerFaceColor',color(wIdx,:));
%                 plot3(robots(wIdx).pos_hist(1,:), robots(wIdx).pos_hist(2,:), robots(wIdx).pos_hist(3,:),'r','LineWidth',1);
%                 
%             elseif wIdx == 2
%                 plot3(robots(wIdx).position(1), robots(wIdx).position(2), robots(wIdx).position(3),'o','MarkerSize',8,...
%                     'MarkerEdgeColor','k','MarkerFaceColor',color(wIdx,:));
%                 plot3(robots(wIdx).pos_hist(1,:), robots(wIdx).pos_hist(2,:), robots(wIdx).pos_hist(3,:),'b','LineWidth',1);
%                 
%             elseif wIdx == 3
%                 plot3(robots(wIdx).position(1), robots(wIdx).position(2), robots(wIdx).position(3),'o','MarkerSize',8,...
%                     'MarkerEdgeColor','k','MarkerFaceColor',color(wIdx,:));
%                 plot3(robots(wIdx).pos_hist(1,:), robots(wIdx).pos_hist(2,:), robots(wIdx).pos_hist(3,:),'g','LineWidth',1);
%                 
%             elseif wIdx == 4
                plot3(robots(wIdx).position(1), robots(wIdx).position(2), robots(wIdx).position(3),'o','MarkerSize',8,...
                    'MarkerEdgeColor','k','MarkerFaceColor',color_1(wIdx,:));
                plot3(robots(wIdx).pos_hist(1,:), robots(wIdx).pos_hist(2,:), robots(wIdx).pos_hist(3,:),'Color',color_1(wIdx,:),'LineWidth',1);
%             end
            
%              plot3([real_landmark(1) real_landmark2(1)], [real_landmark(2) real_landmark2(2)], [real_landmark(3) real_landmark2(3)],'c','LineWidth',1.5)
%              frame_mod(MT,'r',1);
             
             %           line([robot_u.position(1)+robots(wIdx).read_distance(round)*cos(robots(wIdx).read_angle(round)), robot_u.position(1)+robots(wIdx).read_distance(round)*cos(robots(wIdx).read_angle(round))+2*cos(robots(wIdx).position(3))], ...
             %                [robot_u.position(2)+robots(wIdx).read_distance(round)*sin(robots(wIdx).read_angle(round)), robot_u.position(2)+robots(wIdx).read_distance(round)*sin(robots(wIdx).read_angle(round))+2*sin(robots(wIdx).position(3))],'Color','b');
             %
%              title('Movement according to motion-controller');
%              title('Movimiento Acción de control');
             set(gca,'FontSize',11); set(gcf,'Color','White');
             %           axis([-10, 10, -10, 10, -10, 10]);
%              axis([-10, 8, -1, 9, 0, 7]);
             xlabel('x(m)')
             ylabel('y(m)')
             zlabel('z(m)')
             %           rotate3d
%              axis square
             
%              axis([-3, 3, -2, 4, 1, 5]);
             axis equal
             view(45,45)
             grid on
             
             
         end
                   
      end
      
%    pause        




      
      
      %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %     % Move the actual robot and other robots
      %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       robot_u.position = moveParticle(robot_u.position, movement_command, movement_variance);
%       robot_u.pos_hist = [robot_u.pos_hist, robot_u.position];
      

% k_1(:,wIdx) = robots(wIdx).position;

      % Move the W robots
      for wIdx = 1:number_of_robots    
          % if state = mobile
          if robots(wIdx).state == 1
              % Posicion x,y, rotacion de los robots W
              dtx = 0;
              dty = 0.15;
              dtz = 0.15;
%               dtz = 0;
              psi = 0.06;
%               psi = 0;
              
%               dtx = rand * 0.2;
%               dty = rand * 0.2;
%               psi = rand * 0.2;
              
              if round < 20
                  sent = -1;
              elseif round >20 && round < 40
                  sent = 1;
              elseif round > 40 && round < 80
                  sent = -1;
              else
                  sent = 1;
              end
              
%               if round == 12
%                   results
%               elseif round == 12*2
%                   results
%               elseif round == 12*3
%                   results    
%               end
%               sent = 1;

%               Different Commands Movements
%               switch wIdx
%                   case 1
%                       movement_command2 = [dtx; dty; dtz; 0; 0; psi];
% %                       movement_command2 = [dtx; dty; 0; 0; 0; psi];
% %                       movement_command2 = [dtx; -dty; 0; 0; 0; sent*rand ];
%                   case 2
%                       movement_command2 = [dtx; dty; dtz; 0; 0; psi];
% %                       movement_command2 = [dtx; dty; 0; 0; 0; psi];
% %                       movement_command2 = [-dtx; -dty; 0; 0; 0; -sent*rand];
%                   case 3
%                       movement_command2 = [dtx; dty; dtz; 0; 0; psi];
% %                       movement_command2 = [dtx; dty; 0; 0; 0; psi];
% %                       movement_command2 = [dtx; dty; 0; 0; 0; -sent*rand];
%                   case 4
%                       movement_command2 = [dtx; dty; dtz; 0; 0; psi];
% %                       movement_command2 = [dtx; dty; 0; 0; 0; psi];                      
% %                       movement_command2 = [-dtx; dty; 0; 0; 0;  sent*rand];
%               end
              
              movement_command2 = [dtx; dty; dtz; 0; 0; psi];
              
              [new_pos, new_rot] = moveParticle_3d( ...
                  robots(wIdx).position, robots(wIdx).MT(1:3,1:3), movement_command2, movement_variance);
              robots(wIdx).position = new_pos;
              robots(wIdx).MT(1:3,1:3) = new_rot;
              robots(wIdx).MT(1:3,4) = new_pos;
              % Rotacion de las particulas
              %         robots(wIdx).position(3) = robot_u.position(3);
              robots(wIdx).pos_hist = [robots(wIdx).pos_hist ,robots(wIdx).position];                        
          end
          robots(wIdx).pos_hist_all = [robots(wIdx).pos_hist_all ,robots(wIdx).position];  
      end
      
      
      
      pause(0.05)

  end
  
%%
%%%% RUN RESULTS & Save
axis square
% results;
results_v20;
  
% end

% end




%% 
%%%% SHOW SAVED RESULTS 
% show_saved_data;











