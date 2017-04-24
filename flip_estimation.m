%%%% Estimaciòn K - 1 

%%%%%%%%%%%%%%%%%%%%%%
          % EStimation before round K to check consistency and avoid flip
          % ambiguity
          %%%%%%%%%%%%%%%%%%%%%
          for wIdx = 1:number_of_robots
              
              %----------------------------------------
              % 3: broadcast puk?1|uk, ?uk?1|uk, ?uk?1
              %----------------------------------------
%               robots(wIdx).odometry =    get_Odometry_pos( ...
%                   (robots(wIdx).position - robots(wIdx).pos_hist(:,1)), (k_1(:,wIdx) - robots(wIdx).pos_hist(:,1)), odometry_noise);
%               
%               od_x = robots(wIdx).odometry(1)
%               od_y = robots(wIdx).odometry(2)
%               od_z = robots(wIdx).odometry(3)
              
              %----------------------------------------
              % 3: broadcast puk?1|uk, ?uk?1|uk, ?uk?1
              %----------------------------------------
              robots(wIdx).odometry =    get_Odometry_pos( ...
                  (robots(wIdx).position - robots(wIdx).pos_hist(:,1)), (k_1(:,wIdx) - robots(wIdx).pos_hist(:,1)), odometry_noise);
              
              od_x_n = robots(wIdx).odometry(1);
              od_y_n = robots(wIdx).odometry(2);
              od_z_n = robots(wIdx).odometry(3);
              
              
              
              % Get measurement of altitud
              dreal_z_k_1 = k_1(3,wIdx);
              dreal_z_k   = k_1(3,wIdx) + robots(wIdx).odometry(3);
              
              % Traveled Distance with odometry
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
              for uIdx = 1: number_of_robots
                  %               Take a real (noisy) measurement (DISTANCE)
                  %               from the robot to the landmark
                  %----------------------------------------
                  % 4: receive pwk?1|wk, ?wk?1|wk, ?uk?1 for w ? Nuk
                  %----------------------------------------
                  % if state = mobile then
                  if (robots(wIdx).state == 1) && (wIdx ~= uIdx)
                      
                      delta_z_k_1 = k_1(3,uIdx) - dreal_z_k_1 ;
                      delta_z_k   = k_1(3,uIdx) - dreal_z_k ;
                      
                      
                      
                      % Measurement Distance and angle K_1
                      real_landmark = k_1(:,uIdx);
                      %                       [z_real, G] = get_Distance(k_1(:,wIdx), real_landmark, [0;0]);
                      [z_real,a_x,a_y,a_z] = get_Distance_3d(k_1(:,wIdx), real_landmark, [0;0]);
                      k_1_sinruido = z_real;
                      
%                       [z_real,a_x,a_y,a_z] = get_Distance_3d(k_1(:,wIdx), real_landmark, measurement_variance);
                      robots(wIdx).distk_1  = z_real;
                      k_1_conruido = z_real;
                      %                       robots(wIdx).anglek_1 = z_real(2);
                      
                      real_landmark2 = robots(uIdx).position;
                      % Measurement Distance and angle K
                      %                       [z_real, G] = get_Distance(robots(wIdx).position, real_landmark2, [0;0]);
                      [z_real,a_x,a_y,a_z] = get_Distance_3d(robots(wIdx).position, real_landmark2, [0;0]);
                      k_sin_ruido = z_real;
                      
                      [z_real,a_x,a_y,a_z] = get_Distance_3d(robots(wIdx).position, real_landmark2, measurement_variance);
                      robots(wIdx).distk    = z_real;
                      k_con_ruido = z_real;
                      %                       aux(uIdx) = z_real(1);
                      %                       robots(wIdx).anglek   = z_real(2);
                       
                       
                      % Vector en x_y
                      v_xy_k   = sqrt(robots(wIdx).distk^2 - delta_z_k^2);
                      v_xy_k_1 = sqrt(robots(wIdx).distk_1^2 - delta_z_k_1^2);
                      
                      
                      
                      % 6: ?uk ? ? ˆwk|uk through Eq. (4-5)
                      % 7: ? ˆwk|uk ? use Eq. (6-7) ?w ? Nuk
                      % 8: use previous state resolve flip in ?uk
                      
                      %
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      % % RELATIVE POSITION AND ORIENTATION ESTIMATION
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      dk = v_xy_k;
                      %                       robots(wIdx).dk = dk;
                      dk_1 = v_xy_k_1;
                      %               dk_1 = robots(wIdx).distk_1;
                      %
                      luk = luk_xy;
                      arg_cseno = (luk^2 + dk^2 - dk_1^2)/(2*luk * dk);
                      
                      robots(wIdx).beta = acos((luk^2 + dk^2 - dk_1^2)/(2*luk * dk));
                      if (luk^2 + dk^2 - dk_1^2)/(2*luk * dk) < -1
                          robots(wIdx).beta = acos(-1);
                      elseif (luk^2 + dk^2 - dk_1^2)/(2*luk * dk) > 1
                          robots(wIdx).beta = acos(1);
                      end
                      beta = robots(wIdx).beta;
%                       robots(wIdx).beta = 0.5;
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
                      
%                       theta_flip(uIdx,:) = theta_wk_uk
%                       robots(wIdx).pos_flip(uIdx,:) = theta_wk_uk;
                      
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      % ESTIMATED RELATIVE POSITION
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      x_estimado = dk*cos(robots(wIdx).theta_wk_uk(1) + robots(wIdx).epsi);
                      y_estimado = dk*sin(robots(wIdx).theta_wk_uk(1) + robots(wIdx).epsi);
                      z_estimado = delta_z_k;
%                         pause
                      %                       pos_s1_x(uIdx) = x_estimado;
                      %
                      x_estimado_2 = dk*cos(robots(wIdx).theta_wk_uk(2)+robots(wIdx).epsi);
                      y_estimado_2 = dk*sin(robots(wIdx).theta_wk_uk(2)+robots(wIdx).epsi);
                      
                      v_estimado_flip =  [x_estimado,y_estimado,x_estimado_2,y_estimado_2];
                      robots(wIdx).pos_flip(uIdx,:) = [x_estimado,y_estimado,x_estimado_2,y_estimado_2];
                      robots(wIdx).theta_flip(uIdx,:) = theta_wk_uk;
                      
                                            
                      x_resolve = 0;
                      y_resolve = 0;
                      
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      % REAL RELATIVE POSITION
                      real_landmark = robots(wIdx).position;
                      real_landmark2 = robots(uIdx).position;
                      x_real = real_landmark2(1) - real_landmark(1);
                      y_real = real_landmark2(2) - real_landmark(2);
                      z_real = real_landmark2(3) - real_landmark(3);
                      %%%%%%%
                      % Distancia Medida
                      x_meas = real_landmark(1) + cos(a_x) * robots(wIdx).distk;
                      y_meas = real_landmark(2) + cos(a_y) * robots(wIdx).distk;
                      z_meas = real_landmark(3) + cos(a_z) * robots(wIdx).distk;
                      % Distancia Medida Local
                      x_meas_loc = cos(a_x) * robots(wIdx).distk;
                      y_meas_loc = cos(a_y) * robots(wIdx).distk;
                      z_meas_loc = cos(a_z) * robots(wIdx).distk;
                      
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      % GRAFICAS DISTANCIA Y VECTORES ESTIMACION
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       graficas_vec_3d;
%                       pause 
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                      
                      
                  end
              end

              if (robots(wIdx).state == 1)
%                   robots(wIdx).EMC_x = [robots(wIdx).EMC_x; EMC_x];
%                   EMC_x_acu = robots(wIdx).EMC_x
%                   robots(wIdx).EMC_y = [robots(wIdx).EMC_y; EMC_y];
%                   robots(wIdx).EMC_z = [robots(wIdx).EMC_z; EMC_z];
                  
%                   robots(wIdx).theta_flip = [robots(wIdx).theta_flip; theta_flip];
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