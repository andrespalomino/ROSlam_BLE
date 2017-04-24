%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Copyright 2010 Randolph Voorhies
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
function [z, angle_alfa, angle_beta,angle_gamma] = get_Distance_3d(pos, landmark_pos, observation_variance)
%   Given a landmark position and a robot position, this method will return a
%   "measurement" z that contains the distance and the angle to the landmark.
%   Gaussian random noise is added to both based on the variances given in the
%   diagonal of the observation_variance matrix. Note that this method is used
%   both to take a "real" measurement in the simulation, as well as to assess 
%   what kind of measurement each of our hypothetical particles would take.
%   This method also computes the Jacobian of the measurent function for use
%   in an extended Kalman filter.



  % Compute the distance from the current position to the landmark, and add
  % some Gaussian noise to make things interesting. Note that we are using
  % a smaller variance in this Gaussian distribution, as the algorithm seems
  % to work better when it underestimates the quality of the sensor. 
  vector_to_landmark = [landmark_pos(1) - pos(1); landmark_pos(2) - pos(2); landmark_pos(3) - pos(3)];
  landmark_distance = norm(vector_to_landmark);
  
  d_real = landmark_distance;

    sigma_meas = 0.1;
  % Compute the angle from the given pos to the landmark
%   angle_xy = atan2(vector_to_landmark(2), vector_to_landmark(1));
%   dist_xy= sqrt((landmark_pos(1) - pos(1))^2 + (landmark_pos(2) - pos(2))^2);   
%   angle_z = atan2(dist_xy,vector_to_landmark(3));
  
  
  %%% Angulos directores COSENOS Simulación 
  angle_alfa = acos(vector_to_landmark(1)/landmark_distance);
  angle_beta = acos(vector_to_landmark(2)/landmark_distance);
  angle_gamma = acos(vector_to_landmark(3)/landmark_distance);
  
  
  %%% Varianza de medida angulo y distancia
  angle_alfa = angle_alfa + normrnd(0, observation_variance(2));
  angle_beta = angle_beta + normrnd(0, observation_variance(2));
  angle_gamma = angle_gamma + normrnd(0, observation_variance(2));
  

%   landmark_distance = landmark_distance + normrnd(0, observation_variance(1));
  landmark_distance = landmark_distance + observation_variance(1)*randn;
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bluetooth Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RSSI dB en funcion de Distancia
Y_dB = -58.455 -10*(3.3687)*log10(d_real) + normrnd(0, observation_variance(1)); 

% Distance mts en funcion de RSSI dB
d_measured = 0.0186*exp(-0.068 * Y_dB); 

landmark_distance = d_measured;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% 
%   % Compute the Jacobian of this measurement function
%   q = landmark_distance^2.0;
%   H = [-(landmark_pos(1) - pos(1))/sqrt(q), -(landmark_pos(2) - pos(2))/sqrt(q), 0.0;
%         (landmark_pos(2) - pos(2))/q,       -(landmark_pos(1) - pos(1))/q,      -1.0;
%         0.0,                                0.0,                               1.0];
% 
%   z = [landmark_distance; 
%        landmark_angle;
%        vector_to_landmark];
   
   
  z = landmark_distance;
end