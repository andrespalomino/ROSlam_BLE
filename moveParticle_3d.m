 % Compute how the robot should move from "pos" given the requested movement
% and some Gaussian random noise using the motion model. This
% method is used to move the simulated robot as well as each of the
% hypothetical particles.

% rot [phi,theta,psi]
% [dist_travel, rot]
% vector (1x4)
function [newpos,Rnew] = moveParticle_3d(pos, R, movement, variance)
 
  travel_dist_x   = normrnd(movement(1), variance(1)* 0.25);
  travel_dist_y   = normrnd(movement(2), variance(1)* 0.25);
  travel_dist_z   = normrnd(movement(3), variance(1)* 0.25);
  phi           = normrnd(movement(4), variance(2)* 0.25);
  theta         = normrnd(movement(5), variance(2)* 0.25);
  psi           = normrnd(movement(6), variance(2)* 0.25);

  v_travel_dist = [travel_dist_x,travel_dist_y,travel_dist_z]';

  % Euler angles for rotatioin XYZ
%   R = R_euler(phi, theta, psi, 0);
  
  
  rot_X = [1 0 0;...
      0 cos(phi) -sin(phi);...
      0 sin(phi) cos(phi)];
  rot_Y = [cos(theta) 0 sin(theta);...
      0 1 0;...
      -sin(theta) 0  cos(theta)];
  rot_Z = [cos(psi) -sin(psi) 0;...
      sin(psi) cos(psi) 0;...
      0 0 1];
  
%   XYZ
   R_delta = rot_Z*rot_X*rot_Y;
    
%   delta = zeros(4,1);
%   delta(1,1) = travel_dist *cos(pos(4)+phi);
%   delta(2,1) = travel_dist *sin(pos(4)+phi);
%   delta(3,1) = travel_dist *sin(pos(4)+rotation2);
 
  Rnew = R*R_delta;

  newpos = pos + Rnew*v_travel_dist;
  
end
