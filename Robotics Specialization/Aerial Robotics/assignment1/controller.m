function [ u ] = pd_controller(~, s, s_des, params)
%PD_CONTROLLER  PD controller for the height
%
%   s: 2x1 vector containing the current state [z; v_z]
%   s_des: 2x1 vector containing desired state [z; v_z]
%   params: robot parameters


% FILL IN YOUR CODE HERE
%error
ep = s_des(1)-s(1);
ev = s_des(2)-s(2);
z_ddot_des = 0;
Kp = 120;
Kv = 20;

u = params.mass*(z_ddot_des + Kp*ep + Kv*ev + params.gravity);

end

