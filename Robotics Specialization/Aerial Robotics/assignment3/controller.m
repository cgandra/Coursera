function [F, M] = controller(t, state, des_state, params)
%CONTROLLER  Controller for the quadrotor
%
%   state: The current state of the robot with the following fields:
%   state.pos = [x; y; z], state.vel = [x_dot; y_dot; z_dot],
%   state.rot = [phi; theta; psi], state.omega = [p; q; r]
%
%   des_state: The desired states are:
%   des_state.pos = [x; y; z], des_state.vel = [x_dot; y_dot; z_dot],
%   des_state.acc = [x_ddot; y_ddot; z_ddot], des_state.yaw,
%   des_state.yawdot
%
%   params: robot parameters

%   Using these current and desired states, you have to compute the desired
%   controls


% =================== Your code goes here ===================

% Thrust
F = 0;

Kp = [100;100;800];
Kv = [40;40;20];

ep = des_state.pos - state.pos;
ev = des_state.vel - state.vel;

cmd_acc = des_state.acc + Kv.*ev + Kp.*ep;
F = params.mass*(params.gravity + cmd_acc(3));
if F<params.minF
    F=params.minF;
end
if F>params.maxF
    F=params.maxF;
end

% Moment
M = zeros(3,1);

Kp_ang = [160;160;160];
Kv_ang = [1;1;1];

phi_des = (cmd_acc(1)*sin(des_state.yaw) - cmd_acc(2)*cos(des_state.yaw))/params.gravity;
theta_des = (cmd_acc(1)*cos(des_state.yaw) + cmd_acc(2)*sin(des_state.yaw))/params.gravity;

rot_des = [phi_des; theta_des; des_state.yaw];
omega_des = [0; 0; des_state.yawdot];

er = rot_des-state.rot;
eo = omega_des-state.omega;
M = Kp_ang.*er + Kv_ang.*eo;

% =================== Your code ends here ===================

end
