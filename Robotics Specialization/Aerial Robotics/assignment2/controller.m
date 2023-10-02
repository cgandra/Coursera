function [ u1, u2 ] = controller(~, state, des_state, params)
%CONTROLLER  Controller for the planar quadrotor
%
%   state: The current state of the robot with the following fields:
%   state.pos = [y; z], state.vel = [y_dot; z_dot], state.rot = [phi],
%   state.omega = [phi_dot]
%
%   des_state: The desired states are:
%   des_state.pos = [y; z], des_state.vel = [y_dot; z_dot], des_state.acc =
%   [y_ddot; z_ddot]
%
%   params: robot parameters

%   Using these current and desired states, you have to compute the desired
%   controls

u1 = 0;
u2 = 0;

% FILL IN YOUR CODE HERE
Kpz = 80;
Kvz = 20;
Kpy = 20;
Kvy = 5;
Kpphi = 1000;
Kvphi = 10;

epy = des_state.pos(1) - state.pos(1);
evy = des_state.vel(1) - state.vel(1);

epz = des_state.pos(2) - state.pos(2);
evz = des_state.vel(2) - state.vel(2);


y_acc = -params.gravity*state.rot(1);
phic = -(des_state.acc(1) + Kvy*evy + Kpy*epy)/params.gravity;
phic_v = -(Kvy*(des_state.acc(1) - y_acc) + Kpy*evy)/params.gravity;
phic_v = 0;

epphi = phic - state.rot(1);
evphi = phic_v - state.omega(1);

u1 = params.mass*(params.gravity + des_state.acc(2) + Kvz*evz + Kpz*epz);
u2 = params.Ixx*(Kvphi*evphi + Kpphi*epphi);

end