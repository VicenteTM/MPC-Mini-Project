addpath(fullfile('..', 'src'));

close all
clear all
clc

Ts = 1/20; 
rocket = Rocket(Ts);

[xs, us] = rocket.trim(); % Compute steady−state for which 0 = f(xs,us)
sys = rocket.linearize(xs, us);

[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us);

% x open loop controller 
H_x = 8;
mpc_x = MpcControl_x(sys_x, Ts, H_x);
x0_x = [0, 0, 0, 0]';
ref_x = -4;
[u_x, T_opt_x, X_opt_x, U_opt_x] = mpc_x.get_u(x0_x, ref_x);
U_opt_x(:,end+1) = nan;
ph = rocket.plotvis_sub(T_opt_x, X_opt_x, U_opt_x, sys_x, xs, us); 

% x closed loop controller
Tf_x = 8;
[T_x, X_sub_x, U_sub_x] = rocket.simulate_f(sys_x, x0_x, Tf_x, @mpc_x.get_u, ref_x);
ph_x = rocket.plotvis_sub(T_x, X_sub_x, U_sub_x, sys_x, xs, us, ref_x);

% y open loop controller
H_y = 8;
mpc_y = MpcControl_y(sys_y, Ts, H_y);
x0_y = [0, 0, 0, 0]';
ref_y = -4;
[u_y, T_opt_y, X_opt_y, U_opt_y] = mpc_y.get_u(x0_y, ref_y);
U_opt_y(:,end+1) = nan;
ph = rocket.plotvis_sub(T_opt_y, X_opt_y, U_opt_y, sys_y, xs, us); 

% y closed loop controller
Tf_y = 8;
[T_y, X_sub_y, U_sub_y] = rocket.simulate_f(sys_y, x0_y, Tf_y, @mpc_y.get_u, ref_y);
ph_y = rocket.plotvis_sub(T_y, X_sub_y, U_sub_y, sys_y, xs, us, ref_y);

% z open loop controller
H_z = 8;
mpc_z = MpcControl_z(sys_z, Ts, H_z);
x0_z = [0, 0]';
ref_z = -4;
[u_z, T_opt_z, X_opt_z, U_opt_z] = mpc_z.get_u(x0_z, ref_z);
U_opt_z = U_opt_z + us(3);
U_opt_z(:,end+1) = nan;
ph = rocket.plotvis_sub(T_opt_z, X_opt_z, U_opt_z, sys_z, xs, us); 

% z closed loop controller
Tf_z = 8;
[T_z, X_sub_z, U_sub_z] = rocket.simulate_f(sys_z, x0_z, Tf_z, @mpc_z.get_u, ref_z);
ph_z = rocket.plotvis_sub(T_z, X_sub_z, U_sub_z, sys_z, xs, us, ref_z);

% roll open loop controller
H_roll = 8;
mpc_roll = MpcControl_roll(sys_roll, Ts, H_roll);
x0_roll = [0, 0]';
ref_roll = deg2rad(35);
[u_roll, T_opt_roll, X_opt_roll, U_opt_roll] = mpc_roll.get_u(x0_roll, deg2rad(35));
U_opt_roll(:,end+1) = nan;
ph = rocket.plotvis_sub(T_opt_roll, X_opt_roll, U_opt_roll, sys_roll, xs, us); 

% roll closed loop controller
Tf_roll = 8;
[T_roll, X_sub_roll, U_sub_roll] = rocket.simulate_f(sys_roll, x0_roll, Tf_roll, @mpc_roll.get_u, ref_roll);
ph_roll = rocket.plotvis_sub(T_roll, X_sub_roll, U_sub_roll, sys_roll, xs, us, ref_roll);


%% TODO: This file should produce all the plots for the deliverable