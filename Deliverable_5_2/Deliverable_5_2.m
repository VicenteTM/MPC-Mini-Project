addpath(fullfile('..', 'src'));

close all
clear all
clc

Ts = 1/20; 
rocket = Rocket(Ts);
H = 6;
[xs, us] = rocket.trim();
sys = rocket.linearize(xs, us);
[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us);

mpc_x = MpcControl_x(sys_x, Ts, H);
mpc_y = MpcControl_y(sys_y, Ts, H);
mpc_zwe = MpcControl_zwe(sys_z, Ts, H);
mpc_roll = MpcControl_roll(sys_roll, Ts, H);

%offset-free tracking controller
mpc_zwe_merged = rocket.merge_lin_controllers(xs, us, mpc_x, mpc_y, mpc_zwe, mpc_roll);

% Setup reference function
Tf = 16.2;
ref = @(t_, x_) ref_TVC(t_);
x0 = zeros(12,1);

rocket.mass = 2.13; % Manipulate mass for simulation
rocket.mass_rate = - 0.27;
rocket.anim_rate = 3;
[T_we, X_we, U_we, Ref_we, Z_we] = rocket.simulate_est_z(x0, Tf, @mpc_zwe_merged.get_u, ref, mpc_zwe, sys_z);

ph_we = rocket.plotvis(T_we, X_we, U_we, Ref_we);


