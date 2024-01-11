addpath(fullfile('..', 'src'));


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
mpc_z = MpcControl_z(sys_z, Ts, H);
mpc_zwe = MpcControl_zwe(sys_z, Ts, H);
mpc_roll = MpcControl_roll(sys_roll, Ts, H);

% To assess the influence of offset-free tracking we compare
% the offset-free tracking controller with the original controller
%original controller
mpc = rocket.merge_lin_controllers(xs, us, mpc_x, mpc_y, mpc_z, mpc_roll);
%offset-free tracking controller
mpc_zwe_merged = rocket.merge_lin_controllers(xs, us, mpc_x, mpc_y, mpc_zwe, mpc_roll);

% Setup reference function
Tf = 30;
ref = @(t_, x_) ref_TVC(t_);
x0 = zeros(12,1);

rocket.mass = 2.13; % Manipulate mass for simulation
[T, X, U, Ref] = rocket.simulate(x0, Tf, @mpc.get_u, ref);
[T_we, X_we, U_we, Ref_we, Z_we] = rocket.simulate_est_z(x0, Tf, @mpc_zwe_merged.get_u, ref, mpc_zwe, sys_z);


%comparison results
figure;
subplot(2,1,1);
plot(T, X(12,:),T, Ref(3,:));
legend('Original controller','Reference')
title('Original controller');
xlabel('T(s)');
ylabel('Z(m)');

subplot(2,1,2);
plot(T_we, X_we(12,:),T_we, Ref_we(3,:));
legend('Offset-free tracking controller','Reference')
title('Offset-free tracking controller');
xlabel('T(s)');
ylabel('Zwe(m)');

% %simulation plots
rocket.anim_rate = 3;
ph_we = rocket.plotvis(T_we, X_we, U_we, Ref_we);
ph = rocket.plotvis(T, X, U, Ref);




