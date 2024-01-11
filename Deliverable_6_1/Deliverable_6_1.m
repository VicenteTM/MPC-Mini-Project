clc; 
clear; 
close all;

%% Setup
Ts = 1/20;
rocket = Rocket(Ts);

H = 4; % Horizon length in seconds
nmpc = NmpcControl(rocket, H);
x0 = zeros(12,1);

%% Open loop
ref = [2 2 2 deg2rad(40)]';
[u, T_opt, X_opt, U_opt] = nmpc.get_u(x0, ref);
U_opt(:,end+1) = nan;
ph = rocket.plotvis(T_opt, X_opt, U_opt, ref);
ph.fig.Name = 'Nonlin. MPC in nonlinear open loop';
%% Simulation: MPC reference with default maximum roll = 15 deg
Tf = 30;
ref = @(t_, x_) ref_TVC(t_);
[T, X, U, Ref] = rocket.simulate(x0, Tf, @nmpc.get_u, ref);
ph = rocket.plotvis(T, X, U, Ref);
ph.fig.Name = 'Nonlin. MPC in nonlinear simulation, max roll = 15°';

%% Simulation: MPC reference with specified maximum roll = 50 deg
roll_max = deg2rad(50);
Tf = 30;
ref = @(t_, x_) ref_TVC(t_,roll_max);
[T, X, U, Ref] = rocket.simulate(x0, Tf, @nmpc.get_u, ref);
ph = rocket.plotvis(T, X, U, Ref);
ph.fig.Name = 'Nonlin. MPC in nonlinear simulation, max roll = 50';