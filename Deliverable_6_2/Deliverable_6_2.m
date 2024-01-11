clc; 
clear; 
close all;

%% Setup
Ts = 1/40;
rocket = Rocket(Ts);
H = 4; % Horizon length in seconds
x0 = zeros(12,1);
ref = [0.5, 0, 1, deg2rad(65)]';
Tf = 2.5;
rocket.mass = 1.75;


%% Small delay
Td = 8;
nmpc = NmpcControl(rocket, H,Td);

rocket.delay = Td; 
[T, X, U, Ref] = rocket.simulate(x0, Tf, @nmpc.get_u, ref);
ph = rocket.plotvis(T, X, U, Ref);
% ph.fig.Name = sprintf('Delay steps= %d, Delay Time= %ds', Td, Td*Ts); write Delay time with 1 decimal
ph.fig.Name = sprintf('Stable - Delay steps= %d, Delay Time= %ds', Td, Td*Ts);

%% Lost of performance delay
Td = 16;
nmpc = NmpcControl(rocket, H,Td);

rocket.delay = Td; 
[T, X, U, Ref] = rocket.simulate(x0, Tf, @nmpc.get_u, ref);
ph = rocket.plotvis(T, X, U, Ref);
ph.fig.Name = sprintf('Lost performance - Delay steps= %d, Delay Time= %ds', Td, Td*Ts);

%% Unstable delay
Td = 28;
nmpc = NmpcControl(rocket, H,Td);

rocket.delay = Td;
[T, X, U, Ref] = rocket.simulate(x0, Tf, @nmpc.get_u, ref);
ph = rocket.plotvis(T, X, U, Ref);
ph.fig.Name = sprintf('unstable - Delay steps= %d, Delay Time= %ds', Td, Td*Ts);

%% Uncompensated delay
Td = 8;
lost_delay = 4;
expected_delay = Td-lost_delay;
nmpc = NmpcControl(rocket, H,expected_delay);

rocket.delay = Td;
[T, X, U, Ref] = rocket.simulate(x0, Tf, @nmpc.get_u, ref);
ph = rocket.plotvis(T, X, U, Ref);
% ph.fig.Name = sprintf('Delay steps= %d, Delay Time= %ds', Td, Td*Ts); write Delay time with 1 decimal
ph.fig.Name = sprintf('Uncompensated delay - Delay steps=%d, expected delay steps=%d',Td,expected_delay);

