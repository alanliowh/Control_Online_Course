%Part_7c: peak shaving steady-state
clc; close all; clear all;
%% model complexity
mode = 'omega+tower';
% "omega": turbine model with (1) drive-train and (2) wind speed
%     x = [omega;wsp] :: omega [rad/s], wsp [m/s]; u = [Qe,theta] :: Qe [Nm], theta [deg]
% "tower+omega": turbine model with (1) drive-train, (2) tower and (3) wind speed
%    x = [omega;v; dot_xt ; xt ; ] ::  omega [rad/s], v [m/s], dot_xt [m/s], xt [m];  u = [Qe,theta] :: Qe [Nm], theta [deg]

%% choose wind speed here

wind_no = 3;

% 0: Part 1 wind speeds.
% 1: for step wind speed, time = [0,1200]
% 2: for stochastic wind speed, mean wind speed: 8  m/s
% 3: for stochastic wind speed, mean wind speed: 12 m/s
% 4: for stochastic wind speed, mean wind speed: 15 m/s
% 5: for stochastic wind speed, mean wind speed: 18 m/s
% 6: Part 4 step from 10 m/s to 12 m/s
% 7: Part 5 step from 11 m/s to 25 m/s
% 8: Part 7 step from 4 m/s to 15 m/s
sim.Tend = 550; % 1200s for step

%% controller parameters

controller.type = 'CL'; % 'OL' : open-loop, 'CL': closed-loop, 'PI': PI Study Region 3
%%% --- Open-loop parameters
controller.OpenLoop_Torque = 10e6/1.005; %[Nm]
controller.OpenLoop_Pitch = 6.78; %[deg]

%%%% ---- Closed-loop parameters
controller.Kopt = 9.8182e+06 ; 
controller.Kp25 = 6.4949e7  ;
controller.Ki25 =1.4575e7; 
controller.Kp3 = 1; % rad/(rad/s)
controller.Ki3 =0.2 ; % rad/(rad)
controller.KK1 = 0; % deg
controller.KK2 = 0; % deg^2
controller.TorqueCtrlRatio = 1; % constant 15 constant power =1 ;constant torque =0;

%% simulation script
main_script;

gen_plot;
AEP = sum(u(1,100/turbine.dt:end).*x(1,100/turbine.dt:end))*turbine.dt*1/(60*60*1e6); %J -> MWh
disp(['AEP_baseline = ',num2str(AEP,5),'MWh'])

%% simulate your design
controller.rel_sp_open_Qg = x;
main_script;
gen_plot;
AEP = sum(u(1,100/turbine.dt:end).*x(1,100/turbine.dt:end))*turbine.dt*1/(60*60*1e6); %J -> MWh
disp(['AEP_TorqLimit = ',num2str(AEP,5),'MWh'])