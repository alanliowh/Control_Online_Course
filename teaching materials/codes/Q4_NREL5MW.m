%% Q4
clear all; close all;

addpath WT_Data lib WindFiles;

%% model complexity
mode = 'omega+tower';
% "omega": turbine model with (1) drive-train and (2) wind speed
%     x = [omega;wsp] :: omega [rad/s], wsp [m/s]; u = [Qe,theta] :: Qe [Nm], theta [deg]
% "tower+omega": turbine model with (1) drive-train, (2) tower and (3) wind speed
%    x = [omega;v; dot_xt ; xt ; ] ::  omega [rad/s], v [m/s], dot_xt [m/s], xt [m];  u = [Qe,theta] :: Qe [Nm], theta [deg]

%% choose wind speed here

wind_no = 8;

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
controller.OpenLoop_Torque = 1; %[Nm]
controller.OpenLoop_Pitch = 1; %[deg]

%%%% ---- Closed-loop parameters
controller.Kopt = 0.18e7;
controller.Kp25 = 0  ;
controller.Ki25 = 0 ; 
controller.Kp3 = 0; % rad/(rad/s)
controller.Ki3 =0 ; % rad/(rad)
controller.KK1 = 0; % deg
controller.KK2 = 0; % deg^2
controller.TorqueCtrlRatio = 1; % constant 15 constant power =1 ;constant torque =0;

%% simulation script
main_script_NREL5MW;

%% gen plot
gen_plot;

AEP = sum(u(1,100/turbine.dt:end).*x(1,100/turbine.dt:end))*turbine.dt*1/(60*60*1e6); %J -> MWh
disp(['AEP_baseline = ',num2str(AEP,5),'MWh'])
