%% Part 2: PI exercise
clear all; close all;

addpath WT_Data lib WindFiles;

%% model complexity
mode = 'omega';
% "omega": turbine model with (1) drive-train and (2) wind speed
%     x = [omega;wsp] :: omega [rad/s], wsp [m/s]; u = [Qe,theta] :: Qe [Nm], theta [deg]
% "tower+omega": turbine model with (1) drive-train, (2) tower and (3) wind speed
%    x = [omega;v; dot_xt ; xt ; ] ::  omega [rad/s], v [m/s], dot_xt [m/s], xt [m];  u = [Qe,theta] :: Qe [Nm], theta [deg]

%% choose wind speed here


wind_no = 0;
% 0: Part 1 wind speeds.
% 1: for step wind speed, time = [0,1200]
% 2: for stochastic wind speed, mean wind speed: 8  m/s
% 3: for stochastic wind speed, mean wind speed: 12 m/s
% 4: for stochastic wind speed, mean wind speed: 15 m/s
% 5: for stochastic wind speed, mean wind speed: 18 m/s

sim.Tend = 400; % 1200s for step

%% controller parameters

controller.type = 'PI'; % 'OL' : open-loop, 'CL': closed-loop, 'PI': PI Study Region 3
%%% --- Open-loop parameters
controller.OpenLoop_Torque = 10e6/1.005; %[Nm]
controller.OpenLoop_Pitch = 6.78; %[deg]

%%%% ---- Closed-loop parameters
controller.Kopt = 0; 
controller.Kp25 = 0  ;
controller.Ki25 = 0 ; 
controller.Kp3 = 10; % rad/(rad/s)
controller.Ki3 =2 ; % rad/(rad)
controller.KK1 = 11.4; % deg
controller.KK2 = 402.9; % deg^2
controller.TorqueCtrlRatio = 0; % constant 15 constant power =1 ;constant torque =0;

%% simulation script
main_script;

%% gen plot

 
figure(1)
subplot(2,2,1); 
plot(t,x(1,:),'linewidth',2); hold on; plot(t,ones(1,length(t))*controller.ratedOmega,'--','linewidth',2);
xlabel('Time [s]');
ylabel('omega [rad/s]');
grid on;

subplot(2,2,2);
plot(t,x(1,:).*u(1,:)/1e6,'linewidth',2); hold on;
xlabel('Time [s]');
ylabel('Power [MW]');
grid on;


subplot(2,2,3);
plot(t(5:end),u(2,5:end),'linewidth',2); hold on;
xlabel('Time [s]');
ylabel('Pitch [deg]');
grid on;

subplot(2,2,4);
plot(t,wsp,'linewidth',2); hold on;
xlabel('Time [s]');
ylabel('Wind Speed [m/s]');
grid on;


