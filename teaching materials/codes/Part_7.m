%% Part 4: Region 2.5 exercise: DTUWEC vs NREL Controller
clear all; close all;

addpath WT_Data lib WindFiles;

%% model complexity
mode = 'omega+tower';
% "omega": turbine model with (1) drive-train and (2) wind speed
%     x = [omega;wsp] :: omega [rad/s], wsp [m/s]; u = [Qe,theta] :: Qe [Nm], theta [deg]
% "tower+omega": turbine model with (1) drive-train, (2) tower and (3) wind speed
%    x = [omega;v; dot_xt ; xt ; ] ::  omega [rad/s], v [m/s], dot_xt [m/s], xt [m];  u = [Qe,theta] :: Qe [Nm], theta [deg]

%% choose wind speed here


wind_no = 1;
% 0: Part 1 wind speeds.
% 1: for step wind speed, time = [0,1200]
% 2: for stochastic wind speed, mean wind speed: 8  m/s
% 3: for stochastic wind speed, mean wind speed: 12 m/s
% 4: for stochastic wind speed, mean wind speed: 15 m/s
% 5: for stochastic wind speed, mean wind speed: 18 m/s

sim.Tend = 550; % 1200s for step

%% controller parameters

controller.type = 'CL'; % 'OL' : open-loop, 'CL': closed-loop, 'PI': PI Study Region 3
%%% --- Open-loop parameters
controller.OpenLoop_Torque = 10e6/1.005; %[Nm]
controller.OpenLoop_Pitch = 6.78; %[deg]

%%%% ---- Closed-loop parameters
controller.Kopt = 9.8182e+06 ; 
controller.Kp25 = 6.4949e7  ;
controller.Ki25 = 1.4575e7; 
controller.Kp3 = 1; % rad/(rad/s)
controller.Ki3 =0.2 ; % rad/(rad)
controller.KK1 = 11.4; % deg
controller.KK2 = 402.9; % deg^2
controller.TorqueCtrlRatio = 0; % constant 15 constant power =1 ;constant torque =0;

%% simulation script
main_script;

%% gen plot

pitch = u(2,:);
genTorq = u(1,:);
omega = x(1,:);
v = x(2,:);
power = u(1,:).*x(1,:);
if (size(x,1)>2), x_t = x(4,:); end

TimeStart = 300;
TimeEnd  = 450;

figure(1)
subplot(2,2,1); 
plot(t,x(1,:),'linewidth',2); hold on; %plot(t,ones(1,length(t))*controller.ratedOmega,'k--','linewidth',2);
xlabel('Time [s]');
ylabel('omega [rad/s]');
grid on;
%xlim([300 450])

subplot(2,2,2);
plot(t,wsp,'linewidth',2); hold on;
xlabel('Time [s]');
ylabel('Wind Speed [m/s]');
grid on;
%xlim([300 450])
lambda = turbine.r*x(1,:)./wsp;

subplot(2,2,3);
plot(t,lambda,'linewidth',2); hold on;
xlabel('Time [s]');
ylabel('Tip-speed ratio [-]');
grid on;
%xlim([300 450])

subplot(2,2,4);
plot(t,turbine.Cp(lambda,u(2,:)),'linewidth',2); hold on;
xlabel('Time [s]');
ylabel('Cp [-]');
grid on;
%xlim([300 450])



figure(2);
subplot(2,2,1);
plot(t,u(1,:).*x(1,:),'linewidth',2); hold on;
xlabel('Time [s]');
ylabel('Power [W]');
grid on;
AEP = sum(u(1,100/turbine.dt:end).*x(1,100/turbine.dt:end))*turbine.dt*1/(60*60*1e6); %J -> MWh
%title(['AEP = ',num2str(AEP,2),'MWh']);
disp(['AEP_for_DTU = ',num2str(AEP,2),'MWh'])
subplot(2,2,2);
plot(t,u(1,:),'linewidth',2); hold on;
xlabel('Time [s]');
ylabel('GenTorq [Nm]');
grid on;

subplot(2,2,3);
plot(t,u(2,:),'linewidth',2); hold on;
xlabel('Time [s]');
ylabel('Pitch [deg]');
grid on;

if size(x,1)>2 
subplot(2,2,4);
plot(t,x(4,:),'linewidth',2); hold on;
xlabel('Time [s]');
ylabel('Tower fore-aft [m]');
grid on;
end

%% running NREL controller;
clear x u strCR c
main_script_Part7_NREL;

%% re-do the plot

 
figure(1)
subplot(2,2,1); 
plot(t,x(1,:),'linewidth',2); hold on; plot(t,ones(1,length(t))*controller.ratedOmega,'k-','linewidth',2);
xlabel('Time [s]');
ylabel('omega [rad/s]');
grid on;
%xlim([300 450])

subplot(2,2,2);
plot(t,wsp,'linewidth',2); hold on;
xlabel('Time [s]');
ylabel('Wind Speed [m/s]');
grid on;
%xlim([300 450])
lambda = turbine.r*x(1,:)./wsp;

subplot(2,2,3);
plot(t,lambda,'linewidth',2); hold on;
xlabel('Time [s]');
ylabel('Tip-speed ratio [-]');
grid on;
%xlim([300 450])

subplot(2,2,4);
plot(t,turbine.Cp(lambda,u(2,:)),'linewidth',2); hold on;
xlabel('Time [s]');
ylabel('Cp [-]');
grid on;
%xlim([300 450])
legend('DTU','NREL')


figure(2);
subplot(2,2,1);
plot(t,u(1,:).*x(1,:),'linewidth',2); hold on;
xlabel('Time [s]');
ylabel('Power [W]');
grid on;
AEP = sum(u(1,100/turbine.dt:end).*x(1,100/turbine.dt:end))*turbine.dt*1/(60*60*1e6); %J -> MWh
disp(['AEP_for_NREL = ',num2str(AEP,2),'MWh'])

subplot(2,2,2);
plot(t,u(1,:),'linewidth',2); hold on;
xlabel('Time [s]');
ylabel('GenTorq [Nm]');
grid on;

subplot(2,2,3);
plot(t,u(2,:),'linewidth',2); hold on;
xlabel('Time [s]');
ylabel('Pitch [deg]');
grid on;
legend('DTU','NREL')
if size(x,1)>2 
subplot(2,2,4);
plot(t,x(4,:),'linewidth',2); hold on;
xlabel('Time [s]');
ylabel('Tower fore-aft [m]');
grid on;
end
