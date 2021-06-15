%% Q3
clear all; close all;

addpath WT_Data lib WindFiles;

%% model complexity
mode = 'omega+tower';
% "omega": turbine model with (1) drive-train and (2) wind speed
%     x = [omega;wsp] :: omega [rad/s], wsp [m/s]; u = [Qe,theta] :: Qe [Nm], theta [deg]
% "tower+omega": turbine model with (1) drive-train, (2) tower and (3) wind speed
%    x = [omega;v; dot_xt ; xt ; ] ::  omega [rad/s], v [m/s], dot_xt [m/s], xt [m];  u = [Qe,theta] :: Qe [Nm], theta [deg]

%% choose wind speed here

wind_no = 5;

% 0: Part 1 wind speeds.
% 1: for step wind speed, time = [0,1200]
% 2: for stochastic wind speed, mean wind speed: 8  m/s
% 3: for stochastic wind speed, mean wind speed: 12 m/s
% 4: for stochastic wind speed, mean wind speed: 15 m/s
% 5: for stochastic wind speed, mean wind speed: 18 m/s
% 6: Part 4 step from 10 m/s to 12 m/s
% 7: Part 5 step from 11 m/s to 25 m/s
sim.Tend = 550; % 1200s for step

%% controller parameters

controller.type = 'CL'; % 'OL' : open-loop, 'CL': closed-loop, 'PI': PI Study Region 3
%%% --- Open-loop parameters
controller.OpenLoop_Torque = 1; %[Nm]
controller.OpenLoop_Pitch = 1; %[deg]

%%%% ---- Closed-loop parameters
controller.Kopt = 0.185446E+05;
controller.Kp25 = 0   ;
controller.Ki25 = 0 ; 
controller.Kp3 = 0; % rad/(rad/s)
controller.Ki3 =0 ; % rad/(rad)
controller.KK1 = 0; % deg
controller.KK2 = 0; % deg^2
controller.TorqueCtrlRatio = 1; % constant 15 constant power =1 ;constant torque =0;

%% simulation script
main_script_DTU850k;

%% gen plot
gen_plot;
%% procedure of calculating dQdtheta from a Cp surface
turbine.wsp_ss = wsp(149/sim.dt:50/sim.dt:end);
turbine.pitch_ss = pitch(149/sim.dt:50/sim.dt:end);
turbine.tsr_ss = controller.ratedOmega*turbine.r./turbine.wsp_ss;

%%% plot Cp surface
figure(3)
contourf(turbine.pitchList,turbine.tsrList,max(turbine.CpTable,0));hold on;
plot(turbine.pitch_ss,turbine.tsr_ss,'rx','markersize',10);
xlabel('\theta [deg]');
ylabel('\lambda [-]')

%%% gen dQdtheta
epsilon = 0.1;
for i = 1:length(turbine.pitch_ss)
   Q = 1/2*turbine.rho*turbine.r^2*pi.*turbine.wsp_ss(i)^3.*turbine.Cp(turbine.tsr_ss(i),turbine.pitch_ss(i))/controller.ratedOmega;
   Qep = 1/2*turbine.rho*turbine.r^2*pi.*turbine.wsp_ss(i)^3.*turbine.Cp(turbine.tsr_ss(i),turbine.pitch_ss(i)-epsilon)/controller.ratedOmega; 
   turbine.dQdtheta_ss(i) = (Q-Qep)./epsilon;
end

%%% calculate dQdtheta0, KK1 and KK2
[dQdtheta0,KK1,KK2] = find_KK1_KK2(turbine.pitch_ss,turbine.dQdtheta_ss);

%%% plot dQdtheta vs pitch
figure(4)
plot(turbine.pitch_ss,turbine.dQdtheta_ss,'x','linewidth',2,'markersize',8); hold on;       % plot data
est_dQdtheta = dQdtheta0.*(1+turbine.pitch_ss/KK1+ turbine.pitch_ss.^2/KK2);                % calculate the fitted dQdtheta
plot([0,turbine.pitch_ss],[dQdtheta0,est_dQdtheta],'o','linewidth',2,'markersize',8);       % plot fitted dQdtheta
legend('Original','Fitted');
ylabel('dQdtheta [Nm/deg]'); xlabel('\theta [deg]');
grid on

