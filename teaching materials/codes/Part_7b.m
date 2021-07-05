%Part_7b: design your own peak shaving strategy

%% model complexity
mode = 'omega';
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
sim.Tend = 650; % 1200s for step

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
pitch_wpdata = importdata('WT_Data\wpdata_DTU10MW_Part7.100');
controller.pitch_wpdata = @(v) interp1(pitch_wpdata.data(:,1),pitch_wpdata.data(:,2),v);

%% simulation script
main_script;

pitch = u(2,:);
genTorq = u(1,:);
omega = x(1,:);
wsp = x(2,:);
%% procedure of calculating dQdtheta from a Cp surface
turbine.wsp_ss = wsp(149/sim.dt:50/sim.dt:end);
turbine.pitch_ss = pitch(149/sim.dt:50/sim.dt:end);
turbine.omega_ss = omega(149/sim.dt:50/sim.dt:end);
turbine.tsr_ss = turbine.omega_ss*turbine.r./turbine.wsp_ss;
turbine.Ct_ss = turbine.Ct(turbine.tsr_ss,turbine.pitch_ss);
turbine.thrust_ss = 1/2*turbine.rho*turbine.r^2*pi.*turbine.wsp_ss.^2.*turbine.Ct_ss;

turbine.power_ss = 1/2*turbine.rho*turbine.r^2*pi.*turbine.wsp_ss.^3.*turbine.Cp(turbine.tsr_ss,turbine.pitch_ss);
subplot(2,2,1);
plot(turbine.wsp_ss,turbine.power_ss,'x-'); hold on; xlabel('wsp [m/s]'); ylabel('Power [W]');
subplot(2,2,2);
plot(turbine.wsp_ss,turbine.omega_ss,'x-'); hold on; xlabel('wsp [m/s]');  ylabel('Rotor Speed [rad/s]');
subplot(2,2,3);
plot(turbine.wsp_ss,turbine.thrust_ss,'x-'); hold on; xlabel('wsp [m/s]'); ylabel('Thrust [N]');
subplot(2,2,4);
plot(turbine.wsp_ss,turbine.pitch_ss,'x-'); hold on; xlabel('wsp [m/s]'); ylabel('Pitch [deg]');


