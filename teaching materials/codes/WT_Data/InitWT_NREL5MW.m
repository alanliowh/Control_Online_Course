% InitWT_10MW() returns variables (turbine) and (controller) based on NREL
% 5MW RWT
function [turbine,controller] = InitWT_NREL5MW(turbine,controller)

%%% Cp surface

load CP_NREL;
cptable.lam =  turbine.tsrList';
cptable.th = turbine.pitchList;
cptable.cp = turbine.CpTable';
cptable.ct = turbine.CtTable';

%load cptable; 
[Cp_tab2,lambda_tab2,pitch_tab2]=Cpt_curve2(cptable,'Cp');
turbine.Cp_tab2 = Cp_tab2;% offset makes it correct
[Ct_tab2,lambda_tab2,pitch_tab2]=Cpt_curve2(cptable,'Ct');
turbine.Ct_tab2 = Ct_tab2; % 0.04 offset makes it correct

%[Cm_tab2,lambda_tab2,pitch_tab2]=Cpt_curve2(cptable,'Cm');
%turbine.Cm_tab2 = Cm_tab2; % 0.04 offset makes it correct
turbine.lambda_tab2 = lambda_tab2;
turbine.pitch_tab2 = pitch_tab2;



% tower
turbine.omega_t = 0.32*2*pi; 
turbine.zeta_t = 0.2;
turbine.m_t = 449.9e3;

A = [-2*turbine.zeta_t*turbine.omega_t -turbine.omega_t^2;1 0];
B = [1/turbine.m_t;0];
turbine.dt = 0.01; % s
sys_tower = c2d(ss(A,B,[],[]),turbine.dt);
turbine.Ad_tower = sys_tower.A;
turbine.Bd_tower = sys_tower.B;



turbine.J= 4.6519e+07;
turbine.r = 63;
turbine.n_gear = 1;
turbine.rho = 1.225;
turbine.neff = 0.944; 
turbine.Prated = 5e6;
turbine.vrated = 11.4;
turbine.dt = 0.01; % s
turbine.Cpmax = max(max(turbine.Cp_tab2));
[x,y] = find(turbine.Cp_tab2==turbine.Cpmax);
turbine.lambda_opt = turbine.lambda_tab2(x);
turbine.pitch_opt = turbine.pitch_tab2(y);

controller.Kp3Pow = 0.4e-8; % rad/W
controller.Ki3Pow = 0.4e-8;% rad/(Ws)
%controller.switchPitch = 2; % deg
controller.minOmega = 0.7435; % rad/s
controller.ratedOmega = 1.2671; %rad/s
controller.TorqMax = 4.6e6;
controller.maxPitch = 90*pi/180;
controller.minPitch = 0;
controller.minth = 2 ;% deg
controller.filt_wind.tau = 2*2*pi/controller.ratedOmega; % constant 36: Time constant for wind speed low pass filter for minimum pitch [1/1P]
controller.filt_pitch.tau = 3*2*pi/controller.ratedOmega; % constant 37: Time constant for low pass filter for gain-scheduling [1/1P]
controller.filt_omega.f0 = 0.4; % constant 8: Frequency of generator speed filter [Hz]
controller.filt_omega.zeta = 0.7; %constant 9: Damping ratio of speed filter [-]
controller.pitch_velmax = 10; % constant 7: maximum pitch velocity [deg/s]
pitch_wpdata = importdata('WT_Data\wpdata_NREL5MW.100');
controller.pitch_wpdata = @(v) interp1(pitch_wpdata.data(:,1),pitch_wpdata.data(:,2),v);
%%% Do not use. Super slow.
[tsrGrid, pitchGrid]= meshgrid(turbine.lambda_tab2,turbine.pitch_tab2);
 turbine.Cp = @(tsr,theta) interp2(tsrGrid,pitchGrid,turbine.Cp_tab2',tsr,theta);
 turbine.Ct = @(tsr,theta) interp2(tsrGrid,pitchGrid,turbine.Ct_tab2',tsr,theta);
 

% turbine.dQdtheta = @(omega,U,theta) interp3(UGrid,omegaGrid,pitchGrid,turbine.dQdthetaTable,U,omega,theta,'spline');
% turbine.dQdv = @(omega,U,theta) interp3(UGrid,omegaGrid,pitchGrid,turbine.dQdvTable,U,omega,theta,'spline');
% turbine.dQdomega = @(omega,U,theta) interp3(UGrid,omegaGrid,pitchGrid,turbine.dQdomegaTable,U,omega,theta,'spline');
% 
% 


