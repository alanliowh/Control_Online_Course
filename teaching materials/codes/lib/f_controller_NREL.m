% f_controller DTU Basic controller.
%   [elecTorq,pitch] = f_controller(omega,v,turbine,controller) returns the
%   turbine input, generator torque (elecTorq) in Nm and pitch angle (pitch) in deg
%   based upon the measurement rotor speed (omega)

%
%   Written by Alan Wai Hou Lio (wali@dtu.dk)
%   Date 20.12.11

%%% Parameters for DTU10MW
% turbine.J =  1.616464739783273e+08 ;  % inertia of the rotor
% turbine.r = 89.15;
% turbine.n_gear = 1;
% turbine.rho = 1.225;
% turbine.neff = 0.94; % ignored the efficiency. We can offset it later
% turbine.Prated = 10e6;
% turbine.vrated = 11.4;
% turbine.dt = 0.01; % s
%
% controller.deratingMode = 0; controller.dr = 50; % \%
% controller.Kopt = 13013100;
% controller.Kp25 = 0.683456e8  ;
% controller.Ki25 = 0.153367e8 ;
% controller.Kp3 = 1.06713  ; % rad/(rad/s)
% controller.Ki3 = 0.242445 ; % rad/(rad)
% controller.Kp3Pow = 0.4e-8; % rad/W
% controller.Ki3Pow = 0.4e-8;% rad/(Ws)
% controller.KK1 = 11.4; % deg
% controller.KK2 = 402.9; % deg^2
% controller.switchPitch = 2; % deg
% controller.minOmega = 0.628; % rad/s
% controller.ratedOmega = 1.005; %rad/s
% controller.TorqMax = 15.6e6;
% controller.maxPitch = deg2rad(90);
% controller.minPitch = 0;
% controller.minth = 2 ;% deg
%controller.filt_wind.tau = 2; % constant 36: Time constant for wind speed low pass filter for minimum pitch [1/1P]
%controller.filt_pitch.tau = 1; % constant 37: Time constant for low pass filter for gain-scheduling [1/1P]
% controller.TorqueCtrlRatio = 1;%  constant power = 1 ; constant torque = 0
% controller.filt_omega.f0 = 0.4; % constant 8: Frequency of generator speed filter [Hz]
% controller.filt_omega.zeta = 0.7; % constant 9: Damping ratio of speed filter [-]
% controller.pitch_velmax = 10; % constant 7: maximum pitch velocity [deg/s]
% controller.optimal_lambda = 7.8;%  constant 49: 	Optimal tip speed ratio [-]

function [elecTorq,pitch,out] = f_controller(omega,turbine,controller,v,in)
%% assignment of local parameters
minOmega = controller.minOmega;
ratedOmega = controller.ratedOmega;
dt = turbine.dt;
TorqMax = controller.TorqMax;

Kopt = controller.Kopt;
Kp25 = controller.Kp25;
Ki25 = controller.Ki25;
Kp3 = controller.Kp3;
Ki3 = controller.Ki3;
KK1 = controller.KK1;
KK2 = controller.KK2;
maxPitch = controller.maxPitch;
Prated = turbine.Prated;

% persistent Kiterm2 Kiterm3 region  vf vprev  pitchfilt pitchprev PID_pit_var filt_pitch filt_omega filt_wind PID_gen_var
% if isempty(region)
%     Kiterm2 = 0; Kiterm3 = 0; region = 2;  vf = 0; vprev = 0; pitchfilt = 0; pitchprev = 0;
%     PID_pit_var.Kpro = controller.Kp3; PID_pit_var.Kint = controller.Ki3; PID_pit_var.Kdif = 0;
%     PID_pit_var.velmax = deg2rad(controller.pitch_velmax);
%     PID_gen_var.Kpro = controller.Kp25; PID_gen_var.Kint = controller.Ki25; PID_gen_var.Kdif = 0;
% end

Kiterm2 = in.Kiterm2 ;  Kiterm3 = in.Kiterm3; region = in.region;
vf = in.vf;  vprev = in.vprev ;  pitchfilt = in.pitchfilt ;
pitchprev = in.pitchprev ; 
filt_wind = in.filt_wind;  filt_omega = in.filt_omega; filt_pitch = in.filt_pitch;
PID_pit_var.Kpro =   in.PID_pit_var.Kpro ;PID_pit_var.Kint = in.PID_pit_var.Kint; PID_pit_var.Kdif = in.PID_pit_var.Kdif ;
PID_pit_var.velmax  =    in.PID_pit_var.velmax ;
PID_gen_var.Kpro =    in.PID_gen_var.Kpro; PID_gen_var.Kint =  in.PID_gen_var.Kint ; PID_gen_var.Kdif  = in.PID_gen_var.Kdif ;



if nargin <=3
    minPitch = controller.minPitch;
else
    [vf,filt_wind] = f_lowpass1orderfilt(v,filt_wind,dt);
    minPitch = deg2rad(GetOptPitch(vf));
end




%%
%%%%%%%%%%%%%%%%%%%%%%%%
% torque controller
%%%%%%%%%%%%%%%%%%%%%%%%%
if omega > 0.5*(minOmega+ratedOmega)
    omega_ref = ratedOmega;
else
    omega_ref = minOmega;
end


[omegaf,filt_omega] = f_lowpass2orderfilt(omega, filt_omega, dt);
omega = omegaf;

err = omega- omega_ref;

GenTorqueMin_full = Prated/ratedOmega*(1-controller.TorqueCtrlRatio)+ Prated/omega*controller.TorqueCtrlRatio;
% constraint on the torque output
rel_sp_open_Qg = 0.95; % constant 35
GenSpeed_min1 = minOmega;
GenSpeed_min2 = minOmega/rel_sp_open_Qg;
GenSpeed_max1 = (2*rel_sp_open_Qg - 1)*ratedOmega;
GenSpeed_max2 = rel_sp_open_Qg*ratedOmega;


if region == 2
    switchVar = f_switch(rad2deg(minPitch),controller.minth + rad2deg(minPitch),pitchprev);
elseif region ==3
    switchVar = 1;
end

% if omega < minOmega
%     minTorq = -1e7;
% end

%[elecTorq,Kiterm2] = f_PID_fast(err,Kiterm2,Kp25,Ki25,dt,maxTorq,minTorq);
if omega<=minOmega
    elecTorq = 0;
elseif omega > minOmega && omega < GenSpeed_min2
    elecTorq = (Kopt*GenSpeed_min2^2-0)/(GenSpeed_min2-minOmega)*(omega-minOmega);
elseif omega>=GenSpeed_min2 && omega<GenSpeed_max1
    elecTorq = Kopt*omega^2;
elseif omega>= GenSpeed_max1 && omega < ratedOmega
    elecTorq = Kopt*GenSpeed_max1^2 + (Prated/ratedOmega - Kopt*GenSpeed_max1^2)/(ratedOmega-GenSpeed_max1)*(omega-GenSpeed_max1);
elseif omega>=ratedOmega
    elecTorq=GenTorqueMin_full;
end
%PID_gen_var.outmin = minTorq; PID_gen_var.outmax = maxTorq; 
%[elecTorq,PID_gen_var] = f_PID(dt,1,PID_gen_var,err);

%%
%%%%%%%%%%%%%%%
% pitch controller % pitch (rad)
%%%%%%%%%%%%%%%%

minPitch = deg2rad(controller.pitch_wpdata(vf));

errPitch = omega - ratedOmega;
%%% gain-scheduling
[pitchfilt,filt_pitch] = f_lowpass1orderfilt( pitchprev,filt_pitch,dt);
pitchfilt = deg2rad(pitchfilt); % pitchfilt in rad
if KK1 >0
    invKK1 = pitchfilt/deg2rad(KK1);
else
    invKK1 =0;
end
if KK2 >0
    invKK2 = pitchfilt^2/deg2rad(deg2rad(KK2));
else
      invKK2 = 0;
end

k_gain = 1/(1+invKK1 + invKK2);
%k_gain = (1+pitchfilt/deg2rad(KK1) + pitchfilt^2/deg2rad(deg2rad(KK2)));

[pitch,Kiterm3] = f_PID_fast(errPitch,Kiterm3,Kp3,Ki3,dt,maxPitch,minPitch,k_gain);
%PID_pit_var.outmin = minPitch; PID_pit_var.outmax = maxPitch;
%[pitch,PID_pit_var] = f_PID(dt,k_gain,PID_pit_var,errPitch);
pitch = rad2deg(pitch);
pitchprev = pitch;
%% defining region
if pitch   < 0.01 + rad2deg(minPitch)
    region = 2;
elseif pitch > controller.minth + rad2deg(minPitch)
    region = 3;
end



out.region = region;
out.minPitch = minPitch;
out.pitchfilt = pitchfilt;
out.switchVar = switchVar;
out.kgain = k_gain;

out.Kiterm2 = Kiterm2 ;  out.Kiterm3 = Kiterm3; out.region = region;
out.vf = vf;  out.vprev = vprev ;  out.pitchfilt = pitchfilt ;
out.pitchprev = pitchprev ;
out.filt_wind = filt_wind;  out.filt_omega = filt_omega; out.filt_pitch = filt_pitch;

out.PID_pit_var.Kpro =   PID_pit_var.Kpro ; out.PID_pit_var.Kint = PID_pit_var.Kint; out.PID_pit_var.Kdif = PID_pit_var.Kdif ;
out.PID_pit_var.velmax  =    PID_pit_var.velmax ;
out.PID_gen_var.Kpro =  PID_gen_var.Kpro; out.PID_gen_var.Kint =  PID_gen_var.Kint ; out.PID_gen_var.Kdif  = PID_gen_var.Kdif ;


end


