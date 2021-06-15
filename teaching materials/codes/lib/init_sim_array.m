%%% preparing simulation arrays

function [x,u, in ] = init_sim_array(mode,wsp,controller)

switch mode
    case 'omega'
        % x = [omega;v] :: omega [rad/s], v [m/s]
        % u = [Qe,theta] :: Qe [Nm], theta [deg]
        
        x = [0;wsp(1)]*ones(1,length(wsp)); % x = [omega;v;dot_xt;xt]
        u = [0;0]*ones(1,length(wsp)); % u = [Qe;theta], Qe [Nm], theta [deg]
    case 'omega+tower'
        % x = [omega;v; dot_xt ; xt ; ] ::  omega [rad/s], v [m/s], dot_xt [m/s], xt [m]
        % u = [Qe,theta] :: Qe [Nm], theta [deg]
        
        x = [0;wsp(1);0;0]*ones(1,length(wsp));
        u = [0;0]*ones(1,length(wsp)); % u = [Qe;theta], Qe [Nm], theta [deg]
end
switch controller.type
    case 'OL'
        x = [1.05;wsp(1)]*ones(1,length(wsp)); % x = [omega;v;dot_xt;xt]
        u = [controller.OpenLoop_Torque;controller.OpenLoop_Pitch ]*ones(1,length(wsp)); % u = [Qe;theta], Qe [Nm], theta [deg]
    case 'PI'
        x = [1.05;wsp(1)]*ones(1,length(wsp)); % x = [omega;v;dot_xt;xt]
        u = [controller.OpenLoop_Torque;controller.OpenLoop_Pitch ]*ones(1,length(wsp)); % u = [Qe;theta], Qe [Nm], theta [deg]
end
%% initialise controller
    in.Kiterm2 = 0; in.Kiterm3 = 0; in.region = 2;  in.vf = 0; in.vprev = 0; in.pitchfilt = 0; in.pitchprev = 0;
    in.PID_pit_var.Kpro = controller.Kp3; in.PID_pit_var.Kint = controller.Ki3; in.PID_pit_var.Kdif = 0;
    in.PID_pit_var.velmax = deg2rad(controller.pitch_velmax);
    in.PID_gen_var.Kpro = controller.Kp25; in.PID_gen_var.Kint = controller.Ki25; in.PID_gen_var.Kdif = 0;
    in.filt_pitch.tau = controller.filt_pitch.tau/controller.ratedOmega; % error in DTUWEC. We keep the wrong value.
    in.filt_wind.tau = controller.filt_wind.tau/controller.ratedOmega;% error in DTUWEC. We keep the wrong value.
    in.filt_omega.zeta = controller.filt_omega.zeta; in.filt_omega.f0 = controller.filt_omega.f0;
    in.region = 2;
end
