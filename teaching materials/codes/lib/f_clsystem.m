% closed-loop system

function [x_out,u_out,out] = f_clsystem(x,u,turbine,controller,mode)



%%% initialise
persistent in
if isempty(in)
    in.Kiterm2 = 0; in.Kiterm3 = 0; in.region = 2;  in.vf = 0; in.vprev = 0; in.pitchfilt = 0; in.pitchprev = 0;
    in.PID_pit_var.Kpro = controller.Kp3; in.PID_pit_var.Kint = controller.Ki3; in.PID_pit_var.Kdif = 0;
    in.PID_pit_var.velmax = deg2rad(controller.pitch_velmax);
    in.PID_gen_var.Kpro = controller.Kp25; in.PID_gen_var.Kint = controller.Ki25; in.PID_gen_var.Kdif = 0;
    in.filt_pitch.tau = controller.filt_pitch.tau/controller.ratedOmega; % error in DTUWEC. We keep the wrong value.
    in.filt_wind.tau = controller.filt_wind.tau/controller.ratedOmega;% error in DTUWEC. We keep the wrong value.
    in.filt_omega.zeta = controller.filt_omega.zeta; in.filt_omega.f0 = controller.filt_omega.f0;
end

x_out = f_OL(x,u,turbine,mode);
switch controller.type
    case 'CL'
        [u_out(1),u_out(2),out] = f_controller(x_out(1),turbine,controller,x_out(2),in);
    case 'OL'
        u_out(1) = controller.OpenLoop_Torque;
        u_out(2) = controller.OpenLoop_Pitch;
        out = 0;
    case 'PI'
        [u_out(1),u_out(2),out] = f_controller(x_out(1),turbine,controller,x_out(2),in);

end
    in = out;

% additional output
% maxTorq = out.maxTorq;
% minTorq = out.minTorq;
% region = out.region;
% minPitch = out.minPitch;

