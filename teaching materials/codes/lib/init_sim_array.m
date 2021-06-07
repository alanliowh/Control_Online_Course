%%% preparing simulation arrays

function [x,u] = init_sim_array(mode,wsp,controller)

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
end
