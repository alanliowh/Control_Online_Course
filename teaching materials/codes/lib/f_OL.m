% Open-loop system

function [x_out] = f_OL(x,u,turbine,mode)

switch mode
    case 'omega'
        % x = [omega;v] :: omega [rad/s], v [m/s]
        % u = [Qe,theta] :: Qe [Nm], theta [deg]
        [x(1),~]=f_omega(x(2),u(2),x(1),u(1),turbine);
        xTemp(2) = x(2);
        x_out = x;
    case 'omega+tower'
        % x = [omega;v; dot_xt ; xt ; ] ::  omega [rad/s], v [m/s], dot_xt [m/s], xt [m]
        % u = [Qe,theta] :: Qe [Nm], theta [deg]
        [xTemp(1),~]=f_omega(x(2)-x(3),u(2),x(1),u(1),turbine);
        [xTemp(3),xTemp(4)] = f_tower(x(3),x(4),x(2)-x(3),u(2),x(1),turbine);
        xTemp(2) = x(2);
        x_out = xTemp';
    otherwise
        warning('*** Check which model you want. For example, ***');
        warning('*** "omega": turbine model with (1) drive-train and (2) wind speed  ***');
        warning('*** "omega+tower": turbine model with (1) drive-train, (2) tower and (3) wind speed  ***')
end
%system.x = x;

%x_out(2) = f_wind(x(2),turbine);
%[u_out(1),u_out(2),out] = f_controller(x_out(1),turbine,controller,x_out(2));

