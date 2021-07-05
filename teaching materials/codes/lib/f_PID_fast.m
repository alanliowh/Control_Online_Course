function [out,Kiterm] = f_PID_fast(err,Kiterm,Kp0,Ki0,dt,outmax,outmin,k_gain)



%%% gain-scheduling
if nargin > 7
    Kp = Kp0*k_gain;
    Ki = Ki0*k_gain;
else
    Kp = Kp0;
    Ki = Ki0;
end



Kiterm = Kiterm + Ki*err*dt;
Kpterm = Kp*err;
out = Kpterm+ Kiterm;


if nargin>5
    % limiting on the control action
    if out > outmax
        Kiterm = (outmax - Kpterm );
        out = outmax;
    end
    if out < outmin
        Kiterm = (outmin -Kpterm);
        out = outmin;
    end
end

end