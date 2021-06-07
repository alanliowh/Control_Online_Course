function v_out = f_wind(v_in,turbine);


persistent vf
if isempty(vf)
    vf = 0;
end

dt = turbine.dt;
fpass = 5;% Hz
alpha = exp(-dt*fpass);
vf = alpha * vf + (1-alpha) * v_in; 


v_out = vf;