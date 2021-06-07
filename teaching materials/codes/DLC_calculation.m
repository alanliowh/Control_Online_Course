%%% run simulation

clear all; close all;
addpath WT_Data lib WindFiles;

%% User input
mode = 'omega';
% "omega": turbine model with (1) drive-train and (2) wind speed
%     x = [omega;wsp] :: omega [rad/s], wsp [m/s]; u = [Qe,theta] :: Qe [Nm], theta [deg]
% "tower+omega": turbine model with (1) drive-train, (2) tower and (3) wind speed
%    x = [omega;v; dot_xt ; xt ; ] ::  omega [rad/s], v [m/s], dot_xt [m/s], xt [m];  u = [Qe,theta] :: Qe [Nm], theta [deg]

wind_file_name{1} = 'WindFiles/Kaimal_8ms.hh';
wind_file_name{2} = 'WindFiles/Kaimal_12ms.hh';
wind_file_name{3} = 'WindFiles/Kaimal_15ms.hh';


for case_indx = 1:length(wind_file_name)
%% Initialise parameters and wind speed files

turbine = []; controller = [];
[turbine,controller] = InitWT_DTU10(turbine,controller);

%%% simulation time
sim.dt = 0.1;
sim.Tend = 800;
t = [0:sim.dt:sim.Tend];
wsp = MakeWSP(wind_file_name{case_indx},t); % read wind data



%% Simulution

[x,u] = init_sim_array(mode,wsp);

tic
for i = 1:length(t)-1
    x(2,i) = wsp(i); % update the wind speed
    [x(:,i+1),u(:,i+1),out] = f_clsystem(x(:,i),u(:,i),turbine,controller,mode);
end
toc

res{case_indx}.x = x;
res{case_indx}.u = u;

end

Postprocessing(res,turbine);