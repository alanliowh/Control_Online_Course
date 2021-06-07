%%% run simulation


%% Initialise parameters and wind speed files
turbine = [];
[turbine,controller] = InitWT_DTU10(turbine,controller);

%%% simulation time
sim.dt = 0.01;
t = [0:sim.dt:sim.Tend];
[wsp,wsp_file_name] = MakeWSP(wind_no,t); % read wind data



%% Simulution

[x,u] = init_sim_array(mode,wsp,controller);

textprogressbar(['Wind Speed file: ',wsp_file_name]);
tic
for i = 1:length(t)-1
    x(2,i) = wsp(i); % update the wind speed
    [x(:,i+1),u(:,i+1),out] = f_clsystem(x(:,i),u(:,i),turbine,controller,mode);
    if rem(i,100) == 0
    textprogressbar(100*i/(length(t)-1));
    end
end
toc



%% gen plot
% "omega": turbine model with (1) drive-train and (2) wind speed
%     x = [omega;wsp] :: omega [rad/s], wsp [m/s]; u = [Qe,theta] :: Qe [Nm], theta [deg]
% "tower+omega": turbine model with (1) drive-train, (2) tower and (3) wind speed
%    x = [omega;v; dot_xt ; xt ; ] ::  omega [rad/s], v [m/s], dot_xt [m/s], xt [m];  u = [Qe,theta] :: Qe [Nm], theta [deg]
% 
% figure(1)
% subplot(2,1,1);
% plot(t,x(1,:)); hold on; plot(t,ones(1,length(t))*controller.ratedOmega,'--','linewidth',1.5);
% xlabel('Time [s]');
% ylabel('omega [rad/s]');
% 
% subplot(2,1,2);
% plot(t,x(2,:));hold on;
% xlabel('Time [s]');
% ylabel('v [m/s]');
% 
% figure(2)
% subplot(2,1,1);
% plot(t,u(1,:));hold on;
% xlabel('Time [s]');
% ylabel('Qe [Nm]');
% 
% subplot(2,1,2);
% plot(t,u(2,:));hold on;
% xlabel('Time [s]');
% ylabel('theta [deg]');
% 
% if size(x,1) >2
%     figure(3)
%     subplot(2,1,1);
%     plot(t,x(3,:));hold on;
%     xlabel('Time [s]');
%     ylabel('dot_xt [m/s]');
%     
%     subplot(2,1,2);
%     plot(t,x(4,:));hold on;
%     xlabel('Time [s]');
%     ylabel('xt [m]');
%     
% end

