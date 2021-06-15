pitch = u(2,:);
genTorq = u(1,:);
omega = x(1,:);
v = x(2,:);
power = u(1,:).*x(1,:);
if (size(x,1)>2), x_t = x(4,:); end

TimeStart = 0 ;
TimeEnd  = 450;

figure(1)
subplot(2,2,1); 
plot(t,x(1,:),'linewidth',2); hold on; plot(t,ones(1,length(t))*controller.ratedOmega,'k--','linewidth',2);
xlabel('Time [s]');
ylabel('omega [rad/s]');
grid on;
xlim([TimeStart TimeEnd])

subplot(2,2,2);
plot(t,wsp,'linewidth',2); hold on;
xlabel('Time [s]');
ylabel('Wind Speed [m/s]');
grid on;
xlim([TimeStart TimeEnd])
lambda = turbine.r*x(1,:)./wsp;

subplot(2,2,3);
plot(t,lambda,'linewidth',2); hold on;
xlabel('Time [s]');
ylabel('Tip-speed ratio [-]');
grid on;
xlim([TimeStart TimeEnd])

subplot(2,2,4);
plot(t,turbine.Cp(lambda,u(2,:)),'linewidth',2); hold on;
xlabel('Time [s]');
ylabel('Cp [-]');
grid on;
xlim([TimeStart TimeEnd])



figure(2);
subplot(2,2,1);
plot(t,u(1,:).*x(1,:),'linewidth',2); hold on;
xlabel('Time [s]');
ylabel('Power [W]');
grid on;
AEP = sum(u(1,100/turbine.dt:end).*x(1,100/turbine.dt:end))*turbine.dt*1/(60*60*1e6); %J -> MWh
%title(['AEP = ',num2str(AEP,2),'MWh']);
xlim([TimeStart TimeEnd])

subplot(2,2,2);
plot(t,u(1,:),'linewidth',2); hold on;
xlabel('Time [s]');
ylabel('GenTorq [Nm]');
grid on;
xlim([TimeStart TimeEnd])

subplot(2,2,3);
plot(t,u(2,:),'linewidth',2); hold on;
xlabel('Time [s]');
ylabel('Pitch [deg]');
grid on;
xlim([TimeStart TimeEnd])

if size(x,1)>2 
subplot(2,2,4);
plot(t,x(4,:),'linewidth',2); hold on;
xlabel('Time [s]');
ylabel('Tower fore-aft [m]');
grid on;
xlim([TimeStart TimeEnd])

end