%% 1-Hz equivalent fatigue load calculation
% https://toolbox.pages.windenergy.dtu.dk/WindEnergyToolbox/fatigue_tools/fatigue_nb.html


function DEL1Hz = calc_DEL(signal,m, dt);
% dt: sampling time

n = length(signal);
Neq = n*dt; % in s

c = rainflow(signal);
sum_c = 0;
for ci = 1:length(c)
    sum_c  = sum_c + (c(ci,1)*c(ci,2).^m);
end
DEL1Hz = (sum_c / Neq)^(1/m);
