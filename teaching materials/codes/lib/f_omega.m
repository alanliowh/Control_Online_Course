function [omega,tau_a]=f_omega(v,pitch,omega,elec_torq,turbine)

% --- assign local variable
r = turbine.r;
rho = turbine.rho;
J = turbine.J;
neff = turbine.neff;
Cp_tab2 = turbine.Cp_tab2;
%Cm_tab2 = turbine.Cm_tab2;
%pitch_tab2 = round(turbine.pitch_tab2,1);
%lambda_tab2 = round(turbine.lambda_tab2,1);
pitch_tab2 = turbine.pitch_tab2;
lambda_tab2 = turbine.lambda_tab2;
dt = turbine.dt;

% --- --- find pitch_index
pitch_index_temp = find(pitch>=pitch_tab2);
pitch_index = max(pitch_index_temp(end),1);

lambda = omega*r/v;
lambda_index_temp = find(lambda>=lambda_tab2);
lambda_index = max(lambda_index_temp(end),1);
% 
%   % --- --- ensure lambda is not empty
%         if (lambda>max(lambda_tab2))
%         lambda = max(lambda_tab2);
%         warning('omg, the lambda is at max')
%         else if (lambda<min(lambda_tab2))
%             lambda = min(lambda_tab2);
%              warning('omg, the lambda is at MIN')
%              disp(['[v,lambda,omega]=[',num2str(v),',',num2str(lambda),',',num2str(omega),']']); 
%             end
%         end
%         lambda_index = find(lambda_tab2==round(lambda,1));
%         
%    %%% ensure Cp > 0 if lambda > 0
%    if Cp_tab2(lambda_index,pitch_index) <= 0
%        Cp_tab2(lambda_index,pitch_index) =0.04;  
%    end
%    
        
 %%% optimal      
%tau_a = 1/2*rho*r^2*pi*v.^3*(Cp_tab2(lambda_index,pitch_index)-0.0082)/omega; % calucalate the aerodynamic torque
%omega = omega + dt/(J)* (tau_a- 1.0035*elec_torq/0.94-2.2e3*omega); 
%%%%%%%%%%%%%%%%%%%%%%%%
% 
omega = max(omega,0.001); % ensure tau =/= Inf.
tau_a = 1/2*rho*r^2*pi*v.^3*Cp_tab2(lambda_index,pitch_index)/omega; % calucalate the aerodynamic torque
%tau_a = 1/2*rho*r^2*pi*v.^3*turbine.Cp(lambda,pitch)/omega; % calucalate the aerodynamic torque
%tau_a = 1/2*rho*r^2*pi*v.^3*(turbine.Cp(lambda,pitch))/omega; % calucalate the aerodynamic torque

omega = omega + dt/(J)*(tau_a-elec_torq/turbine.neff); %-dt*0.00023*omega; % 4e4   1.003*elec`_toq6.5e3*omega  


end   