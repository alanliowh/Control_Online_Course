% discrete-time function to calculate the fore-aft tower-top displacement
function [dotp,p]=f_tower(dotp,p,v,pitch,omega,turbine)
% p_vec =  [dotp,p]

% define fore-aft tower-top dynmiacs parameters
p_vect = [dotp;p];

m_t = turbine.m_t;

zeta_t = turbine.zeta_t;
omega_t = turbine.omega_t;


% --- assign local variable
r = turbine.r;
rho = turbine.rho;
%J = turbine.J;
%neff = turbine.neff;
Ct_tab2 = turbine.Ct_tab2;
pitch_tab2 = round(turbine.pitch_tab2,1);
lambda_tab2 = round(turbine.lambda_tab2,1);

% --- --- find pitch_index
pitch_index_temp = find(pitch>=pitch_tab2);
pitch_index = max(pitch_index_temp(end),1);

lambda = omega*r/v;
lambda_index_temp = find(lambda>=lambda_tab2);
lambda_index = max(lambda_index_temp(end),1);


% calucalate the aerodynamic thrust        
F_t = 1/2*rho*r^2*pi*v.^2*Ct_tab2(lambda_index,pitch_index); 

% x = [\dot{p},{p}]
% A = [-d_t/m_t, -k_t/m_t;1 0]
% B = [1/m_t;0]
% Discretize

% A = [-2*zeta_t*omega_t -omega_t^2;1 0];
% B = [1/m_t;0];
% dt = turbine.dt;


Ad = turbine.Ad_tower;
Bd = turbine.Bd_tower;
p_vect = Ad *p_vect +Bd * F_t;
dotp = p_vect(1);
p = p_vect(2);

end   
