% Cpt_curve2 Building a higher resolution Cp curve
%   [cpt_out,lam_out,th_out] = Cpt_curve2(cptable_in) takes cptable_in from BEM
%   code as an input cptable(theta,lambda) and returns a higher resolution output cpt_out as
%   cp(lambda,theta).
%
%   [cpt_out,lam_out,th_out] = Cpt_curve2(cptable_in,n) returns a higher
%   resolution table Ct if n = 'Ct'.
%
%   Written by Alan Wai Hou Lio (09-2019)
%   log:
%   version (19.09.30)
%   - the input format is cp(theta,lambda), output format is
%   cp(lambda,theta)
%   version (19.09.20)
%   - simplify the process of removing the negative part of Cp using max().
%   - the code uses the .mat files that are consistent to the BEM code.

function [Cpt_tab2,lambda_tab2,pitch_tab2]=Cpt_curve2(cptable,n)

pitch_tab = [-10;cptable.th];
lambda_tab = [-2;cptable.lam];
[cptable.pitch_grid,cptable.lambda_grid] = meshgrid(pitch_tab,lambda_tab);

if nargin<2
   n = 'Cp';
end
   

switch n
    case 'Cp'
    cptable.cp = cptable.cp';
    Cpt_tab = cptable.cp;
    case 'Ct'
    cptable.ct = cptable.ct';
    Cpt_tab = cptable.ct;
end
%%% ensuring Cp surface is not zero when lambda is 0.
Cpt_tab = [zeros(length(lambda_tab)-1,1),Cpt_tab];
Cpt_tab = [zeros(1,length(pitch_tab));Cpt_tab];

%% build a higher resolution Cp curve
pitch_tab2 = [-10:0.01:30]; % user input
lambda_tab2 = [0:0.1:22.5]; % user input


[pitch_grid2,lambda_grid2] = meshgrid(pitch_tab2,lambda_tab2);
pitch_tab2 = pitch_grid2(1,:)';
lambda_tab2 = lambda_grid2(:,1);
switch n
    case 'Cp'
Cpt_tab2 = interp2(cptable.pitch_grid,cptable.lambda_grid,Cpt_tab,pitch_grid2,lambda_grid2);
    case 'Ct'
Cpt_tab2 = interp2(cptable.pitch_grid,cptable.lambda_grid,Cpt_tab,pitch_grid2,lambda_grid2);
end

Cpt_tab2 = max(Cpt_tab2,0);


end  


