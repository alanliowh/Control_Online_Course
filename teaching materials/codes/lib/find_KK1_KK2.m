%% find dQdtheta0, KK1 and KK2

function [dQdtheta0,KK1,KK2]=find_KK1_KK2(pitch_ss,dQdtheta)
% dQdtheta = dQdthata|0 (1 + pitch/KK1 + pitch^2 /KK2 )
% dQtheta = dQdthata|0 + alpha * pitch + beta * pitch^2
%
% Perform least sqaure
% Ax = y   where y = dQdtheta,  A = [1 pitch pitch^2], x = [dQdetha|0 alpha beta]
y = dQdtheta';
A = [pitch_ss'.^0 pitch_ss'.^1 pitch_ss'.^2];
x = A\y; % 

dQdtheta0 = x(1);
alpha = x(2);
beta = x(3);

KK1 = dQdtheta0/alpha;
KK2 = dQdtheta0/beta;