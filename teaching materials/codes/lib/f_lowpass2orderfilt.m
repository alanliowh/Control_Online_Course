% Second-order low-pass filter
%  [y,filt] = f_lowpass2orderfilt(x,dt,filt) performs a second-order
%  low-pass filter on signal input x and returns a filtered value y. The
%  filter cut-off frequency is set as filter.f0 in Hz and damping ratio as
%  filter.zeta. Sampling period is denoted as dt. 
%
%  Example:
%  clear;
%  x = cumsum(randn(100,1));
%  filt.f0 = 0.4; % [Hz]
%  filt.zeta = 0.7;
%  dt = 0.1;
%  for i = 1:length(x)
%   [y(i),filt] = f_lowpass2orderfilt(x(i),filt,dt);
%  end
%  plot(x); hold on; plot(y);
%
%  Written by Alan Wai Hou Lio (wali@dtu.dk) 
%  Date: 17.12.2020

function [y,filt] = f_lowpass2orderfilt(x,filt,dt)


if ~isfield(filt,'x1')
    filt.x1 = x;
    filt.x2 = x;
    filt.x1_old = filt.x1;
    filt.x2_old = filt.x2;
    filt.y1 = x;
    filt.y2 = x;
    filt.y1_old = filt.y1;
    filt.y2_old = filt.y2;
    y = x;
    
else
    filt.x1_old = filt.x1;
    filt.x2_old = filt.x2;
    filt.y1_old = filt.y1;
    filt.y2_old = filt.y2;
    
    f0 = filt.f0; % [Hz]
    zeta = filt.zeta; % [-]
    denom = 3 + 6* zeta*pi*f0*dt + 4*pi^2*f0^2*dt^2;
    a1 = (6-4*pi^2*f0^2*dt^2)/denom;
    a2 = (-3 +6*zeta*pi*f0*dt-4*pi^2*f0^2*dt^2)/denom;
    b0 = 4*pi^2*f0^2*dt^2/denom;
    b1 = b0;
    b2 = b0;
    y = a1*filt.y1_old + a2*filt.y2_old + b0*x + b1*filt.x1_old + b2*filt.x2_old;
end
filt.x2 = filt.x1_old;
filt.x1 = x;
filt.y2 = filt.y1_old;
filt.y1 = y;

filt.dydt  = (y - filt.y1_old)/dt;
end
    
    
