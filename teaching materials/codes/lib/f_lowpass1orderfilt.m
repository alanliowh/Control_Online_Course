% First-order low-pass filter
%  [y,filt] = f_lowpass1orderfilt(x, filt, dt) performs real-time 1st order
%  low-pass filtering. The input signal x is filterd and returns an output
%  y. The time-constant of the filter is set as filt.tau. Sampling period
%  is dt.
%
%  Example:
%
%  clear;
%  x = cumsum(randn(100,1));
%  fpass = 0.4; %[Hz]
%  filt.tau = 1/(2*pi*fpass);
%  dt = 0.1;
%  for i = 1:length(x)
%      [y(i),filt] = f_lowpass1orderfilt(x(i),filt,dt);
%  end
%  plot(x); hold on; plot(y);
%
%  Reference:
%  http://techteach.no/simview/lowpass_filter/doc/filter_algorithm.pdf
%
%  Written by Alan Wai Hou Lio (wali@dtu.dk) 
%  Date: 17-12-2020

function [y,filt] = f_lowpass1orderfilt(x, filt, dt)

if ~isfield(filt,'x1_old')
    filt.x1_old = x;
    filt.y1_old = x;
    y = x;
else
    filt.x1_old = filt.x1;
    filt.y1_old = filt.y1;

    tau = filt.tau;
    a1 = (2*tau - dt) / (2*tau +dt);
    b0 = dt/(2*tau +dt);
    b1 = b0;
    y = a1*filt.y1_old + b0*x + b1* filt.x1_old;
end
filt.x1 = x;
filt.y1 = y;

end