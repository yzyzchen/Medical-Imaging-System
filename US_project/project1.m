HW%	c			speed of sound (mm/usec)
%	dx			transducer element spacing (mm)
%
% Get needed variables
%
clear all 
clc
load data23;
f0 = 5;		% MHz
fs = 20;	% MHz
c = 1.54;	% mm/us
% 
% Output Array Parameters
%
% Determine the array spacing dx in mm
dx = 2*c/f0;

deltat=1/fs;
[ntime, nelem] = size(Data);		% # time samples, # array elements
disp(sprintf('f0=%g MHz, deltat=%g usec, dx=%g mm', f0, deltat, dx))
disp(sprintf('# of Time Samples=%g,  # of Array Elements=%g',ntime,nelem))

%
% --> QUESTION a. <--
% Make a "wavefield" plot of the raw Data
% Comment this out while you debug other parts of your program
%
showimage3(Data, 4)  % plotting the transpose here...time along x, tranducer x along y.
disp 'hit key', pause
