function F = msd_caged_2Dfun(x,xdata)

% fit function for caged or confined diffusion in 2D
%
% Kusumi et al. 200x

offset = x(1);
D0     = x(2);
L      = x(3);

F = offset+(L^2/3*(1-exp(-12*D0*xdata/(L^2))));