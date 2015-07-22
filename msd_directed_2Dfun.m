function F = msd_directed_2Dfun(x,xdata)

% fit function for directed diffusion
 
offset = x(1);
D      = x(2);
v      = x(3);

F = offset + 4*D*xdata + v^2*(xdata.^2);


