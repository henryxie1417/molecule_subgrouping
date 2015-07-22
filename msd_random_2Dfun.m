function F = msd_random_2Dfun(x,xdata)

% fit function for random, free, or brownian difusion

offset = x(1);
D      = x(2);

F = offset + 4*D*xdata;

