degree = 3;
nsub = [30 10];
vol_frac = 0.49;
rmin = 2;
f = 1e-5;
omega = 2*pi*f;
method = 'MMA';
change_max = 0.01;
max_iter = 1000;

dynamic_cantilever(degree, nsub, vol_frac, rmin, omega, method, change_max, max_iter)