function [x_new, low, upp] = updateMMA(x, xold1, xold2, loop, compliance, ...
    restriction, dc, d2c, dv, vol_frac, low, upp)
m = 1;
n = length(x(:));
xmin = zeros(n,1);
xmax = ones(n,1);
a0 = 1;
a = zeros(m,1);
c_MMA = 10000*ones(m,1);
d = zeros(m,1);
f0val = compliance;
fval = restriction;
df0dx = dc(:);
dfdx = dv;

[x_new, ~, ~, ~, ~, ~, ~, ~, ~, low, upp] = ...
mmasub(m, n, loop, x(:), xmin, xmax, xold1, xold2, ...
f0val, df0dx, d2c, fval, dfdx, 0*dfdx,low,upp,a0,a,c_MMA,d);

end