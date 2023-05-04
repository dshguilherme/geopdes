%  This is the file three.m
%  which calculates function values and derivatives
%  for the "three bar truss problem" defined in the
%  file threemain.m.
%
function [f0val,df0dx,df0dx2,fval,dfdx,dfdx2] = three(x);
%
%  Written in May 1999 by
%  Krister Svanberg <krille@math.kth.se>
%  Department of Mathematics
%  SE-10044 Stockholm, Sweden.
%
e = [1 1 1]';
f0val = e'*x;
df0dx = e;
df0dx2 = 0*e;

D = diag(x);
sq2 = 1/sqrt(2);
R = [ 1 sq2 0
      0 sq2 1 ];

p1 =  [1 0]';
p2 =  [1 1]';
p3 =  [0 1]';

K = R*D*R';
u1 = K\p1;
u2 = K\p2;
u3 = K\p3;

compl1 = p1'*u1 - 1;
compl2 = p2'*u2 - 1;
compl3 = p3'*u3 - 1;

fval = [compl1 compl2 compl3]';

rtu1 = R'*u1;
rtu2 = R'*u2;
rtu3 = R'*u3;

dcompl1 = -rtu1.*rtu1;
dcompl2 = -rtu2.*rtu2;
dcompl3 = -rtu3.*rtu3;

dfdx = [dcompl1 dcompl2 dcompl3]';

% We do not want to calculate second derivatives,
% so we simply set them to zero

dfdx2 = 0*dfdx;
