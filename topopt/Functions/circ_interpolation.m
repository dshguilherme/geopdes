function circ = circ_interpolation(p)
assert(p>1 && p<=5,'p must be between 2 and 5');
lookup = [20, 15, 10, 9]; % 1e-4 error
n = lookup(p-1) +1;

T = linspace(0,1,n+1);
U = zeros(1,n+p+4);
U(1:p+1) = 0;
U(end-p:end) = 1;
for i=p+2:n+3
    sum = 0;
    for j=i:i+p-1
        sum = sum +T(j-p-1);
    end
    U(i) = sum/p;
end
M = zeros(n+3);
Q = zeros(n+3,2);
for t=1+2:n-1+2
    u = T(t-1);
    s = findspan(n+p,p,u,U);
    N = basisfun(s,u,p,U);
    M(t,t-1:t-1+p) = N;
    theta = 2*pi*u;
    Q(t,:) = [cos(theta), sin(theta)];
end
M(1,1) = 1;
M(2,1:2) = [-1 1];
M(end-1,end-1:end) = [-1 1];
M(end,end) = 1;
gamma_0 = (U(p+2)-U(1))/p;
D0 = 2*pi/(n+1);
Dn = 2*pi/(n+1);
gamma_n = (U(end)-U(end-p-1))/p;
Q(1,:) = [cos(0), sin(0)];
Q(end,:) = Q(1,:);
Q(2,:) = [gamma_0*D0, gamma_0*D0];
Q(end-1,:) = [gamma_n*Dn, gamma_n*Dn];


Px = M\Q(:,1);
Py = M\Q(:,2);


ctrl = [Px'; Py'; zeros(size(Px')); ones(size(Px'))];
circ = nrbmak(ctrl,U);
% figure
% nrbctrlplot(circ);
% figure
% nrbkntplot(circ);
end
% nrbplot(circ,100)
