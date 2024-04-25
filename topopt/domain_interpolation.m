function interp = domain_interpolation(n,domain)
T = linspace(0,1,n);
    p = 2;
    % define U
    U = zeros(1,n+p+4);
    U(1:p+1) = 0;
    U(end-p:end) = 1;
    idx = 0;
    for i=p+2:length(U)-p-1
        sum = 0;
        for j=0:p-1
            frac = (idx+j)/n;
            sum = sum+frac;
        end
        U(i) = sum/p;
        idx = idx+1;
    end
    M = zeros(n+3);
    Q = zeros(n+3,3);
    for t=3:n+1
        u = (t-2)/n;
        s = findspan(n+3,p,u,U);
        N = basisfun(s,u,p,U);
        M(t,t-1:t-1+p) = N;
    % Sample points from the domain
        pt = nrbeval(domain,u);
        Q(t,:) = pt';
    end
     M(1:2,1:2) = [1 0; -1 1];
    M(end-1:end,end-1:end) = [-1 1; 0 1];
    gamma_0 = (U(p+2)-U(1))/p;
    D0 = 2*pi/n;
    Dn = 2*pi/n;
    gamma_n = (U(end)-U(end-p-1))/p;
    pt = nrbeval(domain,0);
    Q(1,:) = pt';
    pt = nrbeval(domain,1);
    Q(end,:) = pt';
    Q(2,:) = [gamma_0*D0, gamma_0*D0, 0];
    Q(end-1,:) = [gamma_n*Dn, gamma_n*Dn, 0];


    Px = M\Q(:,1);
    Py = M\Q(:,2);
    
    ctrl = [Px'; Py'; zeros(size(Px')); ones(size(Px'))];
    interp = nrbmak(ctrl,U);
end