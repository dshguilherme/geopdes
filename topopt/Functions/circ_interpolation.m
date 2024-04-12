function circ = circ_interpolation(n)
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

    % U = [0, 0, 0, 0.041667, 0.125, 0.208333, 0.291667, 0.375, 0.458333, 0.541667, 0.625, 0.70833, 0.791667, 0.8775, 0.958333, 1, 1, 1];
    M = zeros(n+3);
    Q = zeros(n+3,2);
    for t=3:n+1
        u = (t-2)/n;
        s = findspan(length(U)-p,p,u,U);
        N = basisfun(s,u,p,U);
        M(t,t-1:t-1+p) = N;
    % Sample points from the circle
        theta = 2*pi*u;
        Q(t,:) = [cos(theta), sin(theta)];
    end
    M(1:2,1:2) = [1 0; -1 1];
    M(end-1:end,end-1:end) = [-1 1; 0 1];
    gamma_0 = (U(p+2)-U(1))/p;
    D0 = 2*pi/n;
    Dn = 2*pi/n;
    gamma_n = (U(end)-U(end-p-1))/p;
    Q(1,:) = [cos(0), sin(0)];
    Q(end,:) = Q(1,:);
    Q(2,:) = [gamma_0*D0, gamma_0*D0];
    Q(end-1,:) = [gamma_n*Dn, gamma_n*Dn];


    Px = M\Q(:,1);
    Py = M\Q(:,2);


    ctrl = [Px'; Py'; zeros(size(Px')); ones(size(Px'))];
    circ = nrbmak(ctrl,U);
end
% nrbplot(circ,100)
