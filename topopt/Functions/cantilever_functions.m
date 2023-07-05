function [f0val, fval] = cantilever_functions(xval)
%% Filtering
x = xval;
load('init.mat');
x = apply_x_filter(filter_options, x);

%% Assembly
u = zeros(sp.ndof,1);
Ks = zeros(sp.ndof);
if ~strcmp(objective_function,"compliance") % if the objective is AIP or mixed
    M = Ks;
    for e=1:msh.nel
        k_e = (YOUNG_MIN +(x(e)^3)*(YOUNG - YOUNG_MIN))*squeeze(Ke(e,:,:));
        if x(e) > 0.1
            m_e = (RHO_MIN +x(e)*(RHO -RHO_MIN))*squeeze(Me(e,:,:));
        elseif x(e) <= 0.1
            m_e = (RHO_MIN +(x(e)^9)*(RHO -RHO_MIN))*squeeze(Me(e,:,:));
        end
        idx = lm(e,:)';
        Ks(idx,idx) = Ks(idx,idx) +k_e;
        M(idx,idx) = M(idx,idx) +m_e;
    end
    Ks = sparse(Ks); M = sparse(M);
    C = alpha_*M +beta_*Ks;
    Kd = Ks +1j*omega*C -omega*omega*M;
else % if the objective is compliance
    for e=1:msh.nel
    k_e = (YOUNG_MIN +(x(e)^3)*(YOUNG - YOUNG_MIN))*squeeze(Ke(e,:,:));
    idx = lm(e,:)';
    Ks(idx,idx) = Ks(idx,idx) +k_e;
    end
    Ks = sparse(Ks);
end


%% Solve problem
switch objective_function
    case "compliance"
        u(dr_dofs) = 0;
        F(free_dofs) = F(free_dofs) -Ks(free_dofs, dr_dofs)*u(dr_dofs);
        u(free_dofs) = Ks(free_dofs,free_dofs)\F(free_dofs);
    case "AIP"
        u(dr_dofs) = 0;
        F(free_dofs) = F(free_dofs) -Kd(free_dofs, dr_dofs)*u(dr_dofs);
        u(free_dofs) = Kd(free_dofs,free_dofs)\F(free_dofs);
    otherwise
        u(dr_dofs) = 0;
        us = u;
        Fs = F;
        
        F(free_dofs) = F(free_dofs) -Kd(free_dofs, dr_dofs)*u(dr_dofs);
        Fs(free_dofs) = Fs(free_dofs) -Ks(free_dofs, dr_dofs)*us(dr_dofs);
        
        u(free_dofs) = Kd(free_dofs,free_dofs)\F(free_dofs);
        us(free_dofs) = Ks(free_dofs,free_dofs)\Fs(free_dofs);
end

%% Objective Functions and Constraints
switch objective_function
    case "compliance"
        f0val = F'*u; % Static Compliance
    case "scaled compliance"
        Cs = Fs'*us;
        f0val = 100*Cs/Cs0;
    case "dB compliance"
        Cs = Fs'*us;
        Cs_db = 100 +10*log10(Cs);
        Cs0_db = 100 +10*log10(Cs0);
        f0val = 100*Cs_db/Cs0_db;
    case "AIP"
        W = 0.5*omega*omega*real((u')*C*u); % Active input power
        f0val = W;
    case "dB AIP"
        W = 0.5*omega*omega*real((u')*C*u); % Active input power
        W_db = 100 +10*log10(W);
        W0_db = 100 +10*log10(W0);
        f0val = 100*(W_db/W0_db);
    case "mixed"
        Cs = Fs'*us;
%         Cs_db = 100 +10*log10(Cs);
%         Cs0_db = 100 +10*log10(Cs0);
%         Cs_scaled = 100*Cs_db/Cs0_db;
        Cs_scaled = 100*Cs/Cs0;
        
        W = 0.5*omega*omega*real((u')*C*u); % Active input power
        W_db = 100 +10*log10(W);
        W0_db = 100 +10*log10(W0);
        W_scaled = 100*W_db/W0_db;
        
        f0val = neta*W_scaled +(1-neta)*Cs_scaled;
    case "Cs0"
        f0val = Fs'*us; % Static Compliance
    case "W0"
        W = 0.5*omega*omega*real((u')*C*u); % Active input power
        f0val = 100 +10*log10(W); % dB scale    
    case "History"
        W = 0.5*omega*omega*real((u')*C*u);
        W_db = 100 +10*log10(W);
        W0_db = 100 +10*log10(W0);
        W_scaled = 100*W_db/W0_db;
        Cs = Fs'*us;
        f0val = [W_scaled, 100*Cs/Cs0];

end
fval = sum(x.*Ve)-vol_frac*sum(Ve); % Volume Constraint

end