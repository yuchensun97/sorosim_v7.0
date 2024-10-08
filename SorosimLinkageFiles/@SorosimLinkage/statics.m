function q = statics(Tr, qu0, u_xi, u_rho)

    n_sact = Tr.n_sact;
    n_ract = Tr.n_ract;

    if n_sact == 0
        u_xi = 0;
    else
        sz = size(u_xi);
        if sz(1)~= n_sact || sz(2)~= 1
            error('u_xi must match the size of cable actuators');
        end
    end

    if n_ract == 0
        u_rho = 0;
    else
        sz = size(u_rho);
        if sz(1) ~= n_ract || sz(2)~=1
            error('u_rho must match the size of radial actuators');
        end
    end
    
    Func    = @(qu)equilibrium(Tr,qu,u_xi,u_rho);
    options = optimoptions('fsolve','Algorithm','trust-region-dogleg','Display','iter','MaxFunctionEvaluations',2e10);
    tic
    q      = fsolve(Func,qu0,options);
    toc
end
