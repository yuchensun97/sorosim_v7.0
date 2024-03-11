function q = statics(Tr, qu0, u_xi, u_rho)

    n_sact = Tr.n_sact;

    if n_sact == 0
        u_xi = 0;
    else
        sz = size(u_xi);
        if sz(1)~= n_sact || sz(2)~= 1
            error('u_xi must match the size of cable actuators');
        end
        
        if any(u_xi >0)
            error('Cable can only apply negative force');
        end
    end
    
    Func    = @(qu)equilibrium(Tr,qu,u_xi,u_rho);
    options = optimoptions('fsolve','Algorithm','trust-region-dogleg','Display','iter','MaxFunctionEvaluations',2e10);%,'Jacobian','on'); 'trust-region-dogleg' (default), 'trust-region', and 'levenberg-marquardt'.
    tic
    q      = fsolve(Func,qu0,options);
    toc
end
