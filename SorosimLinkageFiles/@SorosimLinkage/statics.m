function q = statics(Tr, qu0)
    Func    = @(qu)equilibrium(Tr,qu,0,0);
    options = optimoptions('fsolve','Algorithm','trust-region-dogleg','Display','iter','MaxFunctionEvaluations',2e10);%,'Jacobian','on'); 'trust-region-dogleg' (default), 'trust-region', and 'levenberg-marquardt'.
    tic
    q      = fsolve(Func,qu0,options);
    toc
end