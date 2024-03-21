%Function for the dynamic simulation of the linkage
function [t, qqd] = dynamics(Tr, qqd0, odetype, dt, tmax)
    ndof_xi = Tr.ndof_xi;
    ndof_rho = Tr.ndof_rho;
    odetype = 'ode45';
    dt = 0.01;
    tmax = 2;
    if nargin==1||isempty(qqd0)
        q0 = zeros(ndof_xi+ndof_rho,1);
        qd0 = zeros(ndof_xi+ndof_rho,1);
        qqd0 = [q0;qd0];
    end

    uqt_xi = 0;
    uqt_rho = 0;

    switch odetype
        case 'ode45'
            options = odeset('RelTol',1e-3);
            tic
            [t, qqd] = ode45(@(t,qqd)Tr.derivatives(t, qqd,uqt_xi, uqt_rho),0:dt:tmax,qqd0,options);
            toc
        case 'ode23'
            options = odeset('RelTol',1e-3);
            tic
            [t, qqd] = ode23(@(t,qqd)Tr.derivatives(t, qqd,uqt_xi, uqt_rho),0:dt:tmax,qqd0,options);
            toc
        case 'ode113'
            options = odeset('RelTol',1e-3,'AbsTol',1e-6);
            tic
            [t, qqd] = ode113(@(t,qqd)Tr.derivatives(t, qqd,uqt_xi, uqt_rho),0:dt:tmax,qqd0,options);
            toc
        case 'ode23s'
            options = odeset('RelTol',1e-3,'AbsTol',1e-6);
            tic
            [t, qqd] = ode23s(@(t,qqd)Tr.derivatives(t, qqd,uqt_xi, uqt_rho),0:dt:tmax,qqd0,options);
            toc
        case 'ode15s'
            options = odeset('RelTol',1e-3,'AbsTol',1e-6);
            tic
            [t, qqd] = ode15s(@(t,qqd)Tr.derivatives(t, qqd,uqt_xi, uqt_rho),0:dt:tmax,qqd0,options);
            toc
        case 'ode1'
            tic
            qqd = ode1(@(t,qqd)Tr.derivatives(t,qqd,uqt_xi,uqt_rho),0:dt:tmax, qqd0);
            toc
            t=0:dt:tmax;
        case 'ode2'
            tic
            qqd = ode1(@(t,qqd)Tr.derivatives(t,qqd,uqt_xi,uqt_rho),0:dt:tmax, qqd0);
            toc
            t=0:dt:tmax;
        otherwise
            error("Error: Wrong ode type")
    end
end
