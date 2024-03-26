%Function for the dynamic simulation of the linkage
%TODO: provide default parameters using input parser.
function [t, qqd] = dynamics(Tr, qqd0, uqt_xi, uqt_rho, odetype, dt, tmax)
    ndof_xi = Tr.ndof_xi;
    ndof_rho = Tr.ndof_rho;
    n_sact = Tr.n_sact;
    n_ract = Tr.n_ract;

    if nargin==1 || isempty(qqd0)
        q0 = zeros(ndof_xi+ndof_rho,1);
        qd0 = zeros(ndof_xi+ndof_rho,1);
        qqd0 = [q0;qd0];
        odetype = 'ode15s';
        n_sact = 0;
        n_ract = 0;
        dt = 0.01;
        tmax = 10;
    end

    if n_sact == 0
        uqt_xi = @(t)0;
    else
        if ~isa(uqt_xi, 'function_handle')
            error('uqt_xi should be a function handle');
        end
        u_xi = uqt_xi(0);
        sz = size(u_xi);
        if sz(1)~=n_sact || sz(2)~=1
            error('uqt_xi must match the size of cable actuator');
        end
    end

    if n_ract == 0
        uqt_rho = @(t)0;
    else
        if ~isa(uqt_rho, 'function_handle')
            error('uqt_rho should be a function handle');
        end
        u_rho = uqt_rho(0);
        sz = size(u_rho);
        if sz(1)~=n_ract || sz(2)~=1
            error('uqt_rho must match the size of cable actuator');
        end
    end

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
