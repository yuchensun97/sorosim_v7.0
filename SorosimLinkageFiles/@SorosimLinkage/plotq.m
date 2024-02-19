function fh = plotq(Tr, q_xi, q_rho)
    if nargin==1
        q_xi = zeros(Tr.ndof_xi, 1);
        q_rho = zeros(Tr.ndof_rho, 1);
    end

    if isrow(q_xi)
        q_xi = q_xi';
    end

    if isrow(q_rho)
        q_rho = q_rho';
    end

    PlottingParameters = Tr.PlotParameters;

    if PlottingParameters.ClosePrevious
        close all
    end

    %Plot options
    fh = figure(1);
    fh.Units = 'normalized';
    fh.OuterPosition = [0 0 1 1];

    set(gca, 'CameraPosition', PlottingParameters.CameraPosition,...
       'CameraTarget', PlottingParameters.CameraTarget,...
       'CameraUpVector', PlottingParameters.CameraUpVector,...
       'FontSize', 18)

    if PlottingParameters.Light
        camlight(PlottingParameters.Az_light, PlottingParameters.El_light)
    end

    axis equal
    grid on
    hold on
    xlabel('S (m)')
    ylabel('X1 (m)')
    zlabel('X2 (m)')

    axis([PlottingParameters.X_lim PlottingParameters.Y_lim PlottingParameters.Z_lim])

    %%Fwd Kinematics

    % joint
    dof_xi_joint = Tr.Twists(1).dof_xi;
    dof_rho_joint = Tr.Twists(1).dof_rho;
    q_xi_joint = q_xi(1:dof_xi_joint);
    q_rho_joint = q_rho(1:dof_rho_joint);
    B_xi_joint = Tr.Twists(1).B_xi;
    B_rho_joint = Tr.Twists(1).B_rho;
    xi_star_joint = Tr.Twists(1).xi_star;
    rho_star_joint = Tr.Twists(1).rho_star;

    g_here = Tr.g_base;
    rho_here = Tr.rho_base;

    if dof_xi_joint == 0
        g_joint = eye(4);
    else
        xi_joint = B_xi_joint*q_xi_joint + xi_star_joint;
        g_joint = variable_expmap_ga(xi_joint);
    end

    g_here = g_here*g_joint;

    if dof_rho_joint == 0
        rho_joint = 1;
    else
        rho_joint = B_rho_joint*q_rho_joint + rho_star_joint;
    end

    rho_here = rho_joint;

    n_r = Tr.Link.n_r;
    n_l = Tr.Link.n_l;
    color = Tr.Link.color;

    % soft body starts here
    q_xi = q_xi(dof_xi_joint+1:end);
    q_rho = q_rho(dof_rho_joint+1:end);
    dof_xi = Tr.Twists(2).dof_xi;
    dof_rho = Tr.Twists(2).dof_rho;
    xi_starfn = Tr.Twists(2).xi_starfn;
    rho_starfn = Tr.Twists(2).rho_starfn;
    Bh_xi = Tr.Twists(2).Bh_xi;
    Bh_rho = Tr.Twists(2).Bh_rho;
    B_xi_dof = Tr.Twists(2).B_xi_dof;
    B_xi_odr = Tr.Twists(2).B_xi_odr;
    B_rho_dof = Tr.Twists(2).B_rho_dof;
    B_rho_odr = Tr.Twists(2).B_rho_odr;
    L = Tr.Link.L;

    Xs = linspace(0, 1, n_l);
    H = Xs(2) - Xs(1);

    Z = (1/2) * H; % Zanna quadrature coefficients

    Xpatch = zeros(n_r, (n_r-1)*(n_l-2)+2);
    Ypatch = zeros(n_r, (n_r-1)*(n_l-2)+2);
    Zpatch = zeros(n_r, (n_r-1)*(n_l-2)+2);
    i_patch = 1;

    r_fn = Tr.Link.r_fn;
    r = r_fn(0);
    theta = linspace(0, 2*pi, n_r);
    x = zeros(1, n_r);
    y = rho_here*r*cos(theta);
    z = rho_here*r*sin(theta);
    pos = [x; y; z;ones(1, n_r)];

    pos_here = g_here*pos;
    x_here = pos_here(1, :);
    y_here = pos_here(2, :);
    z_here = pos_here(3, :);
    plot3(x_here, y_here, z_here)
    hold on

    Xpatch(:, i_patch) = x_here';
    Ypatch(:, i_patch) = y_here';
    Zpatch(:, i_patch) = z_here';
    i_patch = i_patch + 1;

    x_pre = x_here;
    y_pre = y_here;
    z_pre = z_here;

    Lscale = L;

    for ii=1:n_l-1
        r = r_fn(Xs(ii+1));
        theta = linspace(0, 2*pi, n_r);
        if dof_rho~=0
            rho_here = Bh_rho(Xs(ii+1), B_rho_dof, B_rho_odr)*q_rho + rho_starfn(Xs(ii+1));
        else
            rho_here = 1;
        end
        x = zeros(1, n_r);
        y = r*cos(theta);
        z = r*sin(theta);
        pos = [x; y; z;ones(1, n_r)];

        X = Xs(ii);
        X_Z = X + Z;

        xi_Zhere = xi_starfn(X_Z);
        xi_Zhere(1:3) = xi_Zhere(1:3) * L;

        if dof_xi~=0
            xi_Zhere = Bh_xi(X_Z, B_xi_dof, B_xi_odr)*q_xi + xi_Zhere;
        end
        Gamma_here = H * xi_Zhere;
        Gamma_here(4:6) = Gamma_here(4:6) * Lscale;
        gh = variable_expmap_g(Gamma_here);
        g_here = g_here*gh;
        
        g_rho = [rho_here * g_here(1:3,1:3), g_here(1:3,4);
                0 0 0 1];
        pos_here = g_rho*pos;
        x_here = pos_here(1, :);
        y_here = pos_here(2, :);
        z_here = pos_here(3, :);

        % for debug only
        plot3(x_here, y_here, z_here)
        hold on
        % plot soft link here
        for jj=1:n_r-1
            Xpatch(1:5,i_patch)   = [x_pre(jj) x_here(jj) x_here(jj+1) x_pre(jj+1) x_pre(jj)]';
            Xpatch(6:end,i_patch) = x_pre(jj)*ones(n_r-5,1);
            Ypatch(1:5,i_patch)   = [y_pre(jj) y_here(jj) y_here(jj+1) y_pre(jj+1) y_pre(jj)]';
            Ypatch(6:end,i_patch) = y_pre(jj)*ones(n_r-5,1);
            Zpatch(1:5,i_patch)   = [z_pre(jj) z_here(jj) z_here(jj+1) z_pre(jj+1) z_pre(jj)]';
            Zpatch(6:end,i_patch) = z_pre(jj)*ones(n_r-5,1);
            i_patch = i_patch+1;
        end
        x_pre = x_here;
        y_pre = y_here;
        z_pre = z_here;
    end

    Xpatch(:, i_patch) = x_here';
    Ypatch(:, i_patch) = y_here';
    Zpatch(:, i_patch) = z_here';

    patch(Xpatch, Ypatch, Zpatch, color, 'EdgeColor', 'none');
    hold off
    % drawnow

end
