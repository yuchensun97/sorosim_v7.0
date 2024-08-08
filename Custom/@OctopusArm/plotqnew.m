function plotqnew(Tr, t, qqd, stamps, filename)
    % plot the arm posture at different time steps
    close all

    if nargin<=4
        filename = 'snapshots';
    end

    % get the number of time steps
    PlottingParameters = Tr.PlotParameters;
    FrameRate = PlottingParameters.FrameRateValue;

    tmax = max(t);
    lt = length(0:1/FrameRate:tmax);
    ptip = NaN(3, lt);

    if PlottingParameters.ClosePrevious
        close all
    end
    fh = figure(1);
    fh.Units = 'normalized';
    fh.OuterPosition = [0 0 1 1];

    set(gca, 'CameraPosition', PlottingParameters.CameraPosition,...
        'CameraTarget', PlottingParameters.CameraTarget,...
        'CameraUpVector', PlottingParameters.CameraUpVector,...
        'FontSize', 28, 'FontName', 'Times New Roman');
 
     if PlottingParameters.Light
         camlight(PlottingParameters.Az_light, PlottingParameters.El_light)
     end
     view(30, 40);
     axis equal
     grid on
     hold on
     xlabel('$S$ (m)', 'Interpreter', 'latex', 'FontSize', 28);
     ylabel('$X_1$ (m)', 'Interpreter', 'latex', 'FontSize', 28);
     zlabel('$X_2$ (m)', 'Interpreter', 'latex', 'FontSize', 28);

     axis([PlottingParameters.X_lim PlottingParameters.Y_lim PlottingParameters.Z_lim])

     n_r = Tr.Link.n_r;
     n_l = Tr.Link.n_l;
 
     dof_xi = Tr.Twists(2).dof_xi;
     dof_rho = Tr.Twists(2).dof_rho;
     rho_starfn = Tr.Twists(2).rho_starfn;
     Bh_xi = Tr.Twists(2).Bh_xi;
     Bh_rho = Tr.Twists(2).Bh_rho;
     B_xi_dof = Tr.Twists(2).B_xi_dof;
     B_xi_odr = Tr.Twists(2).B_xi_odr;
     B_rho_dof = Tr.Twists(2).B_rho_dof;
     B_rho_odr = Tr.Twists(2).B_rho_odr;
     L = Tr.Link.L;
     color = Tr.Link.color;
 
     Xs = linspace(0, 1, n_l);
     H = Xs(2) - Xs(1);
 
     Z = (1/2) * H; % Zanna quadrature coefficients
 
     r_fn = Tr.Link.r_fn;

     k = 1; % index for the time stamps
     it = 1; % index for the tip position

    for tt=0:1/FrameRate:tmax
        % delete(findobj('type', 'patch'));
        qqdtt = interp1(t, qqd, tt);
        q_xi = qqdtt(:, 1:Tr.ndof_xi)';
        q_rho = qqdtt(:, Tr.ndof_xi+1:Tr.ndof_rho+Tr.ndof_xi)';

        Xpatch = zeros(n_r, (n_r-1)*(n_l-2)+2);
        Ypatch = zeros(n_r, (n_r-1)*(n_l-2)+2);
        Zpatch = zeros(n_r, (n_r-1)*(n_l-2)+2);
        i_patch = 1;

        r = r_fn(0);
        theta = linspace(0, 2*pi, n_r);
        if dof_rho ~= 0
            rho_here = Bh_rho(Xs(1), B_rho_dof, B_rho_odr)*q_rho + rho_starfn(Xs(1));
        else
            rho_here = 1;
        end
        x = zeros(1, n_r);
        y = rho_here*r*cos(theta);
        z = rho_here*r*sin(theta);
        pos = [x; y; z; ones(1, n_r)];

        g_here = eye(4);

        pos_here = g_here*pos;
        x_here = pos_here(1, :);
        y_here = pos_here(2, :);
        z_here = pos_here(3, :);

        Xpatch(:, i_patch) = x_here';
        Ypatch(:, i_patch) = y_here';
        Zpatch(:, i_patch) = z_here';
        i_patch = i_patch + 1;

        x_pre = x_here;
        y_pre = y_here;
        z_pre = z_here;

        for ii=1:n_l-1
            r = r_fn(Xs(ii+1));
            theta = linspace(0, 2*pi, n_r);
            if dof_rho ~= 0
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

            xi_Zhere = [0 0 0 1 0 0]';
            if dof_xi ~= 0
                xi_Zhere = Bh_xi(X_Z, B_xi_dof, B_xi_odr)*q_xi + xi_Zhere;
            end

            Gamma_here = H * L * xi_Zhere;
            gh = variable_expmap_g(Gamma_here);
            g_here = g_here*gh;

            g_rho = [rho_here * g_here(1:3, 1:3), g_here(1:3, 4);
                     0 0 0 1];
            pos_here = g_rho*pos;
            x_here = pos_here(1, :);
            y_here = pos_here(2, :);
            z_here = pos_here(3, :);

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
        ptip(:, it) = g_here(1:3, 4);
        it = it+1;

        Xpatch(:, i_patch) = x_here';
        Ypatch(:, i_patch) = y_here';
        Zpatch(:, i_patch) = z_here';

        if k <= length(stamps) && abs(tt-stamps(k))<1/(FrameRate*2)
            alpha = 0.3;
            if k== length(stamps)
                alpha = 1;
            end
            patch(Xpatch, Ypatch, Zpatch, color, 'EdgeColor', 'none', 'FaceAlpha', alpha);
            k = k + 1;
        end
    end

    % plot the tip trajectory
    plot3(ptip(1, :), ptip(2, :), ptip(3, :), 'r', 'LineWidth', 2)

    % save the figure
    exportgraphics(fh, filename, 'ContentType', 'vector')
end
