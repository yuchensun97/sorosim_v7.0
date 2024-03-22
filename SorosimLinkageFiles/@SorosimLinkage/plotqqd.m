function plotqqd(Tr, t, qqd)
    close all

    PlottingParameters = Tr.PlotParameters;

    tic
    tmax = max(t);
    v = VideoWriter('./Dynamics');
    FrameRate = PlottingParameters.FrameRateValue;
    v.FrameRate = FrameRate;
    open(v);

    hfig      = figure('units','normalized','outerposition',[0 0 1 1]);
    hfig.WindowState = 'maximized';
    set(gca,'CameraPosition',PlottingParameters.CameraPosition,...
            'CameraTarget',PlottingParameters.CameraTarget,...
            'CameraUpVector',PlottingParameters.CameraUpVector,...
            'FontSize',24);
        
    if PlottingParameters.Light
        camlight(PlottingParameters.Az_light,PlottingParameters.El_light)
    end
    %     view(45,45);
    axis equal
    grid on
    hold on
    xlabel('X (m)')
    ylabel('Y (m)')
    zlabel('Z (m)') 
    axis([PlottingParameters.X_lim PlottingParameters.Y_lim PlottingParameters.Z_lim]);
    drawnow

    n_r = Tr.Link.n_r;
    n_l = Tr.Link.n_l;

    qqdtt = interp1(t, qqd, tt);
    q_xi = qqdtt(:,1:Tr.ndof_xi)';
    q_rho = qqdtt(:, Tr.ndof_xi+1:Tr.ndof_rho+Tr.ndof_xi)';

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

    for tt = 0:1/FrameRate:tmax
        delete(findobj('type', 'patch'));
        title(strcat('t= ', num2str(tt)));

        for ii=1:n_l-1
        end

    end

end
