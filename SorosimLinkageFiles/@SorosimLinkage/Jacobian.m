function [J_xi J_rho] = Jacobian(Tr, q_xi, q_rho)

    if isrow(q_xi)
        q_xi = q_xi';
    end

    if isrow(q_rho)
        q_rho = q_rho';
    end

    nsig = Tr.nsig;
    ndof_xi = Tr.ndof_xi;
    ndof_rho = Tr.ndof_rho;

    g_here = Tr.g_base; % g at base
    rho_here = Tr.rho_base; % rho at base

    J_xi = zeros(6*nsig, ndof_xi);
    J_rho = zeros(nsig, ndof_rho);

    J_xi_here = zeros(6, ndof_xi);
    J_rho_here = zeros(1, ndof_rho);

    % joint
    dof_xi_joint = Tr.Twists(1).dof_xi;
    dof_rho_joint = Tr.Twists(1).dof_rho;
    B_xi_joint = Tr.Twists(1).B_xi;
    B_rho_joint = Tr.Twists(1).B_rho;
    q_xi_joint = q_xi(1:dof_xi_joint);
    q_rho_joint = q_rho(1:dof_rho_joint);
    xi_star_joint = Tr.Twists(1).xi_star;
    rho_star_joint = Tr.Twists(1).rho_star;

    if dof_xi_joint == 0
        g_joint = eye(4);
        TgB_joint = zeros(6, ndof_xi);
    else
        xi = B_xi_joint*q_xi_joint + xi_star_joint;
        [g_joint, Tg] = variable_expmap_gTg(xi);
        TgB_joint = zeros(6, ndof_xi);
        TgB_joint(:, 1:dof_xi_joint) = Tg*B_xi_joint;
    end

    if dof_rho_joint == 0
        rho_joint = 1;
        J_rho_joint = zeros(1, ndof_rho);
    else
        rho_joint = B_rho_joint*q_rho_joint + rho_star_joint;
        J_rho_joint = B_rho_joint;
    end

    g_here = g_here*g_joint;
    rho_here = rho_joint;
    
    J_xi_here = dinamico_Adjoint(ginv(g_joint))*...
                (TgB_joint + J_xi_here);
    J_rho_here = J_rho_joint;

end
