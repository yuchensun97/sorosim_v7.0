function [J_xi, J_rho] = Jacobian(Tr, q_xi)
%Jacobian: iterative Jacobian of the soft link
%   Detailed explanation goes here
%   Tr      : SorosimLinkage class
%   q_xi    : general coordinates for strains
%   returns:
%   J_xi    : 6nip x ndof_xi Jacobian for strain variables
%   J_rho   : nip x ndof_rho Jacobian for inflation ratio

    if isrow(q_xi)
        q_xi = q_xi';
    end

    nsig = Tr.nsig;
    ndof_xi = Tr.ndof_xi;
    ndof_rho = Tr.ndof_rho;

    g_here = Tr.g_base; % g at base

    J_xi = zeros(6*nsig, ndof_xi);
    J_rho = zeros(nsig, ndof_rho);

    J_xi_here = zeros(6, ndof_xi);

    % joint
    dof_xi_joint = Tr.Twists(1).dof_xi;
    B_xi_joint = Tr.Twists(1).B_xi;
    q_xi_joint = q_xi(1:dof_xi_joint);
    xi_star_joint = Tr.Twists(1).xi_star;

    dof_rho_joint = Tr.Twists(1).dof_rho;
    B_rho_joint = Tr.Twists(1).B_rho;

    if dof_xi_joint == 0
        g_joint = eye(4);
        TgB_joint = zeros(6, ndof_xi);
    else
        xi = B_xi_joint*q_xi_joint + xi_star_joint;
        [g_joint, Tg] = variable_expmap_gTg(xi);
        TgB_joint = zeros(6, ndof_xi);
        TgB_joint(:, 1:dof_xi_joint) = Tg*B_xi_joint;
    end

    if dof_rho_joint > 0
        J_rho(:, 1:dof_rho_joint) = B_rho_joint;
    end

    g_here = g_here*g_joint;
    
    J_xi_here = dinamico_Adjoint(ginv(g_joint))*...
                (TgB_joint + J_xi_here);

    ndof_xi = ndof_xi - dof_xi_joint;

    % soft body
    xi_star = Tr.Twists(2).xi_star; % xi_star at initial pose

    ndof_xi = ndof_xi - dof_xi_joint;

    q_xi = q_xi(dof_xi_joint+1:end);

    Xs = Tr.Twists(2).Xs;
    Lscale = Tr.Link.L;

    if Tr.Z_order==4
        B_Z1 = Tr.Twists(2).B_Z1_xi;
        B_Z2 = Tr.Twists(2).B_Z2_xi;
    else    % Z_order==2
        B_Z = Tr.Twists(2).B_Z_xi;
    end

    B_rho = Tr.Twists(2).B_rho;

    J_xi(1:6, :) = J_xi_here;

    % for xi_star = [0 0 0 1 0 0]' is correct
    % need to be modified later
    for ii=2:nsig
        H = Lscale*(Xs(ii) - Xs(ii-1));
        % xi
        if Tr.Z_order == 4
            xi_Z1here = xi_star(6*(ii-2)+1:6*(ii-1), 2);
            xi_Z2here = xi_star(6*(ii-2)+1:6*(ii-1), 3);
            
            B_Z1here = B_Z1(6*(ii-2)+1:6*(ii-1), :);
            B_Z2here = B_Z2(6*(ii-2)+1:6*(ii-1), :);
            if ndof_xi > 0
                xi_Z1here = B_Z1here*q_xi + xi_Z1here;
                xi_Z2here = B_Z2here*q_xi + xi_Z2here;
            end
            ad_xi_Z1here = dinamico_adj(xi_Z1here);
            BGamma_here = (H/2)*(B_Z1here + B_Z2here)+...
                          ((sqrt(3)*H^2)/12)*(ad_xi_Z1here*B_Z2here-...
                          dinamico_adj(xi_Z2here)*B_Z1here);
            Gamma_here = (H/2)*(xi_Z1here + xi_Z2here)+...
                         ((sqrt(3)*H^2)/12)*ad_xi_Z1here*xi_Z2here;
        else    % Z_order == 2
            xi_Zhere = xi_star(6*(ii-2)+1:6*(ii-1), 4);
            B_Zhere = B_Z(6*(ii-2)+1:6*(ii-1), :);
            if ndof_xi > 0
                xi_Zhere = B_Zhere*q_xi + xi_Zhere;
            end
            BGamma_here = H*B_Zhere;
            Gamma_here = H*xi_Zhere;
        end
        [gh, TGamma_here] = variable_expmap_gTg(Gamma_here);
        TBGamma_here = zeros(6, ndof_xi+dof_xi_joint);
        TBGamma_here(:, dof_xi_joint+1:end) = TGamma_here*BGamma_here;

        % updating g, Jacobian_xi
        g_here = g_here*gh;
        J_xi_here = dinamico_Adjoint(ginv(gh))*(J_xi_here + TBGamma_here);
        J_heret = J_xi_here;
        J_heret(4:6, :) = J_heret(4:6, :)*Lscale;
        J_xi(6*(ii-1)+1:6*ii, :) = J_heret;
    end
    J_rho(:, 1+dof_rho_joint:end) = B_rho;
end
