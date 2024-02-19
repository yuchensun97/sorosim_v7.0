function [g, rho] = FwdKinematics(Tr, q_xi, q_rho)
%FwdKinematics: numericallly map general coordinates to Euclidean space
%   Detailed explanation goes here
%   Tr      : SorosimTwists class
%   q_xi    : general coordinates for strains
%   q_rho   : general coordinates for inflation ratio
%   returns:
%   g       : 4nip x 4 homogeneous transformation matrix
%   rho     : nip x 1 inflation ratio

    % initialization
    nsig = Tr.nsig;
    ndof_xi = Tr.ndof_xi;
    ndof_rho = Tr.ndof_rho;

    g_here = Tr.g_base; % g at base
    rho_here = Tr.rho_base; % rho at base

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
    else
        xi = B_xi_joint*q_xi_joint + xi_star_joint;
        g_joint = variable_expmap_g(xi);
    end

    if dof_rho_joint == 0
        rho_joint = 1;
    else
        rho_joint = B_rho_joint*q_rho_joint + rho_star_joint;
    end

    g_here = g_here*g_joint;
    rho_here = rho_joint;

    % soft body
    xi_star = Tr.Twists(2).xi_star; % xi_star at initial pose
    rho_star = Tr.Twists(2).rho_star; % rho_star at initial pose

    ndof_xi = ndof_xi - dof_xi_joint;
    ndof_rho = ndof_rho - dof_rho_joint;

    q_xi = q_xi(dof_xi_joint+1:end);
    q_rho = q_rho(dof_rho_joint+1:end);

    g = zeros(4*nsig, 4);
    rho = zeros(nsig, 1);
    g(1:4, :) = g_here;
    rho(1) = rho_here;
    Xs = Tr.Twists(2).Xs;
    Lscale = Tr.Link.L;

    if Tr.Z_order==4
        B_Z1 = Tr.Twists(2).B_Z1_xi;
        B_Z2 = Tr.Twists(2).B_Z2_xi;
    else    % Z_order==2
        B_Z = Tr.Twists(2).B_Z_xi;
    end

    B_rho = Tr.Twists(2).B_rho;

    for ii = 2:nsig
        H = Xs(ii) - Xs(ii-1);

        %% update g
        if Tr.Z_order==4
            xi_Z1here = xi_star(6*(ii-2)+1:6*(ii-1), 2);
            xi_Z2here = xi_star(6*(ii-2)+1:6*(ii-1), 3);
            xi_Z1here(1:3) = xi_Z1here(1:3)*Lscale; % why scaling
            xi_Z2here(1:3) = xi_Z2here(1:3)*Lscale;
            if ndof_xi > 0
                B_Z1here = B_Z1(6*(ii-2)+1:6*(ii-1), :);
                B_Z2here = B_Z2(6*(ii-2)+1:6*(ii-1), :);
                xi_Z1here = B_Z1here*q_xi + xi_Z1here;
                xi_Z2here = B_Z2here*q_xi + xi_Z2here;
            end
            ad_xi_Z1here = dinamico_adj(xi_Z1here);
            Gamma_here = (H/2)*(xi_Z1here + xi_Z2here) +...
                         ((sqrt(3)*H^2)/12)*ad_xi_Z1here*xi_Z2here;
        else % Z_order==2
            xi_Zhere = xi_star(6*(ii-2)+1:6*(ii-1), 4);
            xi_Zhere(1:3) = xi_Zhere(1:3)*Lscale; %scaling
            if ndof_xi > 0
                B_Zhere = B_Z(6*(ii-2)+1:6*(ii-1), :);
                xi_Zhere = B_Zhere*q_xi + xi_Zhere;
            end
            Gamma_here = H*xi_Zhere;
        end
        Gamma_here(4:6) = Gamma_here(4:6)*Lscale;
        gh = variable_expmap_g(Gamma_here);

        g_here = g_here*gh;
        g(4*(ii-1)+1:4*ii, :) = g_here;

        %% update rho
        if ndof_rho > 0
            rho_here = rho_star(ii) + B_rho(ii,:)*q_rho;
        else % if ndof_rho ==0, then rho is always 1
            rho_here = 1;
        end
        
        rho(ii) = rho_here;
    end

end %eof
