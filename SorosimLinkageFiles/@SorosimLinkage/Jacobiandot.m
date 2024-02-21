% Function that calculates the time derivative of the Jacobian w.r.t. strain
function Jd = Jacobiandot(Tr, q_xi, qd_xi)
%Jacobiandot: time derivative of the Jacobian of the soft link
%   Detailed explanation goes here
%   Tr      : SorosimLinkage class
%   q_xi    : general coordinates for strains
%   qd_xi   : time derivative of general coordinates for strains
%   returns:
%   Jd   : 6nip x ndof_xi  time derivative Jacobian for strain variables
    if isrow(q_xi)
        q_xi = q_xi';
    end

    if isrow(qd_xi)
        qd_xi = qd_xi';
    end

    ndof = Tr.ndof_xi;
    g_here = Tr.g_base;
    Jd_here = zeros(6, ndof);
    eta_here = zeros(6, 1);
    nsig = Tr.nsig;

    Jd = zeros(6*nsig, ndof);


    % joint
    dof_joint = Tr.Twists(1).dof_xi;
    q_here = q_xi(1:dof_joint);
    qd_here = qd_xi(1:dof_joint);
    B_here = Tr.Twists(1).B_xi;
    xi_star = Tr.Twists(1).xi_star;

    if dof_joint == 0 %fixed joint
        g_joint = eye(4);
        TgB_here = zeros(6, ndof);
        TgBd_here = zeros(6, ndof);
    else
        xi = B_here*q_here + xi_star;
        xid = B_here*qd_here;
        [g_joint, Tg, Tgd] = variable_expmap_gTg_mex(xi, xid);

        TgB_here = zeros(6, ndof);
        TgB_here(:, 1:dof_joint) = Tg*B_here;
        TgBd_here = zeros(6, ndof);
        TgBd_here(:, 1:dof_joint) = dinamico_adj(eta_here)*Tg*B_here+Tgd*B_here;
    end

    g_here = g_here*g_joint;
    Ad_g_joint_inv = dinamico_Adjoint(ginv(g_joint));
    Jd_here = Ad_g_joint_inv*(Jd_here + TgBd_here);
    eta_here = Ad_g_joint_inv*(eta_here + TgB_here(:, 1:dof_joint)*qd_here);

    % soft body
    dof_here = ndof - dof_joint;
    q_here = q_xi(dof_joint+1:end);
    qd_here = qd_xi(dof_joint+1:end);
    xi_star = Tr.Twists(2).xi_star;
    L = Tr.Link.L;
    gi = Tr.Link.gi;
    Xs = Tr.Twists(2).Xs;
    nip = Tr.Twists(2).nip;

    g_here = g_here*gi;
    Ad_gi_inv = dinamico_Adjoint(ginv(gi));
    Jd_here = Ad_gi_inv*Jd_here;
    eta_here = Ad_gi_inv*eta_here;

    Jd(1:6, :) = Jd_here;

    %scaling
    Lscale = L;

    for ii = 2:nip
        H = L*(Xs(ii) - Xs(ii-1));

        if Tr.Z_order==4
            xi_Z1here = xi_star(6*(ii-2)+1:6*(ii-1), 2);
            xi_Z2here = xi_star(6*(ii-2)+1:6*(ii-1), 3);

            B_Z1here = Tr.Twists(2).B_Z1_xi(6*(ii-2)+1:6*(ii-1), :);
            B_Z2here = Tr.Twists(2).B_Z2_xi(6*(ii-2)+1:6*(ii-1), :);

            if dof_here > 0
                xi_Z1here = B_Z1here*q_here + xi_Z1here;
                xi_Z2here = B_Z2here*q_here + xi_Z2here;

                xid_Z1here = B_Z1here*qd_here;

                ad_xi_Z1here = dinamico_adj(xi_Z1here);

                BGamma_here = (H/2)*(B_Z1here+B_Z2here)+...
                              ((sqrt(3)*H^2)/12)*(ad_xi_Z1here*B_Z2here-dinamico_adj(xi_Z2here)*B_Z1here);
                Gammadd_Z4_dq_here = ((sqrt(3)*H^2)/6)*dinamico_adj(xid_Z1here)*B_Z2here;
                Gammad_here = BGamma_here*qd_here;
            else
                ad_xi_Z1here = dinamico_adj(xi_Z1here);
                BGamma_here = (H/2)*(B_Z1here+B_Z2here)+...
                              ((sqrt(3)*H^2)/12)*(ad_xi_Z1here*B_Z2here-dinamico_adj(xi_Z2here)*B_Z1here);
                Gammadd_Z4_dq_here = zeros(6, ndof);
                Gammad_here = zeros(6, 1);
            end
            Gamma_here = (H/2)*(xi_Z1here+xi_Z2here)+...
                         ((sqrt(3)*H^2)/12)*ad_xi_Z1here*xi_Z2here;
        else % order 2
            xi_Zhere = xi_star(6*(ii-2)+1:6*(ii-1), 4);

            B_Zhere = Tr.Twists(2).B_Z_xi(6*(ii-2)+1:6*(ii-1), :);

            if dof_here > 0
                xi_Zhere = B_Zhere*q_here + xi_Zhere;
                BGamma_here = H*B_Zhere;
                Gammad_here = BGamma_here*qd_here;
            else
                BGamma_here = H*B_Zhere;
                Gammad_here = zeros(6, 1);
            end
            Gamma_here = H * xi_Zhere;
        end

        [gh, TGamma_here, TGammad_here] = variable_expmap_gTgTgd_mex(Gamma_here, Gammad_here);
        TBGamma_here = zeros(6, ndof+dof_joint);
        TBGamma_here(:, dof_joint+1:end) = TGamma_here*BGamma_here;
        TBGammad_here = zeros(6, ndof+dof_joint);
        TBGammad_here(:, dof_joint+1:end) = dinamico_adj(eta_here)*TBGamma_here(:, dof_joint+1:end)+TGammad_here*BGamma_here;

        if Tr.Z_order==4
            TBGammad_here(:, dof_joint+1:end) = TBGammad_here(:, dof_joint+1:end) + TGamma_here*Gammadd_Z4_dq_here;
        end

        %updating g, Jacobian, Jacobian_dot and eta
        g_here = g_here*gh;
        Ad_gh_inv = dinamico_Adjoint(ginv(gh));
        Jd_here = Ad_gh_inv*(Jd_here + TBGammad_here);
        eta_here = Ad_gh_inv*(eta_here + TBGamma_here(:, dof_joint+1:end)*qd_here);

        Jd_heret = Jd_here;
        Jd(6*(ii-1)+1:6*ii, :) = Jd_heret;
    end

end
