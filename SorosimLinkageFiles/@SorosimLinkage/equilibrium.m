function err = equilibrium(Tr, qu, uqt_xi, uqt_rho) %unscaled
    % compute the time derivatives of q_xi and q_rho
    % t is the time
    % qu = [q_xi, q_rho]
    % uqt_xi = TODO: fill me in future
    % uqt_rho = TODO: fill me in future
    % returns:
    % ydot = [qd_xi, qd_rho, qdd_xi, qdd_rho]

    ndof_xi = Tr.ndof_xi;
    ndof_rho = Tr.ndof_rho;
    q_xi = qu(1:ndof_xi);
    q_rho = qu(ndof_xi+1:ndof_rho+ndof_xi);

    nsig = Tr.nsig;

    F_xi = zeros(ndof_xi, 1);

    g = zeros(4*nsig, 4);
    J_xi = zeros(6*nsig, ndof_xi);

    rho = zeros(nsig, 1);
    J_rho = Tr.Twists(2).B_rho;
    rho_star = Tr.Twists(2).rho_star;

    g_here = Tr.g_base;
    J_here_xi = zeros(6, ndof_xi);
    
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
        xid = B_xi_joint * qd_xi_joint;
        [g_joint, Tg] = variable_expmap_gTg(xi, xid);
        TgB_joint = zeros(6, ndof_xi);
        TgB_joint(:, 1:dof_xi_joint) = Tg*B_xi_joint;
    end

    if dof_rho_joint == 0
        rho_joint = 1;
    else
        rho_joint = B_rho_joint*q_rho_joint+rho_star_joint;
    end

    g_here = g_here*g_joint;
    rho(1) = rho_joint;

    G = Tr.G; % gravity
    Ad_g_joint_inv = dinamico_Adjoint(ginv(g_joint));
    J_here_xi = Ad_g_joint_inv*...
                (TgB_joint + J_here_xi);

    % soft body
    ndof_xi = ndof_xi - dof_xi_joint;
    ndof_rho = ndof_rho - dof_rho_joint;
    q_xi = q_xi(dof_xi_joint+1:end);
    q_rho = q_rho(dof_rho_joint+1:end);

    gi = Tr.Link.gi;
    g_here = g_here*gi;
    Ad_gi_inv = dinamico_Adjoint(ginv(gi));
    J_here_xi = Ad_gi_inv *J_here_xi;

    J_here_rho = J_rho(1, :);
    
    %at X=0
    g(1:4, :) = g_here;
    J_xi(1:6, :) = J_here_xi;

    rho(1) = J_here_rho * q_rho + rho_star(1);

    q_here_xi = q_xi;

    ld = Tr.Link.L;
    Ms = Tr.Twists(2).Ms;
    Xs = Tr.Twists(2).Xs;
    Ws = Tr.Twists(2).Ws;
    nip = Tr.Twists(2).nip;

    for ii = 2:nip
        H = (Xs(ii) - Xs(ii-1)) *ld;

        xi_Z1here = [0 0 0 1 0 0]'; %modify later to include xi_star
        xi_Z2here = [0 0 0 1 0 0]';
        
        B_Z1here = Tr.Twists(2).B_Z1_xi(6*(ii-2)+1:6*(ii-1),:);
        B_Z2here = Tr.Twists(2).B_Z2_xi(6*(ii-2)+1:6*(ii-1),:);

        xi_Z1here = B_Z1here * q_here_xi + xi_Z1here;
        xi_Z2here = B_Z2here * q_here_xi + xi_Z2here;

        ad_xi_Z1here = dinamico_adj(xi_Z1here);

        BGamma_here = (H/2)*(B_Z1here+B_Z2here)+...
                        ((sqrt(3)*H^2)/12)*(ad_xi_Z1here*B_Z2here-dinamico_adj(xi_Z2here)*B_Z1here);

        Gamma_here = (H/2)*(xi_Z1here+xi_Z2here)+...
                        ((sqrt(3)*H^2)/12)*ad_xi_Z1here*xi_Z2here;

        [gh, TGamma_here] = variable_expmap_gTg(Gamma_here);
        TBGamma_here = zeros(6, ndof_xi+dof_xi_joint);
        TBGamma_here(:, dof_xi_joint+1:end) = TGamma_here*BGamma_here;

        %updating g, Jacobian, Jacobian_dot, eta
        g_here = g_here*gh;
        Ad_gh_inv = dinamico_Adjoint(ginv(gh));
        J_here_xi = Ad_gh_inv*(J_here_xi + TBGamma_here);

        g((ii-1)*4+1:ii*4, :) = g_here;
        J_xi((ii-1)*6+1:ii*6, :) = J_here_xi;

        %updating rho, Jacobian, Jacobian_prime, rho_dot
        if ndof_rho > 0
            J_rho_here = J_rho(dof_rho_joint+ii,:);
            rho_here = rho_star(ii) + J_rho_here*q_rho; %J is Phi
        else
            rho_here = 1;
        end

        rho(ii) = rho_here;

        %integrals evaluation
        if Ws(ii)>0
            W_here = Ws(ii);
            Ms_here = Ms(6*(ii-1)+1:6*ii,:);

            if Tr.Gravity
                F_xi = F_xi + ld*W_here*J_here_xi'*Ms_here*dinamico_Adjoint(ginv(g_here))*G;
            end
        end
    end

    Bq_xi = 0;
    Bq_rho = 0;

    if Tr.PointForce
        % hard code
        % map point force to local frame
        Fp_vec = [0 0 0 -500 0 0]';
        g_here = g(end-3:end,:);
        g_here(1:3,4) = zeros(3,1);
        Ad_g_here_inv = dinamico_Adjoint(ginv(g_here));
        Fp_vec = Ad_g_here_inv*Fp_vec;
        F_xi = F_xi + J_xi(end-5:end,:)'*Fp_vec;
    end

    K_xi = Tr.K_xi;
    K_xi_bar = Tr.K_xi_bar;
    E_xi = K_xi*q_xi+K_xi_bar*q_rho-Bq_xi*uqt_xi-F_xi;

    K_rho = Tr.K_rho_part;
    K_rho_bar = Tr.K_rho_bar;
    E_rho = K_rho*q_rho+K_rho_bar*q_xi-Bq_rho*uqt_rho;
    err = [E_xi;E_rho]*1e7;
end
    