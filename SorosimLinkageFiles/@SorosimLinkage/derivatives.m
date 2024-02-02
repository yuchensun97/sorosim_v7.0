function [ydot_xi, ydot_rho] = derivatives(Tr, t, qqd_xi, qqd_rho, uqt_xi, uqt_rho)
% compute the time derivatives of q_xi and q_rho
% t is the time
% qqd_xi = [q_xi, qdot_xi]
% qqd_rho = [q_rho, qdot_rho]
% uqt_xi = TODO: fill me in future
% uqt_rho = TODO: fill me in future
% returns:
% ydot_xi = [qdot_xi, qddot_xi]
% ydot_rho = [qdot_rho, qddot_rho]

    persistent tlast
    if t==0
        tlast=cputime;
    end
    if cputime-tlast>0.5
        tlast=cputime;
        disp(['t = ', num2str(t)]);
    end

    density = Tr.Link.Rho0;

    ndof_xi = Tr.ndof_xi;
    q_xi = qqd_xi(1:ndof_xi);
    qd_xi = qqd_xi(ndof_xi+1:end);
    ndof_rho = Tr.ndof_rho;
    q_rho = qqd_rho(1:ndof_rho);
    qd_rho = qqd_rho(ndof_rho+1:end);

    nsig = Tr.nsig;
    M_xi = zeros(ndof_xi, ndof_xi);
    C_xi = zeros(ndof_xi, ndof_xi);
    F_xi = zeros(ndof_xi, 1);

    K_rho = zeros(ndof_rho, ndof_rho);
    F_rho = zeros(ndof_rho, 1);

    g = zeros(4*nsig, 4);
    J_xi = zeros(6*nsig, ndof_xi);
    eta = zeros(6*nsig, 1);
    Jd_xi = zeros(6*nsig, ndof_xi);

    rho = zeros(nsig, 1);
    J_rho = Tr.Twists(2).B_rho;
    rhod = zeros(nsig, ndof_rho);
    Jp_rho = Tr.Jacobianprime();
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
    qd_xi_joint = qd_xi(1:dof_xi_joint);
    qd_rho_joint = qd_rho(1:dof_rho_joint);
    xi_star_joint = Tr.Twists(1).xi_star;
    rho_star_joint = Tr.Twists(1).rho_star;
    eta_joint = zeros(6,1);
    Jd_joint_xi = zeros(6,1);

    if dof_xi_joint == 0
        g_joint = eye(4);
        TgB_joint = zeros(6, ndof_xi);
        TgBd_joint = zeros(6, ndof_xi);
    else
        xi = B_xi_joint*q_xi_joint + xi_star_joint;
        xid = B_xi_joint * qd_xi_joint;
        [g_joint, Tg, Tgd] = variable_expmap_gTg(xi, xid);
        TgB_joint = zeros(6, ndof_xi);
        TgB_joint(:, 1:dof_xi_joint) = Tg*B_xi_joint;
        TgBd_joint = zeros(6, ndof_xi);
        TgBd_joint(:, 1:dof_xi_joint) = dinamico_adj(eta_joint)*Tg*B_xi_joint+...
                                        + Tgd*B_xi_joint;
    end

    if dof_rho_joint == 0
        rho_joint = 1;
        rhod_joint = 0;
    else
        rho_joint = B_rho_joint*q_rho_joint+rho_star_joint;
        rhod_joint = B_rho_joint*qd_rho_joint;
    end

    g_here = g_here*g_joint;
    rho(1) = rho_joint;
    rhod(1) = rhod_joint;

    G = Tr.G; % gravity
    Ad_g_joint_inv = dinamico_Adjoint(ginv(g_joint));
    J_here_xi = Ad_g_joint_inv*...
                (TgB_joint + J_here_xi);
    Jd_here_xi = Ad_g_joint_inv*(Jd_joint_xi + TgBd_joint);
    eta_here = Ad_g_joint_inv*(eta_joint + TgB_joint(:, 1:dof_xi_joint)*qd_xi_joint);

    % soft body
    ndof_xi = ndof_xi - dof_xi_joint;
    q_xi = q_xi(dof_xi_joint+1:end);
    qd_xi = qd_xi(dof_xi_joint+1:end);
    q_rho = q_rho(dof_xi_joint+1:end);
    qd_rho = qd_rho(dof_xi_joint+1:end);

    gi = Tr.Link.gi;
    g_here = g_here*gi;
    Ad_gi_inv = dinamico_Adjoint(ginv(gi));
    J_here_xi = Ad_gi_inv *J_here_xi;
    Jd_here_xi = Ad_gi_inv * Jd_here_xi;
    eta_here = Ad_gi_inv*eta_here;
    
    J_here_rho = J_rho(1, :);
    
    g(1:4, :) = g_here;
    J_xi(1:6, :) = J_here_xi;
    Jd_xi(1:6, :) = Jd_here_xi;
    eta(1:6) = eta_here;

    rho(1) = J_here_rho * q_rho + rho_star(1);
    rhod(1) = J_here_rho * qd_rho;

    q_here_xi = q_xi;
    qd_here_xi = qd_xi;

    ld = Tr.Link.L;
    Ms = Tr.Twists(2).Ms;
    Xs = Tr.Twists(2).Xs;
    Ws = Tr.Twists(2).Ws;
    nip = Tr.Twists(2).nip;

    for ii = 2:nip
        H = (Xs(ii) - Xs(ii-1)) *ld;

        xi_Z1here = [0 0 0 1 0 0]';
        xi_Z2here = [0 0 0 1 0 0]';
        
        B_Z1here = Tr.Twists(2).B_Z1(6*(ii-2)+1:6*(ii-1),:);
        B_Z2here = Tr.Twists(2).B_Z2(6*(ii-2)+1:6*(ii-1),:);

        xi_Z1here = B_Z1here * q_here_xi + xi_Z1here;
        xi_Z2here = B_Z2here * q_here_xi + xi_Z2here;

        xid_Z1here = B_Z1here * qd_here_xi;
        ad_xi_Z1here = dinamico_adj(xi_Z1here);

        BGamma_here = (H/2)*(B_Z1here+B_Z2here)+...
                      ((sqrt(3)*H^2)/12)*(ad_xi_Z1here*B_Z2here-dinamico_adj(xi_Z2here)*B_Z1here);

        Gammadd_Z4_dq_here = ((sqrt(3)*H^2)/6)*dinamico_adj(xid_Z1here)*B_Z2here;
        Gammad_here = BGamma_here * qd_here;

        Gamma_here = (H/2)*(xi_Z1here+xi_Z2here)+...
                        ((sqrt(3)*H^2)/12)*ad_xi_Z1here*xi_Z2here;

        [gh, TGamma_here, TGammad_here] = variable_expmap_gTg_mex(Gamma_here, Gammad_here);
        TBGamma_here = zeros(6, ndof_xi+dof_xi_joint);
        TBGamma_here(:, dof_xi_joint+1:end) = TGamma_here*BGamma_here;
        TBGammad_here = dinamico_adj(eta_here)*TBGamma_here+TGammad_here*BGamma_here+...
                        +TGamma_here * Gammadd_Z4_dq_here;

        %updating g, Jacobian, Jacobian_dot, eta
        g_here = g_here*gh;
        Ad_gh_inv = dinamico_Adjoint(ginv(gh));
        J_here_xi = Ad_gh_inv*(J_here_xi + TBGamma_here);
        Jd_here_xi = Ad_gh_inv*(Jd_here_xi+TBGammad_here);
        eta_here = Ad_gh_inv*(eta_here+TBGamma_here*qd_here_xi);

        g((ii-1)*4+1:ii*4, :) = g_here;
        J_xi((ii-1)*6+1:ii*6, :) = J_here_xi;
        Jd_xi((ii-1)*6+1:ii*6, :) = Jd_here_xi;
        eta((ii-1)*6+1:ii*6, :) = eta_here;

        %updating rho, Jacobian, Jacobian_prime, rho_dot
        if ndof_rho > 0
            J_rho_here = J_rho(dof_joint_rho+ii,:);
            rho_here = rho_star(ii) + J_rho_here*q_rho;
            rhod_here = J_rho(ii,:)*qd_rho;
        else
            rho_here = 1;
            rhod_here = 0;
            J_rho_here = zeros(1, dof_joint_rho+ndof_rho);
        end

        rho(ii) = rho_here;
        rhod(ii) = rhod_here;
        rho_star_here = rho_star(ii);

        %integrals evaluation
        if Ws(ii)>0
            W_here = Ws(ii);
            Ms_here = Ms(6*(ii-1)+1:6*ii,:);
            I11 = Ms_here(2, 2);
            I22 = Ms_here(3, 3);
            omega_1 = eta_here(2);
            omega_2 = eta_here(3);
            omega_3 = eta_here(1);
            MI_here = density*(I11*(omega_2^2+omega_3^2)+I22*(omega_1^2+omega_3^2));
            Ms_here(1:3,1:3) = rho_here^4 * Ms_here(1:3, 1:3);
            Ms_here(4:6,4:6) = rho_here^2 * Ms_here(4:6, 4:6);
            Ms_dot_here = Ms_here;
            Ms_dot_here(1:3, 1:3) = 4*rho_here^3 * Ms_dot_here(1:3, 1:3);
            Ms_dot_here(4:6, 4:6) = 2*rho_here * Ms_dot_here(4:6, 4:6);
            Ms_dot_here = rhod_here * Ms_dot_here;

            if Tr.Gravity
                F_xi = F_xi + ld*W_here*J_here_xi'*Ms_here*dinamico_Adjoint(ginv(g_here))*G;
            end
            Qtemp = ld*W_here*J_here_xi'*Ms_here;
            M_xi = M_xi + Qtemp*J_here_xi;
            C_xi = C_xi + Qtemp*Jd_here +...
                   ld*W_here*J_here_xi'*Ms_dot_here*J_here_xi+...
                   ld*W_here*J_here_xi'*dinamico_coadj(eta_here)*Ms_here*J_here_xi;
            K_rho_temp = K_rho - ld*W_here*J_rho_here'*MI_here*J_rho_here;% last term of K_rho
            F_rho = F_rho + ld*W_here*J_rho_here'*MI_here*rho_star_here;
        end
    end

    Bq_xi = 0;
    u_xi = 0;
    Bq_rho = 0;
    u_rho = 0;

    if Tr.Damped
        D_xi = Tr.D;
        D_xi_bar = Tr.D_xi_bar;
    else
        D_xi = 0;
        D_xi_bar = 0;
    end

    K_xi = Tr.K_xi;
    K_xi_bar = Tr.K_xi_bar;
    qdd_xi = M_xi\(Bq_xi*u_xi+F_xi-K_xi*q_xi-(C_xi+D_xi)*qd_xi-...
            K_xi_bar*q_rho - D_xi_bar*qd_rho);
    ydot_xi = [qd_xi;qdd_xi];

    if ndof_rho+dof_rho_joint == 0
        ydot_rho = zeros(ndof_rho+dof_rho_joint,1);
        return
    end

end
