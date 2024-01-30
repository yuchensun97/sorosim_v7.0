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

    M_rho = zeros(ndof_rho, ndof_rho);
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
    B_xi_joint = Tr.Twists(1).B_xi;
    q_xi_joint = q_xi(1:dof_xi_joint);
    xi_star_joint = Tr.Twists(1).xi_star;

    if dof_xi_joint == 0
        g_joint = eye(4);
        TgB_joint = zeros(6, ndof_xi);
    else
        xi = B_xi_joint*q_xi_joint + xi_star_joint;
        [g_joint, Tg] = variable_expmap_gTg(xi);
        TgB_joint = zeros(6, ndof_xi);
        TgB_joint(:, 1:dof_xi_joint) = Tg*B_xi_joint;
    end

    g_here = g_here*g_joint;

    G = Tr.G; % gravity

    J_here_xi = dinamico_Adjoint(ginv(g_joint))*...
                (TgB_joint + J_here_xi);
    Jd_here_xi = zeros(6, ndof_xi);
    eta_here = zeros(6, 1);

    % soft body
    ndof_xi = ndof_xi - dof_xi_joint;
    q_xi = q_xi(dof_xi_joint+1:end);
    
    % TODO: update J(0) of the body
    J_here_rho = J_rho(1, :);
    
    g(1:4, :) = g_here;
    J_xi(1:6, :) = J_here_xi;
    Jd_xi(1:6, :) = Jd_here_xi;
    eta(1:6) = eta_here;

    rho(1) = J_here_rho * q_rho + rho_star;
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
            rho_here = rho_star(ii) + J_rho(ii,:)*q_rho;
        else
            rho_here = 1;
        end

        rho(ii) = rho_here;
        J_here_rho = J_rho(ii, :);

        if ndof_rho>0
            rhod_here = J_here_rho * qd_rho;
        else
            rhod_here = 0;
        end

        rhod(ii) = rhod_here;

        %integrals evaluation
        if Ws(ii)>0
            W_here = Ws(ii);
            Ms_here = Ms(6*(ii-1)+1:6*ii,:);
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

            % TODO: add M terms related to rho
        end
    end

    % TODO: currently ignore point force
    Bq_xi = 0;
    u_xi = 0;
    Bq_rho = 0;
    u_rho = 0;

end
