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
    qdot_xi = qqd_xi(ndof_xi+1:end);
    ndof_rho = Tr.ndof_rho;
    q_rho = qqd_rho(1:ndof_rho);
    qdot_rho = qqd_rho(ndof_rho+1:end);

    [~, rho] = Tr.FwdKinematics(q_xi, q_rho);
     

end
