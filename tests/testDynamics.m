clc;
clear;
rng(5);

% initialization
B_xi = [1 1 1 1 1 1;
        0 0 0 0 0 0]';

B_rho = [1 1];

L = createLinkage(B_xi, B_rho);
ndof_xi = L.ndof_xi;
ndof_rho = L.ndof_rho;
q = L.statics(zeros(ndof_rho+ndof_xi,1), 0, 0);
qqd0 = [q; zeros(ndof_rho+ndof_xi,1)];
[t, qqd] = L.dynamics(qqd0);
% ndof_xi = L.ndof_xi;
% ndof_rho = L.ndof_rho;
% q_xi = qqd(:,1:ndof_xi);
% q_rho = qqd(:,ndof_xi+1:ndof_xi+ndof_rho);

function L = createLinkage(B_xi, B_rho)
    S = SorosimLink();
    S.L = 0.5;
    S.basisType = 'mixed';
    S.B_xi = B_xi;
    S.B_rho = B_rho;
    L = SorosimLinkage(S, Damped=false, Gravity=true);
end
