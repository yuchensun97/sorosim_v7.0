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

q = L.statics(zeros(ndof_rho+ndof_xi,1));
q_xi = q(1:ndof_xi,:);
q_rho = q(ndof_xi+1:end,:);
[g,rho]=L.FwdKinematics(q_xi,q_rho);
f = L.plotq(q_xi, q_rho);

function L = createLinkage(B_xi, B_rho)
    S = SorosimLink();
    S.basisType = 'robin';
    S.L = 0.5;
    S.B_xi = B_xi;
    S.B_rho = B_rho;
    L = SorosimLinkage(S);
end
