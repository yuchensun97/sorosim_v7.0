clc;
clear;
rng(5);

% initialization
B_xi = [1 1 1 1 1 1;
        0 0 0 0 0 0]';

B_rho = [0 1];

L = createLinkage(B_xi, B_rho);
ndof_xi = L.ndof_xi;
ndof_rho = L.ndof_rho;
qu = L.statics(zeros(ndof_rho+ndof_xi,1));
q_xi = qu(1:ndof_xi);
q_rho = qu(ndof_xi+1:end);

function L = createLinkage(B_xi, B_rho)
    S = SorosimLink();
    S.B_xi = B_xi;
    S.B_rho = B_rho;
    L = SorosimLinkage(S);
end
