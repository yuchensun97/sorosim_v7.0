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
% qu = [0.584 0.985 0.4569 0.46289 0.2168 0.8543]';
% err = L.equilibrium(qu, 0, 0);

q = L.statics(zeros(ndof_rho+ndof_xi,1));

function L = createLinkage(B_xi, B_rho)
    S = SorosimLink();
    S.L = 1;
    S.B_xi = B_xi;
    S.B_rho = B_rho;
    L = SorosimLinkage(S);
end
