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
[t, qqd] = L.dynamics();

function L = createLinkage(B_xi, B_rho)
    S = SorosimLink();
    S.L = 0.5;
    S.basisType = 'mixed';
    S.B_xi = B_xi;
    S.B_rho = B_rho;
    L = SorosimLinkage(S, Damped=false, Gravity=true);
end
