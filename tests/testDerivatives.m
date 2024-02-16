clc;
clear;
rng(5);

% initialization
B_xi = [1 1 1 1 1 1;
        0 1 1 2 1 0]';

B_rho = [1 1];

L = createLinkage(B_xi, B_rho);
ndof_xi = L.ndof_xi;
ndof_rho = L.ndof_rho;
qqd = zeros(2*(ndof_rho+ndof_xi),1);
dqdt = L.derivatives(0,qqd,0,0);

function L = createLinkage(B_xi, B_rho)
    S = SorosimLink();
    S.B_xi = B_xi;
    S.B_rho = B_rho;
    L = SorosimLinkage(S);
end
