clc;
clear;
rng(5);

% initialization
B_xi = [1 1 1 1 1 1;
        2 2 2 1 0 0]';

B_rho = [0 1];

L = createLinkage(B_xi, B_rho);
q = [0.23, 0.347, 0.899, 0.176, 0.898, 0.768, 0.556, 0.237, 0.23, 0.89, 0.78, 0.177, 0.278]';
err = L.equilibrium(q, 0, 0);

function L = createLinkage(B_xi, B_rho)
    S = SorosimLink();
    S.basisType = 'mixed';
    S.L = 0.5;
    S.B_xi = B_xi;
    S.B_rho = B_rho;
    L = SorosimLinkage(S, Gravity=true, PointForce=true);
end
