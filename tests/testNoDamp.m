clc;
clear;

% initialization
B_xi = [1 1 1 1 1 1;
        2 2 2 1 0 0]';

B_rho = [1 1];

L = createLinkage(B_xi, B_rho);
ndof_xi = L.ndof_xi;
ndof_rho = L.ndof_rho;
[t, qqd] = L.dynamics();
q_xi = qqd(:, 1:ndof_xi);
q_rho = qqd(:, ndof_xi+1: ndof_rho+ndof_xi);
figure;
plot(t, q_rho(:,1));
xlabel('time');
ylabel('q rho main');
exportgraphics(gcf, './figures/nodamp.pdf','ContentType','vector');

function L = createLinkage(B_xi, B_rho)
    S = SorosimLink();
    S.L = 0.5;
    S.basisType = 'mixed';
    S.B_xi = B_xi;
    S.B_rho = B_rho;
    L = SorosimLinkage(S, Damped=false, Gravity=true);
end
