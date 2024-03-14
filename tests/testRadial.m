clc;
clear;
rng(5);

% initialization
B_xi = [1 1 1 1 1 1;
        2 2 2 1 0 0]';

B_rho_base = [0 1];
B_rho = [1 1];

%% baseline vs w/ inflation
L = createLinkage(B_xi, B_rho);
ndof_xi = L.ndof_xi;
ndof_rho = L.ndof_rho;
u_rho = -[0.12 0.35]'*1e5;

q = L.statics(zeros(ndof_rho+ndof_xi,1), 0, u_rho);
q_xi = q(1:ndof_xi,:);
q_rho = q(ndof_xi+1:end,:);
[g, rho] = L.FwdKinematics(q_xi, q_rho);
f = L.plotq(q_xi, q_rho);

function L = createLinkage(B_xi, B_rho)
    S = SorosimLink();
    S.basisType = 'mixed';
    S.L = 0.5;
    S.B_xi = B_xi;
    S.B_rho = B_rho;

    act1 = Radial(0.2, 0.5);
    act2 = Radial(0.6, 0.8);

    Radials = RadialActuation(act1, act2);

    L = SorosimLinkage(S, Gravity=false,...
                          ActuationR=true,...
                          RadialActuator=Radials);
end
