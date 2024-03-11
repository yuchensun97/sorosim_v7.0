clc;
clear;
rng(5);

% initialization
B_xi = [1 1 1 1 1 1;
        2 2 2 1 0 0]';

B_rho = [1 1];

L = createLinkage(B_xi, B_rho);
ndof_xi = L.ndof_xi;
ndof_rho = L.ndof_rho;
u_xi = -[125 125 125 125]';

q = L.statics(zeros(ndof_rho+ndof_xi,1), u_xi, 0);
q_xi = q(1:ndof_xi,:);
q_rho = q(ndof_xi+1:end,:);
[g,rho]=L.FwdKinematics(q_xi,q_rho);
f = L.plotq(q_xi, q_rho);

function L = createLinkage(B_xi, B_rho)
    S = SorosimLink();
    S.basisType = 'dirichlet';
    S.L = 0.5;
    S.B_xi = B_xi;
    S.B_rho = B_rho;

    act1_y = @(X)0.018;
    act1_z = @(X)0;
    act1 = Cable(act1_y, act1_z);

    act2_y = @(X)0;
    act2_z = @(X)0.018;
    act2 = Cable(act2_y, act2_z);

    act3_y = @(X)-0.018;
    act3_z = @(X)0;
    act3 = Cable(act3_y, act3_z);

    act4_y = @(X)0;
    act4_z = @(X)-0.018;
    act4 = Cable(act4_y, act4_z);

    Cables = CableActuation(act1, act2, act3, act4);

    L = SorosimLinkage(S, ActuationL=true,...
                       CableActuator=Cables);
end
