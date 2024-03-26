clc;
clear;

% initialization
B_xi = [1 1 1 1 1 1;
        2 2 2 1 0 0]';

B_rho = [1 1];

%% baseline vs w/ inflation
L = createLinkage(B_xi, B_rho);
ndof_xi = L.ndof_xi;
ndof_rho = L.ndof_rho;
uqt_xi = @(t)-[50 50 50 50]';
u_xi = uqt_xi(0);

uqt_rho = @(t)0;

q = L.statics(zeros(ndof_rho+ndof_xi,1), u_xi, 0);
qs_xi = q(1:ndof_xi,:);
qs_rho = q(ndof_xi+1:end,:);
qqd0 = zeros(2*(ndof_rho+ndof_xi), 1);
[t, qqd] = L.dynamics(qqd0, uqt_xi, uqt_rho, 'ode15s', 0.01, 2);

%% play video
L.plotqqd(t, qqd, 'cable');


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

    L = SorosimLinkage(S, Damped=true,...
                          ActuationL=true,...
                          CableActuator=Cables);
end
