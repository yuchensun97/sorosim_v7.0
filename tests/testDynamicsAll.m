clc;
clear;

% initialization
B_xi = [1 1 1 1 1 1;
        2 2 2 1 0 0]';

B_rho = [1 2];

%% baseline vs w/ inflation
L = createLinkage(B_xi, B_rho);
ndof_xi = L.ndof_xi;
ndof_rho = L.ndof_rho;
uqt_rho = @(t)-[1.25 0.5]'*1e5; % Pa
uqt_xi = @(t)-[15 6 8 5]'; % N

u_xi = uqt_xi(0);
u_rho = uqt_rho(0);

q = L.statics(zeros(ndof_rho+ndof_xi,1), u_xi, u_rho);
qs_xi = q(1:ndof_xi,:);
qs_rho = q(ndof_xi+1:end,:);
qqd0 = zeros(2*(ndof_rho+ndof_xi), 1);
[t, qqd] = L.dynamics(qqd0, uqt_xi, uqt_rho, 'ode15s', 0.01, 15);

%% play video
L.plotqqd(t, qqd, 'actuators');

function L = createLinkage(B_xi, B_rho)
    S = SorosimLink();
    S.basisType = 'mixed';
    S.L = 0.5;
    S.B_xi = B_xi;
    S.B_rho = B_rho;

    act1_y = @(X)0.018*sin(X);
    act1_z = @(X)0.018*cos(X);
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

    act1_r = Radial(0.2, 0.4);
    act2_r = Radial(0.7, 0.9);

    Radials = RadialActuation(act1_r, act2_r);

    L = SorosimLinkage(S, Damped= true, Gravity=true,...
                          ActuationL=true,...
                          ActuationR=true,...
                          CableActuator=Cables,...
                          RadialActuator=Radials);
end
