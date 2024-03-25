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
u_rho = -[0 0]';
u_xi = -[0 0 0 15]';

q = L.statics(zeros(ndof_rho+ndof_xi,1), u_xi, u_rho);
q_xi = q(1:ndof_xi,:);
q_rho = q(ndof_xi+1:end,:);
[g, rho] = L.FwdKinematics(q_xi, q_rho);
figure;
f = L.plotq(q_xi, q_rho);

figure;
Xs = L.Twists(2).Xs;
plot(Xs, rho);
xlabel('Xs');
ylabel('rho');
grid on

function L = createLinkage(B_xi, B_rho)
    S = SorosimLink();
    S.basisType = 'mixed';
    S.L = 0.5;
    S.B_xi = B_xi;
    S.B_rho = B_rho;

    act1_y = @(X)0.018*sin(X);
    act1_z = @(X)0.018*cos(X);
%     act1_y = @(X)0.018;
%     act1_z = @(X)0;
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

    L = SorosimLinkage(S, Gravity=true,...
                          ActuationL=true,...
                          ActuationR=true,...
                          CableActuator=Cables,...
                          RadialActuator=Radials);
end
