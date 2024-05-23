clc;
clear;
close all;

%%
% load the data
load('./Custom/results/reaching.mat');
OctopusLink = SorosimLink('Octopus.json');
Octopus = OctopusArm(OctopusLink);
ndof_xi = Octopus.ndof_xi;
ndof_rho = Octopus.ndof_rho;
L = Octopus.Link.L;

q = qqd(:, 1:ndof_xi+ndof_rho);
qd = qqd(:, ndof_xi+ndof_rho+1:end);

Bh_xi = Octopus.Twists(2).Bh_xi;
Bh_rho = Octopus.Twists(2).Bh_rho;
B_xi_dof = Octopus.Twists(2).B_xi_dof;
B_rho_dof = Octopus.Twists(2).B_rho_dof;
B_xi_odr = Octopus.Twists(2).B_xi_odr;
B_rho_odr = Octopus.Twists(2).B_rho_odr;
xi_star = [0 0 0 1 0 0]';
rho_star = 1;

%% compute the arm length and bend point velocity
dXs = 0.01;
Z = 0.5 * dXs;
Xs = 0:dXs:1;

arm_length = [];
bp_vel = [];


for i = 1:length(t)
    tt = t(i);
    q_xi = q(i, 1:ndof_xi)';

    curr_length = 0;
    max_curve = 0;
    J_xi_here = zeros(6, ndof_xi);
    J_xi_bp = zeros(6, ndof_xi);
    for j = 1:length(Xs)
        xx = Xs(j);
        Bh_ = Bh_xi(xx+Z, B_xi_dof, B_xi_odr);
        xi_ = Bh_* q_xi + xi_star;
        Gamma_ = dXs * L * xi_;
        BGamma_ = dXs * L * Bh_;
        [gh, TGamma_] = variable_expmap_gTg(Gamma_);
        TBGamma_ = TGamma_*BGamma_;
        J_xi_here = dinamico_Adjoint(ginv(gh))*(J_xi_here + TBGamma_);
        
        k = xi_(2);
        if abs(k) > max_curve
            max_curve = abs(k);
            J_xi_bp = J_xi_here;
        end

        ds = sqrt(xi_(4).^2 + xi_(6).^2) * dXs * L;
        curr_length = curr_length + ds;
    end

    arm_length = [arm_length; curr_length];

    % compute bend point velocity
    qd_xi = qd(i, 1:ndof_xi)';
    eta_curr = J_xi_bp * qd_xi;
    
    vel_curr = norm(eta_curr(4:end));
    bp_vel = [bp_vel; vel_curr];
end

%% plot the arm length and bend point velocity
figure(1);
plot(t, arm_length);
grid on;
xlabel('Time (s)');
ylabel('Arm Length (m)');
title('Arm Length vs Time');

if ~exist('./figures', 'dir')
    mkdir('./figures');
end
exportgraphics(gcf, './figures/arm_length.pdf','ContentType','vector');

%%
figure(2);
plot(t, bp_vel);

grid on;
xlabel('Time (s)');
ylabel('Bend Point Velocity (m/s)');
title('Bend Point Velocity vs Time');
exportgraphics(gcf, './figures/bp_velocity.pdf','ContentType','vector');
