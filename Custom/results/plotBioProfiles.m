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
bp_pos = [];

for i = 1:length(t)
    tt = t(i);
    q_xi = q(i, 1:ndof_xi)';

    curr_length = 0;
    max_curve = 0;
    J_xi_here = zeros(6, ndof_xi);
    J_xi_bp = zeros(6, ndof_xi);
    g_here = eye(4);
    g_bp = zeros(4, 4);
    for j = 1:length(Xs)
        xx = Xs(j);
        Bh_ = Bh_xi(xx+Z, B_xi_dof, B_xi_odr);
        xi_ = Bh_* q_xi + xi_star;
        Gamma_ = dXs * L * xi_;
        BGamma_ = dXs * L * Bh_;
        [gh, TGamma_] = variable_expmap_gTg(Gamma_);
        g_here = g_here * gh;
        TBGamma_ = TGamma_*BGamma_;
        J_xi_here = dinamico_Adjoint(ginv(gh))*(J_xi_here + TBGamma_);
        
        k = xi_(2);
        if abs(k) < 40 && abs(k) > max_curve && xx<0.95
            max_curve = abs(k);
            J_xi_bp = J_xi_here;
            g_bp = g_here;
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

    % get bend point position
    bp_pos_now = g_bp([1, 3], 4);
    bp_pos = [bp_pos bp_pos_now];
end

%% plot the arm length and bend point velocity
figure(1);
fl = fit(t, 100 * arm_length, 'poly5');
l_fit = fl(t);
plot(t, 100 *arm_length, '-*');
set(gca,'FontSize',28, 'FontName', 'Times New Roman')
hold on
plot(t, l_fit, '-', 'LineWidth', 4);
hold off

grid on;
legend('sim', 'fit', 'Location','southeast');
xlabel('$t$ (s)', 'Interpreter','latex', 'FontSize',28);
ylabel('$L(t)$ (cm)', 'Interpreter','latex', 'FontSize',28);
% title('Arm Length vs Time');

if ~exist('./figures', 'dir')
    mkdir('./figures');
end
exportgraphics(gcf, './figures/arm_length.pdf','ContentType','vector');

%%
figure(2);
fv = fit(t, 100 * bp_vel, 'poly5');
bp_vel_fit = fv(t);
plot(t, 100* bp_vel, '-*');
set(gca,'FontSize',28, 'FontName', 'Times New Roman');
hold on
plot(t, bp_vel_fit, '-', 'LineWidth', 4);
hold off

grid on;
legend('sim', 'fit');
xlabel('$t$ (s)', 'Interpreter','latex', 'FontSize',28);
ylabel('$u_b(t)$ (cm/s)', 'Interpreter','latex', 'FontSize',28);
% title('Bend Point Velocity vs Time');
exportgraphics(gcf, './figures/bp_velocity.pdf','ContentType','vector');

%%
figure(3);
x = 100 * bp_pos(1,:);
y = 100 * bp_pos(2, :);
f = fit(x', y', 'poly5');
x_fit = linspace(min(x), max(x), 1000)';
y_fit = f(x_fit);
bp_pos_fit = [x_fit y_fit]';
d_bp = sum(vecnorm(bp_pos_fit(:, 2:end) - bp_pos_fit(:, 1:end-1)));
plot(x, y, '-*');
set(gca,'FontSize',28, 'FontName', 'Times New Roman')
hold on
plot( x', f(x'), '-', 'LineWidth',4);
hold off

grid on;
legend('sim', 'fit');
xlabel('$X$ (cm)', 'Interpreter','latex', 'FontSize',28);
ylabel('$Y$ (cm)', 'Interpreter','latex', 'FontSize',28);
axis equal
exportgraphics(gcf, './figures/bp_pos.pdf','ContentType','vector');
