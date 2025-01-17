clc;
clear;
close all;

%%
% Load the data
load('./Custom/results/reaching.mat');
OctopusLink = SorosimLink('Octopus.json');
Octopus = OctopusArm(OctopusLink);
ndof_xi = Octopus.ndof_xi;
ndof_rho = Octopus.ndof_rho;
L = Octopus.Link.L;
Lam = Octopus.Link.Lamda;
G = Octopus.Link.G;
Eta = Octopus.Link.Eta;

q = qqd(1:10:end, 1:ndof_xi+ndof_rho);
qd = qqd(1:10:end, ndof_xi+ndof_rho+1:end);
t = t(1:10:end);

Bh_xi = Octopus.Twists(2).Bh_xi;
Bh_rho = Octopus.Twists(2).Bh_rho;
Bh_rho_prime = Octopus.Twists(2).Bh_rho_prime;
B_xi_dof = Octopus.Twists(2).B_xi_dof;
B_rho_dof = Octopus.Twists(2).B_rho_dof;
B_xi_odr = Octopus.Twists(2).B_xi_odr;
B_rho_odr = Octopus.Twists(2).B_rho_odr;
xi_star = [0 0 0 1 0 0]';
rho_star = 1;
r_fn = Octopus.Link.r_fn;

Xs = 0:0.01:1;
vol = [];

%% plot \nu_2 and \rho
n = length(t);

f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;
alphaVals = linspace(0.05, 1, n); 

for i=1:n
    tt = t(i);
    q_xi = q(i, 1:ndof_xi)';
    q_rho = q(i, ndof_xi+1:end)';
    qd_rho = qd(i, ndof_xi+1:end)';

    nu2 = [];
    rho = [];
    r = [];
    r0 = [];
    nu3 = [];
    qq = [];
    Q = [];
    for xx=Xs
        xi_ = Bh_xi(xx, B_xi_dof, B_xi_odr)*q_xi + xi_star;
        rho_ = Bh_rho(xx, B_rho_dof, B_rho_odr)*q_rho + rho_star;
        rhod_ = Bh_rho(xx, B_rho_dof, B_rho_odr)*qd_rho;
        rhop_ = Bh_rho_prime(xx, B_rho_dof, B_rho_odr)*q_rho;
        rhopd_ = Bh_rho_prime(xx, B_rho_dof, B_rho_odr)*qd_rho;
        r0_ = r_fn(xx);
        I_ = 0.5 * pi * r0_^4;
        A_ = pi * r0_^2;
        r_ = rho_ * r0_;
        Q_ = G * I_ * rhop_ + Eta * I_ * rhopd_;
        qq_ = 4 * (Lam + G) * A_ * (rho_ - 1) + ...
             2 * Lam * A_ * (xi_(4) - 1) + ...
             4 * Eta * A_ * rhod_;
        r0 = [r0; r0_];
        nu2 = [nu2; xi_(2)];
        nu3 = [nu3; xi_(4)];
        rho = [rho; rho_];
        r = [r; r_];
        Q = [Q; Q_];
        qq = [qq; qq_];
    end
    
    curr_vol = sum(r.^2 .* nu3)/sum(r0.^2) - 1;
    vol = [vol; curr_vol];

    figure(f1);
    plot(Xs, -nu2, 'LineWidth', 2, 'Color', [0, 0, 1, alphaVals(i)]);
    hold on;
    
    figure(f2);
    plot(Xs, rho, 'LineWidth', 2, 'Color', [0, 0, 1, alphaVals(i)]);
    hold on;

    figure(f3);
    plot(Xs, Q, 'LineWidth', 2, 'Color', [0, 0, 1, alphaVals(i)]);
    hold on;

    figure(f4);
    plot(Xs, qq, 'LineWidth', 2, 'Color', [0, 0, 1, alphaVals(i)]);
    hold on;

end

figure(f1);
grid on;
% legend(legendInfo);
set(gca,'FontSize',28, 'FontName', 'Times New Roman');
xlabel('$X=s/L$', 'Interpreter','latex', 'FontSize', 28);
ylabel('$\kappa_2$', 'Interpreter','latex','FontSize',28);

if ~exist('./figures', 'dir')
    mkdir('./figures');
end
exportgraphics(gcf, './figures/nu_over_time.pdf','ContentType','vector');

figure(f2);
grid on;
% legend(legendInfo);
set(gca,'FontSize',28, 'FontName', 'Times New Roman');
xlabel('$X=s/L$', 'Interpreter','latex', 'FontSize', 28);
ylabel('$\rho$', 'Interpreter','latex','FontSize',28);
exportgraphics(gcf, './figures/rho_over_time.pdf','ContentType','vector');

figure(f3);
grid on;
% legend(legendInfo);
set(gca,'FontSize',28, 'FontName', 'Times New Roman');
xlabel('$X=s/L$', 'Interpreter','latex', 'FontSize', 28);
ylabel('$Q$', 'Interpreter','latex','FontSize',28);
exportgraphics(gcf, './figures/Q_over_time.pdf','ContentType','vector');

figure(f4);
grid on;
% legend(legendInfo);
set(gca,'FontSize',28, 'FontName', 'Times New Roman');
xlabel('$X=s/L$', 'Interpreter','latex', 'FontSize', 28);
ylabel('$q$', 'Interpreter','latex','FontSize',28);
exportgraphics(gcf, './figures/qq_over_time.pdf','ContentType','vector');

% figure;
% plot(t, 100 * vol);
% grid on;
% xlabel('time (s)');
% ylabel('Volume Change (%)');
% exportgraphics(gcf, './figures/vol_over_time.pdf','ContentType','vector');