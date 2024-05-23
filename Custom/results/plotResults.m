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

q = qqd(1:10:end, 1:ndof_xi+ndof_rho);
t = t(1:10:end);

Bh_xi = Octopus.Twists(2).Bh_xi;
Bh_rho = Octopus.Twists(2).Bh_rho;
B_xi_dof = Octopus.Twists(2).B_xi_dof;
B_rho_dof = Octopus.Twists(2).B_rho_dof;
B_xi_odr = Octopus.Twists(2).B_xi_odr;
B_rho_odr = Octopus.Twists(2).B_rho_odr;
xi_star = [0 0 0 1 0 0]';
rho_star = 1;

Xs = 0:0.01:1;

%% plot \nu_2 and \rho
n = length(t);

f1 = figure;
f2 = figure;
for i=1:n
    tt = t(i);
    q_xi = q(i, 1:ndof_xi)';
    q_rho = q(i, ndof_xi+1:end)';

    nu2 = [];
    rho = [];
    for xx=Xs
        xi_ = Bh_xi(xx, B_xi_dof, B_xi_odr)*q_xi + xi_star;
        rho_ = Bh_rho(xx, B_rho_dof, B_rho_odr)*q_rho + rho_star;
        nu2 = [nu2; xi_(2)];
        rho = [rho; rho_];
    end
    legendInfo{i} = ['t = ' num2str(tt)];
    figure(f1);
    plot(Xs, nu2);
    hold on;
    
    figure(f2);
    plot(Xs, rho);
    hold on;
    
end

figure(f1);
grid on;
legend(legendInfo);
xlabel('Xs');
ylabel('\nu_2');

if ~exist('./figures', 'dir')
    mkdir('./figures');
end
exportgraphics(gcf, './figures/nu_over_time.pdf','ContentType','vector');

figure(f2);
grid on;
legend(legendInfo);
xlabel('Xs');
ylabel('\rho');
exportgraphics(gcf, './figures/rho_over_time.pdf','ContentType','vector');
