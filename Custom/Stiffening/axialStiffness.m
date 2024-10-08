clc;
clear;
close all;
%% create Link
OctopusLink = SorosimLink('Stiffness.json');
LOM = createLOM(OctopusLink);
TM = createTM();
Octopus = SorosimLinkage(OctopusLink, Damped=true, ...
                         Gravity=false, Water=true, PointForce=false, ...
                         ActuationL=true, ActuationR=true, ...
                         CableActuator=LOM, RadialActuator=TM);
ndof_xi = Octopus.ndof_xi;
ndof_rho = Octopus.ndof_rho;

%% initializations
Xs = Octopus.Twists(2).Xs;
nip = Octopus.Twists(2).nip;
n_sact = LOM.get_n_sact();
L = Octopus.Link.L;
E = Octopus.Link.E;
G = Octopus.Link.G;
nu = Octopus.Link.Poi;
A = pi * (Octopus.Link.r_base)^2;
Ke = 0.01 * E * A/L; % passive axial stiffness, Neumann BC, N/cm

Fmax = 1.25;

% TM
Pmax = 8e3; % maximum boundary stress, Pa
u_rho_init = 1;
stiff_len = [];
stiff_force = [];

[UTM, ULM] = meshgrid(0:0.2:1, 0:0.2:1);
UTM = Pmax * UTM;
ULM = Fmax * ULM;

DeltaL = arrayfun(@(lm, tm)getShortening(Octopus,ndof_xi,ndof_rho,n_sact,lm,tm), ULM, UTM);

%% plot stiffness diagram
font_size = 26;
figure(1)
AxialForce = ULM * 4;
interp = scatteredInterpolant(100 * DeltaL(:), AxialForce(:), UTM(:));
pcolor(100 * DeltaL, AxialForce, UTM);
hold on;
shading interp;
colorbarHandle = colorbar;
title(colorbarHandle, '$P$ (Pa)', 'Interpreter', 'latex', 'FontSize', font_size);

% plot passive stiffness
d = 0:0.1:3;
Fe = Ke .* d; % passive F-L curve
h0 = plot(d, Fe, 'r-', 'LineWidth', 3);
hold on;

% plot active stiffness
F1 = 2 * Ke.*d; % 2*Ke
validIdx1 = (F1 >= 0 & F1 <= 5);
tm1 = interp(d(validIdx1), F1(validIdx1));
h1 = plot(d, F1, 'g-','LineWidth', 2);

F2 = Fe + (Ke * 2 .* d.^2);
validIdx2 = (F2 >= 0 & F2 <= 5);
tm2 = interp(d(validIdx2), F2(validIdx2));
h2 = plot(d, F2, 'm-', 'LineWidth', 2);

F3 = linspace(0, 5, length(d));
validIdx3 = (F3 >= 0 & F3 <= 5);
d3 = zeros(1,length(d));
tm3 = interp(d3(validIdx3), F3(validIdx3));
h3 = plot(d3, F3, 'k-', 'LineWidth', 3);

% customize the plot
grid on;
% set(gca,'DataAspectRatio',[5 1 1])
xlim([0, 3]);
ylim([0, 5]);
% title('Contour Plot of Axial Stiffness');
set(gca,'FontSize',font_size, 'FontName', 'Times New Roman',...
    'OuterPosition', [0, 0, 1, 0.97]);
xlabel('$\Delta L$(cm)', 'Interpreter','latex', 'FontSize',font_size);
ylabel('$F$ (N)', 'Interpreter', 'latex', 'FontSize', font_size);
legend([h0, h1, h2, h3], ...
        {'$k_e$', '$2k_e$', '$k_a$', '$k_d$'},...
        'Interpreter', 'latex', 'FontSize', font_size, 'Location', 'southeast');
if ~exist('./figures', 'dir')
    mkdir('./figures');
end
exportgraphics(gcf, './figures/AxialStiff.pdf','ContentType','vector');

%% plot F-TM load diagram
figure(2)
f0 = plot(0:1:5, zeros(1, 6), 'r-', 'LineWidth', 2);
hold on
f1 = plot(F1(validIdx1), tm1, 'g-', 'LineWidth', 2);
hold on
f2 = plot(F2(validIdx2), tm2, 'm-', 'LineWidth', 2);
hold on
f3 = plot(F3(validIdx3), tm3, 'k-', 'LineWidth', 2);
lm0 = 0:0.5:5;
tm0 = lm0/nu;
tm0 = tm0/(2*A);
f4 = plot(lm0, tm0, 'b--', 'LineWidth', 2);
grid on
% title('LM vs TM loads of Tunable Stiffness');
set(gca,'FontSize',font_size, 'FontName', 'Times New Roman');
xlabel('$F$ (N)', 'Interpreter','latex','FontSize',font_size);
ylabel('$P$ (Pa)', 'Interpreter','latex','FontSize',font_size);
legend([f0, f1, f2, f3, f4], ...
        {'$k_e$','$2k_e$', '$k_a$', '$k_d$', '$^a k_d$'},...
        'Interpreter', 'latex','FontSize', font_size, 'Location', 'northwest');
exportgraphics(gcf, './figures/AxialStiffFTM.pdf','ContentType','vector');

%% usefull functions
function LOM = createLOM(OctopusLink)
    % create longitudinal and oblique muscles
    % input
    %   OctopusLink: SorosimLink object
    % output
    %   LM: CableActuation objects
    rt = OctopusLink.r_tip;
    rb = OctopusLink.r_base;

    % define longitudinal muscles on SCALED domain
    LM1_y = @(X)0;
    LM1_z = @(X)0.8*(rb-(rb-rt)*X);


    LM2_y = @(X)0.8*(rb-(rb-rt)*X);
    LM2_z = @(X)0;

    LM3_y = @(X)0;
    LM3_z = @(X)-0.8*(rb-(rb-rt)*X);

    LM4_y = @(X)-0.8*(rb-(rb-rt)*X);
    LM4_z = @(X)0;

    LM1 = Cable(LM1_y, LM1_z);
    LM2 = Cable(LM2_y, LM2_z);
    LM3 = Cable(LM3_y, LM3_z);
    LM4 = Cable(LM4_y, LM4_z);

    LOM = CableActuation(LM1, LM2, LM3, LM4);
end

function TM = createTM()
    act = Radial(0, 1);
    TM = RadialActuation(act);
end

function deltaL = getShortening(Octopus,ndof_xi,ndof_rho, n_sact, lm, tm)
    q0 = zeros(ndof_xi+ndof_rho, 1);
    u_xi = -ones(n_sact, 1) * lm;
    u_rho = -tm;
    q = Octopus.statics(q0, u_xi, u_rho);
    q_xi = q(1:ndof_xi,:);
    q_rho = q(ndof_xi+1:end, :);
    [g, rho] = Octopus.FwdKinematics(q_xi, q_rho);
    deltaL = Octopus.Link.L - g(end-3, 4);
end
