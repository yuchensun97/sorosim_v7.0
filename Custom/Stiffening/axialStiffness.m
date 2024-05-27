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

%% initializations, see testCable.m
Xs = Octopus.Twists(2).Xs;
nip = Octopus.Twists(2).nip;
n_sact = LOM.get_n_sact();
L = Octopus.Link.L;

u_xi_init = ones(n_sact, 1);

% TM
Pmax = 20e3; % maximum boundary stress, Pa
u_rho_init = 1;
stiff_len = [];
stiff_force = [];

i = 1;
for tm = 0:0.2:1
    P = -Pmax*tm;
    u_rho = P * u_rho_init;

    shorten = [];
    force = [];
    for lm=0:0.1:1
        u_xi = -lm * u_xi_init;
        q0 = zeros(ndof_xi+ndof_rho, 1);
        q = Octopus.statics(q0, u_xi, u_rho);
        q_xi = q(1:ndof_xi,:);
        q_rho = q(ndof_xi+1:end, :);
        [g, ~] = Octopus.FwdKinematics(q_xi, q_rho);
        deltaL = L - g(end-3, 4);
        shorten = [shorten; deltaL];
        force = [force; lm*4];
    end
    stiff_len = [stiff_len, shorten];
    stiff_force = [stiff_force, force];

    stiff_info{i} = ['P = ' num2str(-P)];
    i = i + 1;
end

%% plot
figure(1)
for j=1:i-1
    plot(stiff_len(:, j), stiff_force(:, j));
    hold on
end
grid on;
legend(stiff_info);
xlabel('\delta L (m)');
ylabel('F (N)');
title('Stiffening effect');
if ~exist('./figures', 'dir')
    mkdir('./figures');
end
exportgraphics(gcf, './figures/AxialStiff.pdf','ContentType','vector');

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

    % TODO: create oblique muscles in future

end

function TM = createTM()
    act = Radial(0, 1);
    TM = RadialActuation(act);
end
