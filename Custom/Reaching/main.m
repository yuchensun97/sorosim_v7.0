clc;
clear;
close all;
%% create Link
OctopusLink = SorosimLink('Octopus.json');
LOM = createLOM(OctopusLink);
TM = createTM();
Octopus = OctopusArm(OctopusLink, Damped=true, ...
                                  Gravity=false, Water=true, PointForce=false, ...
                                  ActuationL=true, ActuationR=true, ...
                                  CableActuator=LOM, RadialActuator=TM);
ndof_xi = Octopus.ndof_xi;
ndof_rho = Octopus.ndof_rho;

%% reaching 
Fmax = 0.8;
Fmin = 0.05;
fstart = 0.1;
fend = 0.8; % cable force end at fend
Tp = 2;
Xs = Octopus.Twists(2).Xs;
nip = Octopus.Twists(2).nip;
n_sact = LOM.get_n_sact();

% LM
uqt_xi = cell(n_sact, 1);
uqt_xi{1} = @(t)LMrelease(t, Xs, Fmax, 0, Tp, fstart, fend);
for i = 2:n_sact
    uqt_xi{i} = @(t)LMcontract(t, Xs, Fmin, Tp, fstart, fend);
end

u_xi = zeros(nip, n_sact);
u_xi(:, 1) = uqt_xi{1}(0);

% TM
Pmax = 18e2; % maximum boundary stress, Pa
uqt_rho = @(t)TMcontract(t, Xs, Pmax, fend, Tp);

%% statics
% starts from bending position
q0 = zeros(ndof_xi+ndof_rho, 1);
qb = Octopus.statics(q0, u_xi, zeros(nip, 1));
qb_xi = qb(1:ndof_xi,:);
qb_rho = qb(ndof_xi+1:end, :);
% fb = Octopus.plotq(qb_xi, qb_rho);

Bh_xi = Octopus.Twists(2).Bh_xi;
B_xi_dof = Octopus.Twists(2).B_xi_dof;
B_xi_ord = Octopus.Twists(2).B_xi_odr;
xi_star = [0 0 0 1 0 0]';

%% dynamics
dt = 0.01;
tmax = 3;

qqd_r = [qb; zeros(ndof_xi+ndof_rho,1)];
[t, qqd] = Octopus.dynamics(qqd_r, uqt_xi, uqt_rho, 'ode15s', dt, tmax);
save("./Custom/results/reaching.mat", "t", "qqd");
Octopus.plotqqd(t, qqd, 'Octopus_reaching');

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
