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

%% fetching
Xs = Octopus.Twists(2).Xs;
nip = Octopus.Twists(2).nip;
n_sact = LOM.get_n_sact();

% TODO: test statics
uqt_xi = cell(n_sact, 1);
uqt_xi{1} = @(t)LMrelease(t, Xs, 0.25, 0.02, 2.5, 0.3, 0.9);
for i = 2:4
    uqt_xi{i} = @(t)LMcontract(t, Xs, 0.218, 0.05, 2.5, 0.183, 0.65);
end

u_xi = zeros(nip, n_sact);
for i = 1:n_sact
    u_xi(:, i) = uqt_xi{i}(0);
end

% TM
Pmax = 8e2; % maximum boundary stress, Pa
uqt_rho = @(t)TMcontract(t, Xs, Pmax, 3, 0.28, 0.7);
uqt_rho = @(t)zeros(nip, 1);
u_rho = uqt_rho(0);

q0 = zeros(ndof_xi+ndof_rho, 1);
qb = Octopus.statics(q0, u_xi, u_rho);
qb_xi = qb(1:ndof_xi,:);
qb_rho = qb(ndof_xi+1:end, :);
fb = Octopus.plotq(qb_xi, qb_rho);
[g, rho] = Octopus.FwdKinematics(qb_xi, qb_rho);

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

    % TODO: setup oblique muscles
    OM1_x = @(X)0.8*(rb-(rb-rt)*X) * cos(- 12 * pi * X);
    OM1_y = @(X)0.8*(rb-(rb-rt)*X) * sin(- 12 * pi * X);

    OM2_x = @(X)0.8*(rb-(rb-rt)*X) * cos(pi + 12 * pi * X);
    OM2_y = @(X)0.8*(rb-(rb-rt)*X) * sin(pi + 12 * pi * X);

    OM1 = Cable(OM1_x, OM1_y);
    OM2 = Cable(OM2_x, OM2_y);

    LOM = CableActuation(LM1, LM2, LM3, LM4, OM1, OM2);
end

function TM = createTM()
    act = Radial(0, 1);
    TM = RadialActuation(act);
end
