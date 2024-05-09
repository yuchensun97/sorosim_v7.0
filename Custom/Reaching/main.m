clc;
clear;

%% create Link
OctopusLink = SorosimLink('Octopus.json');
LOM = createLOM(OctopusLink);
Octopus = OctopusArm(OctopusLink, Damped=true, ...
                                  Gravity=false, PointForce=false, ...
                                  ActuationL=true, ActuationR=false, ...
                                  CableActuator=LOM);
ndof_xi = Octopus.ndof_xi;
ndof_rho = Octopus.ndof_rho;

%% assign actuation load
Fmax = 30;
Xs = Octopus.Twists(2).Xs;
nip = Octopus.Twists(2).nip;
n_sact = LOM.get_n_sact();
u_xi = zeros(nip, n_sact);
u_xi(:, 1) = -Fmax * (ones(nip, 1) - Xs/0.7);
u_xi(Xs>0.7, 1)= 0;

%% statics
q0 = zeros(ndof_xi+ndof_rho, 1);
q_static = Octopus.statics(q0, u_xi, 0);
q_xi = q_static(1:ndof_xi,:);
q_rho = q_static(ndof_xi+1:end,:);
f = Octopus.plotq(q_xi, q_rho);

%% TODO: dynamics

%% usefull functions
function LOM = createLOM(OctopusLink)
    % create longitudinal and oblique muscles
    % input
    %   OctopusLink: SorosimLink object
    % output
    %   LM: CableActuation objects
    rt = OctopusLink.r_tip;
    rb = OctopusLink.r_base;
    L = OctopusLink.L;

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

