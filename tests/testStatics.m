clc;
clear;
rng(5);

% initialization
B_xi = [1 1 1 1 1 1;
        2 2 2 1 0 0]';

B_rho = [1 1];

L = createLinkage(B_xi, B_rho);
ndof_xi = L.ndof_xi;
ndof_rho = L.ndof_rho;

q = L.statics(zeros(ndof_rho+ndof_xi,1));
q_xi = q(1:ndof_xi,:);
q_rho = q(ndof_xi+1:end,:);
[g,rho]=L.FwdKinematics(q_xi,q_rho);
f = L.plotq(q_xi, q_rho);

function L = createLinkage(B_xi, B_rho)
    S = SorosimLink();
    S.basisType = 'mixed';
    S.L = 0.5;
    S.B_xi = B_xi;
    S.B_rho = B_rho;
    Fp_loc = [0.2, 0.5, 1]';
    LocalForce = [false, false, false]';
    Fp_vec_1 = @(t)[0 0 0 0 0 -5]';
    Fp_vec_2 = @(t)[0 0 0 -2 -3 5]';
    Fp_vec_3 = @(t)[0 0 0 0 5 5]';
    Fp_vec = {Fp_vec_1; Fp_vec_2; Fp_vec_3};
%       Fp_loc = [1];
%       LocalForce = [false];
%       Fp_vec_h = @(t)[0 0 0 -15 0 0]';
%       Fp_vec = {Fp_vec_h};
    L = SorosimLinkage(S, PointForce=true,...
                          Fp_loc = Fp_loc,...
                          LocalForce = LocalForce,...
                          Fp_vec = Fp_vec);
end
