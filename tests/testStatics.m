clc;
clear;
rng(5);

% initialization
B_xi = [1 1 1 1 1 1;
        0 0 0 0 0 0]';

B_rho = [1 1];

L = createLinkage(B_xi, B_rho);
ndof_xi = L.ndof_xi;
ndof_rho = L.ndof_rho;
% qu = [0.584 0.985 0.4569 0.46289 0.2168 0.8543]';
% err = L.equilibrium(qu, 0, 0);

q = L.statics(zeros(ndof_rho+ndof_xi,1));
q_xi = q(1:ndof_xi,:);
q_rho = q(ndof_xi+1:end,:);
[g,rho]=L.FwdKinematics(q_xi,q_rho);
% f = L.plotq(q_xi, q_rho);

% compute volume
V = 0;
np = length(rho)-1;
for ii=1:np
    rb = 0.02*rho(ii);
    rt = 0.02*rho(ii+1);
    h = g(4*ii+1,4)-g(4*(ii-1)+1,4);
    V = V+tCone(rb, rt,h);
end

function L = createLinkage(B_xi, B_rho)
    S = SorosimLink();
    S.basisType = 'hermite full';
    S.L = 0.5;
    S.B_xi = B_xi;
    S.B_rho = B_rho;
    L = SorosimLinkage(S);
end

function V = tCone(rb, rt, h)
    V = pi * h * (rb^2+rb*rt+rb^2)/3;
end
