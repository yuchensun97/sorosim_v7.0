clc;
clear;

% initialization
B_xi = [1 1 1 1 1 1;
        0 1 1 2 1 0]';

B_rho = [1 2];

L = createLinkage(B_xi, B_rho);
[t, qqd] = L.dynamics();
q_xi = qqd(:,1);
q_rho = qqd(:,2);

%fwd kinemtics per .5 sec
figure;
j = 1;
for ii=1:50:length(t)
    t_now = t(ii);
    q_xi_now = q_xi(ii);
    q_rho_now = q_rho(ii);
    fh = L.plotq(q_xi_now, q_rho_now);
    legendinfo{j} = ['t = ' num2str(t_now)];
    j=j+1;
end
legend(legendInfo);

function L = createLinkage(B_xi, B_rho)
    S = SorosimLink();
    S.B_xi = B_xi;
    S.B_rho = B_rho;
    L = SorosimLinkage(S);
end
