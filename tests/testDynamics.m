clc;
clear;
rng(5);

% initialization
B_xi = [1 1 1 1 1 1;
        0 0 0 0 0 0]';

B_rho = [1 1];

L = createLinkage(B_xi, B_rho);
% qqd(1:5) = qqd(1:5)/L.Link.L;
% qqd(12:16) = qqd(12:16)/L.Link.L;
% dqdt = L.derivatives(0,qqd,0,0);
[t, qqd] = L.dynamics();
ndof_xi = L.ndof_xi;
ndof_rho = L.ndof_rho;
q_xi = qqd(:,1:ndof_xi);
q_rho = qqd(:,ndof_xi+1:ndof_xi+ndof_rho);

%fwd kinemtics per .5 sec
% figure;
% j = 1;
% for ii=1:50:length(t)
%     t_now = t(ii);
%     q_xi_now = q_xi(ii,:);
%     q_rho_now = q_rho(ii,:);
%     fh = L.plotq(q_xi_now, q_rho_now);
%     legendinfo{j} = ['t = ' num2str(t_now)];
%     j=j+1;
% end
% legend(legendinfo);

function L = createLinkage(B_xi, B_rho)
    S = SorosimLink();
    S.B_xi = B_xi;
    S.B_rho = B_rho;
    L = SorosimLinkage(S);
end
