clc;
clear;
rng(23);

% initialization
B_xi = cell(8);
B_xi{1} = [1 1 1 1 1 1; 
           0 0 0 0 0 0]'; 
B_xi{2} = [1 1 1 1 1 1; 
           1 1 1 1 1 1]';
B_xi{3} = [1 1 1 1 1 1;
           2 2 2 2 2 2]';
B_xi{4} = [1 1 1 1 1 1;
           0 1 1 2 1 0]';
for i = 5:7
    B_xi{i} = B_xi{1};
end
B_xi{8} = B_xi{4};

B_rho = cell(8);
for i=1:4
    B_rho{i} = [0 0];
end
B_rho{5} = [1 0];
B_rho{6} = [1 1];
B_rho{7} = [1 2];
B_rho{8} = [1 2];

for i=1:4
    S = createLink(B_xi{i}, B_rho{i});
    [L, q_xi, q_rho] = createLinkage(S);
    f = L.plotq(q_xi, q_rho);
%     filename = sprintf('./figures/plotq_%d.pdf', i);
%     saveas(gcf, filename);
    close(f);
end

function S = createLink(B_xi, B_rho)
    S = SorosimLink();
    S.B_xi = B_xi;
    S.B_rho = B_rho;
end

function [L, q_xi, q_rho] = createLinkage(S)
    L = SorosimLinkage(S);
    dof_xi = L.Twists(2).dof_xi;
    dof_rho = L.Twists(2).dof_rho;
    q_xi = rand(dof_xi, 1);
    q_rho = rand(dof_rho, 1);
end
