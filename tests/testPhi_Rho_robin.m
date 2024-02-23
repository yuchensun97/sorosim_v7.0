clc;
clear;

X = linspace(0, 1, 50); % X varies from 0 to 1
Bdof = 1;
Bodr = 2;
% Remove the line that assigns an empty array to 'results'

results = Phi_Rho_Hermitian_robin(X(1), Bdof, Bodr);
for i = 2:length(X)
    result = Phi_Rho_Hermitian_robin(X(i), Bdof, Bodr);
    results = [results; result];
end

% plot results for each colounm
figure;
for i = 1:size(results, 2)
    plot(X, results(:, i));
    hold on;
end
grid on;
% add legend to each line
for i = 1:size(results, 2)
    legendInfo{i} = ['n = ' num2str(i)];
end
legend(legendInfo);

% add title and labels
title('Phi Rho robin');
xlabel('X');
ylabel('Phi Rho');

% save figure as .pdf in ${workspaceDir}/figures
% if the folder does not exist, create it
if ~exist('./figures', 'dir')
    mkdir('./figures');
end
exportgraphics(gcf, './figures/Phi_Rho_robin.pdf','ContentType','vector');

