% Specify the file path
clc;
clear;

filePath = './Custom/results/vol.csv';

% Read the CSV file
data = readmatrix(filePath);
t = data(:, 1);
noTM = 100 * data(:, 2);
classic = 100 * data(:, 3);
TM = 100 * data(:, 4);

figure;
f1 = plot(t, classic, 'LineWidth', 2);
hold on
f2 = plot(t, noTM, 'LineWidth', 2);
hold on
f3 = plot(t, TM, 'LineWidth', 2);
hold off
grid on
set(gca,'FontSize',28, 'FontName', 'Times New Roman');
legend('classic', 'extended', 'extended+TM', 'Location', 'best', 'FontSize', 24);
xlabel('$t$ (s)', 'Interpreter','latex');
ylabel('$\Delta V$ (\%)', 'Interpreter','latex');

exportgraphics(gcf, './figures/vol_change.pdf','ContentType','vector');

% Display the data
% disp(data);