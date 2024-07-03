clc;
clear;
close;
%% starts
% initialization
t = 0:0.25:3.5;
Xs = 0:0.01:1;
Tp = 2;
LM = [];
LMrest = [];
TM = [];

% get values
for ts=t
    uLM = LMrelease(ts, Xs', 0.2, 0.02, 2.5, 0.3, 0.9);
    oLM = LMcontract(ts, Xs', 0.218, 0.05, 2.5, 0.183, 0.65);
    uTM = TMcontract(ts, Xs', 8e2, 3, 0.28, 0.7);

    LM = [LM uLM];
    LMrest = [LMrest oLM];
    TM = [TM uTM];
end

%% plot
font_size = 20;
figure(1)
for i=1:length(t)
    curr_LM = -LM(:, i);
    % LMinfo{i} = ['s = ' num2str(Xs(i))];
    plot(Xs, curr_LM, 'b-','LineWidth',2);
    hold on
end
grid on;
% legend(LMinfo);
ylim([0, 0.25]);
xlabel('$X=s/L$', 'Interpreter','latex','FontSize', font_size);
ylabel('$u_1$', 'Interpreter','latex', 'FontSize', font_size);
% title('LM 1 release propagation');

if ~exist('./figures', 'dir')
    mkdir('./figures');
end
exportgraphics(gcf, './figures/LMrelease.pdf','ContentType','vector');

figure(2)
for i=1:length(t)
    curr_LMrest = -LMrest(:, i);
    % LMinfo{i} = ['s = ' num2str(Xs(i))];
    plot(Xs, curr_LMrest,'b-', 'LineWidth',2);
    hold on
end
grid on;
% legend(LMinfo);
xlabel('$X=s/L$', 'Interpreter','latex', 'FontSize', font_size);
ylabel('$u_2$', 'Interpreter','latex', 'FontSize', font_size);
% title('LM contract propagation');
exportgraphics(gcf, './figures/LMcontract.pdf','ContentType','vector');

figure(3)
for i=1:length(t)
    curr_TM = -TM(:, i);
    % TMinfo{j} = ['s = ' num2str(Xs(j))];
    plot(Xs, curr_TM, 'b-', 'LineWidth', 2);
    hold on
end
grid on;
% legend(TMinfo);
xlabel('$X=s/L$', 'Interpreter','latex', 'FontSize', font_size);
ylabel('$u_3$', 'Interpreter','latex', 'FontSize', font_size);
% title('TM contract propagation');
exportgraphics(gcf, './figures/TMcontract.pdf','ContentType','vector');
