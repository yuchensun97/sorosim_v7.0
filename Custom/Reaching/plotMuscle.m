clc;
clear;
close;
%% starts
% initialization
t = 0:0.01:8;
Xs = 0:0.2:1;
Fmax = 2.5;
Pmax = 16e3;
fend = 1;
bp_s = 0.5;
bp_e = 1;
Tp = 3;
LM = [];
TM = [];

% get values
for ts=t
    uLM = LMrelease(ts, Xs, Fmax, Tp, bp_s, bp_e);
    uTM = TMcontract(ts, Xs, Pmax, fend);

    LM = [LM uLM];
    TM = [TM uTM];
end

%% plot
figure(1)
for i=1:length(Xs)
    curr_LM = -LM(i, :);
    LMinfo{i} = ['s = ' num2str(Xs(i))];
    plot(t, curr_LM);
    hold on
end
grid on;
legend(LMinfo);
xlabel('time');
ylabel('F');
title('LM release propagation');

if ~exist('./figures', 'dir')
    mkdir('./figures');
end
exportgraphics(gcf, './figures/LMrelease.pdf','ContentType','vector');

figure(2)
for j=1:length(Xs)
    curr_TM = -TM(j, :);
    TMinfo{j} = ['s = ' num2str(Xs(j))];
    plot(t, curr_TM);
    hold on
end
grid on;
legend(TMinfo);
xlabel('time');
ylabel('P');
title('TM contract propagation');
exportgraphics(gcf, './figures/TMcontract.pdf','ContentType','vector');
