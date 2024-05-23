clc;
clear;
close;
%% starts
% initialization
t = 0:0.01:3;
Xs = 0:0.1:1;
Fmax = 4;
Fmin = 0.20;
Pmax = 20e3;
fend = 1;
bp_s = 0.2;
bp_e = 0.8;
Tp = 2;
LM = [];
TM = [];

% get values
for ts=t
    uLM = LMrelease(ts, Xs', Fmax, Fmin, Tp, bp_s, bp_e);
    uTM = TMcontract(ts, Xs', Pmax, fend, Tp);

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

%% TODO: plot LM and TM over space, label the time
