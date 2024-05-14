clc;
close;
%% starts
% initialization
t = 0:0.01:8;
Xs = 0:0.25:1;
Fmax = 3.5;
Pmax = 16e3;
fend = 0.7;
LM = [];
TM = [];

% get values
for ts=t
    uLM = LMrelease(ts, Xs, Fmax, fend);
    uTM = TMcontract(ts, Xs, Fmax, fend);

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
title('LM release propangation');

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
title('TM contract propangation');
exportgraphics(gcf, './figures/TMcontract.pdf','ContentType','vector');
