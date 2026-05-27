clearvars

% Load target data
load('~/ChIP_data_for_fitting/ChIPseq_NC14_maxes_3.mat') % me3 max ChIP signal
maxes_3 = maxes;
load('~/ChIP_data_for_fitting/ChIPseq_NC14_mins_1.mat') % me1 min ChIP signal
mins_1 = mins;

% Load V sweep data
load('~/V_sweep/V_sweep_final_19-Mar-2026 12:06:47/V_vals.mat')
load('~/V_sweep/V_sweep_final_19-Mar-2026 12:06:47/simArrays.mat')

%% Compare simulation results with target data

stepsPerMin = 4; % time steps of 15 seconds like in lundkvist
tMit = [[0:10]*9 100 112 132 192 257 307].*stepsPerMin;

% Specify timing of the measurements
measure14E = tMit(14)+15*stepsPerMin; % NC14E
measure14M = tMit(14)+35*stepsPerMin; % NC14M
measure14L = tMit(14)+60*stepsPerMin; % NC14L


data_ratios_3 = maxes_3./maxes_3(end);
data_ratios_1 = mins_1./mins_1(1);

ratios_all_3 = zeros(length(simArrays),3); % first col: nc14E, second col: nc14M, third col: nc14L
diffs_all_3 = zeros(length(simArrays),3);
ratios_all_1 = zeros(length(simArrays),3); % first col: nc14E, second col: nc14M, third col: nc14L
diffs_all_1 = zeros(length(simArrays),3);
for i = 1:length(simArrays)

    simArray = simArrays{i};
    threeOverT = sum(simArray==3,2:4);
    oneOverT = sum(simArray==1,2:4);

    ratios_3 = threeOverT([measure14E, measure14M, measure14L])./threeOverT(measure14L);
    ratios_1 = oneOverT([measure14E, measure14M, measure14L])./oneOverT(measure14E);

    ratios_all_3(i,:) = ratios_3;
    diffs_all_3(i,:) = ratios_3-data_ratios_3;
    ratios_all_1(i,:) = ratios_1;
    diffs_all_1(i,:) = ratios_1-data_ratios_1;

end

%% Plot sweep results based on me3 (Fig. SM3 B of Modeling Supplement)
close all

mean_diffs_3 = mean(abs(diffs_all_3),2);

f = figure;
plot(V_vals, mean_diffs_3,'linewidth',1.5,'color',[250, 40, 7]./255)
ylim([0 1])
xlim([0.4e-03 1.2e-03])
xticks([0.4:0.2:1.2].*10^-3)
ax = gca;
ax.FontSize = 6;
ax.Box = 'off';
ax.LineWidth = 0.5;
xlabel('{\it V}','fontsize',8)
ylabel('< |Modeled - measured relative me3| >','fontsize',8)
title('Parameter sweep results','fontsize',8)
set(f,'units','centimeters','position',[27.9753 23.9889 5.9619 5.3269])

% exportgraphics(f,'...path to figure/sweep_results_me3.pdf',...
% 'Resolution',300,'BackgroundColor','white')

Vs = V_vals(21:end);
Vs(find(mean_diffs_3(21:end)==min(mean_diffs_3(21:end))))*10000

%% Plot sweep results based on me1 (Fig. SM3 A of Modeling Supplement)
close all

mean_diffs_1 = mean(abs(diffs_all_1),2);

f = figure;
plot(V_vals, mean_diffs_1,'linewidth',1.5,'color',[250, 40, 7]./255)
ylim([0 1])
xlim([0.4e-03 1.2e-03])
xticks([0.4:0.2:1.2].*10^-3)
ax = gca;
ax.FontSize = 6;
ax.Box = 'off';
ax.LineWidth = 0.5;
xlabel('{\it V}','fontsize',8)
ylabel('< |Modeled - measured relative me1| >','fontsize',8)
title('Parameter sweep results','fontsize',8)
set(f,'units','centimeters','position',[27.9753 23.9889 5.9619 5.3269])

% exportgraphics(f,'...path to figure/sweep_results_me1.pdf',...
% 'Resolution',300,'BackgroundColor','white')

Vs = V_vals(21:end);
Vs(find(mean_diffs_1(21:end)==min(mean_diffs_1(21:end))))*10000

