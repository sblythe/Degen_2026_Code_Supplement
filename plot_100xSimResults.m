
% Load maternal and paternal simulation results (that were saved by running simulate_K27me_final.m)
load('.../22-Mar-2026 18:40:31_mat/simArray.mat')
mat_sim = simArray;

load('.../22-Mar-2026 18:41:58_pat/simArray.mat')
pat_sim = simArray;

% Load parameters used to run simulations (should be the same for mat/pat simulations)
load('.../22-Mar-2026 18:40:31_mat/params.mat')

%% Plot and save simulation results (Fig. 2C & C' of Degen et al., 2026)
close all

plotDynamics(mat_sim,'\it Maternal chromosomes','maternal.pdf')
plotDynamics(pat_sim,'\it Paternal chromosomes','paternal.pdf')


%% Count methylation groups & plot mat vs. pat contributions
close all

for i = 1:size(mat_sim,1)

    mat = mat_sim(i,:,:,:);
    pat = pat_sim(i,:,:,:);

    zeros_mat(i) = sum(mat(:)==0);
    ones_mat(i) = sum(mat(:)==1);
    twos_mat(i) = sum(mat(:)==2);
    threes_mat(i) = sum(mat(:)==3);

    zeros_pat(i) = sum(pat(:)==0);
    ones_pat(i) = sum(pat(:)==1);
    twos_pat(i) = sum(pat(:)==2);
    threes_pat(i) = sum(pat(:)==3);
    
end

% Maternal vs. paternal contributions to H3K27me1 (Fig SM7 A of Modeling Supplement of Degen, 2026)
ones_all = ones_mat + ones_pat;
f = figure;
plot(ones_mat./ones_all,'linewidth',1.5,'color',[181, 16, 7]./255)
hold on
plot(ones_pat./ones_all,'linewidth',1.5,'color',[0.7 0.7 0.7])
xlim([-6.5*params.stepsPerMin params.nUpdates])
xticks(params.tMit)
minTicks = params.tMit./params.stepsPerMin;
xticklabels({minTicks(1), '', minTicks(3), '', minTicks(5), '', minTicks(7), '', minTicks(9),...
    '', minTicks(11), '', minTicks(13), minTicks(14), minTicks(15), minTicks(16), minTicks(17)})
ylim([-0.1 1.1])
yticks(0:0.2:1)
ax = gca;
ax.FontSize = 6;
ax.Box = 'off';
ax.LineWidth = 0.5;
xlabel('Developmental time (min)','fontsize',8)
ylabel('Fraction of H3K27me1 groups','fontsize',8)
title('Maternal vs. Paternal H3K27me1','fontsize',8)
set(f,'units','centimeters','position',[35 24.5181 5.85 5])

% exportgraphics(f,'path to figure',...
%         'Resolution',300,'BackgroundColor','white')

% Maternal vs. paternal contributions to H3K27me2 (Fig SM7 B in Modeling Supplement of Degen, 2026)
twos_all = twos_mat + twos_pat;
f = figure;
plot(twos_mat./twos_all,'linewidth',1.5,'color',[181, 16, 7]./255)
hold on
plot(twos_pat./twos_all,'linewidth',1.5,'color',[0.7 0.7 0.7])
xlim([-6.5*params.stepsPerMin params.nUpdates])
xticks(params.tMit)
minTicks = params.tMit./params.stepsPerMin;
xticklabels({minTicks(1), '', minTicks(3), '', minTicks(5), '', minTicks(7), '', minTicks(9),...
    '', minTicks(11), '', minTicks(13), minTicks(14), minTicks(15), minTicks(16), minTicks(17)})
ylim([-0.1 1.1])
yticks(0:0.2:1)
ax = gca;
ax.FontSize = 6;
ax.Box = 'off';
ax.LineWidth = 0.5;
xlabel('Developmental time (min)','fontsize',8)
ylabel('Fraction of H3K27me2 groups','fontsize',8)
title('Maternal vs. Paternal H3K27me2','fontsize',8)
set(f,'units','centimeters','position',[35 24.5181 5.85 5])

% exportgraphics(f,'path to figure',...
%         'Resolution',300,'BackgroundColor','white')


% Maternal vs. paternal contributions to H3K27me3 (Fig 2D & Fig SM7 C of Degen, 2026)
threes_all = threes_mat + threes_pat;
f = figure;
plot(threes_mat./threes_all,'linewidth',1.5,'color',[181, 16, 7]./255)
hold on
plot(threes_pat./threes_all,'linewidth',1.5,'color',[0.7 0.7 0.7])
xlim([-6.5*params.stepsPerMin params.nUpdates])
xticks(params.tMit)
minTicks = params.tMit./params.stepsPerMin;
xticklabels({minTicks(1), '', minTicks(3), '', minTicks(5), '', minTicks(7), '', minTicks(9),...
    '', minTicks(11), '', minTicks(13), minTicks(14), minTicks(15), minTicks(16), minTicks(17)})
ylim([-0.1 1.1])
yticks(0:0.2:1)
ax = gca;
ax.FontSize = 6;
ax.Box = 'off';
ax.LineWidth = 0.5;
xlabel('Developmental time (min)','fontsize',8)
ylabel('Fraction of H3K27me3 groups','fontsize',8)
title('Maternal vs. Paternal H3K27me3','fontsize',8)
set(f,'units','centimeters','position',[35 24.5181 5.85 5])


% exportgraphics(f,'path to figure',...
%         'Resolution',300,'BackgroundColor','white')



%% Plotting Function

function plotDynamics(simArray, titleString, filename)

    NClengths = [9*ones(1,10) 10 12 20 60 65 50];
    numSims = 100; % 100 for the paper
    nNucleosomes = 300;
    stepsPerMin = 4; % time steps of 15 seconds like in lundkvist
    nUpdates = sum(NClengths)*stepsPerMin; 
    tMit = [[0:10]*9 100 112 132 192 257 307].*stepsPerMin;

    % Find the fractions across all simulations
    fracZeroOverT = sum(simArray==0,2:4)./(nNucleosomes*numSims*2);
    fracOneOverT = sum(simArray==1,2:4)./(nNucleosomes*numSims*2);
    fracTwoOverT = sum(simArray==2,2:4)./(nNucleosomes*numSims*2);
    fracThreeOverT = sum(simArray==3,2:4)./(nNucleosomes*numSims*2);
    
    f = figure;
    plot(fracZeroOverT,'linewidth',1.5,'color',[0.75 0.75 0.75])
    hold on
    plot(fracOneOverT,'linewidth',1.5,'color',[0 0.8 0.8])
    hold on
    plot(fracTwoOverT,'linewidth',1.5,'color',[4, 120, 209]./255)
    hold on
    plot(fracThreeOverT,'linewidth',1.5,'color',[2, 41, 66]./255)
%     hold on
%     plot(fracZeroOverT+fracOneOverT+fracTwoOverT+fracThreeOverT)
    
    xlim([-6.5*stepsPerMin nUpdates])
    xticks(tMit)
    minTicks = tMit./stepsPerMin;
    xticklabels({minTicks(1), '', minTicks(3), '', minTicks(5), '', minTicks(7), '', minTicks(9),...
        '', minTicks(11), '', minTicks(13), minTicks(14), minTicks(15), minTicks(16), minTicks(17)})
    ylim([-0.1 1])
    yticks(0:0.2:1)
    
    
    ax = gca;
    ax.FontSize = 6;
    ax.Box = 'off';
    ax.LineWidth = 0.5;
    ylabel('Fraction of histones','fontsize',8)
    xlabel('Developmental time (min)','fontsize',8)
    title(titleString,'fontsize',8)
    
    % l = legend('K27me0','K27me1','K27me3');
    % l.FontSize = 14;
    % l.Location = 'Northwest';
    % l.Box = 'off';
    
    set(f,'units','centimeters','position',[35,25,11.25,4.5])
%     exportgraphics(f,strcat('path to figure',filename),...
%         'Resolution',300,'BackgroundColor','white')

end






