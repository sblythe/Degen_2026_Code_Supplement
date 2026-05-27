
clearvars

% load ez and pho concentration dynamics
load('.../Ez_Pho_concentrations/ez_dynamics.mat')
load('.../Ez_Pho_concentrations/pho_dynamics.mat')

% load ez and pho nM estimates at NC14
load('.../Ez_Pho_concentrations/EzPho_conc_nM.mat')

% define parameters
params = struct;
params.nNucleosomes = 300;
params.stepsPerMin = 4;
params.deltaT = 1; % 1 time step = 15 sec
params.NClengths = [9*ones(1,10) 10 12 20 60 65 50];
params.nUpdates = sum(params.NClengths)*params.stepsPerMin; 
params.Ke = 10; % dissociation constant of e(z) binding
params.Kp = 10; % dissociation constant of pho binding
params.tMit = [[0:10]*9 100 112 132 192 257 307].*params.stepsPerMin; % The timepoints of mitosis
params.tRep = params.tMit+4*params.stepsPerMin; % replication happens 4 min following mitosis
params.allo_3 = 11;
params.allo_2 = 3;
params.ez = ez;
params.pho = pho;
params.constantPho_nM = conc_nM.Pho;
params.V = 7.3500e-04; % V = 7.3500e-04 for the revised paper, based on V_sweep_final_19-Mar-2026 12:06:47
params.numSims = 10; % n = 10 for the paper
params.hitRun = 0;

%% Run IC sweep

IC_me1 = 0; % 0 to recreate Fig. SM1A, 0.25 to recreate Fig. SM1B
IC_me3 = 0:0.05:(1-IC_me1);
allOutputs = zeros(1228,4,length(IC_me3));
for i = 1:length(IC_me3)

 if IC_me1 == 0
    params.fiberIC = [2*ones(2,round(params.nNucleosomes*(1-IC_me3(i)))) 3*ones(2,round(params.nNucleosomes*IC_me3(i)))]; % maternal IC
 else
     params.fiberIC = [ones(2,round(params.nNucleosomes*IC_me1)) 2*ones(2,round(params.nNucleosomes*(1-IC_me1-IC_me3(i)))) 3*ones(2,round(params.nNucleosomes*IC_me3(i)))]; % maternal IC
 end
 [simArray rxnsPerTimestep] = simulateK27(ez,pho,params);
%      plotSimulation(simArray,params,'maternal, dynamic pho')

fracZeroOverT = sum(simArray==0,2:4)./(params.nNucleosomes*params.numSims*2);
fracOneOverT = sum(simArray==1,2:4)./(params.nNucleosomes*params.numSims*2);
fracTwoOverT = sum(simArray==2,2:4)./(params.nNucleosomes*params.numSims*2);
fracThreeOverT = sum(simArray==3,2:4)./(params.nNucleosomes*params.numSims*2);

fractions = [fracZeroOverT, fracOneOverT, fracTwoOverT, fracThreeOverT];
allOutputs(:,:,i) = fractions;
end


%% Plot IC sweep results
close all

cmap = turbo(size(allOutputs,3));

measure13  = params.tMit(13)+10*params.stepsPerMin;
measure14E = params.tMit(14)+15*params.stepsPerMin; % NC14E
measure14M = params.tMit(14)+35*params.stepsPerMin; % NC14M
measure14L = params.tMit(14)+60*params.stepsPerMin; % NC14L

me3 = squeeze(allOutputs(:,4,:));
me2 = squeeze(allOutputs(:,3,:));

f = figure;
subplot(2,1,1)
for i = 1:size(me2,2)
    plot(me2(:,i),'linewidth',1,'color',cmap(i,:))
    hold on
end
c = colorbar;
colormap(gca, 'turbo');
% ylabel(c,'me3 IC (fraction of histones)','fontsize',6);
c.Ticks = [0 0.5 1];
c.FontSize = 6;
c.Position = [0.92 0.5922 0.025 0.3216];
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
ylabel('Fraction of histones','fontsize',8)
xlabel('Developmental time (min)','fontsize',8)
title('H3K27me2','fontsize',10)

subplot(2,1,2)
for i = 1:size(me2,2)
    plot(me3(:,i),'linewidth',1,'color',cmap(i,:))
    hold on
end
c = colorbar;
colormap(gca, 'turbo');
% ylabel(c,'me3 IC (fraction of histones)','fontsize',6);
c.Ticks = [0 0.5 1];
c.FontSize = 6;
c.Position = [0.92 0.1098 0.025 0.3294];
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
ylabel('Fraction of histones','fontsize',8)
xlabel('Developmental time (min)','fontsize',8)
title('H3K27me3','fontsize',10)

set(f,'units','centimeters','position',[35,25,12,9])
% exportgraphics(f,'/Volumes/BlytheLab_Files/Imaging/Ellie/Polycomb_Modeling/Natalie_paper/SupplementalFigures/newSuppFig/ICsweep_newModel_noMe1.pdf',...
%     'Resolution',300,'BackgroundColor','white')
