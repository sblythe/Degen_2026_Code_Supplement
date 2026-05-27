
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
params.ez = ez;
params.pho = pho;
params.constantPho_nM = conc_nM.Pho;
params.V = 2*7.3500e-04; % V = 7.3500e-04 for the revised paper, based on V_sweep_final_19-Mar-2026 12:06:47
params.numSims = 10; % n = 10 for the paper
params.allo_3 = 11; % baseline allo_3 parameter
params.allo_2 = 3; % baseline allo_2 parameter
params.hitRun = 0;

%% Perform allosteric parameter sweep

allo_2 = 1:5;
multiply = 1:10;
allOutputs = zeros(1228,4,length(allo_2),length(multiply));
for i = 1:length(allo_2)

    for j = 1:length(multiply)

        params.fiberIC = [2*ones(2,params.nNucleosomes*0.2) 3*ones(2,params.nNucleosomes*0.8)]; % maternal IC
        params.allo_2 = allo_2(i); % dissociation constant of e(z) binding
        params.allo_3 = multiply(j)*allo_2(i); % dissociation constant of pho binding
        [simArray rxnsPerTimestep] = simulateK27(ez,pho,params);
        
        fracZeroOverT = sum(simArray==0,2:4)./(params.nNucleosomes*params.numSims*2);
        fracOneOverT = sum(simArray==1,2:4)./(params.nNucleosomes*params.numSims*2);
        fracTwoOverT = sum(simArray==2,2:4)./(params.nNucleosomes*params.numSims*2);
        fracThreeOverT = sum(simArray==3,2:4)./(params.nNucleosomes*params.numSims*2);
        
        fractions = [fracZeroOverT, fracOneOverT, fracTwoOverT, fracThreeOverT];
        allOutputs(:,:,i,j) = fractions;
    end
end


%% Make heatmaps of me2 and me3 at NC8 and NC14M (Fig. SM4 of Modeling Supplement)
close all

me2 = squeeze(allOutputs(:,3,:,:));
me2_8 = squeeze(me2(params.tMit(8)+4*params.stepsPerMin,:,:));
me2_14M = squeeze(me2(params.tMit(14)+35*params.stepsPerMin,:,:));

me3 = squeeze(allOutputs(:,4,:,:));
me3_8 = squeeze(me3(params.tMit(8)+4*params.stepsPerMin,:,:));
me3_14M = squeeze(me3(params.tMit(14)+35*params.stepsPerMin,:,:));

f = figure;
s1 = subplot(2,1,1);
imagesc(flip(me2_8',1))
ax = gca;
ax.FontSize = 6;
ax.Box = 'off';
ax.LineWidth = 0.5;
yticks(1:2:9)
yticklabels(flip(2:2:10))
ylabel('allo_3/allo_2','fontsize',8)
xlabel('allo_2','fontsize',8)
title('H3K27me2 at NC8','fontsize',10)
c = colorbar;
caxis([0 1])
c.Label.String = 'Fraction of Histones';
c.Label.FontSize = 8;

s2 = subplot(2,1,2);
imagesc(flip(me2_14M',1))
ax = gca;
ax.FontSize = 6;
ax.Box = 'off';
ax.LineWidth = 0.5;
yticks(1:2:9)
yticklabels(flip(2:2:10))
ylabel('allo_3/allo_2','fontsize',8)
xlabel('allo_2','fontsize',8)
title('H3K27me2 at NC14M','fontsize',10)
c = colorbar;
caxis([0 1])
c.Label.String = 'Fraction of Histones';
c.Label.FontSize = 8;

f.Position = [973 684 184 299];
% exportgraphics(f,'path to figure',...
%     'Resolution',300,'BackgroundColor','white')

f = figure;
s3 = subplot(2,1,1);
imagesc(flip(me3_8',1))
ax = gca;
ax.FontSize = 6;
ax.Box = 'off';
ax.LineWidth = 0.5;
yticks(1:2:9)
yticklabels(flip(2:2:10))
ylabel('allo_3/allo_2','fontsize',8)
xlabel('allo_2','fontsize',8)
title('H3K27me3 at NC8','fontsize',10)
c = colorbar;
caxis([0 1])
c.Label.String = 'Fraction of Histones';
c.Label.FontSize = 8;

s4 = subplot(2,1,2);
imagesc(flip(me3_14M',1))
ax = gca;
ax.FontSize = 6;
ax.Box = 'off';
ax.LineWidth = 0.5;
yticks(1:2:9)
yticklabels(flip(2:2:10))
ylabel('allo_3/allo_2','fontsize',8)
xlabel('allo_2','fontsize',8)
title('H3K27me3 at NC14M','fontsize',10)
c = colorbar;
caxis([0 1])
c.Label.String = 'Fraction of Histones';
c.Label.FontSize = 8;

f.Position = [973 684 184 299];
% exportgraphics(f,'path to figure',...
%     'Resolution',300,'BackgroundColor','white')
