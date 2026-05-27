
clearvars

% load ez and pho concentration dynamics
load('~/Ez_Pho_concentrations/ez_dynamics.mat')
load('~/Ez_Pho_concentrations/pho_dynamics.mat')

% load ez and pho nM estimates at NC14
load('~/Ez_Pho_concentrations/EzPho_conc_nM.mat')

% define parameters
params = struct;
params.nNucleosomes = 300;
params.stepsPerMin = 4;
params.deltaT = 1; % 1 time step = 15 sec
params.NClengths = [9*ones(1,10) 10 12 20 60 65 50];
params.nUpdates = sum(params.NClengths)*params.stepsPerMin; 
params.tMit = [[0:10]*9 100 112 132 192 257 307].*params.stepsPerMin; % The timepoints of mitosis
params.tRep = params.tMit+4*params.stepsPerMin; % replication happens 4 min following mitosis
params.allo_3 = 11;
params.allo_2 = 3;
params.ez = ez;
params.pho = pho;
params.constantPho_nM = conc_nM.Pho;
params.V = 7.3500e-04; % V = 7.3500e-04 for the revised paper, based on V_sweep_final_19-Mar-2026 12:06:47
params.numSims = 10; % n = 100 for the paper
params.hitRun = 0;

% Kds to sweep through
Ke = 5:5:75; 
Kp = 5:5:75;

%% Run Kd parameter sweep

allOutputs = zeros(1228,4,length(Ke),length(Kp));
for i = 1:length(Ke)

    for j = 1:length(Kp)
         % Maternal IC, dynamic Pho
         params.fiberIC = [2*ones(2,params.nNucleosomes*0.2) 3*ones(2,params.nNucleosomes*0.8)]; % maternal IC
         params.Ke = Ke(i); % dissociation constant of e(z) binding
         params.Kp = Kp(j); % dissociation constant of pho binding
         [simArray rxnsPerTimestep] = simulateK27(ez,pho,params);
        
        fracZeroOverT = sum(simArray==0,2:4)./(params.nNucleosomes*params.numSims*2);
        fracOneOverT = sum(simArray==1,2:4)./(params.nNucleosomes*params.numSims*2);
        fracTwoOverT = sum(simArray==2,2:4)./(params.nNucleosomes*params.numSims*2);
        fracThreeOverT = sum(simArray==3,2:4)./(params.nNucleosomes*params.numSims*2);
        
        fractions = [fracZeroOverT, fracOneOverT, fracTwoOverT, fracThreeOverT];
        allOutputs(:,:,i,j) = fractions;
    end
end

%% Make heatmaps of me2 and me3 at NC8 and NC14M (Fig. SM5 of Modeling Supplement)

close all

me2 = squeeze(allOutputs(:,3,:,:));
me2_8 = squeeze(me2(params.tMit(8)+4*params.stepsPerMin,:,:));
me2_14M = squeeze(me2(params.tMit(14)+35*params.stepsPerMin,:,:));

me3 = squeeze(allOutputs(:,4,:,:));
me3_8 = squeeze(me3(params.tMit(8)+4*params.stepsPerMin,:,:));
me3_14M = squeeze(me3(params.tMit(14)+35*params.stepsPerMin,:,:));

f = figure;
subplot(2,1,1)
imagesc(me2_8) % right now removing first row/col bc I tested ke = kp = 0 by accident
xticks(1:15) % 1:8 for short sweep. 1:10 for long sweep
yticks(1:15)
xticklabels(Kp)
yticklabels(Ke)
c = colorbar;
caxis([0 1])
ax = gca;
ax.FontSize = 6;
ax.Box = 'off';
ax.LineWidth = 0.5;
xlabel('K_p (nM)','fontsize',8)
ylabel('K_e (nM)','fontsize',8)
title('H3K27me2 at NC8','fontsize',10)
c = colorbar;
caxis([0 1])
c.Label.String = 'Fraction of Histones';
c.Label.FontSize = 8;

subplot(2,1,2)
imagesc(me2_14M)
xticks(1:15) % 1:8 for short sweep. 1:10 for long sweep
yticks(1:15)
xticklabels(Kp)
yticklabels(Ke)
c = colorbar;
caxis([0 1])
ax = gca;
ax.FontSize = 6;
ax.Box = 'off';
ax.LineWidth = 0.5;
xlabel('K_p (nM)','fontsize',8)
ylabel('K_e (nM)','fontsize',8)
title('H3K27me2 at NC14M','fontsize',10)
c = colorbar;
caxis([0 1])
c.Label.String = 'Fraction of Histones';
c.Label.FontSize = 8;

f.Position = [973 684 184 299];
% exportgraphics(f,'/Volumes/BlytheLab_Files/Imaging/Ellie/Polycomb_Modeling/Natalie_paper/forPaper_final/addressingReviews/Kd_sweep/KdSweep_me2.pdf',...
%     'Resolution',300,'BackgroundColor','white')

f = figure;
subplot(2,1,1)
imagesc(me3_8)
xticks(1:15) % 1:8 for short sweep. 1:10 for long sweep
yticks(1:15)
xticklabels(Kp)
yticklabels(Ke)
c = colorbar;
caxis([0 1])
ax = gca;
ax.FontSize = 6;
ax.Box = 'off';
ax.LineWidth = 0.5;
xlabel('K_p (nM)','fontsize',8)
ylabel('K_e (nM)','fontsize',8)
title('H3K27me3 at NC8','fontsize',10)
c = colorbar;
caxis([0 1])
c.Label.String = 'Fraction of Histones';
c.Label.FontSize = 8;

subplot(2,1,2)
imagesc(me3_14M)
xticks(1:15) % 1:8 for short sweep. 1:10 for long sweep
yticks(1:15)
xticklabels(Kp)
yticklabels(Ke)
c = colorbar;
caxis([0 1])
ax = gca;
ax.FontSize = 6;
ax.Box = 'off';
ax.LineWidth = 0.5;
xlabel('K_p (nM)','fontsize',8)
ylabel('K_e (nM)','fontsize',8)
title('H3K27me3 at NC14M','fontsize',10)
c = colorbar;
caxis([0 1])
c.Label.String = 'Fraction of Histones';
c.Label.FontSize = 8;

f.Position = [973 684 184 299];
% exportgraphics(f,'/Volumes/BlytheLab_Files/Imaging/Ellie/Polycomb_Modeling/Natalie_paper/forPaper_final/addressingReviews/Kd_sweep/KdSweep_me3.pdf',...
%     'Resolution',300,'BackgroundColor','white')
