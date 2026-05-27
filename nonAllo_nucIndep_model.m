 clearvars
 close all

 % load ez and pho concentration dynamics
 load('~/Ez_Pho_concentrations/ez_dynamics.mat')
 load('~/Ez_Pho_concentrations/pho_dynamics.mat')

 % load ez and pho nM estimates at NC14
 load('~/Ez_Pho_concentrations/EzPho_conc_nM.mat')

 %% Simulation
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
 params.allo_2 = 3;
 params.allo_3 = 11;
 params.V = 7.3500e-04; % V = 7.3500e-04 for the revised paper, based on V_sweep_final_19-Mar-2026 12:06:47
 params.numSims = 100; % n = 100 for revised paper
 params.hitRun = 1;

 
%  params.fiberIC = [2*ones(2,params.nNucleosomes*0.2) 3*ones(2,params.nNucleosomes*0.8)]; % maternal IC
 params.fiberIC = zeros(2,params.nNucleosomes); % paternal IC
 [simArray] = simulateK27(ez,pho,params);
 plotSimulation(simArray,params,'nucleation-independent')


function plotSimulation(simArray,params,titleString)

    tMit = params.tMit;
    stepsPerMin = params.stepsPerMin;
    numSims = params.numSims;
    nNucleosomes = params.nNucleosomes;
    nUpdates = params.nUpdates;

    
    % Find the methylation fractions across all simulations
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
    title(titleString,'fontsize',10)
    
    % l = legend('K27me0','K27me1','K27me3');
    % l.FontSize = 14;
    % l.Location = 'Northwest';
    % l.Box = 'off';
    
    set(f,'units','centimeters','position',[35,25,11.25,4.5])
%     exportgraphics(f,'path to figure',...
%         'Resolution',300,'BackgroundColor','white')

end