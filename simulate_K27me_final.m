 clearvars

 % load ez and pho concentration dynamics
 load('.../Ez_Pho_concentrations/ez_dynamics.mat') 
 load('.../Ez_Pho_concentrations/pho_dynamics.mat') 

 % load ez and pho nM estimates at NC14
 load('.../Ez_Pho_concentrations/EzPho_conc_nM.mat') 
 
 %% Run 100x simulations for Fig. 2 of Degen et al., 2026
 close all

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
 params.numSims = 100; % n = 100 for paper

 
 % Maternal IC
 params.hitRun = 0;
 params.fiberIC = [2*ones(2,params.nNucleosomes*0.2) 3*ones(2,params.nNucleosomes*0.8)]; % maternal IC
 [simArray rxnsPerTimestep] = simulateK27(ez,pho,params);
 plotSimulation(simArray,params,'Maternal')

 newDirectory = 'path to folder to save maternal simulation results' ...
     + string(datetime) + '_mat';
 mkdir(newDirectory)
 save(newDirectory+'/params.mat','params')
 save(newDirectory+'/simArray.mat','simArray', '-v7.3')
 save(newDirectory+'/rxnsPerTimestep.mat','rxnsPerTimestep', '-v7.3')

 % Paternal IC
 params.hitRun = 0;
 params.fiberIC = [0*ones(2,params.nNucleosomes)]; % paternal IC
 [simArray rxnsPerTimestep] = simulateK27(ez,pho,params);
 plotSimulation(simArray,params,'Paternal')

 newDirectory = 'path to folder to save paternal simulation results'...
     + string(datetime) + '_pat';
 mkdir(newDirectory)
 save(newDirectory+'/params.mat','params')
 save(newDirectory+'/simArray.mat','simArray', '-v7.3')
 save(newDirectory+'/rxnsPerTimestep.mat','rxnsPerTimestep', '-v7.3')

 %% Parameter sweep: V (using maternal IC)
 close all

 params = struct;
 params.numSims = 10; % For the paper, numSims = 10
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
 params.hitRun = 0;

 sweepDir = 'path to folder to save V sweep results'+string(datetime);
 mkdir(sweepDir)

 V_vals = 0.0004:0.00001:0.0012;
 simArrays = cell(length(V_vals),1);
 for i = 1:length(V_vals)

     params.V = V_vals(i); % for the paper

     % Maternal IC, dynamic Pho
     params.fiberIC = [2*ones(2,params.nNucleosomes*0.2) 3*ones(2,params.nNucleosomes*0.8)]; % maternal IC
     [simArray rxnsPerTimestep] = simulateK27(ez,pho,params);
     plotSimulation(simArray,params,'maternal, V = '+string(V_vals(i)))

     simArrays{i} = simArray;

 end

 save(strcat(sweepDir,'/simArrays.mat'),'simArrays','-v7.3')
 save(strcat(sweepDir,'/V_vals.mat'),'V_vals')
 save(strcat(sweepDir,'/params.mat'),'params')

 %% Functions

function plotSimulation(simArray,params,titleString)

    tMit = params.tMit;
    stepsPerMin = params.stepsPerMin;
    numSims = params.numSims;
    nNucleosomes = params.nNucleosomes;
    nUpdates = params.nUpdates;
    
    % Find the fractions of each methyl group across all simulations
    fracZeroOverT = sum(simArray==0,2:4)./(nNucleosomes*numSims*2);
    fracOneOverT = sum(simArray==1,2:4)./(nNucleosomes*numSims*2);
    fracTwoOverT = sum(simArray==2,2:4)./(nNucleosomes*numSims*2);
    fracThreeOverT = sum(simArray==3,2:4)./(nNucleosomes*numSims*2);
    
    f = figure;
    plot(fracZeroOverT,'linewidth',2,'color',[0.75 0.75 0.75])
    hold on
    plot(fracOneOverT,'linewidth',2,'color',[15, 171, 26]./255)
    hold on
    plot(fracTwoOverT,'linewidth',2,'color','m')
    hold on
    plot(fracThreeOverT,'linewidth',2,'color',[1 0.6 0])
    hold on
    plot(fracZeroOverT+fracOneOverT+fracTwoOverT+fracThreeOverT)
    
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
    
%     l = legend('K27me0','K27me1','K27me3');
%     l.FontSize = 14;
%     l.Location = 'Northwest';
%     l.Box = 'off';
    
    set(f,'units','centimeters','position',[35,25,11.25,4.5])
    % exportgraphics(f,'path to figure')

end
