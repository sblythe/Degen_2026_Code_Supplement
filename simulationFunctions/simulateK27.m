function [simArray rxnsPerTimestep] = simulateK27(ez,pho,params)
 tic 
 rng(1)

 numSims = params.numSims;
 nUpdates = params.nUpdates;
 Ke = params.Ke;
 Kp = params.Kp;
 V = params.V;
 deltaT = params.deltaT;
 tRep = params.tRep;
 nNucleosomes = params.nNucleosomes;

 outputs = cell(numSims,1);
 h = waitbar(0,'Simulating...');

 rxnsPerTimestep = zeros(numSims, sum(params.NClengths));
 for n = 1:numSims
     
     % initialize fiber
     fiber = params.fiberIC;
     fiber(1,:) = fiber(1,randperm(size(fiber,2))); % Randomly columns in first row
     fiber(2,:) = fiber(2,randperm(size(fiber,2))); % Randomly columns in second row
     bothPositions = nan(nUpdates,nNucleosomes,2);
     bothPositions(1,:,1) = fiber(1,:);
     bothPositions(1,:,2) = fiber(2,:);

     t = 1;
     while t < nUpdates % loop through time 

        % Are Pho and E(z) bound to the locus?
        prob_e = ez(t)./(Ke+ez(t));
        prob_p = pho(t)./(Kp+pho(t));

        prob_ep = prob_e * prob_p;

        % Inspect all nucs in random order
        index = randperm(length(fiber));
        rxns = 0; % number of reactions in timestep starts at 0
        for i = 1:nNucleosomes
       
             nuc = fiber(:,index(i));
   
        
             % Allosteric stimulation? Check the me positions at each
             % of the flanking nucleosomes
             if index(i) == 1

                 leftNuc = [];
                 rightNuc = fiber(:,index(i)+1);

             elseif index(i) == length(fiber)

                 rightNuc = [];
                 leftNuc = fiber(:,index(i)-1);

             else

                 leftNuc = fiber(:,index(i)-1);
                 rightNuc = fiber(:,index(i)+1);

             end

             meOrder = randperm(2,2); % Randomly select the me position order 
             for m = 1:2

                 rxn = 0;

                 isDiORTri = [leftNuc, rightNuc] >= 2;
                 if sum(isDiORTri(:)) > 0  % If at least one of the neighbors is dimethylated, there's allostery

                     % If you want different allosteric parameter values
                     % for me3 and me2 (based on Margueron
                     % https://pmc.ncbi.nlm.nih.gov/articles/PMC3772642/#F4)
                     if ~isempty(find([leftNuc, rightNuc]==3,1)) 
                        params.allostery = params.allo_3;
                     elseif ~isempty(find([leftNuc, rightNuc]==2,1)) 
                        params.allostery = params.allo_2;
                     end

                     [meState, rxn] = methylationFunction(fiber(meOrder(m),i),V,params.allostery,prob_e);

                 end

                 if rxn == 0 % If Pho and E(z) are bound, there's nucleation

                     if params.hitRun == 1

                         [meState, rxn] = methylationFunction(fiber(meOrder(m),i),V,1,prob_e);
                     
                     else

                        [meState, rxn] = methylationFunction(fiber(meOrder(m),i),V,1,prob_ep); % NOTE should have prob_ep

                     end

                 end

                 if rxn == 0

                     % If there was no reaction, don't update the me state
                     meState = fiber(meOrder(m),i);

                 end
  
                 fiber(meOrder(m),i) = meState;

                 rxns = rxns + rxn; % count of reactions in the timestep
             end  

            % Demethylate? If at replication time
            checkRep = tRep-t;
            if sum(checkRep == 0) == 1

                r = rand;
                if r <= 0.5 % 50% chance of losing all methylation if DNA replication

                   fiber(:,i) = [0; 0];

                end

            end
        end

        t = t + deltaT;

        rxnsPerTimestep(n,t) = rxns;
          
        bothPositions(t,:,:) = fiber';

     end

     outputs{n} = bothPositions;
    
     waitbar(n/numSims)
 end

 close(h)
 toc

 % reformat cell outputs to a 4D array: time, nucleosomes, position, replicate
 simArray = cat(4,outputs{:});

end
