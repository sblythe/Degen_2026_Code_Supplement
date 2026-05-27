function [meState rxn] = methylationFunction(meState,V,allostery,p_bind)
    
    if meState == 0
    
       k = 13*V*allostery*p_bind; 
       
    elseif meState == 1
       
       k = 4*V*allostery*p_bind;
    
    elseif meState == 2
    
       k = V*allostery*p_bind;
    
    elseif meState >= 3
    
       k = 0;
    
    end
    
    r = rand;
    if r <= k
    
    meState = meState + 1;
    rxn = 1;
    
    else
    
    rxn = 0;
    
    end

end